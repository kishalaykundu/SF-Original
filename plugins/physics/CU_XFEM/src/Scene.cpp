/**
 * @file Scene.cpp
 * @author Kishalay Kundu <kishalay.kundu@gmail.com>
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 * The scene class for the CU_XFEM library. This handles intersection
 * between XFE meshes and blade
 */

#include <ctime>
#include <forward_list>
#include <boost/bind.hpp>

#include "ThreadControl.h"
#include "aabb.h"

#include "Rigid/inc/Mesh.h"

#include "Scene.h"
#include "Submesh.h"

using namespace std;
using namespace boost;
using namespace boost::threadpool;

/*********************** POOLJOB RELATED METHODS ***********************/
namespace SF {
  namespace XFE {

    // default constructor
    PoolJob::PoolJob ()
    : _submesh (NULL), _partitionIndex (0), _bladeCurr (NULL), _bladePrev (NULL), _bladeIndices (NULL)
    {
      for (unsigned int i = 0; i < 2; ++i){
        _bladeNormals [i] = NULL;
      }
    }
    // destructor
    PoolJob::~PoolJob () { }

    // overloaded constructor
    PoolJob::PoolJob (Submesh *sm, unsigned int pIndex, vector <vec> **bCurr, vector <vec> **bPrev, vector <vec> bNormals [2], vector <unsigned int> *bIndices)
    : _submesh (sm), _partitionIndex (pIndex), _bladeCurr (bCurr), _bladePrev (bPrev), _bladeIndices (bIndices)
    {
      for (unsigned int i = 0; i < 2; ++i){
        _bladeNormals [i] = &bNormals [i];
      }
    }

    // gathering method for cells affected by blades
    void
    PoolJob::getAffectedCells ()
    {
      _submesh->getAffectedCells (_partitionIndex, **_bladeCurr, **_bladePrev, *_bladeIndices, _bladeNormals);
    }

    // method to resolve faces
    void
    PoolJob::resolveFaces ()
    {
      _submesh->resolveFaces ();
    }

    // finalizing method for cells
    void
    PoolJob::finalizeCollision ()
    {
      _submesh->finalizeCollision (_partitionIndex, **_bladeCurr, **_bladePrev, *_bladeIndices, _bladeNormals);
    }
  }
}

/*********************** SCENE RELATED METHODS ***********************/
namespace SF {
  namespace XFE {

    const unsigned int TICKS = 5000;
    bool cflag = false;
    clock_t before, before1;
    unsigned int fCounter = 0;
    clock_t p13, p2, p4, tmpc;

    // default constructor
    Scene::Scene ()
    :_bladeWaitIndex (-1), _bladePostIndex (-1), _bladeCurr (&( _bladeVerts [0])), _bladePrev (&(_bladeVerts [1]))
    { }

    // destructor
    Scene::~Scene ()
    {
      _pool.clear ();
    }

    // method to add blade-related variables
    void
    Scene::addBlade (Resource *r)
    {
      SF::RM::Mesh *m = static_cast <SF::RM::Mesh *> (r);

      _bladeSyncControl = &(m->_syncControl);
      _bladeWaitIndex = m->_semIntersectionWaitIndex;
      _bladePostIndex = m->_semIntersectionPostIndex;

      _bladeVerts [0] = *(m->_bladeCurr);
      delete m->_bladeCurr;
      m->_bladeCurr = &(_bladeVerts [0]);
      _bladeCurr = &(_bladeVerts [0]);

      _bladeVerts [1] = *(m->_bladePrev);
      delete m->_bladePrev;
      m->_bladePrev = &(_bladeVerts [1]);
      _bladePrev = &(_bladeVerts [1]);

      updateBladeBounds ();

      _bladeIndices = *(m->_bladeIndices);
      delete m->_bladeIndices;
      m->_bladeIndices = &(_bladeIndices);

      _bladeNormals [0].resize (_bladeIndices.size ()/ 2);
      _bladeNormals [1].resize (_bladeIndices.size ()/ 2);
    }

    // intersection detection and resolution method
    void
    Scene::run ()
    {
      // local variables
      Mesh *m;
      vec e1, e2;
      Submesh *sm = NULL;
      unsigned int i1, i2;
      bool normalComputeFlag;
      forward_list <unsigned int>::iterator iter, biter;

      // push all possible jobs on to a queue
      vector <unsigned int> jobOffsets (_mesh.size());
      vector <boost::shared_ptr <SF::XFE::PoolJob> > collisionJobs;
      for (unsigned int i = 0; i < _mesh.size (); ++i){
        m = _mesh [i].get ();
        jobOffsets [i] = 0;
        for (unsigned int j = 0; j < m->_submesh.size (); ++j){
          sm = m->_submesh [j].get ();
          for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
            collisionJobs.push_back (boost::shared_ptr <SF::XFE::PoolJob> (new PoolJob (sm, k, &_bladeCurr, &_bladePrev, _bladeNormals, &_bladeIndices)));
          }
          if (i < _mesh.size () - 1){
            jobOffsets [i + 1] += sm->_partitions.size ();
          }
        }
      }
      for (unsigned int i = 1; i < jobOffsets.size (); ++i){
        jobOffsets [i] += jobOffsets [i - 1];
      }

      while (true){

        // lock blade
        (*_bladeSyncControl) [_bladeWaitIndex].wait ();

				vector <vec> *tmpp = _bladeCurr;
				_bladeCurr = _bladePrev;
				_bladePrev = tmpp;

        updateBladeBounds ();

        normalComputeFlag = true;

        for (unsigned int i = 0; i < _mesh.size (); ++i){
          // lock mesh
          m = _mesh [i].get ();
          m->_syncControl [m->_semIntersectionWaitIndex].wait ();

          if (_bladeBounds.collide (m->_bbox)){

if (cflag){
  before = clock ();
  if (!fCounter){
    p13 = p2 = p4 = before - before;
  }
  ++fCounter;
}
            // for all submeshes, gather a list of all cells that collide with blade
            for (unsigned int j = 0; j < m->_submesh.size (); ++j){
              sm = m->_submesh [j].get ();

              // push affected cell gathering tasks to threadpool
              if (_bladeBounds.collide (sm->_bbox)){

                // process every partition for new cuts to cells
                for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                  if (_bladeBounds.collide (sm->_partitions [k]._bbox)){

                    // compute blade normals (done once per frame if any submesh partition collides)
                    if (normalComputeFlag){
if (cflag){
  before1 = clock ();
}
                      for (unsigned int l = 0; l < _bladeNormals [0].size (); ++l){
                        i1 = _bladeIndices [2*l];
                        i2 = _bladeIndices [2*l + 1];

                        e1 = (*_bladeCurr) [i2] - (*_bladeCurr) [i1];
                        e2 = (*_bladePrev) [i2] - (*_bladeCurr) [i1];
                        e1.fast_cross (_bladeNormals [0][l], e2);

                        e1 = (*_bladePrev) [i1] - (*_bladePrev) [i2];
                        e2 = (*_bladeCurr) [i1] - (*_bladePrev) [i2];
                        e1.fast_cross (_bladeNormals [1][l], e2);
                      }
                      normalComputeFlag = false;
                    }
if (cflag){
  tmpc = clock () - before1;
  p13 += tmpc;
}
                    if (!fCounter){
                      cflag = true;
                    }
                    // push collision detection and resolution tasks to the threadpool
                    schedule (_pool, bind (&SF::XFE::PoolJob::getAffectedCells, collisionJobs [jobOffsets [i] + j*sm->_partitions.size () +k]));

                  } // end - if (_bladeBounds.collide (sm->_partitions [k]._bbox))
                } // end - for (unsigned int k = 0; k < sm->_partitions.size (); ++i)
              } // end - if (_bladeBounds.collide (sm->_bbox))
            } // end - for (unsigned int j = 0; j < m->_submesh.size (); ++j)

            // synchronize
            _pool.wait ();
if (cflag){
  p2 += clock () - before;
  p2 -= tmpc;
  before1 = clock ();
}
            /** For all submeshes:
            * 1. Shuffle any cells that are in partitions that don't own them
            * 2. Update any face structure that needs to be fixed
            */
            for (unsigned int j = 0; j < m->_submesh.size (); ++j){
              sm = m->_submesh [j].get ();

              bool *shuffleFlags = new bool [sm->_partitions.size ()];

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                shuffleFlags [k] = false;
              }

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                if (!sm->_partitions [k]._cutCells.empty ()){

                  for (iter = sm->_partitions [k]._cutCells.begin (), biter = sm->_partitions [k]._cutCells.before_begin (); iter != sm->_partitions [k]._cutCells.end (); ++iter){
                    if (*iter < sm->_partitions [k]._cellStartIndex){
                      for (int l = k - 1; l >= 0; --l){
                        if (*iter >= sm->_partitions [l]._cellStartIndex){
                          sm->_partitions [l]._cutCells.push_front (*iter);
                          shuffleFlags [l] = true;
                          break;
                        }
                      }
                      if (iter == sm->_partitions [k]._cutCells.begin ()){
                        sm->_partitions [k]._cutCells.pop_front ();
                      } else {
                        sm->_partitions [k]._cutCells.erase_after (biter);
                      }
                      iter = biter;
                    } else if (*iter > sm->_partitions [k]._cellEndIndex){
                      for (unsigned int l = k + 1; l < sm->_partitions.size (); ++l){
                        if (*iter <= sm->_partitions [l]._cellEndIndex){
                          sm->_partitions [l]._cutCells.push_front (*iter);
                          shuffleFlags [l] = true;
                          break;
                        }
                      }
                      if (iter == sm->_partitions [k]._cutCells.begin ()){
                        sm->_partitions [k]._cutCells.pop_front ();
                      } else {
                        sm->_partitions [k]._cutCells.erase_after (biter);
                      }
                      iter = biter;
                    } else {
                      ++biter;
                    }
                  }
                } // end - if (!sm->_partitions [k]._cutCells.empty ())
              } // end - for (unsigned int k = 0; k < sm->_partitions.size (); ++i)

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                if (shuffleFlags [k]){
                  sm->_partitions [k]._cutCells.sort ();
                  sm->_partitions [k]._cutCells.unique ();
                }
              }

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                shuffleFlags [k] = false;
              }

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){

                if (!sm->_partitions [k]._reExaminedCells.empty ()){

                  for (iter = sm->_partitions [k]._reExaminedCells.begin (), biter = sm->_partitions [k]._reExaminedCells.before_begin ();
                    iter != sm->_partitions [k]._reExaminedCells.end (); ++iter){
                    if (*iter < sm->_partitions [k]._cellStartIndex){
                      for (int l = k - 1; l >= 0; --l){
                        if (*iter >= sm->_partitions [l]._cellStartIndex){
                          sm->_partitions [l]._reExaminedCells.push_front (*iter);
                          shuffleFlags [l] = true;
                          break;
                        }
                      }
                      if (iter == sm->_partitions [k]._reExaminedCells.begin ()){
                        sm->_partitions [k]._reExaminedCells.pop_front ();
                      } else {
                        sm->_partitions [k]._reExaminedCells.erase_after (biter);
                      }
                      iter = biter;
                    }
                    else if (*iter > sm->_partitions [k]._cellEndIndex){
                      for (unsigned int l = k + 1; l < sm->_partitions.size (); ++l){
                        if (*iter <= sm->_partitions [l]._cellEndIndex){
                          sm->_partitions [l]._reExaminedCells.push_front (*iter);
                          shuffleFlags [l] = true;
                          break;
                        }
                      }
                      if (iter == sm->_partitions [k]._reExaminedCells.begin ()){
                        sm->_partitions [k]._reExaminedCells.pop_front ();
                      } else {
                        sm->_partitions [k]._reExaminedCells.erase_after (biter);
                      }
                      iter = biter;
                    } else {
                      ++biter;
                    }
                  }
                } // end - if (!sm->_partitions [k]._reExaminedCells.empty ())
              } // end - for (unsigned int k = 0; k < sm->_partitions.size (); ++i)

              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                if (shuffleFlags [k]){
                  sm->_partitions [k]._reExaminedCells.sort ();
                  sm->_partitions [k]._reExaminedCells.unique ();
                }
              }

              delete [] shuffleFlags;

            } // end - for (unsigned int j = 0; j < m->_submesh.size (); ++j)
            for (unsigned int j = 0; j < m->_submesh.size (); ++j){
              sm = m->_submesh [j].get ();
              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                if (!sm->_partitions [k]._cutCells.empty () || !sm->_partitions [k]._reExaminedCells.empty ()){
                  schedule (_pool, bind (&SF::XFE::PoolJob::resolveFaces, collisionJobs [jobOffsets [i] + j*sm->_partitions.size () +k]));
                  break;
                }
              }
            }
            _pool.wait ();

            // rejiggle any vertex that is too near the blade indices
            m->adjustVertices (*_bladeCurr, *_bladePrev, _bladeIndices, _bladeNormals [0], _bladeNormals [1]);
if (cflag){
  p13 += clock () - before1;
  before = clock ();
}
            // generate final triangles
            for (unsigned int j = 0; j < m->_submesh.size (); ++j){
              sm = m->_submesh [j].get ();

              // process every partition for new cuts to cells
              for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
                if (!sm->_partitions [k]._cutCells.empty () || !sm->_partitions [k]._reExaminedCells.empty ()){

                  // push collision detection and resolution tasks to the threadpool
                  schedule (_pool, bind (&SF::XFE::PoolJob::finalizeCollision, collisionJobs [jobOffsets [i] + j*sm->_partitions.size () +k]));

                } // end - (!sm->_partitions [k]._cutCells.empty () || !sm->_partitions [k]._reExaminedCells.empty ())
              } // end - for (unsigned int k = 0; k < sm->_partitions.size (); ++i)
            } // end - for (unsigned int j = 0; j < m->_submesh.size (); ++j)

            // synchronize
            _pool.wait ();
if (cflag){
  p4 += clock () - before;
  if (fCounter >= TICKS){
    cflag = false;
    fprintf (stdout, "\nTimes: 13: %g\t 2: %g\t 4: %g\n", static_cast <double> (p13) * 1000./(TICKS * CLOCKS_PER_SEC),
             static_cast <double> (p2) * 1000./(TICKS * CLOCKS_PER_SEC),
             static_cast <double> (p4) * 1000./(TICKS * CLOCKS_PER_SEC));
		unsigned int cellCount = 0;
		for (unsigned int j = 0; j < m->_submesh.size (); ++j){
    	for (unsigned int k = 0; k < sm->_partitions.size (); ++k){
				for (forward_list <unsigned int>::iterator iter = sm->_partitions [k]._cutCells.begin ();
							iter != sm->_partitions [k]._cutCells.end (); ++iter){
					++cellCount;
				}
				for (forward_list <unsigned int>::iterator iter = sm->_partitions [k]._finishedCells.begin ();
							iter != sm->_partitions [k]._finishedCells.end (); ++iter){
					++cellCount;
				}
			}
		}
		fprintf (stdout, "Total Affected Cells: %u\n", cellCount);
  }
}
          } // end - if (_bladeBounds.collide (m->_bbox))

          // release mesh
          m->_syncControl [m->_semIntersectionPostIndex].post ();

        } // end - for (unsigned int i = 0; i < _mesh.size (); ++i)

        // release blade
        (*_bladeSyncControl) [_bladePostIndex].post ();
      } // end - while (true)
    }

  }
}
