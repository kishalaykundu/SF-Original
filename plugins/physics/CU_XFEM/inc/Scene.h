/**
 * @file Scene.h
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

#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/threadpool.hpp>

#include "Preprocess.h"

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

#include "Common.h"
#include "Partition.h"
#include "Submesh.h"
#include "Mesh.h"

using namespace std;
using namespace boost;
using namespace boost::threadpool;

namespace SF {

    class ThreadControl;
    class aabb;

  namespace XFE {

    class PoolJob {
      public:
        Submesh *_submesh;
        unsigned int _partitionIndex;
        vector <vec> **_bladeCurr;
        vector <vec> **_bladePrev;
        vector <vec> *_bladeNormals [2];
        vector <unsigned int> *_bladeIndices;

      public:
        PoolJob ();
        ~PoolJob ();
        PoolJob (Submesh *sm, unsigned int pIndex, vector <vec> **bCurr, vector <vec> **bPrev, vector <vec> bNormals [2], vector <unsigned int> *bIndices);

        void getAffectedCells ();
        void resolveFaces ();
        void finalizeCollision ();
    };

    class Scene {

      private:
        // threads
        pool _pool;

        // meshes
        vector <boost::shared_ptr <Mesh> > _mesh;

        // cutting tool
        int _bladeWaitIndex;
        int _bladePostIndex;
        ThreadControl *_bladeSyncControl;

        aabb _bladeBounds;
        vector <vec> *_bladeCurr;
        vector <vec> *_bladePrev;
        vector <vec> _bladeVerts [2];
        vector <vec> _bladeNormals [2];

        vector <unsigned int> _bladeIndices;

      public:
        Scene ();
        ~Scene ();

        inline void resizePool (unsigned int n)
        {
          _pool.size_controller ().resize (n);
        }
        inline void addMesh (Resource &r)
        {
          Mesh *m = dynamic_cast <Mesh *> (&r);

          // register certain variables from mesh to submeshes
          for (unsigned int i = 0; i < m->_submesh.size (); ++i){
            m->_submesh [i].get ()->_changeBit = &(m->_faceChangeBits [i]);
            m->_submesh [i].get ()->_vertexInfo = &(m->_vertexInfo);
          }

          _mesh.push_back (boost::shared_ptr <Mesh> (m));
        }
        void addBlade (Resource *r);

        void run ();

      private:
        // private method to update bounding box for blade
        inline void updateBladeBounds ()
        {
          for (int i = 0; i < 3; ++i){
            _bladeBounds._v [0]._v [i] = (*_bladeCurr) [0]._v [i];
          }
          for (int i = 0; i < 3; ++i){
            _bladeBounds._v [1]._v [i] = _bladeBounds._v [0]._v [i];
          }
          for (unsigned int i = 1; i < _bladeCurr->size (); ++i){
            for (int j = 0; j < 3; ++j){
              if (_bladeBounds._v [0]._v [j] > (*_bladeCurr) [i]._v[j]){
                _bladeBounds._v [0]._v [j] = (*_bladeCurr) [i]._v[j];
              } else if (_bladeBounds._v [1]._v [j] < (*_bladeCurr) [i]._v[j]){
                _bladeBounds._v [1]._v [j] = (*_bladeCurr) [i]._v[j];
              }
            }
          }
          for (unsigned int i = 1; i < _bladePrev->size (); ++i){
            for (int j = 0; j < 3; ++j){
              if (_bladeBounds._v [0]._v [j] > (*_bladePrev) [i]._v[j]){
                _bladeBounds._v [0]._v [j] = (*_bladePrev) [i]._v[j];
              } else if (_bladeBounds._v [1]._v [j] < (*_bladePrev) [i]._v[j]){
                _bladeBounds._v [1]._v [j] = (*_bladePrev) [i]._v[j];
              }
            }
          }
        }
    };

  }
}
