/**
 * @file Partition.cpp
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
 * The partition class for the CU_XFEM library. It is a subcontainer
 * of data within the submesh class
 */


#include "Collide/triTriCollide.h"
#include "Collide/lineTriCollide.h"

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "Intersect.h"
#include "Partition.h"

namespace SF {
  namespace XFE {

    static const real CUT_DISTANCE = 0.01;

    // static method to compute barycentric co-ordinates of a point inside a triangle
    inline void calculateBarycentricCoords (vec2 &uv, vec &p, vec &a, vec &b, vec &c)
    {
      vec v0 = c - a;
      vec v1 = b - a;
      vec v2 = p - a;

      real d00 = v0.dot (v0);
      real d01 = v0.dot (v1);
      real d02 = v0.dot (v2);
      real d11 = v1.dot (v1);
      real d12 = v1.dot (v2);

      real id = 1./ (d00 * d11 - d01 * d01);
      uv._v [0] = (d11 * d02 - d01 * d12) * id;
      uv._v [1] = (d00 * d12 - d01 * d02) * id;
    }

    // default constructor
    Partition::Partition ()
    :_cellStartIndex (0), _cellEndIndex (0), _exFaceStartIndex (0), _exFaceEndIndex (0), _inFaceStartIndex (1), _inFaceEndIndex (0),
    _vertInfo (NULL), _tex2D (NULL), _tex3D (NULL), _exVertices (NULL), _exUVCoords (NULL), _ex2DTexCoords (NULL), _exFaceIndices (NULL),
    _inVertices (NULL), _inUVCoords (NULL), _inSurfaceVertexStatus (NULL), _in2DTexCoords (NULL), _in3DTexCoords (NULL), _inFaceIndices (NULL)
    { }

    // destructor
    Partition::~Partition () { }

    // copy constructor
    Partition::Partition (const Partition &p)
    : _bbox (p._bbox),
    _cellStartIndex (p._cellStartIndex), _cellEndIndex (p._cellEndIndex),
    _exFaceStartIndex (p._exFaceStartIndex), _exFaceEndIndex (p._exFaceEndIndex),
    _inFaceStartIndex (p._inFaceStartIndex), _inFaceEndIndex (p._inFaceEndIndex),
    _cutCells (p._cutCells), _reExaminedCells (p._reExaminedCells), _finishedCells (p._finishedCells),
    _collidingVertices (p._collidingVertices), _cuts (p._cuts), _vertInfo (p._vertInfo), _tex2D (p._tex2D), _tex3D (p._tex3D),
    _exVertices (p._exVertices), _exUVCoords (p._exUVCoords), _ex2DTexCoords (p._ex2DTexCoords), _exFaceIndices (p._exFaceIndices),
    _inVertices (p._inVertices), _inUVCoords (p._inUVCoords), _inSurfaceVertexStatus (p._inSurfaceVertexStatus),
    _in2DTexCoords (p._in2DTexCoords), _in3DTexCoords (p._in3DTexCoords), _inFaceIndices (p._inFaceIndices),
    _inEmptyVertices (p._inEmptyVertices), _inEmptyFaces (p._inEmptyFaces), _exEmptyVertices (p._exEmptyVertices), _exEmptyFaces (p._exEmptyFaces)
    { }

    // assignment operator
    Partition &
    Partition::operator = (const Partition &p)
    {
      _bbox = p._bbox;
      _cellStartIndex = p._cellStartIndex;
      _cellEndIndex = p._cellEndIndex;
      _exFaceStartIndex = p._exFaceStartIndex;
      _exFaceEndIndex = p._exFaceEndIndex;
      _inFaceStartIndex = p._inFaceStartIndex;
      _inFaceEndIndex = p._inFaceEndIndex;
      _cutCells = p._cutCells;
      _reExaminedCells = p._reExaminedCells;
      _finishedCells = p._finishedCells;
      _collidingVertices = p._collidingVertices;
      _cuts = p._cuts;
      _vertInfo = p._vertInfo;
      _tex2D = p._tex2D;
      _tex3D = p._tex3D;
      _exVertices = p._exVertices;
      _exUVCoords = p._exUVCoords;
      _ex2DTexCoords = p._ex2DTexCoords;
      _exFaceIndices = p._exFaceIndices;
      _inVertices = p._inVertices;
      _inUVCoords = p._inUVCoords;
      _inSurfaceVertexStatus = p._inSurfaceVertexStatus;
      _in2DTexCoords = p._in2DTexCoords;
      _in3DTexCoords = p._in3DTexCoords;
      _inFaceIndices = p._inFaceIndices;
      _inEmptyVertices = p._inEmptyVertices;
      _inEmptyFaces = p._inEmptyFaces;
      _exEmptyVertices = p._exEmptyVertices;
      _exEmptyFaces = p._exEmptyFaces;

      return *this;
    }

    // method to get affected cells
    void
    Partition::gatherAffectedCells (unsigned int sIndex, vector <Vertex> &vertexInfo, vector <vec> &verts, vector <unsigned int> &indices, vector <Face> &faces,
                                    vector <unsigned int> &iindices, vector <Face> &ifaces, vector <Edge> &edges, vector <Cell> &cells,
                                    vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      // look for new surface cuts
      bool collideFlag = false;
      vec e1, e2, normal;

      // examine external faces
      for (unsigned int i = 3*_exFaceStartIndex; i <= 3*_exFaceEndIndex; i += 3){

        // only examine non-degenerate triangles
        if (faces [i/ 3]._owner < UINT_MAX){

          collideFlag = false;
          e1 = verts [indices [i + 1]] - verts [indices [i]];
          e2 = verts [indices [i + 2]] - verts [indices [i]];
          e1.fast_ncross (normal, e2);

          // test for collisions with external triangles
          for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
            if (triTriCollide ((*(bladeNormals [0])) [j], bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]],
                               normal, verts [indices [i]], verts [indices [i + 1]], verts [indices [i + 2]], e1)){
              _cutCells.push_front (faces [i/ 3]._owner);
              collideFlag = true;
              break;
            }
          }
          if (!collideFlag){
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (triTriCollide ((*(bladeNormals [1])) [j], bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]],
                                 normal, verts [indices [i]], verts [indices [i + 1]], verts [indices [i + 2]], e1)){
                _cutCells.push_front (faces [i/ 3]._owner);
                collideFlag = true;
                break;
              }
            }
          }

          if (collideFlag){
            for (unsigned int j = 0; j < 3; ++j){
              indices [i + j] = 0;
            }
          }

        } // end - if (_faces [i]._owner < UINT_MAX)
      } // end - for (unsigned int i = 3*_exFaceStartIndex; i <= 3*_exFaceEndIndex; i += 3)

      // examine inside faces
      for (unsigned int i = 3*_inFaceStartIndex; i <= 3*_inFaceEndIndex; i += 3){

        // only examine non-degenerate triangles
        if (ifaces [i/ 3]._owner < UINT_MAX){

          collideFlag = false;
          e1 = verts [iindices [i + 1]] - verts [iindices [i]];
          e2 = verts [iindices [i + 2]] - verts [iindices [i]];
          e1.fast_ncross (normal, e2);

          // test for collisions with internal triangles
          for (unsigned int j = 0; j < bladeIndices.size ()/ 2; ++j){
            if (triTriCollide ((*(bladeNormals [0])) [j], bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]],
                               normal, verts [iindices [i]], verts [iindices [i + 1]], verts [iindices [i + 2]], e1)){
              _cutCells.push_front (ifaces [i/ 3]._owner);
              ifaces [i/ 3]._owner = UINT_MAX;
              collideFlag = true;
              break;
            }
          }
          if (!collideFlag){
            for (unsigned int j = 0; j < bladeIndices.size ()/ 2; ++j){
              if (triTriCollide ((*(bladeNormals [1])) [j], bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]],
                                 normal, verts [iindices [i]], verts [iindices [i + 1]], verts [iindices [i + 2]], e1)){
                _cutCells.push_front (ifaces [i/ 3]._owner);
                ifaces [i/ 3]._owner = UINT_MAX;
                collideFlag = true;
                break;
              }
            }
          }
        } // end - if (ifaces [i]._owner < UINT_MAX)
      } // end - for (unsigned int i = 3*_inFaceStartIndex; i <= 3*_inFaceEndIndex; i += 3)

      // quick check to return
      if (_cutCells.empty () && _reExaminedCells.empty ()){
        return;
      }

      // check for vertex collisions for tagged tetrahedra
      forward_list <unsigned int>::iterator biter = _cutCells.before_begin ();

      // resolve newly tagged tetrahedra (list in partition)
      bool reshuffleFlag = false;
      forward_list <unsigned int>::iterator iter;
      forward_list <unsigned int>::iterator lastElement = _cutCells.before_begin ();
      for (iter = _cutCells.begin (); iter != _cutCells.end (); ++iter){
        ++lastElement;
      }
      for (iter = _cutCells.begin (); iter != _cutCells.end (); ++iter){
        if (!cells [*iter].testCellExamFlag ()){
          cellBladeCollide (sIndex, lastElement, vertexInfo, verts, edges, cells, cells [*iter], bladeCurr, bladePrev, bladeIndices, bladeNormals);
          reshuffleFlag |= cells [*iter].testAnyCollisionFlag ();
        }
      }

      _cutCells.sort ();
      _cutCells.unique ();

      // early return if no reshuffling to re-examination list required
      if (!reshuffleFlag){
        return;
      }

      _collidingVertices.sort ();
      _collidingVertices.unique ();

      // push any cell who need to re-examined later to its respective list
      biter = _cutCells.before_begin ();
      for (iter = _cutCells.begin (); iter != _cutCells.end (); ++iter){
        if (cells [*iter].testAnyCollisionFlag ()){
          _reExaminedCells.push_front (*iter);
          iter = biter;
          _cutCells.erase_after (biter);
        } else {
          ++biter;
        }
      }

    }

    // method to do the finalising of all cells
    void
    Partition::finalizeCollision (vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells,
                                  vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      // check edges cells that are in the re-examination queue
      if (!_reExaminedCells.empty ()){
        resolveReExaminedCells (verts, edges, cells, bladeCurr, bladePrev, bladeIndices, bladeNormals);
      }

      // examine each cut cell and finalize triangulation for display
      unsigned int index;
      forward_list <unsigned int>::iterator iter, biter = _cutCells.before_begin ();
      for (iter = _cutCells.begin (); iter != _cutCells.end (); ++iter){

        index = *iter;

        // reset all the flags for next round
        for (unsigned int i = 0; i < 4; ++i){
          (* _vertInfo) [cells [index]._index [i]].reset ();
        }
        for (unsigned int i = 0; i < 6; ++i){
          edges [cells [index]._edgeIndex [i]].reset ();
        }
        cells [index].reset ();

        // update cut-info structure for unfinished cells
        if (cells [index]._cutIndex < 0){
          cells [index]._cutIndex = _cuts.size ();
          _cuts.push_back (Cut ());
        }

        formFaces (cells [index], _cuts [cells [index]._cutIndex], edges, verts, bladeCurr, bladePrev, bladeIndices, bladeNormals);

        // remove completely cut cell from cut-cell list and put it in finished-cell list
        if (cells [index].testCellFinalizeFlag ()){
          _finishedCells.push_front (index);
          iter = biter;
          _cutCells.erase_after (iter);
        } else {
          ++biter;
        }
      }

      // populate vertex information for all finished cells (unfinished cells are already done)
      vec *cellVerts [4];
      real uv [3];
      unsigned int *vertIndexArray, *uvIndexArray;
      for (iter = _finishedCells.begin (); iter != _finishedCells.end (); ++iter){

        index = *iter;
        for (unsigned int i = 0; i < 4; ++i){
          cellVerts [i] = &(verts [cells [index]._index [i]]);
        }

        index = cells [index]._cutIndex;
        vertIndexArray = _cuts [index]._inVertices;
        uvIndexArray = _cuts [index]._inUVCoords;
        for (unsigned int i = 0; i < _cuts [index]._numInVertices; ++i){
          memcpy (uv, (* _inUVCoords) [uvIndexArray [i]]._v, 3 * sizeof (real));
          (* _inVertices ) [vertIndexArray [i]] = *(cellVerts [0]) * uv[0] + *(cellVerts [1]) *uv [1] + *(cellVerts [2]) * uv [2] +
            *(cellVerts [3]) * (1. - uv [0] - uv [1] - uv [2]);
        }
      }

    }

    // private method to resolve cells earmarked for re-examination
    void
    Partition::resolveReExaminedCells (vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells,
                                       vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      for (forward_list <unsigned int>::iterator iter = _reExaminedCells.begin (); iter != _reExaminedCells.end (); ++iter){
        for (unsigned int i = 0; i < 6; ++i){
          edges [cells [*iter]._edgeIndex [i]].reset ();
        }
      }
      unsigned int index, cIndex;
      while (!_reExaminedCells.empty ()){
        cIndex = _reExaminedCells.front ();

        for (unsigned int i = 0; i < 6; ++i){
          index = cells [cIndex]._edgeIndex [i];
          switch (i){
            case 0: // edge 0-1
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [1]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [1]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
            case 1: // edge 0-2
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [2]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [2]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
            case 2: // edge 0-3
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [3]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [0]], verts [cells [cIndex]._index [3]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
            case 3: // edge 1-2
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [1]], verts [cells [cIndex]._index [2]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [1]], verts [cells [cIndex]._index [2]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
            case 4: // edge 1-3
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [1]], verts [cells [cIndex]._index [3]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [1]], verts [cells [cIndex]._index [3]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
            case 5: // edge 2-3
              if (edges [index].testCollisionFlag ()){
                break;
              }
              edges [index].setCollisionFlag ();
              for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
                if (lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [2]], verts [cells [cIndex]._index [3]],
                                    bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                    lineTriCollide (edges [index]._u, verts [cells [cIndex]._index [2]], verts [cells [cIndex]._index [3]],
                                    bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                  break;
                }
              }
              break;
          }
        }
        _cutCells.push_front (cIndex);
        _reExaminedCells.pop_front ();
      }
    }

    // test for and respond to intersection between blade and cell
    void
    Partition::cellBladeCollide (unsigned int sIndex, forward_list <unsigned int>::iterator &lastElement, vector <Vertex> &vInfo,
                                 vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells, Cell &cell,
                                 vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      cell.setCellExamFlag ();

      unsigned int index;

      // test for vertex collisions
      for (unsigned int i = 0; i < 4; ++i){
        index = cell._index [i];
        if (!vInfo [index].testCollisionFlag ()){
          vInfo [index].setCollisionFlag ();
          for (unsigned int j = 0; j < bladeIndices.size ()/ 2; ++j){
            if (pointInTriangle (verts [index], bladeCurr [bladeIndices [2*j]],
                                 bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*bladeNormals [0]) [j]) ||
                pointInTriangle (verts [index], bladePrev [bladeIndices [2*j + 1]],
                                 bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*bladeNormals [1]) [j])){
              _collidingVertices.push_front (index);
              cell.setVertexCollisionFlag (i);
              for (unsigned int k = 0; k < vInfo [index]._numSubmeshes; ++k){
                if (vInfo [cell._index [i]]._owners [k][0] == sIndex){
                  for (unsigned int l = 2; l < vInfo [index]._owners [k][1]; ++l){
                    if (!cells [vInfo [index]._owners [k][l]].testCellExamFlag ()){
                      _cutCells.insert_after (lastElement, vInfo [index]._owners [k][l]);
                      cells [vInfo [index]._owners [k][l]].setThisVertexCollisionFlag (index);
                      ++lastElement;
                    }
                  }
                  break;
                }
              }
              break;
            }
          }
        }
      }

      // test edges for collisions
      real eu1 = 0., eu2 = 0.;
      for (unsigned int i = 0; i < 6; ++i){
        eu1 = 0.;
        eu2 = 0.;
        index = cell._edgeIndex [i];

        switch (i){

          case 0: // edge 0-1
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [0]], verts [cell._index [1]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [0]], verts [cell._index [1]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [0]);
                    _collidingVertices.push_front (cell._index [1]);
                  }
                } else if (cell._index [0] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;

          case 1: // edge 0-2
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [0]], verts [cell._index [2]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [0]], verts [cell._index [2]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [0]);
                    _collidingVertices.push_front (cell._index [2]);
                  }
                } else if (cell._index [0] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;

          case 2: // edge 0-3
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [0]], verts [cell._index [3]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [0]], verts [cell._index [3]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [0]);
                    _collidingVertices.push_front (cell._index [3]);
                  }
                } else if (cell._index [0] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;

          case 3: // edge 1-2
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [1]], verts [cell._index [2]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [1]], verts [cell._index [2]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [1]);
                    _collidingVertices.push_front (cell._index [2]);
                  }
                } else if (cell._index [1] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;

          case 4: // edge 1-3
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [1]], verts [cell._index [3]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [1]], verts [cell._index [3]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [1]);
                    _collidingVertices.push_front (cell._index [3]);
                  }
                } else if (cell._index [1] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;

          case 5: // edge 2-3
            if (edges [index].testCollisionFlag () || edges [index]._u > 0.){
              break;
            }
            edges [index].setCollisionFlag ();
            for (unsigned int j = 0; j < bladeNormals [0]->size (); ++j){
              if (lineTriCollide (eu1, verts [cell._index [2]], verts [cell._index [3]],
                                  bladeCurr [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j + 1]], (*( bladeNormals [0]))[j]) ||
                  lineTriCollide (eu2, verts [cell._index [2]], verts [cell._index [3]],
                                  bladePrev [bladeIndices [2*j + 1]], bladePrev [bladeIndices [2*j]], bladeCurr [bladeIndices [2*j]], (*( bladeNormals [1]))[j])){
                if (eu1 > 0.){
                  edges [index]._u = eu1;
                } else {
                  edges [index]._u = eu2;
                }
                for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                  if (!cells [edges [index]._owner [k]].testCellExamFlag ()){
                    _cutCells.insert_after (lastElement, edges [index]._owner [k]);
                    ++lastElement;
                  }
                }
                if (edges [index]._u > 1.){
                  for (unsigned int k = 0; k < edges [index]._numOwners; ++k){
                    cells [edges [index]._owner [k]].setEdgeCollisionFlag ();
                    _collidingVertices.push_front (cell._index [2]);
                    _collidingVertices.push_front (cell._index [3]);
                  }
                } else if (cell._index [2] != edges [index]._firstVertex){
                  edges [index]._u = 1. - edges [index]._u;
                }
                break;
              }
            }
            break;
        }
      } // end - for (unsigned int i = 0; i < 6; ++i)
    }

    /** Case-based face formulation algorithm. Edge enumeration:
    * 0 : 0-1
    * 1 : 0-2
    * 2 : 0-3
    * 3 : 1-2
    * 4 : 1-3
    * 5 : 2-3
    */
    void
    Partition::formFaces (Cell &cell, Cut &cut, vector <Edge> &edges, vector <vec> &verts,
                          vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      unsigned char choice = 0x00;
      for (unsigned int i = 0; i < 6; ++i){
        if (edges [cell._edgeIndex [i]]._u > 0.){
          choice |= 0x01 << i;
        }
      }

      switch (choice){

        case 0x00: // no edges cut
        break;

        case 0x01: // edge 0 cut
        {
          real u = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2), 0, 1, 2, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x02: // edge 1 cut
        {
          real u = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (0), 0, 2, 3, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x04: // edge 2 cut
        {
          real u = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1), 0, 3, 1, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x08: // edge 3 cut
        {
          real u = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3), 1, 2, 0, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x10: // edge 4 cut
        {
          real u = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (2), 1, 3, 2, 0, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x20: // edge 5 cut
        {
          real u =  edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u = 1. - u;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performOneEdgeCut (u, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3), 2, 3, 0, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x03: // edges 0, 1 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0),
                             0, 1, 2, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x05: // edges 0, 2 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (2),
                             0, 3, 1, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x09: // edges 0, 3 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (0),
                             1, 2, 0, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x11: // edges 0, 4 cut
        {
          real u0 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2),
                             1, 0, 3, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x06: // edges 1, 2 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1),
                             0, 2, 3, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x0A: // edges 1, 3 cut
        {
          real u0 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (0),
                             2, 0, 1, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x22: // edges 1, 5 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (1),
                             2, 3, 0, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x14: // edges 2, 4 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (2),
                             3, 1, 0, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x24: // edges 2, 5 cut
        {
          real u0 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1),
                             3, 0, 2, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x18: // edges 3, 4 cut
        {
          real u0 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3),
                             1, 3, 2, 0, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x28: // edges 3, 5 cut
        {
          real u0 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3),
                             2, 1, 3, 0, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x30: // edges 4, 5 cut
        {
          real u0 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performTwoEdgeCut (u0, u1, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3),
                             3, 2, 1, 0, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x13: // edges 0, 1, 4 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2),
                                         cell.testExternalFaceFlag (3), 0, 1, 2, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x23: // edges 0, 1, 5 cut
        {
          real u0 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (0),
                                         cell.testExternalFaceFlag (2), 0, 2, 3, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x0D: // edges 0, 2, 3 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0),
                                         cell.testExternalFaceFlag (3), 1, 0, 3, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x25: // edges 0, 2, 5 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1),
                                         cell.testExternalFaceFlag (3), 0, 3, 1, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x29: // edges 0, 3, 5 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3),
                                         cell.testExternalFaceFlag (1), 1, 2, 0, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x31: // edges 0, 4, 5 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3),
                                         cell.testExternalFaceFlag (1), 3, 1, 0, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x0E: // edges 1, 2, 3 cut
        {
          real u0 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1),
                                         cell.testExternalFaceFlag (2), 2, 0, 1, 3, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x16: // edges 1, 2, 4 cut
        {
          real u0 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1),
                                         cell.testExternalFaceFlag (0), 0, 3, 1, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x1A: // edges 1, 3, 4 cut
        {
          real u0 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (0),
                                         cell.testExternalFaceFlag (1), 2, 1, 3, 0, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x32: // edges 1, 4, 5 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3),
                                         cell.testExternalFaceFlag (2), 2, 3, 0, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x1C: // edges 2, 3, 4 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (true, u0, u1, u2, cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3),
                                         cell.testExternalFaceFlag (0), 3, 1, 0, 2, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x2C: // edges 2, 3, 5 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          vec *cellVerts [4] = {&(verts [cell._index [0]]), &(verts [cell._index [1]]), &(verts [cell._index [2]]), &(verts [cell._index [3]])};
          performUnfinishedThreeEdgeCut (false, u0, u1, u2, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3),
                                         cell.testExternalFaceFlag (0), 2, 3, 0, 1, cellVerts, cell, cut, bladeCurr, bladePrev, bladeIndices, bladeNormals);
        }
        break;

        case 0x07: // edges 0, 1, 2 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u2 = 1. - u2;
          }
          performFinishedThreeEdgeCut (u0, u1, u2, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1), 0, 2, 1, 3, cell, cut);
        }
        break;

        case 0x19: // edges 0, 3, 4 cut
        {
          real u0 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [1]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [1]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [1]){
            u2 = 1. - u2;
          }
          performFinishedThreeEdgeCut (u0, u1, u2, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3), 1, 3, 0, 2, cell, cut);
        }
        break;

        case 0x2A: // edges 1, 3, 5 cut
        {
          real u0 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          performFinishedThreeEdgeCut (u0, u1, u2, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1), cell.testExternalFaceFlag (3), 2, 1, 0, 3, cell, cut);
        }
        break;

        case 0x34: // edges 2, 4, 5 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [3]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          performFinishedThreeEdgeCut (u0, u1, u2, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (1), 3, 0, 1, 2, cell, cut);
        }
        break;

        case 0x33: // edges 0, 1, 4, 5 cut
        {
          real u0 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          real u3 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u3 = 1. - u3;
          }
          performFourEdgeCut (u0, u1, u2, u3, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (1),
                              cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (2), 0, 1, 2, 3, cell, cut);
        }
        break;

        case 0x1E: // edges 1, 2, 3, 4 cut
        {
          real u0 = edges [cell._edgeIndex [1]]._u;
          if (edges [cell._edgeIndex [1]]._firstVertex != cell._index [2]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [4]]._u;
          if (edges [cell._edgeIndex [4]]._firstVertex != cell._index [3]){
            u2 = 1. - u2;
          }
          real u3 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [3]){
            u3 = 1. - u3;
          }
          performFourEdgeCut (u0, u1, u2, u3, cell.testExternalFaceFlag (0), cell.testExternalFaceFlag (3),
                              cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (1), 2, 0, 1, 3, cell, cut);
        }
        break;

        case 0x2D: // edges 0, 2, 3, 5 cut
        {
          real u0 = edges [cell._edgeIndex [2]]._u;
          if (edges [cell._edgeIndex [2]]._firstVertex != cell._index [0]){
            u0 = 1. - u0;
          }
          real u1 = edges [cell._edgeIndex [0]]._u;
          if (edges [cell._edgeIndex [0]]._firstVertex != cell._index [0]){
            u1 = 1. - u1;
          }
          real u2 = edges [cell._edgeIndex [3]]._u;
          if (edges [cell._edgeIndex [3]]._firstVertex != cell._index [2]){
            u2 = 1. - u2;
          }
          real u3 = edges [cell._edgeIndex [5]]._u;
          if (edges [cell._edgeIndex [5]]._firstVertex != cell._index [2]){
            u3 = 1. - u3;
          }
          performFourEdgeCut (u0, u1, u2, u3, cell.testExternalFaceFlag (2), cell.testExternalFaceFlag (0),
                              cell.testExternalFaceFlag (3), cell.testExternalFaceFlag (1), 0, 3, 1, 2, cell, cut);
        }
        break;

        default:
          ;//PRINT ("Unknown cut configuration: %x\n", choice);
      }

    }

    /** Private method to perform cut operations for one cut-edge case. Arguments:
    * 0     : u-index of the cut-edge
    * 1,2   : Flags to signify if faces are external (right face first)
    * 3,4,5 : Local vertex indices of right face (3-4 same as edge)
    * 6     : Fourth vertex index (args 3,6,4 constitute the left face)
    * 7     : 4-element array containing pointers to 4 cell vertices
    * 8     : Reference to cell-structure for which current method is being called
    * 9     : Reference to the cut-structure for current cell
    * 10-13 : Blade-related parameters
    */
    void
    Partition::performOneEdgeCut (real u, bool faceFlag0, bool faceFlag1, unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3,
                                  vec *verts [4], Cell &cell, Cut &cut,
                                  vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      // optionally allocate cut variables for internal faces
      bool newflag = false;
      if (cut._numInVertices < 4 || cut._numInFaces < 2){
        newflag = true;

        _inMutex->lock ();
        *_inUpdateFlag = true;
        cut.allocateInternalVariables (4, 2, *_inVertices, *_inSurfaceVertexStatus, *_in2DTexCoords, *_in3DTexCoords, *_inFaceIndices, _inEmptyVertices, _inEmptyFaces);
        _inMutex->unlock ();

        unsigned int inds [4] = {cut._inVertices [0], cut._inVertices [1], cut._inVertices [2], cut._inVertices [3]};

        if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag (v1)){
          (* _inSurfaceVertexStatus) [inds [0]] = 1.;
          (* _inSurfaceVertexStatus) [inds [1]] = 1.;
          (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v0]] * (1. - u) + (* _tex2D) [cell._index [v1]] * u;
          (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
        }
        (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v0]] * (1. - u) + (* _tex3D) [cell._index [v1]] * u;
        (* _in3DTexCoords) [inds [1]] = (* _in3DTexCoords) [inds [0]];

        unsigned int tmpu = cut._inFaces [0];
        (* _inFaceIndices) [tmpu] = inds [0];
        (* _inFaceIndices) [tmpu + 1] = inds [3];
        (* _inFaceIndices) [tmpu + 2] = inds [2];

        tmpu = cut._inFaces [1];
        (* _inFaceIndices) [tmpu] = inds [2];
        (* _inFaceIndices) [tmpu + 1] = inds [3];
        (* _inFaceIndices) [tmpu + 2] = inds [1];
      }

      // get edge-points
      (* _inVertices) [cut._inVertices [0]] = (*(verts [v0])) * (1. - u + CUT_DISTANCE) + (*(verts [v1])) * (u - CUT_DISTANCE);
      (* _inVertices) [cut._inVertices [1]] = (*(verts [v0])) * (1. - u - CUT_DISTANCE) + (*(verts [v1])) * (u + CUT_DISTANCE);

      vec point0 = ((* _inVertices) [cut._inVertices [0]] + (* _inVertices) [cut._inVertices [1]]) * .5;

      // normal for trig 0
      vec edge0 = *(verts [v1]) - *(verts [v0]);
      vec edge1 = *(verts [v2]) - *(verts [v0]);
      vec normal = edge0.cross (edge1);

      real dist1, dist2, max = 0.;
      vec point1, point2, maxpt;

      // get point furthest from edge
      bool collideflag1 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v1]), *(verts [v2]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v1]), *(verts [v2]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      vec2 uv;
      if (collideflag1){
        (* _inVertices) [cut._inVertices [2]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [v0]), *(verts [v1]), *(verts [v2]));
        (* _in3DTexCoords) [cut._inVertices [2]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v1]]* uv._v [1] +
          (* _tex3D) [cell._index [v2]]* (1. - uv._v [0] - uv._v [1]);
      }

      // normal for trig 1
      edge0 = *(verts [v3]) - *(verts [v0]);
      edge1 = *(verts [v1]) - *(verts [v0]);
      edge0.fast_cross (normal, edge1);

      // get point furthest from edge
      max = 0.;
      bool collideflag2 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v3]), *(verts [v1]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v3]), *(verts [v1]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      if (collideflag2){
        (* _inVertices) [cut._inVertices [3]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [v0]), *(verts [v3]), *(verts [v1]));
        (* _in3DTexCoords) [cut._inVertices [3]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v3]]* uv._v [1] +
          (* _tex3D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
      }

      // early return
      if (!(faceFlag0 || faceFlag1)){
        return;
      }

      /** Optionally allocate cut variables for external faces **/
      // both trigs are external
      if (faceFlag0 && faceFlag1){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (8, 8, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [8] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3],
            cut._exVertices [4], cut._exVertices [5], cut._exVertices [6], cut._exVertices [7]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [0];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = inds [1];
          (* _exFaceIndices) [tmpu + 1] = inds [3];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [5];
          (* _exFaceIndices) [tmpu] = inds [3];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = inds [3];
          (* _exFaceIndices) [tmpu + 1] = inds [4];
          (* _exFaceIndices) [tmpu + 2] = inds [7];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = inds [3];
          (* _exFaceIndices) [tmpu + 1] = inds [0];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          (* _inSurfaceVertexStatus) [cut._inVertices [2]] = 1.;
          (* _inSurfaceVertexStatus) [cut._inVertices [3]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [7]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices) [cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices) [cut._exVertices [4]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [5]] = *(verts [v1]);
        (* _exVertices) [cut._exVertices [6]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [7]] = *(verts [v3]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [2]], *(verts [v0]), *(verts [v1]), *(verts [v2]));
          (* _in2DTexCoords) [cut._inVertices [2]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
            (* _tex2D) [cell._index [v2]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [2]];
        }
        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [3]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [3]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [3]] = (* _in2DTexCoords) [cut._inVertices [3]];
        }
      }
      // trig 0 is external
      else if (faceFlag0){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [6] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4], cut._exVertices [5]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [0];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [4];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          (* _inSurfaceVertexStatus) [cut._inVertices [2]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [3]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v2]];
        }

        (* _exVertices)[cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices)[cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices)[cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices)[cut._exVertices [3]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [4]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v2]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [2]], *(verts [v0]), *(verts [v1]), *(verts [v2]));
          (* _in2DTexCoords) [cut._inVertices [2]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
            (* _tex2D) [cell._index [v2]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [2]];
        }
      }
      // trig 1 is external
      else if (faceFlag1){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [6] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4], cut._exVertices [5]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [3];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [1];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [4];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [3];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          (* _inSurfaceVertexStatus) [cut._inVertices [3]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [3]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v3]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v1]];
        }

        (* _exVertices)[cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices)[cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices)[cut._exVertices [2]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices)[cut._exVertices [3]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [4]] = *(verts [v3]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v1]);

        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [3]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [3]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [3]];
        }
      }
    }

    /** Private method to perform cut operations for two cut-edges case. Arguments:
    * 0     : u-index of the right cut-edge
    * 1     : u-index of the left cut-edge
    * 2,3,4 : Flags to signify if faces are external (right, left, center)
    * 5,6,7 : Local vertex index of center face (one with two cut edges)
    * 8     : Fourth vertex index (Right - args 5,7,8 Left - args 5,8,6)
    * 9     : 4-element array containing pointers to 4 cell vertices
    * 10    : Reference to cell-structure for which current method is being called
    * 11    : Reference to the cut-structure for current cell
    * 12-15 : Blade-related parameters
    */
    void
    Partition::performTwoEdgeCut (real u0, real u1, bool faceFlag0, bool faceFlag1, bool faceFlag2,
                                  unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, vec *verts [4], Cell &cell, Cut &cut,
                                  vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      // optionally allocate cut variables for internal faces
      bool newflag = false;
      if (cut._numInVertices < 6 || cut._numInFaces < 4){
        newflag = true;

        bool alreadyflag = cut._numInVertices > 0 ? true : false;

        _inMutex->lock ();
        *_inUpdateFlag = true;
        cut.allocateInternalVariables (6, 4, *_inVertices, *_inSurfaceVertexStatus, *_in2DTexCoords, *_in3DTexCoords, *_inFaceIndices, _inEmptyVertices, _inEmptyFaces);
        _inMutex->unlock ();

        // this means that the cell has one edge cut from previous step
        if (alreadyflag){
          vec tmpv = (* _inVertices) [cut._inVertices [2]] - *(verts [v0]);

          vec nrm1 = (*(verts [v2]) - *(verts [v0])).ncross (*(verts [v3]) - *(verts [v0]));
          vec nrm2 = (*(verts [v3]) - *(verts [v0])).ncross (*(verts [v1]) - *(verts [v0]));

          vec2 uv;
          if (ABS (tmpv.dot (nrm1)) < 100.*EPSILON ){
            tmpv += *(verts [v0]);
            (* _inVertices) [cut._inVertices [4]] = tmpv;
            calculateBarycentricCoords (uv, tmpv, *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in3DTexCoords) [cut._inVertices [4]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v2]]* uv._v [1] +
              (* _tex3D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          }
          else if (ABS (tmpv.dot (nrm2)) < 100.*EPSILON ){
            tmpv += *(verts [v0]);
            (* _inVertices) [cut._inVertices [5]] = tmpv;
            calculateBarycentricCoords (uv, tmpv, *(verts [v0]), *(verts [v3]), *(verts [v1]));
            (* _in3DTexCoords) [cut._inVertices [5]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v3]]* uv._v [1] +
              (* _tex3D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          }

          tmpv = (* _inVertices) [cut._inVertices [3]] - *(verts [v0]);

          if (ABS (tmpv.dot (nrm1)) < 100.*EPSILON ){
            tmpv += *(verts [v0]);
            (* _inVertices) [cut._inVertices [4]] = tmpv;
            calculateBarycentricCoords (uv, tmpv, *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in3DTexCoords) [cut._inVertices [4]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v2]]* uv._v [1] +
              (* _tex3D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          }
          else if (ABS (tmpv.dot (nrm2)) < 100.*EPSILON ){
            tmpv += *(verts [v0]);
            (* _inVertices) [cut._inVertices [5]] = tmpv;
            calculateBarycentricCoords (uv, tmpv, *(verts [v0]), *(verts [v3]), *(verts [v1]));
            (* _in3DTexCoords) [cut._inVertices [5]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v3]]* uv._v [1] +
              (* _tex3D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          }

        }

        unsigned int inds [6] = {cut._inVertices [0], cut._inVertices [1], cut._inVertices [2], cut._inVertices [3], cut._inVertices [4], cut._inVertices [5]};

        if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v2))){
          (* _inSurfaceVertexStatus) [inds [0]] = 1.;
          (* _inSurfaceVertexStatus) [inds [1]] = 1.;
          (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v0]] * (1. - u0) + (* _tex2D) [cell._index [v2]] * u0;
          (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
        }
        if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v1))){
          (* _inSurfaceVertexStatus) [inds [2]] = 1.;
          (* _inSurfaceVertexStatus) [inds [3]] = 1.;
          (* _in2DTexCoords) [inds [2]] = (* _tex2D) [cell._index [v0]] * (1. - u1) + (* _tex2D) [cell._index [v1]] * u1;
          (* _in2DTexCoords) [inds [3]] = (* _in2DTexCoords) [inds [2]];
        }

        (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v0]] * (1. - u0) + (* _tex3D) [cell._index [v2]] * u0;
        (* _in3DTexCoords) [inds [1]] = (* _in3DTexCoords) [inds [0]];
        (* _in3DTexCoords) [inds [2]] = (* _tex3D) [cell._index [v0]] * (1. - u1) + (* _tex3D) [cell._index [v1]] * u1;
        (* _in3DTexCoords) [inds [3]] = (* _in3DTexCoords) [inds [2]];

        unsigned int tmpu = cut._inFaces [0];
        (* _inFaceIndices) [tmpu] = inds [0];
        (* _inFaceIndices) [tmpu + 1] = inds [4];
        (* _inFaceIndices) [tmpu + 2] = inds [5];

        tmpu = cut._inFaces [1];
        (* _inFaceIndices) [tmpu] = inds [0];
        (* _inFaceIndices) [tmpu + 1] = inds [5];
        (* _inFaceIndices) [tmpu + 2] = inds [2];

        tmpu = cut._inFaces [2];
        (* _inFaceIndices) [tmpu] = inds [1];
        (* _inFaceIndices) [tmpu + 1] = inds [4];
        (* _inFaceIndices) [tmpu + 2] = inds [3];

        tmpu = cut._inFaces [3];
        (* _inFaceIndices) [tmpu] = inds [5];
        (* _inFaceIndices) [tmpu + 1] = inds [3];
        (* _inFaceIndices) [tmpu + 2] = inds [4];
      }

      // get edge-points
      (* _inVertices) [cut._inVertices [0]] = (*(verts [v0])) * (1. - u0 + CUT_DISTANCE) + (*(verts [v2])) * (u0 - CUT_DISTANCE);
      (* _inVertices) [cut._inVertices [1]] = (*(verts [v0])) * (1. - u0 - CUT_DISTANCE) + (*(verts [v2])) * (u0 + CUT_DISTANCE);
      (* _inVertices) [cut._inVertices [2]] = (*(verts [v0])) * (1. - u1 + CUT_DISTANCE) + (*(verts [v1])) * (u1 - CUT_DISTANCE);
      (* _inVertices) [cut._inVertices [3]] = (*(verts [v0])) * (1. - u1 - CUT_DISTANCE) + (*(verts [v1])) * (u1 + CUT_DISTANCE);

      vec point0 = ((* _inVertices) [cut._inVertices [0]] + (* _inVertices) [cut._inVertices [1]]) * .5;

      // normal for trig 0
      vec edge0 = *(verts [v2]) - *(verts [v0]);
      vec edge1 = *(verts [v3]) - *(verts [v0]);
      vec normal = edge0.cross (edge1);

      real dist1, dist2, max = 0.;
      vec point1, point2, maxpt;

      // get point furthest from edge
      bool collideflag1 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v2]), *(verts [v3]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v2]), *(verts [v3]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      vec2 uv;
      if (collideflag1){
        (* _inVertices) [cut._inVertices [4]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [v0]), *(verts [v2]), *(verts [v3]));
        (* _in3DTexCoords) [cut._inVertices [4]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v2]]* uv._v [1] +
          (* _tex3D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
      }

      point0 = ((* _inVertices) [cut._inVertices [2]] + (* _inVertices) [cut._inVertices [3]]) * .5;

      // normal for trig 1
      edge0 = *(verts [v3]) - *(verts [v0]);
      edge1 = *(verts [v1]) - *(verts [v0]);
      edge0.fast_cross (normal, edge1);

      // get point furthest from edge
      max = 0.;
      bool collideflag2 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v3]), *(verts [v1]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [v0]), *(verts [v3]), *(verts [v1]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      if (collideflag2){
        (* _inVertices) [cut._inVertices [5]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [v0]), *(verts [v3]), *(verts [v1]));
        (* _in3DTexCoords) [cut._inVertices [5]] = (* _tex3D) [cell._index [v0]]* uv._v [0] + (* _tex3D) [cell._index [v3]]* uv._v [1] +
          (* _tex3D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
      }

      // early return
      if (!(faceFlag0 || faceFlag1 || faceFlag2)){
        return;
      }

      // optionally allocate cut variables for external faces
      if (faceFlag0 && faceFlag1 && faceFlag2){ // all trigs are external

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (10, 11, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
            cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [8];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [8];
          (* _exFaceIndices) [tmpu + 2] = inds [9];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [9];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [5];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [9];
          (* _exFaceIndices) [tmpu + 2] = inds [7];

          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [9];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          tmpu = cut._exFaces [8];
          (* _exFaceIndices) [tmpu] = inds [6];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [9];
          (* _exFaceIndices) [tmpu] = inds [7];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [10];
          (* _exFaceIndices) [tmpu] = inds [7];
          (* _exFaceIndices) [tmpu + 1] = inds [8];
          (* _exFaceIndices) [tmpu + 2] = inds [1];

          (* _inSurfaceVertexStatus) [cut._inVertices [4]] = 1.;
          (* _inSurfaceVertexStatus) [cut._inVertices [5]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [2]] = (* _in2DTexCoords)[cut._inVertices [2]];
          (* _ex2DTexCoords) [inds [3]] = (* _in2DTexCoords)[cut._inVertices [3]];
          (* _ex2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [7]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [8]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [9]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices) [cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices) [cut._exVertices [4]] = (* _inVertices) [cut._inVertices [4]];
        (* _exVertices) [cut._exVertices [5]] = (* _inVertices) [cut._inVertices [5]];
        (* _exVertices) [cut._exVertices [6]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [7]] = *(verts [v1]);
        (* _exVertices) [cut._exVertices [8]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [9]] = *(verts [v3]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [4]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
          (* _in2DTexCoords) [cut._inVertices [4]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
            (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords)[cut._inVertices [4]];
        }

        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [5]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [5]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _in2DTexCoords)[cut._inVertices [5]];
        }
      }
      // trigs 0 and 1 are external
      else if (faceFlag0 && faceFlag1){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (10, 8, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
            cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [8];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [8];
          (* _exFaceIndices) [tmpu + 2] = inds [9];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [9];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [5];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [9];
          (* _exFaceIndices) [tmpu + 2] = inds [7];

          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [9];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          (* _inSurfaceVertexStatus) [cut._inVertices [4]] = 1.;
          (* _inSurfaceVertexStatus) [cut._inVertices [5]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [2]] = (* _in2DTexCoords)[cut._inVertices [2]];
          (* _ex2DTexCoords) [inds [3]] = (* _in2DTexCoords)[cut._inVertices [3]];
          (* _ex2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [7]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [8]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [9]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices) [cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices) [cut._exVertices [4]] = (* _inVertices) [cut._inVertices [4]];
        (* _exVertices) [cut._exVertices [5]] = (* _inVertices) [cut._inVertices [5]];
        (* _exVertices) [cut._exVertices [6]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [7]] = *(verts [v1]);
        (* _exVertices) [cut._exVertices [8]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [9]] = *(verts [v3]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [4]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
          (* _in2DTexCoords) [cut._inVertices [4]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
            (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords)[cut._inVertices [4]];
        }

        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [5]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [5]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _in2DTexCoords)[cut._inVertices [5]];
        }
      }
      // trigs 0 and 2 are external
      else if (faceFlag0 && faceFlag2){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (9, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [9] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3],
            cut._exVertices [4], cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [7];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [8];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [8];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [5];
          (* _exFaceIndices) [tmpu] = inds [6];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = inds [6];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [1];

          (* _inSurfaceVertexStatus) [cut._inVertices [4]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [2]] = (* _in2DTexCoords)[cut._inVertices [2]];
          (* _ex2DTexCoords) [inds [3]] = (* _in2DTexCoords)[cut._inVertices [3]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [7]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [8]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices) [cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices) [cut._exVertices [4]] = (* _inVertices) [cut._inVertices [4]];
        (* _exVertices) [cut._exVertices [5]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [6]] = *(verts [v1]);
        (* _exVertices) [cut._exVertices [7]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [8]] = *(verts [v1]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [4]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
          (* _in2DTexCoords) [cut._inVertices [4]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
            (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords)[cut._inVertices [4]];
        }
      }
      // trigs 1 and 2 are external
      else if (faceFlag1 && faceFlag2){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (9, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [9] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3],
            cut._exVertices [4], cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [6];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [8];
          (* _exFaceIndices) [tmpu + 2] = inds [6];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [8];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [4];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = inds [5];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [5];
          (* _exFaceIndices) [tmpu] = inds [6];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = inds [6];
          (* _exFaceIndices) [tmpu + 1] = inds [7];
          (* _exFaceIndices) [tmpu + 2] = inds [1];

          (* _inSurfaceVertexStatus) [cut._inVertices [5]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [2]] = (* _in2DTexCoords)[cut._inVertices [2]];
          (* _ex2DTexCoords) [inds [3]] = (* _in2DTexCoords)[cut._inVertices [3]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [inds [7]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [8]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices) [cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices) [cut._exVertices [4]] = (* _inVertices) [cut._inVertices [5]];
        (* _exVertices) [cut._exVertices [5]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [6]] = *(verts [v1]);
        (* _exVertices) [cut._exVertices [7]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [8]] = *(verts [v3]);

        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [5]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [5]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords)[cut._inVertices [5]];
        }
      }
      // trig 0 is external
      else if (faceFlag0){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [6] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4], cut._exVertices [5]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [0];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [4];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [3];

          (* _inSurfaceVertexStatus) [cut._inVertices [4]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [1]];
          (* _ex2DTexCoords) [inds [3]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v3]];
        }

        (* _exVertices) [cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices) [cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [4]];
        (* _exVertices) [cut._exVertices [3]] = *(verts [v0]);
        (* _exVertices) [cut._exVertices [4]] = *(verts [v2]);
        (* _exVertices) [cut._exVertices [5]] = *(verts [v3]);

        if (collideflag1){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [4]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
          (* _in2DTexCoords) [cut._inVertices [4]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
            (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords)[cut._inVertices [4]];
        }
      }
      // trig 1 is external
      else if (faceFlag1){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [6] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4], cut._exVertices [5]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = inds [3];
          (* _exFaceIndices) [tmpu + 1] = inds [2];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [5];
          (* _exFaceIndices) [tmpu + 2] = inds [1];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [4];
          (* _exFaceIndices) [tmpu + 2] = inds [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = inds [2];
          (* _exFaceIndices) [tmpu + 1] = inds [3];
          (* _exFaceIndices) [tmpu + 2] = inds [4];

          (* _inSurfaceVertexStatus) [cut._inVertices [5]] = 1.;

          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords)[cut._inVertices [2]];
          (* _ex2DTexCoords) [inds [1]] = (* _in2DTexCoords)[cut._inVertices [3]];
          (* _ex2DTexCoords) [inds [3]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v3]];
          (* _ex2DTexCoords) [inds [5]] = (* _tex2D) [cell._index [v1]];
        }

        (* _exVertices)[cut._exVertices [0]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices)[cut._exVertices [1]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices)[cut._exVertices [2]] = (* _inVertices) [cut._inVertices [5]];
        (* _exVertices)[cut._exVertices [3]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [4]] = *(verts [v3]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v1]);

        if (collideflag2){
          calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [5]], *(verts [v0]), *(verts [v3]), *(verts [v1]));
          (* _in2DTexCoords) [cut._inVertices [5]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v3]]* uv._v [1] +
            (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords)[cut._inVertices [5]];
        }
      }
      // trig 2 is external
      else if (faceFlag2){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int inds [2] = {cut._exVertices [1], cut._exVertices [5]};

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = cut._exVertices [4];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = cut._exVertices [3];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = inds [0];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = inds [0];
          (* _exFaceIndices) [tmpu + 1] = inds [1];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

          (* _ex2DTexCoords) [cut._exVertices [0]] = (* _in2DTexCoords) [cut._inVertices [0]];
          (* _ex2DTexCoords) [inds [0]] = (* _in2DTexCoords) [cut._inVertices [1]];
          (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [2]];
          (* _ex2DTexCoords) [cut._exVertices [3]] = (* _in2DTexCoords) [cut._inVertices [3]];
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [inds [1]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [cut._exVertices [6]] = (* _tex2D) [cell._index [v2]];
        }

        (* _exVertices)[cut._exVertices [0]] = (* _inVertices) [cut._inVertices [0]];
        (* _exVertices)[cut._exVertices [1]] = (* _inVertices) [cut._inVertices [1]];
        (* _exVertices)[cut._exVertices [2]] = (* _inVertices) [cut._inVertices [2]];
        (* _exVertices)[cut._exVertices [3]] = (* _inVertices) [cut._inVertices [3]];
        (* _exVertices)[cut._exVertices [4]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [6]] = *(verts [v2]);
      }

    }


    /** Private method to perform cut operations for three cut-edges cases with no severance. Arguments:
    */
    void
    Partition::performUnfinishedThreeEdgeCut (bool upFlag, real u0, real u1, real u2, bool faceFlag0, bool faceFlag1, bool faceFlag2, bool faceFlag3,
                                              unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, vec *verts [4], Cell &cell, Cut &cut,
                                              vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      // optionally allocate cut variables for internal faces
      bool newflag = false;
      if (cut._numInVertices < 8 || cut._numInFaces < 6){
        newflag = true;

        unsigned int alreadyflag = 0;
        if (cut._numInVertices > 5){
          alreadyflag = 2;
        } else if (cut._numInVertices > 3){
          alreadyflag = 1;
        }

        _inMutex->lock ();
        *_inUpdateFlag = true;
        cut.allocateInternalVariables (8, 6, *_inVertices, *_inSurfaceVertexStatus, *_in2DTexCoords, *_in3DTexCoords, *_inFaceIndices, _inEmptyVertices, _inEmptyFaces);
        _inMutex->unlock ();

        if (alreadyflag){

          unsigned int inds [2];
          if (alreadyflag > 1){
            inds [0] = cut._inVertices [4];
            inds [1] = cut._inVertices [5];
          } else {
            inds [0] = cut._inVertices [2];
            inds [1] = cut._inVertices [3];
          }

          unsigned int i0, i1, i2 = v3;
          if (upFlag){
            i0 = v0;
            i1 = v2;
          } else {
            i0 = v2;
            i1 = v1;
          }
          vec tmpv1 = (* _inVertices) [inds [0]] - *(verts [i0]);
          vec tmpv2 = (* _inVertices) [inds [1]] - *(verts [i0]);

          vec nrm = (*(verts [i1]) - *(verts [i0])).ncross (*(verts [i2]) - *(verts [i0]));

          vec2 uv;
          if (ABS (tmpv1.dot (nrm)) < 100.*EPSILON ){
            tmpv1 += *(verts [i0]);
            (* _inVertices) [cut._inVertices [6]] = tmpv1;
            calculateBarycentricCoords (uv, tmpv1, *(verts [i0]), *(verts [i1]), *(verts [i2]));
            (* _in3DTexCoords) [cut._inVertices [6]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
              (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
          }
          else if (ABS (tmpv2.dot (nrm)) < 100.*EPSILON ){
            tmpv2 += *(verts [v0]);
            (* _inVertices) [cut._inVertices [6]] = tmpv2;
            calculateBarycentricCoords (uv, tmpv2, *(verts [i0]), *(verts [i1]), *(verts [i2]));
            (* _in3DTexCoords) [cut._inVertices [6]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
              (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
          }

          i1 = v2;
          if (upFlag){
            i0 = v3;
            i2 = v1;
          } else {
            i0 = v0;
            i2 = v3;
          }

          nrm = (*(verts [i1]) - *(verts [i0])).ncross (*(verts [i2]) - *(verts [i0]));

          if (ABS (tmpv1.dot (nrm)) < 100.*EPSILON ){
            tmpv1 += *(verts [i0]);
            (* _inVertices) [cut._inVertices [7]] = tmpv1;
            calculateBarycentricCoords (uv, tmpv1, *(verts [i0]), *(verts [i1]), *(verts [i2]));
            (* _in3DTexCoords) [cut._inVertices [7]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
              (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
          }
          else if (ABS (tmpv2.dot (nrm)) < 100.*EPSILON ){
            tmpv2 += *(verts [v0]);
            (* _inVertices) [cut._inVertices [7]] = tmpv2;
            calculateBarycentricCoords (uv, tmpv2, *(verts [i0]), *(verts [i1]), *(verts [i2]));
            (* _in3DTexCoords) [cut._inVertices [7]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
              (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
          }
        }

        unsigned int inds [8] = {cut._inVertices [0], cut._inVertices [1], cut._inVertices [2], cut._inVertices [3],
          cut._inVertices [4], cut._inVertices [5], cut._inVertices [6], cut._inVertices [7]};

        if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag (v1)){
          (* _inSurfaceVertexStatus) [inds [2]] = 1.;
          (* _inSurfaceVertexStatus) [inds [3]] = 1.;
          (* _in2DTexCoords) [inds [2]] = (* _tex2D) [cell._index [v0]] * (1. - u1) + (* _tex2D) [cell._index [v1]] * u1;
          (* _in2DTexCoords) [inds [3]] = (* _in2DTexCoords) [inds [2]];
        }
        (* _in3DTexCoords) [inds [2]] = (* _tex3D) [cell._index [v0]] * (1. - u1) + (* _tex3D) [cell._index [v1]] * u1;
        (* _in3DTexCoords) [inds [3]] = (* _in3DTexCoords) [inds [2]];

        if (upFlag){
          if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag (v2)){
            (* _inSurfaceVertexStatus) [inds [0]] = 1.;
            (* _inSurfaceVertexStatus) [inds [1]] = 1.;
            (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v0]] * (1. - u0) + (* _tex2D) [cell._index [v2]] * u0;
            (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
          }
          (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v0]] * (1. - u0) + (* _tex3D) [cell._index [v2]] * u0;

          if (cell.testExternalVertexFlag (v3) && cell.testExternalVertexFlag (v1)){
            (* _inSurfaceVertexStatus) [inds [4]] = 1.;
            (* _inSurfaceVertexStatus) [inds [5]] = 1.;
            (* _in2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v3]] * (1. - u2) + (* _tex2D) [cell._index [v1]] * u2;
            (* _in2DTexCoords) [inds [5]] = (* _in2DTexCoords) [inds [4]];
          }
          (* _in3DTexCoords) [inds [4]] = (* _tex3D) [cell._index [v3]] * (1. - u2) + (* _tex3D) [cell._index [v1]] * u2;
        }
        else {
          if (cell.testExternalVertexFlag (v2) && cell.testExternalVertexFlag (v1)){
            (* _inSurfaceVertexStatus) [inds [0]] = 1.;
            (* _inSurfaceVertexStatus) [inds [1]] = 1.;
            (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v2]] * (1. - u0) + (* _tex2D) [cell._index [v1]] * u0;
            (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
          }
          (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v2]] * (1. - u0) + (* _tex3D) [cell._index [v1]] * u0;

          if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag (v3)){
            (* _inSurfaceVertexStatus) [inds [4]] = 1.;
            (* _inSurfaceVertexStatus) [inds [5]] = 1.;
            (* _in2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v0]] * (1. - u2) + (* _tex2D) [cell._index [v3]] * u2;
            (* _in2DTexCoords) [inds [5]] = (* _in2DTexCoords) [inds [4]];
          }
          (* _in3DTexCoords) [inds [4]] = (* _tex3D) [cell._index [v0]] * (1. - u2) + (* _tex3D) [cell._index [v3]] * u2;
        }

        (* _in3DTexCoords) [inds [1]] = (* _in3DTexCoords) [inds [0]];
        (* _in3DTexCoords) [inds [5]] = (* _in3DTexCoords) [inds [4]];

        unsigned int tmpu = cut._inFaces [0];
        (* _inFaceIndices) [tmpu] = inds [6];
        (* _inFaceIndices) [tmpu + 1] = inds [0];
        (* _inFaceIndices) [tmpu + 2] = inds [2];

        tmpu = cut._inFaces [1];
        (* _inFaceIndices) [tmpu] = inds [2];
        (* _inFaceIndices) [tmpu + 1] = inds [7];
        (* _inFaceIndices) [tmpu + 2] = inds [6];

        tmpu = cut._inFaces [2];
        (* _inFaceIndices) [tmpu] = inds [2];
        (* _inFaceIndices) [tmpu + 1] = inds [4];
        (* _inFaceIndices) [tmpu + 2] = inds [7];

        tmpu = cut._inFaces [3];
        (* _inFaceIndices) [tmpu] = inds [1];
        (* _inFaceIndices) [tmpu + 1] = inds [3];
        (* _inFaceIndices) [tmpu + 2] = inds [6];

        tmpu = cut._inFaces [4];
        (* _inFaceIndices) [tmpu] = inds [3];
        (* _inFaceIndices) [tmpu + 1] = inds [7];
        (* _inFaceIndices) [tmpu + 2] = inds [6];

        tmpu = cut._inFaces [5];
        (* _inFaceIndices) [tmpu] = inds [3];
        (* _inFaceIndices) [tmpu + 1] = inds [5];
        (* _inFaceIndices) [tmpu + 2] = inds [7];

      } // end - if (cut._numInVertices < 8 || cut._numInFaces < 6)

      // get edge-points
      (* _inVertices) [cut._inVertices [2]] = (*(verts [v0])) * (1. - u1 + CUT_DISTANCE) + (*(verts [v1])) * (u1 - CUT_DISTANCE);
      (* _inVertices) [cut._inVertices [3]] = (*(verts [v0])) * (1. - u1 - CUT_DISTANCE) + (*(verts [v1])) * (u1 + CUT_DISTANCE);
      if (upFlag){
        (* _inVertices) [cut._inVertices [0]] = (*(verts [v0])) * (1. - u0 + CUT_DISTANCE) + (*(verts [v2])) * (u0 - CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [1]] = (*(verts [v0])) * (1. - u0 - CUT_DISTANCE) + (*(verts [v2])) * (u0 + CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [4]] = (*(verts [v3])) * (1. - u2 + CUT_DISTANCE) + (*(verts [v1])) * (u2 - CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [5]] = (*(verts [v3])) * (1. - u2 - CUT_DISTANCE) + (*(verts [v1])) * (u2 + CUT_DISTANCE);
      } else {
        (* _inVertices) [cut._inVertices [0]] = (*(verts [v2])) * (1. - u0 + CUT_DISTANCE) + (*(verts [v1])) * (u0 - CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [1]] = (*(verts [v2])) * (1. - u0 - CUT_DISTANCE) + (*(verts [v1])) * (u0 + CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [4]] = (*(verts [v0])) * (1. - u2 + CUT_DISTANCE) + (*(verts [v3])) * (u2 - CUT_DISTANCE);
        (* _inVertices) [cut._inVertices [5]] = (*(verts [v0])) * (1. - u2 - CUT_DISTANCE) + (*(verts [v3])) * (u2 + CUT_DISTANCE);
      }

      vec point0 = ((* _inVertices) [cut._inVertices [0]] + (* _inVertices) [cut._inVertices [1]]) * .5;

      unsigned int i0, i1, i2 = v3;
      if (upFlag){
        i0 = v0;
        i1 = v2;
      } else {
        i0 = v2;
        i1 = v1;
      }

      // normal for trig 0
      vec edge0 = *(verts [i1]) - *(verts [i0]);
      vec edge1 = *(verts [i2]) - *(verts [i0]);
      vec normal = edge0.cross (edge1);

      real dist1, dist2, max = 0.;
      vec point1, point2, maxpt;

      // get point furthest from edge
      bool collideflag1 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [i0]), *(verts [i1]), *(verts [i2]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [i0]), *(verts [i1]), *(verts [i2]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag1 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      vec2 uv;
      if (collideflag1){
        (* _inVertices) [cut._inVertices [6]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [i0]), *(verts [i1]), *(verts [i2]));
        (* _in3DTexCoords) [cut._inVertices [6]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
          (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
      }

      point0 = ((* _inVertices) [cut._inVertices [4]] + (* _inVertices) [cut._inVertices [5]]) * .5;

      i1 = v2;
      if (upFlag){
        i0 = v3;
        i2 = v1;
      } else {
        i0 = v0;
        i2 = v3;
      }

      // normal for trig 1
      edge0 = *(verts [i1]) - *(verts [i0]);
      edge1 = *(verts [i2]) - *(verts [i0]);
      edge0.fast_cross (normal, edge1);

      // get point furthest from edge
      max = 0.;
      bool collideflag2 = false;
      for (unsigned int i = 0; i < bladeNormals [0]->size (); ++i){
        if (triTriIntersect (normal, *(verts [i0]), *(verts [i1]), *(verts [i2]), (*( bladeNormals [0]))[i],
                             bladeCurr [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i + 1]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
        if (triTriIntersect (normal, *(verts [i0]), *(verts [i1]), *(verts [i2]), (*( bladeNormals [1]))[i],
                             bladePrev [bladeIndices [2*i + 1]], bladePrev [bladeIndices [2*i]], bladeCurr [bladeIndices [2*i]], edge1, point1, point2)){
          collideflag2 = true;
          dist1 = (point1 - point0).length ();
          dist2 = (point2 - point0).length ();
          if (dist1 > dist2){
            if (max < dist1){
              max = dist1;
              maxpt = point1;
            }
          } else {
            if (max < dist2){
              max = dist2;
              maxpt = point2;
            }
          }
        }
      }

      if (collideflag2){
        (* _inVertices) [cut._inVertices [7]] = maxpt;
        calculateBarycentricCoords (uv, maxpt, *(verts [i0]), *(verts [i1]), *(verts [i2]));
        (* _in3DTexCoords) [cut._inVertices [7]] = (* _tex3D) [cell._index [i0]]* uv._v [0] + (* _tex3D) [cell._index [i1]]* uv._v [1] +
          (* _tex3D) [cell._index [i2]]* (1. - uv._v [0] - uv._v [1]);
      }

      // early return
      if (!cell.testAnyExternalFaceFlag ()){
        return;
      }

      // optionally allocate cut variables for external faces
      if (faceFlag1 && faceFlag2){

        if (newflag){

          if (faceFlag0){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (11, 10, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [7];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

              tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];
            }
            else {
              unsigned int tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [7];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

              tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [6]] = 1.;
          }
          else if (faceFlag3){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (11, 10, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

              tmpu = cut._exFaces [7];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];
            }
            else {
              unsigned int tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [7];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [7]] = 1.;
          }
          else {
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();
          }

          if (upFlag){
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [0];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

            tmpu = cut._exFaces [3];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

            tmpu = cut._exFaces [4];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

            tmpu = cut._exFaces [5];
            (* _exFaceIndices) [tmpu] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];
          }
          else {
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [0];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [1];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

            tmpu = cut._exFaces [3];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

            tmpu = cut._exFaces [4];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

            tmpu = cut._exFaces [5];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];
          }

          for (unsigned int i = 0; i < 6; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i]];
          }
          (* _ex2DTexCoords) [cut._exVertices [6]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [cut._exVertices [8]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [cut._exVertices [9]] = (* _tex2D) [cell._index [v3]];

        } // end - if (newflag)

        for (unsigned int i = 0; i < 6; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i]];
        }
        (* _exVertices)[cut._exVertices [6]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [7]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [8]] = *(verts [v2]);
        (* _exVertices)[cut._exVertices [9]] = *(verts [v3]);

        if (faceFlag0 && collideflag1){
          (* _exVertices) [cut._exVertices [10]] = (* _inVertices) [cut._inVertices [6]];

          if (upFlag){
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          } else {
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v2]), *(verts [v1]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v2]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          }
          (* _ex2DTexCoords) [cut._exVertices [10]] = (* _in2DTexCoords) [cut._inVertices [6]];
        }
        else if (faceFlag3 && collideflag2){
          (* _exVertices) [cut._exVertices [10]] = (* _inVertices) [cut._inVertices [7]];

          if (upFlag){
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v3]), *(verts [v2]), *(verts [v1]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v3]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
          } else {
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
          }
          (* _ex2DTexCoords) [cut._exVertices [10]] = (* _in2DTexCoords) [cut._inVertices [7]];
        }

      }
      else if (faceFlag0 && faceFlag3){

        if (newflag){

          if (faceFlag1){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (12, 11, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [11];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [10];
              (* _exFaceIndices) [tmpu] = cut._exVertices [11];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
            }
            else {
              unsigned int tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

              tmpu = cut._exFaces [10];
              (* _exFaceIndices) [tmpu] = cut._exVertices [11];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
            }

            (* _ex2DTexCoords) [cut._exVertices [10]] = (* _in2DTexCoords) [cut._inVertices [2]];
            (* _ex2DTexCoords) [cut._exVertices [11]] = (* _in2DTexCoords) [cut._inVertices [3]];
          }
          else if (faceFlag2){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (12, 11, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [10];
              (* _exFaceIndices) [tmpu] = cut._exVertices [3];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [11];
            }
            else {
              unsigned int tmpu = cut._exFaces [8];
              (* _exFaceIndices) [tmpu] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

              tmpu = cut._exFaces [9];
              (* _exFaceIndices) [tmpu] = cut._exVertices [3];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [11];

              tmpu = cut._exFaces [10];
              (* _exFaceIndices) [tmpu] = cut._exVertices [11];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];
            }

            (* _ex2DTexCoords) [cut._exVertices [10]] = (* _in2DTexCoords) [cut._inVertices [2]];
            (* _ex2DTexCoords) [cut._exVertices [11]] = (* _in2DTexCoords) [cut._inVertices [3]];
          }
          else {
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (10, 8, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();
          }

          if (upFlag){
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

            tmpu = cut._exFaces [3];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

            tmpu = cut._exFaces [4];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

            tmpu = cut._exFaces [5];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

            tmpu = cut._exFaces [6];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

            tmpu = cut._exFaces [7];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];
          }
          else {
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

            tmpu = cut._exFaces [3];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

            tmpu = cut._exFaces [4];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

            tmpu = cut._exFaces [5];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

            tmpu = cut._exFaces [6];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [8];

            tmpu = cut._exFaces [7];
            (* _exFaceIndices) [tmpu] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];
          }

          for (unsigned int i = 0; i < 2; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i]];
          }
          for (unsigned int i = 2; i < 4; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i + 2]];
          }
          (* _ex2DTexCoords) [cut._exVertices [6]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [cut._exVertices [8]] = (* _tex2D) [cell._index [v2]];
          (* _ex2DTexCoords) [cut._exVertices [9]] = (* _tex2D) [cell._index [v3]];

          (* _inSurfaceVertexStatus) [cut._inVertices [6]] = 1.;
          (* _inSurfaceVertexStatus) [cut._inVertices [7]] = 1.;

        } // end - if (newflag)

        for (unsigned int i = 0; i < 2; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i]];
        }
        for (unsigned int i = 2; i < 4; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i + 2]];
        }
        (* _exVertices)[cut._exVertices [6]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [7]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [8]] = *(verts [v2]);
        (* _exVertices)[cut._exVertices [9]] = *(verts [v3]);

        if (collideflag1){
          (* _exVertices) [cut._exVertices [4]] = (* _inVertices) [cut._inVertices [6]];

          if (upFlag){
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
          else {
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v2]), *(verts [v1]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v2]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
        }

        if (collideflag2){
          (* _exVertices) [cut._exVertices [5]] = (* _inVertices) [cut._inVertices [7]];

          if (upFlag){
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v3]), *(verts [v2]), *(verts [v1]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v3]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [5]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
          else {
            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [5]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
        }

        if (faceFlag1 || faceFlag2){
          for (unsigned int i = 2; i < 4; ++i){
            (* _exVertices) [cut._exVertices [i + 8]] = (* _inVertices) [cut._inVertices [i]];
          }
        }

      }
      else if (faceFlag1){

        if (newflag){

          if (faceFlag0){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (9, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];
            }
            else {
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [6]] = 1.;
            (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v3]];
          }
          else if (faceFlag3){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (11, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];
            }
            else {
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [7]] = 1.;
            for (unsigned int i = 4; i < 6; ++i){
              (* _ex2DTexCoords) [cut._exVertices [i + 5]] = (* _in2DTexCoords) [cut._inVertices [i]];
            }
            (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v3]];
          }
          else {
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();
          }

          if (upFlag){
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
          }
          else {
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [0];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
          }

          for (unsigned int i = 0; i < 4; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i]];
          }
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [cut._exVertices [6]] = (* _tex2D) [cell._index [v2]];

        } // end - if (newflag)

        for (unsigned int i = 0; i < 4; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i]];
        }
        (* _exVertices)[cut._exVertices [4]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [6]] = *(verts [v2]);

        if (faceFlag0){
          (* _exVertices)[cut._exVertices [7]] = *(verts [v3]);

          if (collideflag1){
            (* _exVertices) [cut._exVertices [8]] = (* _inVertices) [cut._inVertices [6]];

            if (upFlag){
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            } else {
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v2]), *(verts [v1]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v2]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            }
            (* _ex2DTexCoords) [cut._exVertices [8]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
        }
        else if (faceFlag3){
          (* _exVertices)[cut._exVertices [7]] = *(verts [v3]);

          for (unsigned int i = 4; i < 6; ++i){
            (* _exVertices) [cut._exVertices [i + 5]] = (* _inVertices) [cut._inVertices [i]];
          }

          if (collideflag2){
            (* _exVertices) [cut._exVertices [8]] = (* _inVertices) [cut._inVertices [7]];

            if (upFlag){
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v3]), *(verts [v2]), *(verts [v1]));
              (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v3]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
            } else {
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            }
            (* _ex2DTexCoords) [cut._exVertices [8]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
        }

      }
      else if (faceFlag2){

        if (newflag){

          if (faceFlag3){
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (9, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (upFlag){
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];
            }
            else {
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [7]] = 1.;
            (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v2]];
          }
          else if (faceFlag0){
            *_exUpdateFlag = true;
            _exMutex->lock ();
            cut.allocateExternalVariables (11, 7, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();

            if (!upFlag){
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];
            }
            else {
              unsigned int tmpu = cut._exFaces [3];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

              tmpu = cut._exFaces [4];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

              tmpu = cut._exFaces [5];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

              tmpu = cut._exFaces [6];
              (* _exFaceIndices) [tmpu] = cut._exVertices [8];
              (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
              (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];
            }

            (* _inSurfaceVertexStatus) [cut._inVertices [6]] = 1.;
            for (unsigned int i = 0; i < 2; ++i){
              (* _ex2DTexCoords) [cut._exVertices [i + 9]] = (* _in2DTexCoords) [cut._inVertices [i]];
            }
            (* _ex2DTexCoords) [cut._exVertices [7]] = (* _tex2D) [cell._index [v2]];
          }
          else {
            _exMutex->lock ();
            *_exUpdateFlag = true;
            cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
            _exMutex->unlock ();
          }

          if (upFlag){
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];
          }
          else {
            unsigned int tmpu = cut._exFaces [0];
            (* _exFaceIndices) [tmpu] = cut._exVertices [4];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

            tmpu = cut._exFaces [1];
            (* _exFaceIndices) [tmpu] = cut._exVertices [3];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

            tmpu = cut._exFaces [2];
            (* _exFaceIndices) [tmpu] = cut._exVertices [1];
            (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
            (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];
          }

          for (unsigned int i = 0; i < 4; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i + 2]];
          }
          (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v0]];
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _tex2D) [cell._index [v1]];
          (* _ex2DTexCoords) [cut._exVertices [6]] = (* _tex2D) [cell._index [v3]];

        } // end - if (newflag)

        for (unsigned int i = 0; i < 4; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i + 2]];
        }
        (* _exVertices)[cut._exVertices [4]] = *(verts [v0]);
        (* _exVertices)[cut._exVertices [5]] = *(verts [v1]);
        (* _exVertices)[cut._exVertices [6]] = *(verts [v3]);

        if (faceFlag3){
          (* _exVertices)[cut._exVertices [7]] = *(verts [v2]);

          if (collideflag2){
            (* _exVertices) [cut._exVertices [8]] = (* _inVertices) [cut._inVertices [7]];

            if (upFlag){
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v3]), *(verts [v2]), *(verts [v1]));
              (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v3]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
            } else {
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            }
            (* _ex2DTexCoords) [cut._exVertices [8]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
        }
        else if (faceFlag0){
          (* _exVertices)[cut._exVertices [7]] = *(verts [v2]);

          for (unsigned int i = 0; i < 2; ++i){
            (* _exVertices) [cut._exVertices [i + 9]] = (* _inVertices) [cut._inVertices [i]];
          }

          if (collideflag1){
            (* _exVertices) [cut._exVertices [8]] = (* _inVertices) [cut._inVertices [6]];

            if (upFlag){
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            } else {
              calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v2]), *(verts [v1]), *(verts [v3]));
              (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v2]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
                (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            }
            (* _ex2DTexCoords) [cut._exVertices [8]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
        }

      }
      else if (faceFlag0){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

          (* _inSurfaceVertexStatus) [cut._inVertices [6]] = 1.;

          for (unsigned int i = 0; i < 2; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i]];
          }

          if (upFlag){
            (* _ex2DTexCoords) [cut._exVertices [3]] = (* _tex2D) [cell._index [v0]];
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v2]];
          }
          else {
            (* _ex2DTexCoords) [cut._exVertices [3]] = (* _tex2D) [cell._index [v2]];
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v1]];
          }
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _tex2D) [cell._index [v3]];

        } // end - if (newflag)

        for (unsigned int i = 0; i < 2; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i]];
        }

        (* _exVertices)[cut._exVertices [5]] = *(verts [v3]);

        if (upFlag){
          (* _exVertices)[cut._exVertices [3]] = *(verts [v0]);
          (* _exVertices)[cut._exVertices [4]] = *(verts [v2]);

          if (collideflag1){
            (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [6]];

            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
        }
        else {
          (* _exVertices)[cut._exVertices [3]] = *(verts [v2]);
          (* _exVertices)[cut._exVertices [4]] = *(verts [v1]);

          if (collideflag1){
            (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [6]];

            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [6]], *(verts [v2]), *(verts [v1]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [6]] = (* _tex2D) [cell._index [v2]]* uv._v [0] + (* _tex2D) [cell._index [v1]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [6]];
          }
        }
      }
      else if (faceFlag3){

        if (newflag){
          _exMutex->lock ();
          *_exUpdateFlag = true;
          cut.allocateExternalVariables (6, 4, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          _exMutex->unlock ();

          unsigned int tmpu = cut._exFaces [0];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

          tmpu = cut._exFaces [1];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

          tmpu = cut._exFaces [2];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

          (* _inSurfaceVertexStatus) [cut._inVertices [7]] = 1.;

          for (unsigned int i = 0; i < 2; ++i){
            (* _ex2DTexCoords) [cut._exVertices [i]] = (* _in2DTexCoords) [cut._inVertices [i + 4]];
          }

          if (upFlag){
            (* _ex2DTexCoords) [cut._exVertices [3]] = (* _tex2D) [cell._index [v3]];
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v1]];
          }
          else {
            (* _ex2DTexCoords) [cut._exVertices [3]] = (* _tex2D) [cell._index [v0]];
            (* _ex2DTexCoords) [cut._exVertices [4]] = (* _tex2D) [cell._index [v3]];
          }
          (* _ex2DTexCoords) [cut._exVertices [5]] = (* _tex2D) [cell._index [v2]];

        } // end - if (newflag)

        for (unsigned int i = 0; i < 2; ++i){
          (* _exVertices) [cut._exVertices [i]] = (* _inVertices) [cut._inVertices [i + 4]];
        }
        (* _exVertices)[cut._exVertices [5]] = *(verts [v2]);

        if (upFlag){
          (* _exVertices)[cut._exVertices [3]] = *(verts [v3]);
          (* _exVertices)[cut._exVertices [4]] = *(verts [v1]);

          if (collideflag2){
            (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [7]];

            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v3]), *(verts [v2]), *(verts [v1]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v3]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v1]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
        }
        else {
          (* _exVertices)[cut._exVertices [3]] = *(verts [v0]);
          (* _exVertices)[cut._exVertices [4]] = *(verts [v3]);

          if (collideflag2){
            (* _exVertices) [cut._exVertices [2]] = (* _inVertices) [cut._inVertices [7]];

            calculateBarycentricCoords (uv, (* _inVertices) [cut._inVertices [7]], *(verts [v0]), *(verts [v2]), *(verts [v3]));
            (* _in2DTexCoords) [cut._inVertices [7]] = (* _tex2D) [cell._index [v0]]* uv._v [0] + (* _tex2D) [cell._index [v2]]* uv._v [1] +
              (* _tex2D) [cell._index [v3]]* (1. - uv._v [0] - uv._v [1]);
            (* _ex2DTexCoords) [cut._exVertices [2]] = (* _in2DTexCoords) [cut._inVertices [7]];
          }
        }
      }

    }

    /** Private method to perform cut operations for three cut-edges cases with severance. Arguments:
    * 0,1,2   : u-indices of the three cut-edges
    * 3,4,5   : Flags to signify if faces are external (right, left, center)
    * 6,7,8,9 : Local vertex indices in order (cut edges in order: 6-8, 6-7 and 6-9)
    * 10      : Reference to cell-structure for which current method is being called
    * 11      : Reference to the cut-structure for current cell
    */
    void
    Partition::performFinishedThreeEdgeCut (real u0, real u1, real u2, bool faceFlag0, bool faceFlag1, bool faceFlag2,
                                            unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, Cell &cell, Cut &cut)
    {
      cell.finalize ();

      if (cut._numInVertices > 6){
        cut.deallocate (6, cut._numInVertices, cut._inVertices, _inEmptyVertices);
      }
      if (cut._numInFaces > 2){
        cut.deallocate (2, cut._numInFaces, cut._inFaces, _inEmptyFaces);
      }

      _inMutex->lock ();
      *_inUpdateFlag = true;
      if (cut._numInVertices < 6 || cut._numInFaces < 2){
        cut.allocateInternalVariables (6, 2, *_inVertices, *_inSurfaceVertexStatus, *_in2DTexCoords, *_in3DTexCoords, *_inFaceIndices, _inEmptyVertices, _inEmptyFaces);
      }
      cut.allocateInternalUVCoords (* _inUVCoords);
      _inMutex->unlock ();

      unsigned int inds [6] = {cut._inVertices [0], cut._inVertices [1], cut._inVertices [2], cut._inVertices [3], cut._inVertices [4], cut._inVertices [5]};

      if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v1))){
        (* _inSurfaceVertexStatus) [inds [0]] = 1.;
        (* _inSurfaceVertexStatus) [inds [1]] = 1.;
        (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v0]] * (1. - u0) + (* _tex2D) [cell._index [v1]] * u0;
        (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
      }
      (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v0]] * (1. - u0) + (* _tex3D) [cell._index [v1]] * u0;
      (* _in3DTexCoords) [inds [1]] = (* _in3DTexCoords) [inds [0]];

      if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v2))){
        (* _inSurfaceVertexStatus) [inds [2]] = 1.;
        (* _inSurfaceVertexStatus) [inds [3]] = 1.;
        (* _in2DTexCoords) [inds [2]] = (* _tex2D) [cell._index [v0]] * (1. - u1) + (* _tex2D) [cell._index [v2]] * u1;
        (* _in2DTexCoords) [inds [3]] = (* _in2DTexCoords) [inds [2]];
      }
      (* _in3DTexCoords) [inds [2]] = (* _tex3D) [cell._index [v0]] * (1. - u1) + (* _tex3D) [cell._index [v2]] * u1;
      (* _in3DTexCoords) [inds [3]] = (* _in3DTexCoords) [inds [2]];

      if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v3))){
        (* _inSurfaceVertexStatus) [inds [4]] = 1.;
        (* _inSurfaceVertexStatus) [inds [5]] = 1.;
        (* _in2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v0]] * (1. - u2) + (* _tex2D) [cell._index [v3]] * u2;
        (* _in2DTexCoords) [inds [5]] = (* _in2DTexCoords) [inds [4]];
      }
      (* _in3DTexCoords) [inds [4]] = (* _tex3D) [cell._index [v0]] * (1. - u2) + (* _tex3D) [cell._index [v3]] * u2;
      (* _in3DTexCoords) [inds [5]] = (* _in3DTexCoords) [inds [4]];

      unsigned int tmpu = cut._inFaces [0];
      (* _inFaceIndices) [tmpu] = inds [0];
      (* _inFaceIndices) [tmpu + 1] = inds [4];
      (* _inFaceIndices) [tmpu + 2] = inds [2];

      tmpu = cut._inFaces [1];
      (* _inFaceIndices) [tmpu] = inds [1];
      (* _inFaceIndices) [tmpu + 1] = inds [3];
      (* _inFaceIndices) [tmpu + 2] = inds [5];

      switch (v0){
        case 0:
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (1. - u0 + CUT_DISTANCE, 0., u0 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (1. - u0 - CUT_DISTANCE, 0., u0 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (1. - u1 + CUT_DISTANCE, u1 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (1. - u1 - CUT_DISTANCE, u1 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (1. - u2 + CUT_DISTANCE, 0., 0.);
          (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (1. - u2 - CUT_DISTANCE, 0., 0.);
          break;
        case 1:
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (0., 1. - u0 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (0., 1. - u0 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (u1 - CUT_DISTANCE, 1. - u1 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (u1 + CUT_DISTANCE, 1. - u1 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (0., 1. - u2 + CUT_DISTANCE, u2 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (0., 1. - u2 - CUT_DISTANCE, u2 + CUT_DISTANCE);
          break;
        case 2:
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (0., u0 - CUT_DISTANCE, 1. - u0 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (0., u0 + CUT_DISTANCE, 1. - u0 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (u1 - CUT_DISTANCE, 0., 1. - u1 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (u1 + CUT_DISTANCE, 0., 1. - u1 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (0., 0., 1. - u2 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (0., 0., 1. - u2 - CUT_DISTANCE);
          break;
        case 3:
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (u0 - CUT_DISTANCE, 0., 0.);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (u0 + CUT_DISTANCE, 0., 0.);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (0., u1 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (0., u1 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (0., 0., u2 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (0., 0., u2 + CUT_DISTANCE);
          break;
      }

      // early return;
      if (!(faceFlag0 || faceFlag1 || faceFlag2)){
        return;
      }

      // optionally allocate cut variables for external faces
      if (faceFlag0 && faceFlag1 && faceFlag2){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 10 || cut._numExFaces != 9){
          cut.allocateExternalVariables (10, 9, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 6; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i];
        }

        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [7] = (* _tex2D) [cell._index [v1]];
        (* _ex2DTexCoords) [8] = (* _tex2D) [cell._index [v2]];
        (* _ex2DTexCoords) [9] = (* _tex2D) [cell._index [v3]];

        switch (v0){
          case 0:
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            (* _exUVCoords) [7] = vec3 (0., 0., 1.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            break;
          case 2:
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            (* _exUVCoords) [7] = vec3 (1., 0., 0.);
            (* _exUVCoords) [8] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
            break;
        }

        unsigned int einds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
          cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = einds [0];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [4];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [0];
        (* _exFaceIndices) [tmpu + 2] = einds [4];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [7];
        (* _exFaceIndices) [tmpu + 2] = einds [8];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [3];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [9];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [6];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [5];
        (* _exFaceIndices) [tmpu + 2] = einds [9];

        tmpu = cut._exFaces [7];
        (* _exFaceIndices) [tmpu] = einds [5];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [9];

        tmpu = cut._exFaces [8];
        (* _exFaceIndices) [tmpu] = einds [9];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [8];
      }
      else if (faceFlag0 && faceFlag1){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 10 || cut._numExFaces != 6){
          cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 6; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i];
        }

        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [7] = (* _tex2D) [cell._index [v1]];
        (* _ex2DTexCoords) [8] = (* _tex2D) [cell._index [v2]];
        (* _ex2DTexCoords) [9] = (* _tex2D) [cell._index [v3]];

        switch (v0){
          case 0:
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            (* _exUVCoords) [7] = vec3 (0., 0., 1.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            break;
          case 2:
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            (* _exUVCoords) [7] = vec3 (1., 0., 0.);
            (* _exUVCoords) [8] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
            break;
        }

        unsigned int einds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
          cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = einds [0];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [4];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [7];
        (* _exFaceIndices) [tmpu + 2] = einds [8];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [3];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [9];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [5];
        (* _exFaceIndices) [tmpu + 2] = einds [9];
      }
      else if (faceFlag0 && faceFlag2){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 10 || cut._numExFaces != 6){
          cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 6; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i];
        }

        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [7] = (* _tex2D) [cell._index [v1]];
        (* _ex2DTexCoords) [8] = (* _tex2D) [cell._index [v2]];
        (* _ex2DTexCoords) [9] = (* _tex2D) [cell._index [v3]];

        switch (v0){
          case 0:
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            (* _exUVCoords) [7] = vec3 (0., 0., 1.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            break;
          case 2:
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            (* _exUVCoords) [7] = vec3 (1., 0., 0.);
            (* _exUVCoords) [8] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
            break;
        }

        unsigned int einds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
          cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = einds [0];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [0];
        (* _exFaceIndices) [tmpu + 2] = einds [4];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [7];
        (* _exFaceIndices) [tmpu + 2] = einds [8];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = einds [1];
        (* _exFaceIndices) [tmpu + 1] = einds [3];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = einds [5];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [9];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = einds [9];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [8];
      }
      else if (faceFlag1 && faceFlag2){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 10 || cut._numExFaces != 6){
          cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 6; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i];
        }

        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [7] = (* _tex2D) [cell._index [v1]];
        (* _ex2DTexCoords) [8] = (* _tex2D) [cell._index [v2]];
        (* _ex2DTexCoords) [9] = (* _tex2D) [cell._index [v3]];

        switch (v0){
          case 0:
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            (* _exUVCoords) [7] = vec3 (0., 0., 1.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            break;
          case 2:
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            (* _exUVCoords) [7] = vec3 (1., 0., 0.);
            (* _exUVCoords) [8] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [7] = vec3 (0., 1., 0.);
            (* _exUVCoords) [8] = vec3 (1., 0., 0.);
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
            break;
        }

        unsigned int einds [10] = {cut._exVertices [0], cut._exVertices [1], cut._exVertices [2], cut._exVertices [3], cut._exVertices [4],
          cut._exVertices [5], cut._exVertices [6], cut._exVertices [7], cut._exVertices [8], cut._exVertices [9]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [4];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [6];
        (* _exFaceIndices) [tmpu + 1] = einds [0];
        (* _exFaceIndices) [tmpu + 2] = einds [4];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [9];
        (* _exFaceIndices) [tmpu + 2] = einds [7];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = einds [3];
        (* _exFaceIndices) [tmpu + 1] = einds [5];
        (* _exFaceIndices) [tmpu + 2] = einds [9];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = einds [5];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [9];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = einds [9];
        (* _exFaceIndices) [tmpu + 1] = einds [1];
        (* _exFaceIndices) [tmpu + 2] = einds [8];

      }
      else if (faceFlag0){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 7 || cut._numExFaces != 3){
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i];
        }

        (* _ex2DTexCoords) [4] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [5] = (* _tex2D) [cell._index [v1]];
        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v2]];

        switch (v0){
          case 0:
            (* _exUVCoords) [4] = vec3 (1., 0., 0.);
            (* _exUVCoords) [5] = vec3 (0., 1., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [4] = vec3 (0., 1., 0.);
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            break;
          case 2:
            (* _exUVCoords) [4] = vec3 (0., 0., 1.);
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [5] = vec3 (0., 1., 0.);
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
            break;
        }

        unsigned int einds [2] = {cut._exVertices [1], cut._exVertices [5]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

      }
      else if (faceFlag1){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 7 || cut._numExFaces != 3){
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i + 2];
          (* _ex2DTexCoords) [i] = (* _in2DTexCoords) [i + 2];
        }

        (* _ex2DTexCoords) [4] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [5] = (* _tex2D) [cell._index [v2]];
        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v3]];

        switch (v0){
          case 0:
            (* _exUVCoords) [4] = vec3 (1., 0., 0.);
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            break;
          case 1:
            (* _exUVCoords) [4] = vec3 (0., 1., 0.);
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            break;
          case 2:
            (* _exUVCoords) [4] = vec3 (0., 0., 1.);
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
            break;
          case 3:
            (* _exUVCoords) [5] = vec3 (0., 1., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            break;
        }

        unsigned int einds [2] = {cut._exVertices [1], cut._exVertices [6]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

      }
      else if (faceFlag2){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 7 || cut._numExFaces != 3){
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        (* _ex2DTexCoords) [0] = (* _in2DTexCoords) [0];
        (* _ex2DTexCoords) [1] = (* _in2DTexCoords) [1];
        (* _ex2DTexCoords) [2] = (* _in2DTexCoords) [4];
        (* _ex2DTexCoords) [3] = (* _in2DTexCoords) [5];
        (* _ex2DTexCoords) [4] = (* _tex2D) [cell._index [v0]];
        (* _ex2DTexCoords) [5] = (* _tex2D) [cell._index [v3]];
        (* _ex2DTexCoords) [6] = (* _tex2D) [cell._index [v1]];


        (* _exUVCoords) [0] = (* _inUVCoords) [0];
        (* _exUVCoords) [1] = (* _inUVCoords) [1];
        (* _exUVCoords) [2] = (* _inUVCoords) [4];
        (* _exUVCoords) [3] = (* _inUVCoords) [5];

        switch (v0){
          case 0:
            (* _exUVCoords) [4] = vec3 (1., 0., 0.);
            (* _exUVCoords) [5] = vec3 (0., 0., 1.);
            break;
          case 1:
            (* _exUVCoords) [4] = vec3 (0., 1., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
            break;
          case 2:
            (* _exUVCoords) [4] = vec3 (0., 0., 1.);
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
            break;
          case 3:
            (* _exUVCoords) [5] = vec3 (0., 0., 1.);
            (* _exUVCoords) [6  ] = vec3 (1., 0., 0.);
            break;
        }

        unsigned int einds [2] = {cut._exVertices [1], cut._exVertices [6]};

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
        (* _exFaceIndices) [tmpu + 2] = einds [2];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = einds [0];
        (* _exFaceIndices) [tmpu + 1] = einds [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];
      }
    }

    /** Private method to perform cut operations for four cut-edge cases with severance. Arguments:
    * 0,1,2,3   : u-indices of the four cut-edges
    * 4,5,6,7   : Flags to signify if faces are external (front, back, left, right)
    * 8,9,10,11 : Local vertex indices in order
    * 12        : Reference to cell-structure for which current method is being called
    * 13        : Reference to the cut-structure for current cell
    */
    void
    Partition::performFourEdgeCut (real u0, real u1, real u2, real u3, bool faceFlag0, bool faceFlag1, bool faceFlag2, bool faceFlag3,
                                   unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, Cell &cell, Cut &cut)
    {
      cell.finalize ();

      if (cut._numInVertices > 8){
        cut.deallocate (8, cut._numInVertices, cut._inVertices, _inEmptyVertices);
      }
      if (cut._numInFaces > 4){
        cut.deallocate (4, cut._numInFaces, cut._inFaces, _inEmptyFaces);
      }

      _inMutex->lock ();
      *_inUpdateFlag = true;
      if (cut._numInVertices < 8 || cut._numInFaces < 8){
        cut.allocateInternalVariables (8, 4, *_inVertices, *_inSurfaceVertexStatus, *_in2DTexCoords, *_in3DTexCoords, *_inFaceIndices, _inEmptyVertices, _inEmptyFaces);
      }
      cut.allocateInternalUVCoords (* _inUVCoords);
      _inMutex->unlock ();

      unsigned int inds [8] = {cut._inVertices [0], cut._inVertices [1], cut._inVertices [2], cut._inVertices [3],
        cut._inVertices [4], cut._inVertices [5], cut._inVertices [6], cut._inVertices [7]};

      if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v1))){
        (* _inSurfaceVertexStatus) [inds [0]] = 1.;
        (* _inSurfaceVertexStatus) [inds [1]] = 1.;
        (* _in2DTexCoords) [inds [0]] = (* _tex2D) [cell._index [v0]] * (1. - u0) + (* _tex2D) [cell._index [v1]] * u0;
        (* _in2DTexCoords) [inds [1]] = (* _in2DTexCoords) [inds [0]];
      }
      (* _in3DTexCoords) [inds [0]] = (* _tex3D) [cell._index [v0]] * (1. - u0) + (* _tex3D) [cell._index [v1]] * u0;
      (* _in3DTexCoords) [inds [1]] = (* _in3DTexCoords) [inds [0]];

      if (cell.testExternalVertexFlag (v0) && cell.testExternalVertexFlag ((v2))){
        (* _inSurfaceVertexStatus) [inds [2]] = 1.;
        (* _inSurfaceVertexStatus) [inds [3]] = 1.;
        (* _in2DTexCoords) [inds [2]] = (* _tex2D) [cell._index [v0]] * (1. - u1) + (* _tex2D) [cell._index [v2]] * u1;
        (* _in2DTexCoords) [inds [3]] = (* _in2DTexCoords) [inds [2]];
      }
      (* _in3DTexCoords) [inds [2]] = (* _tex3D) [cell._index [v0]] * (1. - u1) + (* _tex3D) [cell._index [v2]] * u1;
      (* _in3DTexCoords) [inds [3]] = (* _in3DTexCoords) [inds [2]];

      if (cell.testExternalVertexFlag (v3) && cell.testExternalVertexFlag ((v2))){
        (* _inSurfaceVertexStatus) [inds [4]] = 1.;
        (* _inSurfaceVertexStatus) [inds [5]] = 1.;
        (* _in2DTexCoords) [inds [4]] = (* _tex2D) [cell._index [v3]] * (1. - u2) + (* _tex2D) [cell._index [v2]] * u2;
        (* _in2DTexCoords) [inds [5]] = (* _in2DTexCoords) [inds [4]];
      }
      (* _in3DTexCoords) [inds [4]] = (* _tex3D) [cell._index [v3]] * (1. - u2) + (* _tex3D) [cell._index [v2]] * u2;
      (* _in3DTexCoords) [inds [5]] = (* _in3DTexCoords) [inds [4]];

      if (cell.testExternalVertexFlag (v3) && cell.testExternalVertexFlag ((v1))){
        (* _inSurfaceVertexStatus) [inds [6]] = 1.;
        (* _inSurfaceVertexStatus) [inds [7]] = 1.;
        (* _in2DTexCoords) [inds [6]] = (* _tex2D) [cell._index [v3]] * (1. - u3) + (* _tex2D) [cell._index [v1]] * u3;
        (* _in2DTexCoords) [inds [7]] = (* _in2DTexCoords) [inds [6]];
      }
      (* _in3DTexCoords) [inds [6]] = (* _tex3D) [cell._index [v3]] * (1. - u3) + (* _tex3D) [cell._index [v1]] * u3;
      (* _in3DTexCoords) [inds [7]] = (* _in3DTexCoords) [inds [6]];

      unsigned int tmpu = cut._inFaces [0];
      (* _inFaceIndices) [tmpu] = inds [0];
      (* _inFaceIndices) [tmpu + 1] = inds [4];
      (* _inFaceIndices) [tmpu + 2] = inds [2];

      tmpu = cut._inFaces [1];
      (* _inFaceIndices) [tmpu] = inds [0];
      (* _inFaceIndices) [tmpu + 1] = inds [6];
      (* _inFaceIndices) [tmpu + 2] = inds [4];

      tmpu = cut._inFaces [2];
      (* _inFaceIndices) [tmpu] = inds [1];
      (* _inFaceIndices) [tmpu + 1] = inds [3];
      (* _inFaceIndices) [tmpu + 2] = inds [5];

      tmpu = cut._inFaces [3];
      (* _inFaceIndices) [tmpu] = inds [1];
      (* _inFaceIndices) [tmpu + 1] = inds [5];
      (* _inFaceIndices) [tmpu + 2] = inds [7];

      if (v0 == 0){
        if (v1 == 1){
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (1. - u0 + CUT_DISTANCE, u0 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (1. - u0 - CUT_DISTANCE, u0 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (1. - u1 + CUT_DISTANCE, 0., u1 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (1. - u1 - CUT_DISTANCE, 0., u1 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (0., 0., u2 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (0., 0., u2 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [6]] = vec3 (0., u3 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [7]] = vec3 (0., u3 + CUT_DISTANCE, 0.);
        } else {
          (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (1. - u0 + CUT_DISTANCE, 0., 0.);
          (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (1. - u0 - CUT_DISTANCE, 0., 0.);
          (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (1. - u1 + CUT_DISTANCE, u1 - CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (1. - u1 - CUT_DISTANCE, u1 + CUT_DISTANCE, 0.);
          (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (0., u2 - CUT_DISTANCE, 1. - u2 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (0., u2 + CUT_DISTANCE, 1. - u2 - CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [6]] = vec3 (0., 0., 1. - u3 + CUT_DISTANCE);
          (*_inUVCoords) [cut._inUVCoords [7]] = vec3 (0., 0., 1. - u3 - CUT_DISTANCE);
        }
      }
      else {
        (*_inUVCoords) [cut._inUVCoords [0]] = vec3 (u0 - CUT_DISTANCE, 0., 1. - u0 + CUT_DISTANCE);
        (*_inUVCoords) [cut._inUVCoords [1]] = vec3 (u0 + CUT_DISTANCE, 0., 1. - u0 - CUT_DISTANCE);
        (*_inUVCoords) [cut._inUVCoords [2]] = vec3 (0., u1 - CUT_DISTANCE, 1. - u1 + CUT_DISTANCE);
        (*_inUVCoords) [cut._inUVCoords [3]] = vec3 (0., u1 + CUT_DISTANCE, 1. - u1 - CUT_DISTANCE);
        (*_inUVCoords) [cut._inUVCoords [4]] = vec3 (0., u2 - CUT_DISTANCE, 0.);
        (*_inUVCoords) [cut._inUVCoords [5]] = vec3 (0., u2 + CUT_DISTANCE, 0.);
        (*_inUVCoords) [cut._inUVCoords [6]] = vec3 (u3 - CUT_DISTANCE, 0., 0.);
        (*_inUVCoords) [cut._inUVCoords [7]] = vec3 (u3 + CUT_DISTANCE, 0., 0.);
      }

      // early return
      if (!cell.testAnyExternalFaceFlag ()){
        return;
      }

      // optionally allocate cut variables for external faces
      if (faceFlag0 && faceFlag2){

        bool flag1 = false, flag2 = false;
        if (faceFlag1){
          flag1 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 9){
            cut.allocateExternalVariables (12, 9, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else if (faceFlag3){
          flag2 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 9){
            cut.allocateExternalVariables (12, 9, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else {
          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 6){
            cut.allocateExternalVariables (12, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }

        for (unsigned int i = 0; i < 8; ++i){
          (* _exUVCoords) [cut._exUVCoords [i]] = (* _inUVCoords) [cut._inUVCoords [i]];
        }

        if (v0 == 0){
          (* _exUVCoords) [9] = vec3 (1., 0., 0.);
          if (v1 == 1){
            (* _exUVCoords) [10] = vec3 (0., 1., 0.);
            (* _exUVCoords) [11] = vec3 (0., 0., 1.);
          } else {
            (* _exUVCoords) [11] = vec3 (0., 1., 0.);
            (* _exUVCoords) [12] = vec3 (0., 0., 1.);
          }
        }
        else {
          (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          (* _exUVCoords) [10] = vec3 (1., 0., 0.);
          (* _exUVCoords) [11] = vec3 (0., 1., 0.);
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [9];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [11];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = cut._exVertices [12];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = cut._exVertices [5];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = cut._exVertices [5];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

        if (flag1){
          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = cut._exVertices [3];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

          tmpu = cut._exFaces [8];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [12];
        }
        else if (flag2){
          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = cut._exVertices [10];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [12];

          tmpu = cut._exFaces [8];
          (* _exFaceIndices) [tmpu] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [12];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];
        }
      }
      else if (faceFlag1 && faceFlag3){

        bool flag1 = false, flag2 = false;
        if (faceFlag0){
          flag1 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 9){
            cut.allocateExternalVariables (12, 9, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else if (faceFlag2){
          flag2 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 9){
            cut.allocateExternalVariables (12, 9, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else {
          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 12 || cut._numExFaces != 6){
            cut.allocateExternalVariables (12, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }

        for (unsigned int i = 0; i < 8; ++i){
          (* _exUVCoords) [cut._exUVCoords [i]] = (* _inUVCoords) [cut._inUVCoords [i]];
        }

        if (v0 == 0){
          (* _exUVCoords) [9] = vec3 (1., 0., 0.);
          if (v1 == 1){
            (* _exUVCoords) [10] = vec3 (0., 1., 0.);
            (* _exUVCoords) [11] = vec3 (0., 0., 1.);
          } else {
            (* _exUVCoords) [11] = vec3 (0., 1., 0.);
            (* _exUVCoords) [12] = vec3 (0., 0., 1.);
          }
        }
        else {
          (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          (* _exUVCoords) [10] = vec3 (1., 0., 0.);
          (* _exUVCoords) [11] = vec3 (0., 1., 0.);
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [3];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [9];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [9];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [12];

        tmpu = cut._exFaces [3];
        (* _exFaceIndices) [tmpu] = cut._exVertices [10];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];

        tmpu = cut._exFaces [4];
        (* _exFaceIndices) [tmpu] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [12];

        tmpu = cut._exFaces [5];
        (* _exFaceIndices) [tmpu] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [12];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        if (flag1){
          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = cut._exVertices [1];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [11];

          tmpu = cut._exFaces [8];
          (* _exFaceIndices) [tmpu] = cut._exVertices [1];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];
        }
        else if (flag2){
          tmpu = cut._exFaces [6];
          (* _exFaceIndices) [tmpu] = cut._exVertices [12];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [4];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

          tmpu = cut._exFaces [7];
          (* _exFaceIndices) [tmpu] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [11];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [10];

          tmpu = cut._exFaces [8];
          (* _exFaceIndices) [tmpu] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [10];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [7];
        }

      }
      else if (faceFlag0){

        bool flag1 = false, flag2 = false;
        if (faceFlag1){
          flag1 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 10 || cut._numExFaces != 6){
            cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else if (faceFlag3){
          flag2 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 10 || cut._numExFaces != 6){
            cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else {
          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 7 || cut._numExFaces != 3){
            cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }

        for (unsigned int i = 0; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i];
        }
        if (v0 == 0){
          (* _exUVCoords) [4] = vec3 (1., 0., 0.);
          if (v1 == 1){
            (* _exUVCoords) [5] = vec3 (0., 1., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
          } else {
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
          }
        } else {
          (* _exUVCoords) [4] = vec3 (0., 0., 1.);
          (* _exUVCoords) [5] = vec3 (1., 0., 0.);
          (* _exUVCoords) [6] = vec3 (0., 1., 0.);
        }

        if (flag1){
          for (unsigned int i = 4; i < 6; ++i){
            (* _exUVCoords) [i + 3] = (* _inUVCoords) [i];
          }
          if (v1 == 3){
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          }
        }
        else if (flag2){
          for (unsigned int i = 6; i < 8; ++i){
            (* _exUVCoords) [i + 1] = (* _inUVCoords) [i];
          }
          if (v1 == 3){
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          }
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

        if (flag1){
          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [6];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
        }
        else if (flag2){
          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [8];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
        }
      }
      else if (faceFlag2){

        bool flag1 = false, flag2 = false;
        if (faceFlag1){
          flag1 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 10 || cut._numExFaces != 6){
            cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else if (faceFlag3){
          flag2 = true;

          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 10 || cut._numExFaces != 6){
            cut.allocateExternalVariables (10, 6, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }
        else {
          _exMutex->lock ();
          *_exUpdateFlag = true;
          if (cut._numExVertices != 7 || cut._numExFaces != 3){
            cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
          }
          cut.allocateExternalUVCoords (* _exUVCoords);
          _exMutex->unlock ();
        }

        for (unsigned int i = 0; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i + 4];
        }
        if (v0 == 0){
          if (v1 == 1){
            (* _exUVCoords) [5] = vec3 (0., 0., 1.);
            (* _exUVCoords) [6] = vec3 (0., 1., 0.);
          } else {
            (* _exUVCoords) [4] = vec3 (0., 0., 1.);
            (* _exUVCoords) [5] = vec3 (0., 1., 0.);
          }
        } else {
          (* _exUVCoords) [5] = vec3 (0., 1., 0.);
          (* _exUVCoords) [6] = vec3 (1., 0., 0.);
        }

        if (flag1){
          for (unsigned int i = 2; i < 4; ++i){
            (* _exUVCoords) [i + 5] = (* _inUVCoords) [i];
          }
          if (v0 == 0){
            (* _exUVCoords) [9] = vec3 (1., 0., 0.);
          } else {
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          }
        }
        else if (flag2){
          for (unsigned int i = 0; i < 2; ++i){
            (* _exUVCoords) [i + 7] = (* _inUVCoords) [i];
          }
          if (v0 == 0){
            (* _exUVCoords) [9] = vec3 (1., 0., 0.);
          } else {
            (* _exUVCoords) [9] = vec3 (0., 0., 1.);
          }
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [2];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [1];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [6];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

        if (flag2){
          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [6];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [8];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [3];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [9];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [2];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
        }
        else if (flag1){
          tmpu = cut._exFaces [3];
          (* _exFaceIndices) [tmpu] = cut._exVertices [8];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [5];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [7];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [0];

          tmpu = cut._exFaces [4];
          (* _exFaceIndices) [tmpu] = cut._exVertices [9];
          (* _exFaceIndices) [tmpu + 1] = cut._exVertices [0];
          (* _exFaceIndices) [tmpu + 2] = cut._exVertices [4];
        }
      }
      else if (faceFlag1){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 7 || cut._numExFaces != 3){
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i + 2];
        }
        if (v0 == 0){
          if (v1 == 1){
            (* _exUVCoords) [4] = vec3 (0., 0., 1.);
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
          } else {
            (* _exUVCoords) [4] = vec3 (0., 1., 0.);
            (* _exUVCoords) [5] = vec3 (0., 0., 1.);
            (* _exUVCoords) [6] = vec3 (1., 0., 0.);
          }
        }
        else {
          (* _exUVCoords) [4] = vec3 (0., 1., 0.);
          (* _exUVCoords) [6] = vec3 (0., 0., 1.);
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [6];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];
      }
      else if (faceFlag3){

        _exMutex->lock ();
        *_exUpdateFlag = true;
        if (cut._numExVertices != 7 || cut._numExFaces != 3){
          cut.allocateExternalVariables (7, 3, *_exVertices, *_ex2DTexCoords, *_exFaceIndices, _exEmptyVertices, _exEmptyFaces);
        }
        cut.allocateExternalUVCoords (* _exUVCoords);
        _exMutex->unlock ();

        for (unsigned int i = 0; i < 2; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i + 6];
        }
        for (unsigned int i = 2; i < 4; ++i){
          (* _exUVCoords) [i] = (* _inUVCoords) [i - 2];
        }

        if (v0 == 0){
          if (v0 == 1){
            (* _exUVCoords) [4] = vec3 (0., 1., 0.);
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
          } else {
            (* _exUVCoords) [5] = vec3 (1., 0., 0.);
            (* _exUVCoords) [6] = vec3 (0., 0., 1.);
          }
        }
        else {
          (* _exUVCoords) [4] = vec3 (0., 0., 1.);
          (* _exUVCoords) [5] = vec3 (1., 0., 0.);
        }

        unsigned int tmpu = cut._exFaces [0];
        (* _exFaceIndices) [tmpu] = cut._exVertices [4];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [3];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [1];

        tmpu = cut._exFaces [1];
        (* _exFaceIndices) [tmpu] = cut._exVertices [0];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [6];

        tmpu = cut._exFaces [2];
        (* _exFaceIndices) [tmpu] = cut._exVertices [6];
        (* _exFaceIndices) [tmpu + 1] = cut._exVertices [2];
        (* _exFaceIndices) [tmpu + 2] = cut._exVertices [5];
      }
    }
  }
}
