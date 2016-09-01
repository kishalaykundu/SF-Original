/**
 * @file Partition.h
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

#pragma once

#include <stack>
#include <vector>
#include <forward_list>
#include <boost/thread/mutex.hpp>

#include "Preprocess.h"

#include "aabb.h"
#include "vec2.h"

#include "Vertex.h"
#include "Cut.h"

using namespace std;

namespace SF {
  namespace XFE {

	  class Cell;
	  class Face;
	  class Edge;
	  class Vertex;

    class Partition {

    public:
      aabb _bbox;

      // start and end index into the submesh cell arrays
      unsigned int _cellStartIndex, _cellEndIndex;

      // start and end index into the submesh face arrays
      unsigned int _exFaceStartIndex, _exFaceEndIndex;
      unsigned int _inFaceStartIndex, _inFaceEndIndex;

      // indices of cells that have undergone partitial cuts
      forward_list <unsigned int> _cutCells;
      forward_list <unsigned int> _reExaminedCells;
      forward_list <unsigned int> _finishedCells;
      forward_list <unsigned int> _collidingVertices;

      vector <Cut> _cuts;

      vector <Vertex> *_vertInfo;
      vector <vec2> *_tex2D;
      vector <vec3> *_tex3D;

      boost::mutex *_exMutex;
      bool *_exUpdateFlag;
      vector <vec> *_exVertices;
      vector <vec3> *_exUVCoords;
      vector <vec2> *_ex2DTexCoords;
      vector <unsigned int> *_exFaceIndices;

      boost::mutex *_inMutex;
      bool *_inUpdateFlag;
      vector <vec> *_inVertices;
      vector <vec3> *_inUVCoords;
      vector <float> *_inSurfaceVertexStatus;
      vector <vec2> *_in2DTexCoords;
      vector <vec3> *_in3DTexCoords;

      vector <unsigned int> *_inFaceIndices;

      stack <unsigned int> _inEmptyVertices;
      stack <unsigned int> _inEmptyFaces;
      stack <unsigned int> _exEmptyVertices;
      stack <unsigned int> _exEmptyFaces;

      Partition ();
      ~Partition ();

      Partition (const Partition &p);

      Partition & operator = (const Partition &p);

      // method to get affected cells
      void gatherAffectedCells (unsigned int sIndex, vector <Vertex> &vertexInfo, vector <vec> &verts, vector <unsigned int> &indices, vector <Face> &faces,
                                vector <unsigned int> &iindices, vector <Face> &ifaces, vector <Edge> &edges, vector <Cell> &cells,
                                vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void finalizeCollision (vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells,
                              vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

    private:
      void resolveReExaminedCells (vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells,
                                   vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void cellBladeCollide (unsigned int sIndex, forward_list <unsigned int>::iterator &lastElement, vector <Vertex> &vInfo,
                             vector <vec> &verts, vector <Edge> &edges, vector <Cell> &cells, Cell &cell,
                             vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void formFaces (Cell &cell, Cut &cut, vector <Edge> &edges, vector <vec> &verts,
                      vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void performOneEdgeCut (real eu, bool faceFlag0, bool faceFlag1, unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3,
                              vec *verts [4], Cell &cell, Cut &cut,
                              vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void performTwoEdgeCut (real eu0, real eu1, bool faceFlag0, bool faceFlag1, bool faceFlag2,
                              unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, vec *verts [4], Cell &cell, Cut &cut,
                              vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void performUnfinishedThreeEdgeCut (bool upFlag, real eu0, real eu1, real eu2, bool faceFlag0, bool faceFlag1, bool faceFlag2, bool faceFlag3,
                                          unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, vec *verts [4], Cell &cell, Cut &cut,
                                          vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

      void performFinishedThreeEdgeCut (real eu0, real eu1, real eu2, bool faceFlag0, bool faceFlag1, bool faceFlag2,
                                        unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, Cell &cell, Cut &cut);

      void performFourEdgeCut (real u0, real u1, real u2, real u3, bool faceFlag0, bool faceFlag1, bool faceFlag2, bool faceFlag3,
                               unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, Cell &cell, Cut &cut);
    };
  }
}
