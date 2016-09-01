/**
 * @file Submesh.cpp
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
 * The submesh class for the CU_XFEM library.
 */

#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "GL/common.h"
#include "Collide/lineTriCollide.h"

#include "Common.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "Partition.h"
#include "Submesh.h"

namespace SF {
	namespace XFE {

	  //static method to read cell info related files
	  static bool
    readCellFiles (const string &prefix, unsigned int numVerts, vector <Cell> &cells)
    {
      int tmpd [4];

      string file (prefix);
      file.append (".tet.ele");

      FILE* fp = fopen (file.c_str (), "r");
      assert (fp);

      int status = fscanf (fp, "%d\n", &(tmpd [0]));
      assert (status != 0);
      if (tmpd [0] <= 0){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd [0], file.c_str ());
        exit (EXIT_FAILURE);
      }
      unsigned int ncells = static_cast <unsigned int> (tmpd [0]);
      cells.reserve (ncells);

      unsigned int tmpu [4];
      for (unsigned int i = 0; i < ncells; ++i){
        status = fscanf (fp, "%d %d %d %d\n", &(tmpu [0]), &(tmpu [1]), &(tmpu [2]), &(tmpu [3]));
        assert (status != 0);
        assert (tmpu [0] >= 0 && tmpu [0] < numVerts);
        assert (tmpu [1] >= 0 && tmpu [1] < numVerts);
        assert (tmpu [2] >= 0 && tmpu [2] < numVerts);
        assert (tmpu [3] >= 0 && tmpu [3] < numVerts);
        cells.push_back (Cell (tmpu));
      }
      fclose (fp);

      assert (cells.size () == ncells);

      file = prefix;
      file.append (".tet.top");

      fp = fopen (file.c_str (), "r");
      assert (fp);

      status = fscanf (fp, "%d\n", &(tmpd [0]));
      assert (status != 0);
      if (tmpd [0] <= 0){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd [0], file.c_str ());
        exit (EXIT_FAILURE);
      }
      ncells = static_cast <unsigned int> (tmpd [0]);
      if (ncells != cells.size ()){
        PRINT ("fatal error: number of elements %u in %s don't match number of elements %lu in cell-indices files\n", ncells, file.c_str (), cells.size ());
        exit (EXIT_FAILURE);
      }

      for (unsigned int i = 0; i < ncells; ++i){
        status = fscanf (fp, "%d %d %d %d\n", &(tmpd [0]), &(tmpd [1]), &(tmpd [2]), &(tmpd [3]));
        assert (status != 0);
        assert (tmpd [0] < static_cast <int> (ncells) && tmpd [1] < static_cast <int> (ncells) &&
                tmpd [2] < static_cast <int> (ncells) && tmpd [3] < static_cast <int> (ncells));
        cells [i].addNeighbors (tmpd);
      }

      return true;
    }

    // static method to read face related files
    static bool
    readFaceFiles (const string &prefix, unsigned int numCells, unsigned int numFaces, unsigned int numVertices,
                   vector <Face> &ofaces, vector <unsigned int> &iindices, vector <Face> &ifaces)
    {
      int tmpd;

      string file (prefix);

      // read outside face related info file
      file.append (".trio.own");

      FILE *fp = fopen (file.c_str (), "r");
      assert (fp);

      int status = fscanf (fp, "%d\n", &tmpd);
      assert (status);
      if (tmpd <= 0 || 3*tmpd != static_cast <int> (numFaces)){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd, file.c_str ());
        exit (EXIT_FAILURE);
      }

      unsigned int nelems = static_cast <unsigned int> (tmpd);
      ofaces.resize (nelems);

      unsigned int tmpu [3];
      for (unsigned int i = 0; i < nelems; ++i){
        status = fscanf (fp, "%u %u\n", &(tmpu[0]), &(tmpu[1]));
        assert (status && tmpu [0] >= 0 && tmpu [0] < numCells && tmpu [1] < 4);
        ofaces [i]._owner = tmpu [0];
        ofaces [i]._index = static_cast <unsigned char> (tmpu [1]);
      }
      fclose (fp);

      // read inside faces
      file = prefix;
      file.append (".trii.ele");

      fp = fopen (file.c_str (), "r");
      assert (fp);

      status = fscanf (fp, "%d\n", &tmpd);
      assert (status);
      if (tmpd < 0){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd, file.c_str ());
        exit (EXIT_FAILURE);
      }

      if (tmpd){
        nelems = static_cast <unsigned int> (tmpd);
        iindices.resize (3*nelems);

        for (unsigned int i = 0; i < nelems; ++i){
          status = fscanf (fp, "%u %u %u\n", &(tmpu[0]), &(tmpu[1]), &(tmpu[2]));
          assert (status && tmpu [0] < numVertices && tmpu [1] < numVertices && tmpu [2] < numVertices);
          for (unsigned int j = 0; j < 3; ++j){
            iindices [3*i + j] = tmpu [j];
          }
        }
        fclose (fp);

        // read inside face related info file
        file = prefix;
        file.append (".trii.own");

        fp = fopen (file.c_str (), "r");
        assert (fp);

        status = fscanf (fp, "%d\n", &tmpd);
        assert (status);
        if (tmpd <= 0 || tmpd != static_cast <int> (nelems)){
          PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd, file.c_str ());
          exit (EXIT_FAILURE);
        }

        ifaces.resize (nelems);
        for (unsigned int i = 0; i < nelems; ++i){
          status = fscanf (fp, "%u %u\n", &(tmpu[0]), &(tmpu[1]));
          assert (status && tmpu [0] >= 0 && tmpu [0] < numCells && tmpu [1] < 4);
          ifaces [i]._owner = tmpu [0];
          ifaces [i]._index = static_cast <unsigned char> (tmpu [1]);
        }
        fclose (fp);
      } else {
        fclose (fp);
      }

      return true;
    }

    // static method to read edge related files
    static bool
    readEdgeFiles (const string &prefix, unsigned int numVerts, vector <unsigned int> &indices, vector <Edge> &edges)
    {
      int tmpd;

      // read edge element files;
      string file (prefix);
      file.append (".edge.ele");

      FILE *fp = fopen (file.c_str (), "r");
      assert (fp);

      int status = fscanf (fp, "%d\n", &tmpd);
      assert (status);
      if (tmpd <= 0){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd, file.c_str ());
        exit (EXIT_FAILURE);
      }

      unsigned int tmpu [2];
      unsigned int nedges = static_cast <unsigned int> (tmpd);
      indices.resize (2*nedges);
      for (unsigned int i = 0; i < 2*nedges; i += 2){
        status = fscanf (fp, "%u %u\n", &tmpu [0], &tmpu [1]);
        assert (status && tmpu [0] < numVerts && tmpu [1] < numVerts);
        indices [i] = tmpu [0];
        indices [i + 1] = tmpu [1];
      }
      fclose (fp);

      // read edge topology file

      file = prefix;
      file.append (".edge.top");

      fp = fopen (file.c_str (), "r");
      assert (fp);

      status = fscanf (fp, "%d\n", &tmpd);
      assert (status);
      if (tmpd <= 0 || tmpd != static_cast <int> (nedges)){
        PRINT ("fatal error: invalid number of elements %d in %s\n", tmpd, file.c_str ());
        exit (EXIT_FAILURE);
      }

      edges.resize (nedges);
      unsigned int *tmpo;
      for (unsigned int i = 0; i < nedges; ++i){
        status = fscanf (fp, "%u", &tmpu [0]);
        assert (status && tmpu [0]);

        tmpo = new unsigned int [tmpu [0]];
        for (unsigned int j = 0; j < tmpu [0]; ++j){
          status = fscanf (fp, "%u", &tmpo [j]);
        }
        edges [i] = Edge (indices [2*i], tmpu [0], tmpo);
        delete [] tmpo;
      }
      fclose (fp);

      return true;
    }

    // static method to update edge info structure of a cell
    static inline void
    updateCellEdgeInfo (unsigned int eIndex, unsigned int vert1, unsigned int vert2, Cell &c)
    {
      unsigned int ind [2];
      unsigned int pos = 6;
      for (unsigned int i = 0; i < 6; ++i){
        switch (i) {
          case 0: // edge 0-1
            ind [0] = c._index [0];
            ind [1] = c._index [1];
            break;
          case 1: // edge 0-2
            ind [0] = c._index [0];
            ind [1] = c._index [2];
            break;
          case 2: // edge 0-3
            ind [0] = c._index [0];
            ind [1] = c._index [3];
            break;
          case 3: // edge 1-2
            ind [0] = c._index [1];
            ind [1] = c._index [2];
            break;
          case 4: // edge 1-3
            ind [0] = c._index [1];
            ind [1] = c._index [3];
            break;
          case 5: // edge 2-3
            ind [0] = c._index [2];
            ind [1] = c._index [3];
        }
        if ((ind [0] == vert1 && ind [1] == vert2) || (ind [0] == vert2 && ind [1] == vert1)){
          pos = i;
          break;
        }
      }
      assert (pos < 6);
      c._edgeIndex [pos] = eIndex;
    }

    // static method to get partition index of a given cell
    static inline unsigned int
    getCellPartitionIndex (const vector<Partition> &partitions, const Cell &cell, const vector <vec> &vertices)
    {
      unsigned int result = partitions.size ();

      bool inFlag = false;
      for(unsigned int i = 0; i < partitions.size (); ++i){
        for (unsigned int j = 0; j < 4; ++j){
          if (partitions [i]._bbox.collide (vertices [cell._index [j]])){
            result = i;
            inFlag = true;
            break;
          }
        }
        if (inFlag){
          break;
        }
      }
      assert (result < partitions.size ());

      // get partition that collides with maximum number of vertices (race to 2)
      unsigned int counter = 0;
      for (unsigned int i = 0; i < partitions.size (); ++i){
        counter = 0;
        if (partitions [i]._bbox.collide (vertices [cell._index [0]])){
          ++counter;
        }
        if (partitions [i]._bbox.collide (vertices [cell._index [1]])){
          if (counter){
            return i;
          }
          ++counter;
        }
        if (partitions [i]._bbox.collide (vertices [cell._index [2]])){
          if (counter){
            return i;
          }
          ++counter;
        }
        if (partitions [i]._bbox.collide (vertices [cell._index [3]])){
          if (counter){
            return i;
          }
        }
      }

      return result;
    }

    // static method to get partition index of an edge (partition that contains max number of owning cells)
    static inline unsigned int
    getEdgePartitionIndex (const Edge &edge, const vector <unsigned int> &offsets, const vector <unsigned int> &counters)
    {
      unsigned int osize = offsets.size ();
      unsigned int *nCells = new unsigned int [osize];
      for (unsigned int i = 0; i < osize; ++i){
        nCells [i] = 0;
      }
      for (unsigned int i = 0; i < edge._numOwners; ++i){
        for (unsigned int j = 0; j < osize; ++j){
          if (edge._owner [i] < offsets [j] + counters [j]){
            ++nCells [j];
            break;
          }
        }
      }
      unsigned int maxInd = 0, maxNum = nCells [0];
      for (unsigned int i = 1; i < osize; ++i){
        if (maxNum < nCells [i]){
          maxInd = i;
          maxNum = nCells [i];
        }
      }
      delete [] nCells;
      return maxInd;
    }

		// private constructors and operators
		Submesh::Submesh () { }
		Submesh::Submesh (const Submesh &s) { }
		Submesh& Submesh::operator =(const Submesh &s) { return *this; }

		// proper constructor
		Submesh::Submesh (const string &config, const string &prefix, unsigned int index, unsigned int maxSurfaceVertexIndex,
                    vector <Vertex> &vi, const FaceChangeStruct &fc, vector <vec> **verts, vector <vec3> *texCoords, vector <unsigned int> &indices)
		: _myIndex (index), _maxSurfaceVertexIndex (maxSurfaceVertexIndex), _vertexInfo (&vi), _meshVertices (verts), _meshVertexTexCoords (texCoords),
		_meshFaceIndices (&indices), _meshSurfaceVertexTexCoords (vector <vec2> (maxSurfaceVertexIndex + 1))
		{
      // initialize OpenGL related attributes
		  initGLAttribs (config);

      // formulate file prefix
      string fpref (prefix);
      {
        char indexStr [4];
        sprintf (indexStr,"%u", index);
        fpref.append (".");
        fpref.append (indexStr);
      }

      // read cell indices and neighboring info
      readCellFiles (fpref, (*_meshVertices)->size (), _cells);

      // read face related info structures
      readFaceFiles (fpref, _cells.size (), _meshFaceIndices->size (), (**_meshVertices).size (), _faces, _insideFaceIndices, _insideFaces);

      // read edge related structure
      {
        vector <unsigned int> eIndices;
        readEdgeFiles (fpref, (**_meshVertices).size (), eIndices, _edges);
        for (unsigned int i = 0; i < _edges.size (); ++i){
          for (unsigned int j = 0; j < _edges [i]._numOwners; ++j){
            updateCellEdgeInfo (i, eIndices [2*i], eIndices [2*i + 1], _cells [_edges [i]._owner [j]]);
          }
        }
      }

      // set surface vertex flags
      for (unsigned int i = 0; i < _cells.size (); ++i){
        for (unsigned int j = 0; j < 4; ++j){
          if (_cells [i]._index [j] <= _maxSurfaceVertexIndex){
            _cells [i].setExternalVertexFlag (j);
          }
        }
      }
      {
        vector <bool> surfaceFlags;
        surfaceFlags.resize ((**_meshVertices).size (), false);
        for (unsigned int i = 0; i < _insideFaceIndices.size (); ++i){
          surfaceFlags [_insideFaceIndices [i]] = true;
        }
        for (unsigned int i = 0; i < _cells.size (); ++i){
          for (unsigned int j = 0; j < 4; ++j){
            if (surfaceFlags [_cells [i]._index [j]]){
              _cells [i].setInternalVertexFlag (j);
            }
          }
        }
      }

      // initialize partition
      {
        string pStr;
        if (!getConfigParameter (config, "num_submesh_partitions", pStr)){
          PRINT ("fatal error: partition size not specified in %s\n", config.c_str ());
          exit (EXIT_FAILURE);
        }
        assert (!pStr.empty ());
        for (unsigned int i = 0; i < pStr.size (); ++i){
          if (!isdigit (pStr [i])){
            PRINT ("fatal error: partition size %s specified in %s is not a number\n", pStr.c_str (), config.c_str ());
            exit (EXIT_FAILURE);
          }
        }
        unsigned int psize = static_cast <unsigned int> (atoi (pStr.c_str ()));
        _partitions.reserve (psize);
        _partitions.resize (psize);

        // intialize bounding volumes for submesh
        vector <vec> *vecs = *verts;
        vec3 min (vecs->at (_meshFaceIndices-> at (0))._v [0], vecs->at (_meshFaceIndices-> at (0))._v [1], vecs->at (_meshFaceIndices-> at (0))._v [2]);
        vec3 max (min);
        for (unsigned int i = 1; i < _meshFaceIndices->size (); ++i){
          for (unsigned int j = 0; j < 3; ++j){
            if (min._v [j] > vecs->at (_meshFaceIndices->at (i))._v [j]){
              min._v [j] = vecs->at (_meshFaceIndices->at (i))._v [j];
            } else if (max._v [j] < vecs->at (_meshFaceIndices->at (i))._v [j]){
              max._v [j] = vecs->at (_meshFaceIndices->at (i))._v [j];
            }
          }
        }
        for (unsigned int i = 0; i < _insideFaceIndices.size (); ++i){
          for (unsigned int j = 0; j < 3; ++j){
            if (min._v [j] > vecs->at (_insideFaceIndices [i])._v [j]){
              min._v [j] = vecs->at (_insideFaceIndices [i])._v [j];
            } else if (max._v [j] < vecs->at (_insideFaceIndices [i])._v [j]){
              max._v [j] = vecs->at (_insideFaceIndices [i])._v [j];
            }
          }
        }
        _bbox = aabb (min, max);

        // find axes with two smallest extents
        unsigned int maxAxis = 0, minAxis1 = 1, minAxis2 = 2; // default: x- greatest extent
        real diff [3] = {max._v [0] - min._v [0], max._v [1] - min._v [1], max._v [2] - min._v [2]};
        if (diff [1] > diff [0]){
          if (diff [2] > diff [1]){ // z - greatest extent
            maxAxis = 2;
            minAxis2 = 0;
          } else { // y - greatest extent
            maxAxis = 1;
            minAxis1 = 0;
          }
        } else if (diff [2] > diff [0]) { // z - greatest extent
          maxAxis = 2;
          minAxis2 = 0;
        }
        diff [maxAxis] /= static_cast <real> (psize);

        for (unsigned int i = 0; i < _partitions.size (); ++i){
          _partitions [i]._bbox._v [0]._v [maxAxis] = _bbox._v [0]._v [maxAxis] + diff [maxAxis]*i;
          _partitions [i]._bbox._v [0]._v [minAxis1] = _bbox._v [0]._v [minAxis1];
          _partitions [i]._bbox._v [0]._v [minAxis2] = _bbox._v [0]._v [minAxis2];

          _partitions [i]._bbox._v [1]._v [maxAxis] = _bbox._v [0]._v [maxAxis] + diff [maxAxis]*(i + 1);
          _partitions [i]._bbox._v [1]._v [minAxis1] = _bbox._v [1]._v [minAxis1];
          _partitions [i]._bbox._v [1]._v [minAxis2] = _bbox._v [1]._v [minAxis2];

          _partitions [i]._vertInfo = _vertexInfo;
          _partitions [i]._tex2D = &_meshSurfaceVertexTexCoords;
          _partitions [i]._tex3D = _meshVertexTexCoords;

          _partitions [i]._exMutex = &(_exMutex);
          _partitions [i]._exUpdateFlag = &(_exUpdateFlag);
          _partitions [i]._exVertices = &(_exVertices);
          _partitions [i]._exUVCoords = &(_exUVCoords);
          _partitions [i]._ex2DTexCoords = &(_ex2DTexCoords);
          _partitions [i]._exFaceIndices = &(_exFaceIndices);

          _partitions [i]._inMutex = &(_inMutex);
          _partitions [i]._inUpdateFlag = &(_inUpdateFlag);
          _partitions [i]._inVertices = &(_inVertices);
          _partitions [i]._inUVCoords = &(_inUVCoords);
          _partitions [i]._inSurfaceVertexStatus = &(_inSurfaceVertexStatus);
          _partitions [i]._in2DTexCoords = &(_in2DTexCoords);
          _partitions [i]._in3DTexCoords = &(_in3DTexCoords);
          _partitions [i]._inFaceIndices = &(_inFaceIndices);
        }
      }

      // reshuffle elements to align them with partitions
      reshuffleElements (index);

      // update bounds
      updateBounds ();
    }

		// destructor
		Submesh::~Submesh () { }

		// method to update bounds
		void
		Submesh::updateBounds ()
		{
      // update partition bounds
      for (unsigned int i = 0; i < _partitions.size (); ++i){
        for (unsigned int j = _partitions [i]._cellStartIndex; j <= _partitions [i]._cellEndIndex; ++j){
          for (unsigned int k = 0; k < 4; ++k){
            for (unsigned int l = 0; l < 3; ++l){
              if (_partitions [i]._bbox._v [0]. _v[l] > (**_meshVertices) [_cells [j]._index [k]]._v [l]){
                _partitions [i]._bbox._v [0]. _v[l] = (**_meshVertices) [_cells [j]._index [k]]._v [l];
              }
              else if (_partitions [i]._bbox._v [1]. _v[l] < (**_meshVertices) [_cells [j]._index [k]]._v [l]){
                _partitions [i]._bbox._v [1]. _v[l] = (**_meshVertices) [_cells [j]._index [k]]._v [l];
              }
            }
          }
        }
      }
      // use partition bounds to update submesh bounds
      for (unsigned int i = 0; i < _partitions.size (); ++i){
        for (unsigned int j = 0; j < 3; ++j){
          if (_bbox._v [0]._v [j] > _partitions [i]._bbox._v [0]._v [j]){
            _bbox._v [0]._v [j] = _partitions [i]._bbox._v [0]._v [j];
          } else if (_bbox._v [1]._v [j] < _partitions [i]._bbox._v [1]._v [j]){
            _bbox._v [1]._v [j] = _partitions [i]._bbox._v [1]._v [j];
          }
        }
      }
		}

    // method to resolve faces info structures and faceChange structure
    void
    Submesh::resolveFaces ()
    {
      for (unsigned int i = 0; i < _partitions.size (); ++i){
        if (!_partitions [i]._cutCells.empty () || !_partitions [i]._reExaminedCells.empty ()){
          for (unsigned int j = 3*_partitions [i]._exFaceStartIndex; j <= 3*_partitions [i]._exFaceEndIndex; j += 3){
            if (_faces [j/ 3]._owner < UINT_MAX && (*_meshFaceIndices) [j] == 0 && (*_meshFaceIndices) [j + 1] == 0){
              _faces [j/ 3]._owner = UINT_MAX;
              _changeBit->_cbit = true;
              if (_changeBit->_cfrom > j){
                _changeBit->_cfrom = j;
              }
              if (_changeBit->_cto < j){
                _changeBit->_cto = j;
              }
            }
          }
        }
      }
    }

		// method to gather a list of all cells affected by collisions
		void
		Submesh::getAffectedCells (unsigned int pIndex, vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      _partitions [pIndex].gatherAffectedCells (_myIndex, *_vertexInfo, **_meshVertices, *_meshFaceIndices, _faces,
                                                _insideFaceIndices, _insideFaces, _edges, _cells, bladeCurr, bladePrev, bladeIndices, bladeNormals);
    }

    // method to do the finalising of all cells
    void
    Submesh::finalizeCollision (unsigned int pIndex, vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2])
    {
      _partitions [pIndex].finalizeCollision((**_meshVertices), _edges, _cells, bladeCurr, bladePrev, bladeIndices, bladeNormals);
    }

    // shuffle cells and other info-structures so that cells with bordering vertices get pushed to the front and compartmentalized into partitions
    void
    Submesh::reshuffleElements (unsigned int myindex)
    {
      // no reshuffling required for one partition in submesh
      if (_partitions.size () < 2){
        _partitions [0]._cellStartIndex = 0;
        _partitions [0]._cellEndIndex = _cells.size () - 1;
        _partitions [0]._exFaceStartIndex = 0;
        _partitions [0]._exFaceEndIndex = _faces.size () - 1;
        if (!_insideFaceIndices.empty ()){
          _partitions [0]._inFaceStartIndex = 0;
          _partitions [0]._inFaceEndIndex = _insideFaces.size () - 1;
        }
        return;
      }

      // get partition index of cells
      vector <unsigned int> partitionIndex;
      partitionIndex.resize (_cells.size ());
      for (unsigned int i = 0; i < partitionIndex.size (); ++i){
        partitionIndex [i] = getCellPartitionIndex (_partitions, _cells [i], **_meshVertices);
      }

      /******************* RESHUFFLE CELLS RELATED DATA STRUCTURES AMONG SUBMESH PARTITIONS *******************/
      vector <unsigned int> cellOffsets, cellCounters, newIndices;

      // get cell-offsets for each partition
      {
        cellOffsets.resize (_partitions.size (), 0);
        cellCounters.resize (_partitions.size (), 0);

        for (unsigned int i = 0; i < partitionIndex.size (); ++i){
          ++cellOffsets [partitionIndex [i]];
        }
        for (int i = cellOffsets.size () - 1; i > 0; --i){
          cellOffsets [i] = cellOffsets [i - 1];
        }
        cellOffsets [0] = 0;
        for(unsigned int i = 2; i < cellOffsets.size (); ++i){
          cellOffsets [i] += cellOffsets [i - 1];
        }
        // copy offsets to partition data
        for (unsigned int i = 0; i < _partitions.size () - 1; ++i){
          _partitions [i]._cellStartIndex = cellOffsets [i];
          _partitions [i]._cellEndIndex = cellOffsets [i + 1] - 1;
        }
        _partitions [_partitions.size () - 1]._cellStartIndex = cellOffsets [cellOffsets.size () - 1];
        _partitions [_partitions.size () - 1]._cellEndIndex = _cells.size () - 1;

        // initialize new index array
        newIndices.resize (_cells.size (), UINT_MAX);

        // re-index cells with all surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (_cells [i].numExternalVertexBits () > 3){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 3 surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numExternalVertexBits () == 3){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 2 surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numExternalVertexBits () == 2){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 1 surface vertex
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numExternalVertexBits () == 1){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }

        // re-index cells with all internal surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numInternalVertexBits () > 3){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 3 surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numInternalVertexBits () == 3){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 2 surface vertices
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numInternalVertexBits () == 2){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
        // re-index cells with 1 surface vertex
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX && _cells [i].numInternalVertexBits () == 1){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }

        // re-index rest of the cells
        for (unsigned int i = 0; i < _cells.size (); ++i){
          if (newIndices [i] == UINT_MAX){
            newIndices [i] = cellOffsets [partitionIndex [i]] + cellCounters [partitionIndex [i]];
            ++cellCounters [partitionIndex [i]];
          }
        }
      }

      // update cell index information for face-info structures
      for (unsigned int i = 0; i < _faces.size (); ++i){
        _faces [i]._owner = newIndices [_faces [i]._owner];
      }
      for (unsigned int i = 0; i < _insideFaces.size (); ++i){
        _insideFaces [i]._owner = newIndices [_insideFaces [i]._owner];
      }
      for (unsigned int i = 0; i < _cells.size (); ++i){
        for (unsigned int j = 0; j < 4; ++j){
          if (_cells [i]._neighbor [j] >= 0){
            _cells [i]._neighbor [j] = newIndices [_cells [i]._neighbor [j]];
          }
          else if (_cells [i].numExternalVertexBits () >= 3){
            _cells [i].setExternalFaceFlag (j);
          }
          else {
            _cells [i].setInternalFaceFlag (j);
          }
        }
      }

      // update cell index information for vertex info structure
      for (unsigned int i = 0; i < _vertexInfo->size (); ++i){
        for (unsigned int j = 0; j < (* _vertexInfo) [i]._numSubmeshes; ++j){
          if ((* _vertexInfo) [i]._owners [j][0] == myindex){
            for (unsigned int k = 2; k < (* _vertexInfo) [i]._owners [j][1]; ++k){
              (* _vertexInfo) [i]._owners [j][k] = newIndices [(* _vertexInfo) [i]._owners [j][k]];
            }
            break;
          }
        }
      }

      // update edge-info structure to include new cell indices
      for (unsigned int i = 0; i < _edges.size (); ++i){
        for (unsigned int j = 0; j < _edges [i]._numOwners; ++j){
          _edges [i]._owner [j] = newIndices [_edges [i]._owner [j]];
        }
      }

      // shuffle cells
      {
        vector <Cell> tmpCells (_cells.size ());
        for (unsigned int i = 0; i < _cells.size (); ++i){
          memcpy (tmpCells [newIndices [i]]._bitset, _cells [i]._bitset, CELL_BIT_ARRAY_SIZE * sizeof (unsigned char));
          memcpy (tmpCells [newIndices [i]]._index, _cells [i]._index, 4 * sizeof (unsigned int));
          memcpy (tmpCells [newIndices [i]]._neighbor, _cells [i]._neighbor, 4 * sizeof (int));
          memcpy (tmpCells [newIndices [i]]._edgeIndex, _cells [i]._edgeIndex, 6 * sizeof (unsigned int));
        }
        _cells.swap (tmpCells);
      }
      partitionIndex.clear ();
      newIndices.clear ();

      /******************* RESHUFFLE EDGES RELATED DATA STRUCTURES AMONG SUBMESH PARTITIONS *******************/
      {
        vector <unsigned int> eOffsets, eCounters;
        eOffsets.resize (_partitions.size (), 0);
        eCounters.resize (_partitions.size (), 0);

        // get partition index of each edge
        partitionIndex.resize (_edges.size ());
        for (unsigned int i = 0; i < partitionIndex.size (); ++i){
          partitionIndex [i] = getEdgePartitionIndex (_edges [i], cellOffsets, cellCounters);
        }

        // calculate the offset for starting edge in each partition
        for (unsigned int i = 0; i < partitionIndex.size (); ++i){
          ++eOffsets [partitionIndex [i]];
        }
        for (int i = eOffsets.size () - 1; i > 0; --i){
          eOffsets [i] = eOffsets [i - 1];
        }
        eOffsets [0] = 0;
        for (unsigned int i = 2; i < eOffsets.size (); ++i){
          eOffsets [i] += eOffsets [i - 1];
        }

        // calculate new index positions for each edge
        newIndices.resize (_edges.size ());
        for (unsigned int i = 0; i < newIndices.size (); ++i){
          newIndices [i] = eOffsets [partitionIndex [i]] + eCounters [partitionIndex [i]];
          ++eCounters [partitionIndex [i]];
        }

        // update cell information structure
        for (unsigned int i = 0; i < _cells.size (); ++i){
          for (unsigned int j = 0; j < 6; ++j){
            _cells [i]._edgeIndex [j] = newIndices [_cells [i]._edgeIndex [j]];
          }
        }

        // shuffle edges
        vector <Edge> edges (_edges.size ());
        for (unsigned int i = 0; i < _edges.size (); ++i){
          edges [newIndices [i]] = _edges [i];
        }
        _edges.swap (edges);

        partitionIndex.clear ();
        newIndices.clear ();
      }

      /******************* RESHUFFLE EXTERNAL FACE RELATED DATA STRUCTURES AMONG SUBMESH PARTITIONS *******************/

      // get partition index of each external face
      partitionIndex.resize (_faces.size (), 0);
      for (unsigned int i = 0; i < _faces.size (); ++i){
        for (unsigned int j = 0; j < cellOffsets.size (); ++j){
          if (_faces [i]._owner < cellOffsets [j] + cellCounters [j]){
            partitionIndex [i] = j;
            break;
          }
        }
      }

      // set start external face index for each partition
      for (unsigned int i = 0; i < partitionIndex.size (); ++i){
        ++_partitions [partitionIndex [i]]._exFaceStartIndex;
      }
      for (int i = _partitions.size () - 1; i > 0; --i){
        _partitions [i]._exFaceStartIndex = _partitions [i - 1]._exFaceStartIndex;
      }
      _partitions [0]._exFaceStartIndex = 0;
      for (unsigned int i = 2; i < _partitions.size (); ++i){
        _partitions [i]._exFaceStartIndex += _partitions [i - 1]._exFaceStartIndex;
      }

      // calculate new index positions for each face (here we reuse exFaceEndIndex of each partition)
      newIndices.resize (_faces.size ());
      for (unsigned int i = 0; i < partitionIndex.size (); ++i){
        newIndices [i] = _partitions [partitionIndex [i]]._exFaceStartIndex + _partitions [partitionIndex [i]]._exFaceEndIndex;
        ++_partitions [partitionIndex [i]]._exFaceEndIndex;
      }

      // shuffle face indices and info structures
      {
        vector <unsigned int> tmpIndices ((*_meshFaceIndices).size ());
        for (unsigned int i = 0; i < (*_meshFaceIndices).size ()/ 3; ++i){
          for (unsigned int j = 0; j < 3; ++j){
            tmpIndices [3*newIndices [i] + j] = (*_meshFaceIndices) [3*i + j];
          }
        }
        _meshFaceIndices->swap (tmpIndices);

        vector <Face> tmpFaces (_faces.size ());
        for (unsigned int i = 0; i < _faces.size (); ++i){
          tmpFaces [newIndices [i]] = _faces [i];
        }
        _faces.swap (tmpFaces);
      }

      // set end external face index for each partition
      for (unsigned int i = 0; i < _partitions.size () - 1; ++i){
        _partitions [i]._exFaceEndIndex = _partitions [i + 1]._exFaceStartIndex - 1;
      }
      _partitions [_partitions.size () - 1]._exFaceEndIndex = _faces.size () - 1;

      partitionIndex.clear ();
      newIndices.clear ();

      /******************* RESHUFFLE INTERNAL FACE RELATED DATA STRUCTURES AMONG SUBMESH PARTITIONS *******************/
      if (_insideFaceIndices.empty ()){
        return;
      } else {
        for (unsigned int i = 0; i < _partitions.size (); ++i){
          _partitions [i]._inFaceStartIndex = 0;
        }
      }

      // get partition index of each external face
      partitionIndex.resize (_insideFaces.size (), 0);
      for (unsigned int i = 0; i < _insideFaces.size (); ++i){
        for (unsigned int j = 0; j < cellOffsets.size (); ++j){
          if (_insideFaces [i]._owner < cellOffsets [j] + cellCounters [j]){
            partitionIndex [i] = j;
            break;
          }
        }
      }

      // set start external face index for each partition
      for (unsigned int i = 0; i < partitionIndex.size (); ++i){
        ++_partitions [partitionIndex [i]]._inFaceStartIndex;
      }
      for (int i = _partitions.size () - 1; i > 0; --i){
        _partitions [i]._inFaceStartIndex = _partitions [i - 1]._inFaceStartIndex;
      }
      _partitions [0]._inFaceStartIndex = 0;
      for (unsigned int i = 2; i < _partitions.size (); ++i){
        _partitions [i]._inFaceStartIndex += _partitions [i - 1]._inFaceStartIndex;
      }

      // calculate new index positions for each face (here we reuse inFaceEndIndex of each partition)
      newIndices.resize (_insideFaces.size ());
      for (unsigned int i = 0; i < partitionIndex.size (); ++i){
        newIndices [i] = _partitions [partitionIndex [i]]._inFaceStartIndex + _partitions [partitionIndex [i]]._inFaceEndIndex;
        ++_partitions [partitionIndex [i]]._inFaceEndIndex;
      }

      // shuffle face indices and info structures
      {
        vector <unsigned int> tmpIndices (_insideFaceIndices.size ());
        for (unsigned int i = 0; i < _insideFaceIndices.size ()/ 3; ++i){
          for (unsigned int j = 0; j < 3; ++j){
            tmpIndices [3*newIndices [i] + j] = _insideFaceIndices [3*i + j];
          }
        }
        _insideFaceIndices.swap (tmpIndices);

        vector <Face> tmpFaces (_insideFaces.size ());
        for (unsigned int i = 0; i < _insideFaces.size (); ++i){
          tmpFaces [newIndices [i]] = _insideFaces [i];
        }
        _insideFaces.swap (tmpFaces);
      }

      // set end external face index for each partition
      for (unsigned int i = 0; i < _partitions.size () - 1; ++i){
        _partitions [i]._inFaceEndIndex = _partitions [i + 1]._inFaceStartIndex - 1;
      }
      _partitions [_partitions.size () - 1]._inFaceEndIndex = _insideFaces.size () - 1;
    }

		// private method to initialize OpenGL parameters
		void
		Submesh::initGLAttribs (const string &config)
		{
		  _inUpdateFlag = false;
      // initialize internal vertex buffer objects
      string iStr;
      getConfigParameter (config, "cut_internal_vbuffer_size", iStr);
      if (iStr.empty ()){
        PRINT ("error: Could not find cut_internal_vbuffer_size in %s. Aborting\n", config.c_str ());
        exit (EXIT_FAILURE);
      }
      for (unsigned int i = 0; i < iStr.size (); ++i){
        assert (isdigit (iStr [i]));
      }
      unsigned int vbsize = static_cast <unsigned int> (atoi (iStr.c_str ()));

      GLenum error;
      {
        vector <vec> verts;
        verts.reserve (vbsize);
        verts.resize (vbsize, vec::ZERO);
        assert (verts.size () == vbsize);

        glGenBuffers (1, &_glInVertexBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glInVertexBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*vbsize, &(verts [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }
      {
        vector <float> status;
        status.reserve (vbsize);
        status.resize (vbsize, 0.);
        assert (status.size () == vbsize);

        glGenBuffers (1, &_glInVertexStatusBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glInVertexStatusBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, sizeof (float)*vbsize, &(status [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }
      {
        vector <vec2> texcoords;
        texcoords.reserve (vbsize);
        texcoords.resize (vbsize, vec2::ZERO);
        assert (texcoords.size () == vbsize);

        glGenBuffers (1, &_glIn2DTexCoordBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glIn2DTexCoordBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*vbsize, &(texcoords [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }
      {
        vector <vec3> texcoords;
        texcoords.reserve (vbsize);
        texcoords.resize (vbsize, vec3::ZERO);
        assert (texcoords.size () == vbsize);

        glGenBuffers (1, &_glIn3DTexCoordBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glIn3DTexCoordBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, 3*sizeof (real)*vbsize, &(texcoords [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }

      getConfigParameter (config, "cut_internal_ibuffer_size", iStr);
      if (iStr.empty ()){
        PRINT ("error: Could not find cut_internal_ibuffer_size in %s. Aborting\n", config.c_str ());
        exit (EXIT_FAILURE);
      }
      for (unsigned int i = 0; i < iStr.size (); ++i){
        assert (isdigit (iStr [i]));
      }
      unsigned int ibsize = static_cast <unsigned int> (atoi (iStr.c_str ()));
      ibsize *= 3;

      {
        vector <unsigned int> indices;
        indices.reserve (ibsize);
        indices.resize (ibsize, 0);
        assert (indices.size () == ibsize);

        glGenBuffers (1, &_glInIndexBufferId);
        checkGLError (error);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glInIndexBufferId);
        checkGLError (error);
        glBufferData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int)*ibsize, &(indices [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
      }
/*
		  _exUpdateFlag = false;
		  // initialize external vertex buffer objects
      getConfigParameter (config, "cut_external_vbuffer_size", iStr);
      if (iStr.empty ()){
        PRINT ("error: Could not find cut_external_vbuffer_size in %s. Aborting\n", config.c_str ());
        exit (EXIT_FAILURE);
      }
      for (unsigned int i = 0; i < iStr.size (); ++i){
        assert (isdigit (iStr [i]));
      }
      vbsize = static_cast <unsigned int> (atoi (iStr.c_str ()));

      {
        vector <vec> verts;
        verts.reserve (vbsize);
        verts.resize (vbsize, vec::ZERO);
        assert (verts.size () == vbsize);

        glGenBuffers (1, &_glExVertexBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glExVertexBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*vbsize, &(verts [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }
      {
        vector <vec2> texcoords;
        texcoords.reserve (vbsize);
        texcoords.resize (vbsize, vec2::ZERO);
        assert (texcoords.size () == vbsize);

        glGenBuffers (1, &_glEx2DTexCoordBufferId);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glEx2DTexCoordBufferId);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*vbsize, &(texcoords [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
      }

      getConfigParameter (config, "cut_external_ibuffer_size", iStr);
      if (iStr.empty ()){
        PRINT ("error: Could not find cut_external_ibuffer_size in %s. Aborting\n", config.c_str ());
        exit (EXIT_FAILURE);
      }
      for (unsigned int i = 0; i < iStr.size (); ++i){
        assert (isdigit (iStr [i]));
      }
      ibsize = static_cast <unsigned int> (atoi (iStr.c_str ()));
      ibsize *= 3;

      {
        vector <unsigned int> indices;
        indices.reserve (ibsize);
        indices.resize (ibsize, 0);
        assert (indices.size () == ibsize);

        glGenBuffers (1, &_glExIndexBufferId);
        checkGLError (error);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glExIndexBufferId);
        checkGLError (error);
        glBufferData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int)*ibsize, &(indices [0]), GL_DYNAMIC_DRAW);
        checkGLError (error);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
      }
*/
		}
	}
}
