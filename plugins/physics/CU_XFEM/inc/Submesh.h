/**
 * @file Submesh.h
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

#pragma once
#include <boost/thread/mutex.hpp>

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>

#else
#include <GL/gl.h>
#endif
}

#include <vector>
#include <string>
#include <forward_list>
#include <boost/shared_ptr.hpp>

#include "aabb.h"

#include "Common.h"
#include "Vertex.h"

using namespace std;
using namespace boost;

namespace SF {

	namespace XFE {

	  class Cell;
	  class Face;
	  class Edge;
	  class Partition;

		class Submesh {

		public:
			aabb _bbox;
			unsigned int _myIndex;
			FaceChangeStruct *_changeBit;
			vector <Partition> _partitions;

			/**************************** DATA RELATED PARAMETERS ****************************/
			unsigned int _maxSurfaceVertexIndex; // index of the last surface vertex in the vertex array
			vector <Vertex> *_vertexInfo; // pointer to the vertexInfo structure of mesh
			vector <vec> **_meshVertices; // pointer to pointer to current mesh vertices
			vector <vec3> *_meshVertexTexCoords; // pointer 3D texture coordinates for all vertices
			vector <unsigned int> *_meshFaceIndices; // pointer to current mesh face indices

			vector <Face> _faces;

			vector <vec2> _meshSurfaceVertexTexCoords; // 2D texture coordinates for surface vertices (used to look up texture atlas)

      /**
      * Structures to store inside faces that border other submeshes.
      * Face structure contains index of cell that owns it and the face enumeration in the cell
      */
      vector <unsigned int> _insideFaceIndices;
      vector <Face> _insideFaces;

      /**
      * Edge info structure. Contains:
      * 1. Number of cells containing current edge.
      * 2. List of indices of cells containing edge
      * 3. UV coordinate of any cut
      */
      vector <Edge> _edges;

      /**
       * Cell Info structure. Contains:
       * 1. indices to the vertices
       * 2. neighbors contains indices of cells that share faces with current cell
       */
      vector <Cell> _cells;

      /** Information regarding cut manifolds (two sets: 1 for internal trigs, 1 for external trigs)
       * 1. Texture co-ordinates (with respect to individual tetrahaderon) for cut vertices. Used
       *    to generate _vertices in each frame. These contain the information about cut.
       * 2. Vertices store the actual vertices calculated from above (updated everytime texture coords change)
       * 3. Face indices of the manifold cuts.
       * 4. 2D texture coordinates to store texture coordinates into texture-coordinate-texture
       * 5. 3D texture coordinates to store texture coordinates into 3D texture
       * 6. Status flags that denote if a vertex belongs to the original surface or not
       * 7. Index of cells that own the cut triangles
       */
       boost::mutex _exMutex;
       vector <vec> _exVertices;
       vector <vec3> _exUVCoords;
       vector <vec2> _ex2DTexCoords;

       vector <unsigned int> _exFaceIndices;

       boost::mutex _inMutex;
       vector <vec> _inVertices;
       vector <vec3> _inUVCoords;
       vector <float> _inSurfaceVertexStatus;
       vector <vec2> _in2DTexCoords;
       vector <vec3> _in3DTexCoords;

       vector <unsigned int> _inFaceIndices;

			/**************************** OPENGL RELATED PARAMETERS ****************************/
      bool _exUpdateFlag;
			GLuint _glExVertexBufferId;
			GLuint _glEx2DTexCoordBufferId;
			GLuint _glExIndexBufferId;
			GLuint _glExRenderVertexArrayId;

      bool _inUpdateFlag;
			GLuint _glInVertexBufferId;
			GLuint _glInVertexStatusBufferId;
			GLuint _glIn2DTexCoordBufferId;
			GLuint _glIn3DTexCoordBufferId;

			GLuint _glInIndexBufferId;
			GLuint _glInRenderVertexArrayId;

		public:
			Submesh (const string &config, const string &prefix, unsigned int i, unsigned int maxSurfaceVertexIndex,
            vector <Vertex> &vi, const FaceChangeStruct &fc, vector <vec> **verts, vector <vec3> *texCoords, vector <unsigned int> &indices);
			~Submesh ();

			inline void plainDraw ()
			{

			}
			inline void texturedDraw1 ()
			{
        glBindBuffer (GL_ARRAY_BUFFER, _glInVertexBufferId);
        glBufferSubData (GL_ARRAY_BUFFER, 0, SF_VECTOR_SIZE * sizeof (real) * _inVertices.size (), &(_inVertices [0]));

        glBindBuffer (GL_ARRAY_BUFFER, _glIn2DTexCoordBufferId);
        glBufferSubData (GL_ARRAY_BUFFER, 0, 2 * sizeof (real) * _in2DTexCoords.size (), &(_in2DTexCoords [0]));

        glBindBuffer (GL_ARRAY_BUFFER, _glIn3DTexCoordBufferId);
        glBufferSubData (GL_ARRAY_BUFFER, 0, 3 * sizeof (real) * _in3DTexCoords.size (), &(_in3DTexCoords [0]));

        if (_inUpdateFlag){
          glBindBuffer (GL_ARRAY_BUFFER, _glInVertexStatusBufferId);
          glBufferSubData (GL_ARRAY_BUFFER, 0, sizeof (real) * _inSurfaceVertexStatus.size (), &(_inSurfaceVertexStatus [0]));

          glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glInIndexBufferId);
          glBufferSubData (GL_ELEMENT_ARRAY_BUFFER, 0, sizeof (unsigned int) * _inFaceIndices.size (), &(_inFaceIndices [0]));

          _inUpdateFlag = false;
        }

        glBindVertexArray (_glInRenderVertexArrayId);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glInIndexBufferId);
        glDrawElements (GL_TRIANGLES, _inFaceIndices.size (), GL_UNSIGNED_INT, 0);
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
			}

      void updateBounds ();

      void resolveFaces ();
      void getAffectedCells (unsigned int pIndex, vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);
      void finalizeCollision (unsigned int pIndex, vector <vec> &bladeCurr, vector <vec> &bladePrev, vector <unsigned int> &bladeIndices, vector <vec> *bladeNormals [2]);

    private:
      void reshuffleElements (unsigned int index);
      void initGLAttribs (const string &config);

		private:
			Submesh ();
			Submesh (const Submesh &);
			Submesh& operator = (const Submesh &);
		};
	}
}
