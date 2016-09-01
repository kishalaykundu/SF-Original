/**
 * @file texture.h
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
 * Common texture related functions used by SF modules
 */

#pragma once

#include <vector>

#include "Preprocess.h"
#include "vec2.h"
#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif
#include "aabb.h"

using namespace std;

namespace SF {

  typedef struct texture3D_t {

    unsigned int _dimension [3];
    real _aspectRatio [3];
    vector <unsigned char> _rgba;

    inline texture3D_t ()
    {
      for (int i = 0; i < 3; ++i){
        _dimension [i] = 0;
      }
      for (int i = 0; i < 3; ++i){
        _aspectRatio [i] = 1.;
      }
    }
  } Texture3D;

  typedef struct fneighbor_t {
    int _v [3];
    inline fneighbor_t ()
    {
      for (int i = 0; i < 3; ++i) {
        _v [i] = -1;
      }
    }
  } FaceNeighbor;

  typedef struct fedge_t {
    unsigned int _v [2];
    inline fedge_t (const fedge_t &e)
    {
      memcpy (_v, e._v, 2*sizeof (unsigned int));
    }
    inline fedge_t (unsigned int v1, unsigned int v2)
    {
      _v [0] = v1;
      _v [1] = v2;
    }
    inline fedge_t & operator = (const fedge_t &e)
    {
      memcpy (_v, e._v, 2*sizeof (unsigned int));
      return *this;
    }
  } FaceEdge;

  // function to write different format data's to PNG file
  void writeRGBToPng (const char *prefix, int index, int dim, const vector <GLubyte> &rgb);
  void writeRGBAToPng (const char *prefix, int index, int dim, const vector <GLubyte> &rgba);
  void writeFloatToPng (const char *prefix, int index, int dim, const vector <GLfloat> &rgba);

  // function to initialize topology information
  void initTopologyInfo (const vector <unsigned int> &faceIndices, vector <FaceEdge> &edges, vector <FaceNeighbor> &neighbors);

  // function to get one-ring of faces around a face-mesh
  void getFaceRings (unsigned int index, const vector <vector <unsigned int> > &faceIndices, vector <unsigned int> &faces);

  // function to calculate parametric coordinates using Tutte's method (used to calculate texture coordinates)
  void calculateParametricCoordinates (unsigned int numSurfaceVerts, const vector <vec> &vertices, const vector <unsigned int> &indices, vector <vec2> &texCoords);

  // function to scale vertices with aspect ratio
  void scaleVertices (const real *aspect, const vector <vec> &src, const aabb &bv, vector <vec> &dest);

  // function to calculate area-weighted vertex normals of a mesh
  void calculateVertexNormals (const vector <vec> &verts, const vector <vector <unsigned int> > &faces, vector <vec> &normals);

  // function to calculate texture atlases
  void initTextureAtlas (GLuint program, int dim, const vector <vec> &verts, const vector <vec2> &texCoords,
                         const vector <unsigned int> &faces, vector <GLfloat> &rgbaData);

  // ray-trace function to reach texture boundary (returns end-position of rays)
  void raytraceThroughVolumef (int dim, const vector <GLfloat> &coData, const vector <GLfloat> &noData, const Texture3D &texture, vector <GLfloat> &rgbaData);

  // ray-trace function to reach texture boundary (returns trilinearly interpolated texel color at end of ray)
  void raytraceThroughVolumeb (int dim, const vector <GLfloat> &coData, const vector <GLfloat> &noData, const Texture3D &texture, vector <GLubyte> &rgbaData);
}
