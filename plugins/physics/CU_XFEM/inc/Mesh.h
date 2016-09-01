/**
 * @file Mesh.h
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
 * The mesh class for the CU_XFEM library. This is derived class from
 * Resource.
 */
#pragma once

#include <vector>
#include <string>
#include <forward_list>
#include <boost/shared_ptr.hpp>

extern "C" {
#include <GL/glx.h>
#include <GL/glxext.h>
#include <cuda.h>
#include <cuda_gl_interop.h>
#include <cudaGL.h>
}

#include "GL/texture.h"
#include "Resource.h"

#include "Common.h"
#include "Vertex.h"

using namespace std;
using namespace boost;

namespace SF {

  class aabb;
  class Driver;
  class ThreadControl;

  namespace XFE {

    class Submesh;
    class Vertex;

    class Mesh : public Resource {

    public:
      aabb _bbox;

      /************************ THREADCONTROL RELATED PARAMETERS *************************/
      ThreadControl _syncControl;
      int _semPhysicsWaitIndex;
      int _semPhysicsPostIndex;
      int _semCollisionWaitIndex;
      int _semCollisionPostIndex;
      int _semIntersectionWaitIndex;
      int _semIntersectionPostIndex;
      int _semGraphicsWaitIndex;
      int _semGraphicsPostIndex;

      /************************ DATA RELATED PARAMETERS *************************/
      unsigned int _numSurfaceVertices;
      vector <vec> _vertices [2]; // on-host memory surface vertex buffers
      vector <Vertex> _vertexInfo;
      vector <vec3> _texCoords3D;
      vector <vec> *_curr;
      vector <vec> *_prev;

      forward_list <unsigned int> _collidingVertices;

      unsigned int _numCells;

      vector <unsigned int> _numFaces;
      vector <vector <unsigned int> > _faceIndices;

      vector <FaceChangeStruct> _faceChangeBits;
      vector <boost::shared_ptr <Submesh> > _submesh;

      // timing related parameters
      clock_t _present;
      clock_t _past;
      real _deltaT;
      real _deltaTminus1;

      /************************ OPENGL RELATED PARAMETERS *************************/
      bool _glBufferFlag; // flag to switch between two vertex buffers
      bool _glTextureFlag; // flag to denote presence of 3D texture data
      bool _glReprogramFlag; // flag to denote reloading of rendering programs

      GLuint _glNormalFramebufferDimensions [2]; // normal framebuffer dimensions
      GLuint _glNormalFramebufferId; // on-device framebuffer holding render texture
      GLuint _glNormalTexCoordBufferId; // on-device texture coordinate buffer
      GLuint _glNormalTextureId; // on-device 2D render texture for storing normals
      GLuint _glNormalVertexArrayId [2]; // vertex array to hold the buffers

      unsigned int _glNumFaces;
      GLuint _glNormalIndexBufferId; // element array to hold the outside indices for normal calculation

      GLuint _glEnvTextureId; // environment map ID (optional)
      GLuint _gl3DTextureId; // on-device 3D texture

      GLuint _glVertexBufferId [2]; // on-device memory vertex buffers

      // submesh-specific GL entities
      vector <GLuint> _glIndexBufferId; // indexed triangles for each octant
      vector <GLuint> _glTexCoordBufferId; // texture coordinates into texture atlas
      vector <GLuint> _gl2DTextureId; // texture atlas Id containing 3D texture coordinates
      vector <GLuint> _glRenderVertexArrayId; // vertex buffer used by rendering programs (2*number-of-submeshes)

      GLint _glEnvTextureLocation;

      // program variable locations (Rendering program 1)
      GLint _glModelviewMatrixLocation1;
      GLint _glProjectionMatrixLocation1;
      GLint _glNormalTextureLocation1;
      GLint _glColorTextureLocation1;
      GLint _glTexCoordTextureLocation1;
      GLint _glColorLocation1;

      GLint _glLightDirLocation1;
      GLint _glLightAmbLocation1;
      GLint _glLightDiffLocation1;
      GLint _glLightSpecLocation1;
      GLint _glLightExpLocation1;

      GLint _glLightDirLocation2;
      GLint _glLightAmbLocation2;
      GLint _glLightDiffLocation2;
      GLint _glLightSpecLocation2;
      GLint _glLightExpLocation2;

      // program variable locations (Rendering program 2)
      GLint _glModelviewMatrixLocation2;
      GLint _glProjectionMatrixLocation2;
      GLint _glColorTextureLocation2;
      GLint _glTexCoordTextureLocation2;
      GLint _glColorLocation2;

      // variables and references to members of Display
      unsigned int _glNumLights;

      real *_glLightDir1;
      real *_glLightAmb1;
      real *_glLightDiff1;
      real _glLightSpec1;
      real _glLightExp1;

      real *_glLightDir2;
      real *_glLightAmb2;
      real *_glLightDiff2;
      real _glLightSpec2;
      real _glLightExp2;

      real *_glModelview;
      real *_glProjection;
      string *_glslPrefixString;

      // color
      real _glColor [3];

      /*
       * Program 0: calculates vertex normals
       * Program 1: renders the objects
       * Program 2: calculates vertex forces using Hooke's Law
       */
      string _glProgramName [3]; // contains the names for the rendering programs
      GLuint _glProgram [3];

      /************************ GLX RELATED PARAMETERS *************************/
      GLXContext _glContext;
      Display *_glDisplay;
      GLXDrawable _glDrawable;
      GLXContext _cudaGLContext;

      /************************ CUDA RELATED PARAMETERS *************************
      string _cudaPTXFileName;
      string _cudaKernelFuncName;

      CUcontext _cudaContext;
      CUfunction _cudaKernelFunc;
      CUdeviceptr _devVertexBufferPtr [2];
      CUgraphicsResource _cudaVertexBufferId [2]; */

    public:
      Mesh (const string &configFile, Driver &driver);
      ~Mesh ();

      void run (); // run method
      void cleanup (); // cleanup method
      void adjustVertices (vector <vec> &curr, vector <vec> &prev, vector <unsigned int> &inds, vector <vec> &normals1, vector <vec> &normals2);
      bool initGPUPrograms (); // initializes all GPU programs

    private:
      Mesh ();
      Mesh (const Mesh &mesh);
      Mesh & operator = (const Mesh &mesh);

      bool initGLBufferObjects (); // initializes non-texture related GL Buffer objects
      bool initGLTextureObjects (unsigned int scale, const string & atlasShader, const Texture3D &texture);
      void rasterizeCharts (unsigned int atlasScale, const string &shader, const Texture3D &texture, const vector <real> &scales);
    };

  }
}
