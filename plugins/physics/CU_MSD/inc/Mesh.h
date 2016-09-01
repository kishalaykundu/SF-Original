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
 * The mesh class for the CU_MSD library. This is derived class from
 * Resource.
 */
#pragma once

#include <vector>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>

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

using namespace std;
using namespace boost;
using namespace boost::posix_time;

namespace SF {

  class aabb;
  class Driver;
  class ThreadControl;

	namespace MSD {

		class Mesh: public Resource {

    public:
      aabb _bbox;

      /************************ THREADCONTROL RELATED PARAMETERS *************************/
      ThreadControl _syncControl;
      int _semPhysicsWaitIndex;
      int _semPhysicsPostIndex;
      int _semCollisionWaitIndex;
      int _semCollisionPostIndex;
      int _semGraphicsWaitIndex;
      int _semGraphicsPostIndex;

      /************************ DATA RELATED PARAMETERS *************************/
      unsigned int _numSurfaceVertices;
      unsigned int _numTotalVertices;
      vector <vec> _vertices [2]; // on-host memory surface vertex buffers
      vector <vec> *_curr;
      vector <vec> *_prev;

      real *_mass;

			unsigned int _numSprings;

			vector <unsigned int> _numFaces;
			vector <vector <unsigned int> > _faceIndices;

			// time-related parameters
			ptime _past;
			ptime _present;
			time_duration _deltaT0;
			time_duration _deltaT1;

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

			GLuint _glSpringIndexBufferId; // on-gpu-memory buffer object for spring indices

			GLuint _glForceFrameBufferDimensions [2]; // framebuffer dimensions
			GLuint _glForceFrameBufferId; // on-gpu framebuffer holding the render texture
			GLuint _glForceTexCoordBufferId; // on-gpu-memory force texture co-ordinate buffers
			GLuint _glForceTextureId; // on-gpu-memory 2D normal render texture
			GLuint _glForceVertexArrayId [2]; // double buffered force vertex

      GLuint _glEnvTextureId; // environment map ID (optional)

      GLuint _glVertexBufferId [2]; // on-device memory vertex buffers
      GLuint _glRestVertexBufferId; // contains the original positions (used by the force calculation program)

      // partition-specific GL entities
      vector <GLuint> _glIndexBufferId; // indexed triangles for each octant
      vector <GLuint> _glTexCoordBufferId; // texture coordinates into texture atlas
      vector <GLuint> _glTextureId; // texture atlas ID containing color
      vector <GLuint> _glRenderVertexArrayId; // vertex buffer used by rendering programs (2*number-of-submeshes)

      // program variable locations (Rendering program 1)
      GLint _glModelviewMatrixLocation;
      GLint _glProjectionMatrixLocation;
      GLint _glNormalTextureLocation;
      GLint _glColorTextureLocation;
      GLint _glColorLocation;

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

      GLint _glEnvTextureLocation;

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
       * Program 2: calculates vertex forces
       */
      string _glProgramName [3]; // contains the names for the rendering programs
      GLuint _glProgram [3];

      /************************ GLX RELATED PARAMETERS *************************/
      GLXContext _glContext;
      Display *_glDisplay;
      GLXDrawable _glDrawable;
      GLXContext _cuglContext;

      /************************ CUDA (4.1) RELATED PARAMETERS *************************/
      string _cuPTXFileName;
      string _cuKernelFuncName [2];

      CUcontext _cuContext;
      CUfunction _cuKernelFunc [2];

      CUdeviceptr _devVertexBufferPtr [2];
      CUdeviceptr _devForceTexCoordBufferPtr;
      CUdeviceptr _devMassArrayPtr;

      CUarray _devForceTexturePtr;
      CUtexref _devTexRef;

      CUgraphicsResource _cuVertexBufferId [2];
      CUgraphicsResource _cuForceTexCoordBufferId;
      CUgraphicsResource _cuForceTextureId;

    public:
      Mesh (const string &configFile, Driver &driver);
      ~Mesh ();

      void run (); // run method
      void cleanup (); // clean up method for managing CUDA threads
      bool initGPUPrograms (); // initializes all GPU programs

    private:
      Mesh ();
      Mesh (const Mesh &mesh);
      Mesh & operator = (const Mesh &mesh);

      bool initGLForceBufferObjects (const vector <unsigned int> &springs);
      bool initGLBufferObjects (); // initializes non-texture related GL Buffer objects
      bool initGLTextureObjects (unsigned int scale, const string & atlasShader, const Texture3D &texture);
      void rasterizeCharts (unsigned int atlasScale, const string &shader, const Texture3D &texture, const vector <real> &scales, vector <vector <vec2> > &texCoords);
		};
	}
}
