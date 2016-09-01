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
 * The mesh class for the Rigid library.
 */
#pragma once

#include <vector>
#include <string>

#include "GL/common.h"
#include "Resource.h"

#include "Common.h"

namespace SF {

  class aabb;
  class Driver;
  class ThreadControl;

  namespace RM {

    class Mesh : public Resource {

    public:
      aabb _bbox;
      bool _transformFlag; // flag to switch between two vertex buffers

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
      size_t _numSurfaceVertices;
      vector <vec> _vertices [2]; // on-host memory surface vertex buffers
      vector <vec> *_curr;
      vector <vec> *_prev;

      vector <size_t> _numFaces;
      vector <vector <unsigned int> > _faceIndices;

      // vertices and indices for blade edge
      vector <vec> *_bladeCurr;
      vector <vec> *_bladePrev;
      vector <unsigned int> *_bladeIndices;

      /************************ OPENGL RELATED PARAMETERS *************************/
      bool _glBufferFlag; // flag to switch between two vertex buffers
      bool _glReprogramFlag; // flag to denote reloading of rendering programs

      GLuint _glNormalFramebufferDimensions [2]; // normal framebuffer dimensions
      GLuint _glNormalFramebufferId; // on-device framebuffer holding render texture
      GLuint _glNormalTexCoordBufferId; // on-device texture coordinate buffer
      GLuint _glNormalTextureId; // on-device 2D render texture for storing normals
      GLuint _glNormalVertexArrayId [2]; // vertex array to hold the buffers

      GLuint _glEnvTextureId; // environment map ID (optional)

      GLuint _glVertexBufferId [2]; // on-device memory vertex buffers
      GLuint _glIndexBufferId; // indexed triangles for each octant
      GLuint _glRenderVertexArrayId [2]; // vertex buffer used by rendering programs (2*number-of-submeshes)

      GLint _glEnvTextureLocation;

      // program variable locations (Rendering program 1)
      GLint _glModelviewMatrixLocation;
      GLint _glProjectionMatrixLocation;
      GLint _glNormalTextureLocation;
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
       */
      string _glProgramName [2]; // contains the names for the rendering programs
      GLuint _glProgram [2];

    public:
      Mesh (const string &configFile, Driver &driver);

      void run (); // run method
      void move (); // method to move
      bool initGPUPrograms (); // initializes all GPU programs

    private:
      Mesh ();
      Mesh (const Mesh &mesh);
      Mesh & operator = (const Mesh &mesh);

      bool initGLBufferObjects (); // initializes non-texture related GL Buffer objects
    };
  }
}
