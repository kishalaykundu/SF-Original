/**
 * @file Mesh.cpp
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

#include <cmath>

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glx.h>
#include <OpenGL/glxext.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glxext.h>
#endif
}


#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

#include "aabb.h"
#include "GL/common.h"
#include "GL/texture.h"

#include "CUDA/common.h"

#include "ThreadControl.h"
#include "Driver.h"
#include "Display.h"

#include "Common.h"
#include "Cell.h"
#include "Edge.h"
#include "Partition.h"
#include "Submesh.h"
#include "Mesh.h"

extern "C" {
#include <cuda_runtime_api.h>
}

namespace SF {

	namespace XFE {

	  static int GLX_ATTRIBUTE_LIST [] = {GLX_RGBA, None};

	  // static function to reload GPU programs
	  static void reloadPrograms (Resource & r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);
	    mptr->_glReprogramFlag = true;
	  }

    // inline function to draw normals
    static inline void drawNormals (Mesh *mptr)
    {
#ifndef NDEBUG
      GLenum error;
#endif

      // update data
      if (!mptr->_glBufferFlag){
        glBindBuffer (GL_ARRAY_BUFFER, mptr->_glVertexBufferId [0]);
#ifndef NDEBUG
        checkGLError (error);
#endif
        glBufferSubData (GL_ARRAY_BUFFER, 0, sizeof (real)*SF_VECTOR_SIZE*mptr->_numSurfaceVertices, &(mptr->_vertices [0] [0]));
#ifndef NDEBUG
        checkGLError (error);
#endif
      } else {
        glBindBuffer (GL_ARRAY_BUFFER, mptr->_glVertexBufferId [1]);
#ifndef NDEBUG
        checkGLError (error);
#endif
        glBufferSubData (GL_ARRAY_BUFFER, 0, sizeof (real)*SF_VECTOR_SIZE*mptr->_numSurfaceVertices, &(mptr->_vertices [1] [0]));
#ifndef NDEBUG
        checkGLError (error);
#endif
      }
      glBindBuffer (GL_ARRAY_BUFFER, 0);

      glUseProgram (mptr->_glProgram [0]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glClampColor (GL_CLAMP_VERTEX_COLOR, GL_FALSE);
      glClampColor (GL_CLAMP_READ_COLOR, GL_FALSE);
      glClampColor (GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glEnable (GL_BLEND);
      glBlendFunc (GL_ONE, GL_ONE);

      // bind normal frame buffer and render
      glBindFramebuffer (GL_FRAMEBUFFER, mptr->_glNormalFramebufferId);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glPushAttrib (GL_VIEWPORT_BIT);
      glViewport (0, 0, mptr->_glNormalFramebufferDimensions [0], mptr->_glNormalFramebufferDimensions [1]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glDrawBuffer (GL_COLOR_ATTACHMENT0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glClearColor (0, 0, 0, 0);
      glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      if (!mptr->_glBufferFlag){
        glBindVertexArray (mptr->_glNormalVertexArrayId [0]);
      } else {
        glBindVertexArray (mptr->_glNormalVertexArrayId [1]);
      }
#ifndef NDEBUG
      checkGLError (error);
#endif

        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glNormalIndexBufferId);
        glDrawElements (GL_TRIANGLES, mptr->_glNumFaces, GL_UNSIGNED_INT, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
      glFlush ();

      glBindVertexArray (0);
      glDrawBuffer (0);
      glBindFramebuffer (GL_FRAMEBUFFER, 0);
      glPopAttrib ();

      glDisable (GL_BLEND);

      glClampColor (GL_CLAMP_VERTEX_COLOR, GL_TRUE);
      glClampColor (GL_CLAMP_READ_COLOR, GL_TRUE);
      glClampColor (GL_CLAMP_FRAGMENT_COLOR, GL_TRUE);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUseProgram (0);
    }

    // drawing function for non-textured datasets
	  static void plainDraw (Resource &r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);

	    mptr->_syncControl [mptr->_semGraphicsWaitIndex].wait ();

	    // reload program if needed
	    if (mptr->_glReprogramFlag) {
        mptr->initGPUPrograms ();
        mptr->_glReprogramFlag = false;
	    }

      drawNormals (mptr);

#ifndef NDEBUG
      GLenum error;
#endif

      for (unsigned int i = 0; i < mptr->_glIndexBufferId.size (); ++i){
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        if (mptr->_faceChangeBits [i]._cbit){
          glBufferSubData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int) * mptr->_faceChangeBits [i]._cfrom,
                           sizeof (unsigned int) * (3 + mptr->_faceChangeBits [i]._cto - mptr->_faceChangeBits [i]._cfrom),
                           &(mptr->_faceIndices [i] [mptr->_faceChangeBits [i]._cfrom]));
#ifndef NDEBUG
      checkGLError (error);
#endif
          mptr->_faceChangeBits [i].reset ();
        }
      }

      /************************ RENDERING OF EXTERNAL SURFACE HAPPENS HERE ************************/

      glUseProgram (mptr->_glProgram [1]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation1, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation1, 1, false, mptr->_glProjection);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniform3f (mptr->_glColorLocation1, mptr->_glColor [0], mptr->_glColor [1], mptr->_glColor [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      if (mptr->_glNumLights){
        glUniform3f (mptr->_glLightDirLocation1, mptr->_glLightDir1 [0], mptr->_glLightDir1 [1], mptr->_glLightDir1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightAmbLocation1, mptr->_glLightAmb1 [0], mptr->_glLightAmb1 [1], mptr->_glLightAmb1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightDiffLocation1, mptr->_glLightDiff1 [0], mptr->_glLightDiff1 [1], mptr->_glLightDiff1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightSpecLocation1, mptr->_glLightSpec1);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightExpLocation1, mptr->_glLightExp1);
#ifndef NDEBUG
      checkGLError (error);
#endif
      }
      if (mptr->_glNumLights > 1){
        glUniform3f (mptr->_glLightDirLocation2, mptr->_glLightDir2 [0], mptr->_glLightDir2 [1], mptr->_glLightDir2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightAmbLocation2, mptr->_glLightAmb2 [0], mptr->_glLightAmb2 [1], mptr->_glLightAmb2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightDiffLocation2, mptr->_glLightDiff2 [0], mptr->_glLightDiff2 [1], mptr->_glLightDiff2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightSpecLocation2, mptr->_glLightSpec2);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightExpLocation2, mptr->_glLightExp2);
#ifndef NDEBUG
      checkGLError (error);
#endif
      }

      glActiveTexture (GL_TEXTURE0);
      glBindTexture (GL_TEXTURE_2D, mptr->_glNormalTextureId );
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniform1i (mptr->_glNormalTextureLocation1, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      if (mptr->_glEnvTextureId){
        glActiveTexture (GL_TEXTURE1);
        glBindTexture (GL_TEXTURE_CUBE_MAP, mptr->_glEnvTextureId);
#ifndef NDEBUG
        checkGLError (error);
#endif
        glUniform1i (mptr->_glEnvTextureLocation, 1);
#ifndef NDEBUG
        checkGLError (error);
#endif
      }

      unsigned int offset = 0;
      if (mptr->_glBufferFlag){
        offset = mptr->_glIndexBufferId.size ();
      }
      for (unsigned int i = 0; i < mptr->_glIndexBufferId.size (); ++i){
        glBindVertexArray (mptr->_glRenderVertexArrayId [i + offset]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glDrawElements (GL_TRIANGLES, mptr->_numFaces [i], GL_UNSIGNED_INT, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }
      glBindTexture (GL_TEXTURE_2D, 0);

      glBindTexture (GL_TEXTURE_CUBE_MAP, 0);
      glUseProgram (0);

      /************************ RENDERING OF CUT SURFACE HAPPENS HERE *************************/
      glUseProgram (mptr->_glProgram [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation2, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation2, 1, false, mptr->_glProjection);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniform3f (mptr->_glColorLocation2, mptr->_glColor [0], mptr->_glColor [1], mptr->_glColor [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      for (unsigned int i = 0; i < mptr->_submesh.size (); ++i){
        mptr->_submesh [i].get ()->plainDraw ();
      }
      glUseProgram (0);

      mptr->_syncControl [mptr->_semGraphicsPostIndex].post ();
	  }

    // drawing function textured datasets
	  static void texturedDraw (Resource &r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);

	    mptr->_syncControl [mptr->_semGraphicsWaitIndex].wait ();
/******************************************************************************************************
	    glXMakeCurrent (mptr->_glDisplay, mptr->_glDrawable, mptr->_glContext);
/******************************************************************************************************/

	    // reload program if needed
	    if (mptr->_glReprogramFlag) {
        mptr->initGPUPrograms ();
        mptr->_glReprogramFlag = false;
	    }

      drawNormals (mptr);

#ifndef NDEBUG
      GLenum error;
#endif

      for (unsigned int i = 0; i < mptr->_glIndexBufferId.size (); ++i){
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        if (mptr->_faceChangeBits [i]._cbit){
          glBufferSubData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int) * mptr->_faceChangeBits [i]._cfrom,
                           sizeof (unsigned int) * (3 + mptr->_faceChangeBits [i]._cto - mptr->_faceChangeBits [i]._cfrom),
                           &(mptr->_faceIndices [i] [mptr->_faceChangeBits [i]._cfrom]));
#ifndef NDEBUG
      checkGLError (error);
#endif
          mptr->_faceChangeBits [i].reset ();
        }
      }

      /************************ RENDERING OF EXTERNAL SURFACE HAPPENS HERE ************************/
      glUseProgram (mptr->_glProgram [1]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation1, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation1, 1, false, mptr->_glProjection);
#ifndef NDEBUG
      checkGLError (error);
#endif

      if (mptr->_glNumLights){
        glUniform3f (mptr->_glLightDirLocation1, mptr->_glLightDir1 [0], mptr->_glLightDir1 [1], mptr->_glLightDir1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightAmbLocation1, mptr->_glLightAmb1 [0], mptr->_glLightAmb1 [1], mptr->_glLightAmb1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightDiffLocation1, mptr->_glLightDiff1 [0], mptr->_glLightDiff1 [1], mptr->_glLightDiff1 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightSpecLocation1, mptr->_glLightSpec1);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightExpLocation1, mptr->_glLightExp1);
#ifndef NDEBUG
      checkGLError (error);
#endif
      }
      if (mptr->_glNumLights > 1){
        glUniform3f (mptr->_glLightDirLocation2, mptr->_glLightDir2 [0], mptr->_glLightDir2 [1], mptr->_glLightDir2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightAmbLocation2, mptr->_glLightAmb2 [0], mptr->_glLightAmb2 [1], mptr->_glLightAmb2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform3f (mptr->_glLightDiffLocation2, mptr->_glLightDiff2 [0], mptr->_glLightDiff2 [1], mptr->_glLightDiff2 [2]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightSpecLocation2, mptr->_glLightSpec2);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1f (mptr->_glLightExpLocation2, mptr->_glLightExp2);
#ifndef NDEBUG
      checkGLError (error);
#endif
      }

      glActiveTexture (GL_TEXTURE0);
      glBindTexture (GL_TEXTURE_3D, mptr->_gl3DTextureId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniform1i (mptr->_glColorTextureLocation1, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glActiveTexture (GL_TEXTURE1);
      glBindTexture (GL_TEXTURE_2D, mptr->_glNormalTextureId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniform1i (mptr->_glNormalTextureLocation1, 1);
#ifndef NDEBUG
      checkGLError (error);
#endif

      if (mptr->_glEnvTextureId){
        glActiveTexture (GL_TEXTURE2);
        glBindTexture (GL_TEXTURE_CUBE_MAP, mptr->_glEnvTextureId);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1i (mptr->_glEnvTextureLocation, 2);
#ifndef NDEBUG
      checkGLError (error);
#endif
      }

      unsigned int offset = 0;
      if (mptr->_glBufferFlag){
        offset = mptr->_glIndexBufferId.size ();
      }

      for (unsigned int i = 0; i < mptr->_glIndexBufferId.size (); ++i){
        glActiveTexture (GL_TEXTURE3);
        glBindTexture (GL_TEXTURE_2D, mptr->_gl2DTextureId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1i (mptr->_glTexCoordTextureLocation1, 3);
#ifndef NDEBUG
      checkGLError (error);
#endif

        glBindVertexArray (mptr->_glRenderVertexArrayId [i + offset]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glDrawElements (GL_TRIANGLES, mptr->_numFaces [i], GL_UNSIGNED_INT, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      glBindTexture (GL_TEXTURE_2D, 0);
      glBindTexture (GL_TEXTURE_CUBE_MAP, 0);
      glBindTexture (GL_TEXTURE_3D, 0);

      glUseProgram (0);

      /************************ RENDERING OF CUT SURFACE HAPPENS HERE ************************/
      glUseProgram (mptr->_glProgram [2]);

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation2, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation2, 1, false, mptr->_glProjection);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glActiveTexture (GL_TEXTURE0);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glBindTexture (GL_TEXTURE_3D, mptr->_gl3DTextureId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniform1i (mptr->_glColorTextureLocation2, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      Submesh *sm;
      for (unsigned int i = 0; i < mptr->_submesh.size (); ++i){
        sm = mptr->_submesh [i].get ();
        if (!sm->_inVertices.empty ()){
          glActiveTexture (GL_TEXTURE1);
#ifndef NDEBUG
      checkGLError (error);
#endif
          glBindTexture (GL_TEXTURE_2D, mptr->_gl2DTextureId [i]);
#ifndef NDEBUG
          checkGLError (error);
#endif
          glUniform1i (mptr->_glTexCoordTextureLocation2, 1);
#ifndef NDEBUG
          checkGLError (error);
#endif

          sm->texturedDraw1 ();
        }
      }

      glBindTexture (GL_TEXTURE_3D, 0);
      glBindTexture (GL_TEXTURE_2D, 0);

      glUseProgram (0);

	    mptr->_syncControl [mptr->_semGraphicsPostIndex].post ();
	  }

		// protected constructor and assignment functions
		Mesh::Mesh () { }
		Mesh::Mesh (const Mesh &m) { }
		Mesh& Mesh::operator = (const Mesh &m) { return *this; }

		// only legitimate constructor
		Mesh::Mesh (const string &config, Driver &driver)
		: _semPhysicsWaitIndex (-1), _semPhysicsPostIndex (-1),
		  _semIntersectionWaitIndex (-1), _semIntersectionPostIndex (-1),
		  _semGraphicsWaitIndex (-1), _semGraphicsPostIndex (-1),
		  _numSurfaceVertices (0), _curr (&(_vertices [0])), _prev (&(_vertices [1])), _numCells (0),
		  _present (clock ()), _past (clock ()), _deltaT (0.), _deltaTminus1 (0.),
		  _glBufferFlag (false), _glTextureFlag (false), _glReprogramFlag (false),
		  _glNormalFramebufferId (0), _glNormalTexCoordBufferId (0), _glNormalTextureId (0), _glNumFaces (0), _glNormalIndexBufferId (0),
		  _glEnvTextureId (driver._display.get ()->_glEnvTextureId), _gl3DTextureId (0),
		  _glEnvTextureLocation (-1),
		  _glModelviewMatrixLocation1 (-1), _glProjectionMatrixLocation1 (-1), _glNormalTextureLocation1 (-1),
		  _glColorTextureLocation1 (-1), _glTexCoordTextureLocation1 (-1), _glColorLocation1 (-1),
		  _glLightDirLocation1 (-1), _glLightAmbLocation1 (-1), _glLightDiffLocation1 (-1), _glLightSpecLocation1 (-1), _glLightExpLocation1 (-1),
		  _glLightDirLocation2 (-1), _glLightAmbLocation2 (-1), _glLightDiffLocation2 (-1), _glLightSpecLocation2 (-1), _glLightExpLocation2 (-1),
		  _glModelviewMatrixLocation2 (-1), _glProjectionMatrixLocation2 (-1),
		  _glColorTextureLocation2 (-1), _glTexCoordTextureLocation2 (-1), _glColorLocation2 (-1),
		  _glNumLights (driver._display.get ()->_numLights),
		  _glLightDir1 (NULL), _glLightAmb1 (NULL), _glLightDiff1 (NULL), _glLightSpec1 (0.), _glLightExp1 (0.),
		  _glLightDir2 (NULL), _glLightAmb2 (NULL), _glLightDiff2 (NULL), _glLightSpec2 (0.), _glLightExp2 (0.),
		  _glModelview (&(driver._display.get ()->_modelview [0])), _glProjection (&(driver._display.get ()->_projection [0])),
		  _glslPrefixString (&(driver._display.get ()->_glslPrefixString))
		{
			assert (!config.empty ());

			_owner = boost::shared_ptr <string> (new string ("CudaXfem"));

#ifndef SF_NO_PRINT
			unsigned int numTotalFaces = 0;
#endif
			{
				string name;
				if (!getConfigParameter (config, "name", name)){
          PRINT ("fatal error: name not specified in %s\n", config.c_str ());
          exit (EXIT_FAILURE);
				}
				_name = boost::shared_ptr <string> (new string (name));

				string folder;
				if (!getConfigParameter (config, "data_folder", folder)){
          PRINT ("fatal error: data-folder not specified in %s\n", config.c_str ());
          exit (EXIT_FAILURE);
				}
				if (folder.at (folder.size () - 1) != '/'){
					folder.append ("/");
				}

				unsigned int depth = 0;
				{
					string depthStr;
					if (!getConfigParameter (config, "max_depth", depthStr)){
            PRINT ("fatal error: max-depth not specified in %s\n", config.c_str ());
            exit (EXIT_FAILURE);
					}
					for (unsigned int i = 0; i < depthStr.size (); ++i){
						if (!isdigit (depthStr.at (i))){
              PRINT ("fatal error: max-depth %s specified in %s is not a number\n", depthStr.c_str (), config.c_str ());
              exit (EXIT_FAILURE);
						}
					}
					depth = static_cast <unsigned int> (atoi (depthStr.c_str ()));

					folder.append (depthStr);
					folder.append ("/");
				}

				unsigned int numSubmeshes = 1;
				for (unsigned int i = 0; i < depth; ++i){
					numSubmeshes *= 8;
				}

				// generate file prefix
				string prefix (folder);
				prefix.append (name);

				/*************************** READ NODE (VERTEX) AND OWNER INFO FILE ***************************/
				int tmpd  [3]; // size 3 because this is used to read triangular elements later

				string file (prefix);
				file.append (".node");

				FILE* fp = fopen (file.c_str (), "r");
				assert (fp);

				int status = fscanf (fp, "%d\n", &(tmpd [0]));
				assert (status != 0);
				if (tmpd [0] <= 0){
          PRINT ("fatal error: invalid number of vertices \'%d\' in %s\n", tmpd [0], file.c_str ());
          exit (EXIT_FAILURE);
				}

				unsigned int nverts = static_cast <unsigned int> (tmpd [0]);

				_vertices [0].reserve (nverts);
				_vertices [1].reserve (nverts);
				{
					real tmpr  [SF_VECTOR_SIZE];
#ifdef SF_VECTOR4_ENABLED
					tmpr [3] = 1.;
#endif
#ifdef SF_DOUBLE_PRECISION
					status = fscanf (fp, "%lf %lf %lf\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#else
					status = fscanf (fp, "%f %f %f\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#endif
					assert (status != 0);
					_vertices [0].push_back (vec (tmpr));

					vec3 min (vec3 (tmpr [0], tmpr [1], tmpr [2]));
					vec3 max (min);

					for (unsigned int i = 1; i < nverts; ++i){
#ifdef SF_DOUBLE_PRECISION
						status = fscanf (fp, "%lf %lf %lf\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#else
						status = fscanf (fp, "%f %f %f\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#endif
						assert (status != 0);
						_vertices [0].push_back (vec (tmpr));
						for (int j = 0; j < 3; ++j){
							if (min._v [j] > tmpr [j]){
								min._v [j] = tmpr [j];
							}
							else if (max._v [j] < tmpr [j]){
								max._v [j] = tmpr [j];
							}
						}
					}
					for (int j = 0; j < 3; ++j){
            min._v [j] -= .05;
					}
					for (int j = 0; j < 3; ++j){
            max._v [j] += .05;
					}
					_bbox = aabb (min, max);
				}
				fclose (fp);

				_vertices [1] = _vertices [0];

				file.append (".own");
				_vertexInfo.resize (_vertices [0].size ());

				fp = fopen (file.c_str (), "r");
				assert (fp);
				{
				  status = fscanf (fp, "%d\n", &(tmpd [0]));
				  assert (status);
          if (tmpd [0] <= 0){
            PRINT ("fatal error: invalid number of vertices \'%d\' in %s\n", tmpd [0], file.c_str ());
            exit (EXIT_FAILURE);
          }
          assert (nverts == static_cast <unsigned int> (tmpd [0]));

          unsigned int nElems, *elems;
          unsigned int *subNums = new unsigned int [numSubmeshes];
          for (unsigned int i = 0; i < nverts; ++i){

            for (unsigned int j = 0; j < numSubmeshes; ++j){
              subNums [j] = 0;
            }
            status = fscanf (fp, "%u", &nElems);
            assert (status);
            nElems *= 2;
            elems = new unsigned int [nElems];
            for (unsigned int j = 0; j < nElems; j += 2){
              status = fscanf (fp, "%u %u", &(elems [j]), &(elems [j + 1]));
              assert (status);
            }
            for (unsigned int j = 0; j < nElems; j += 2){
              ++subNums [elems [j]];
            }
            for (unsigned int j = 0; j < numSubmeshes; ++j){
              if (subNums [j]){
                _vertexInfo [i].allocateSubmeshSpace (j, subNums [j]);
              }
            }
            for (unsigned int j = 0; j < nElems; j += 2){
              _vertexInfo [i].addOwner (elems [j], elems [j + 1]);
            }
            delete [] elems;
          }
          delete [] subNums;
          fclose (fp);
				}

				/*************************** READ TRIANGULAR ELEMENT FILES AND POPULATE SUBMESHES ***************************/
				_numFaces.resize (numSubmeshes, 0);
				_faceIndices.resize (numSubmeshes, vector <unsigned int> ());

				for (unsigned int i = 0; i < numSubmeshes; ++i){

					// generate file name
					file = prefix;
					file.append (".");

					char indexStr [4];
					sprintf (indexStr, "%u", i);
					file.append (indexStr);
					file.append (".trio.ele");

					fp = fopen (file.c_str (), "r");
					assert (fp);

					status = fscanf (fp, "%d\n", &(tmpd [0]));
					assert (status != 0);
          if (tmpd [0] < 0){
            PRINT ("fatal error: invalid number of elements \'%d\' in %s\n", tmpd [0], file.c_str ());
            exit (EXIT_FAILURE);
          }

					_numFaces [i] = static_cast <unsigned int> (tmpd [0]);

					if (_numFaces [i]){
						_faceIndices [i].reserve (3*_numFaces [i]);
#ifndef NDEBUG
            int vertSize = static_cast <int> (_vertices [0].size ());
#endif

						for (unsigned int j = 0; j < _numFaces [i]; ++j){
							status = fscanf (fp, "%d %d %d\n", &(tmpd [0]), &(tmpd [1]), &(tmpd [2]));
							assert (status != 0);
							assert (tmpd [0] >= 0 && tmpd [0] < vertSize);
							assert (tmpd [1] >= 0 && tmpd [1] < vertSize);
							assert (tmpd [2] >= 0 && tmpd [2] < vertSize);
							for (int k = 0; k < 3; ++k){
								_faceIndices [i].push_back (static_cast <unsigned int> (tmpd [k]));
								if (_numSurfaceVertices < static_cast <unsigned int> (tmpd [k])){
									_numSurfaceVertices = static_cast <unsigned int> (tmpd [k]);
								}
							}
						}

	#ifndef SF_NO_PRINT
						numTotalFaces += _numFaces.at (i);
	#endif
						_numFaces [i] *= 3;
					}
					fclose (fp);
				} // end - for (unsigned int i = 0; i < numSubmeshes; ++i)

				++_numSurfaceVertices; // this is done because before this _numSurfaceVertices contains biggest index of triangles

				// initialize submeshes and needed substructures
				_faceChangeBits.resize (numSubmeshes, FaceChangeStruct ());

				_submesh.reserve (numSubmeshes);
				for (unsigned int i = 0; i < numSubmeshes; ++i){
					_submesh.push_back (boost::shared_ptr <Submesh> (new Submesh (config, prefix, i, _numSurfaceVertices - 1, _vertexInfo,
                                                                   _faceChangeBits [i], &_curr, &_texCoords3D, _faceIndices [i])));
				}
				assert (_submesh.size () == numSubmeshes);
			}

      /*************************** INITIALIZE THREAD CONTROL PARAMETERS ***************************/
      {
        string mStr;
        getConfigParameter (config, "num_mutexes", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }

        int numMutex = atoi (mStr.c_str ());
        string msp ("mutex_startval"), msv;
        for (int i = 0; i < numMutex; ++i){
          mStr = msp;
          char indexStr [4];
          sprintf (indexStr, "%d", i + 1);
          mStr.append (indexStr);

          getConfigParameter (config, mStr.c_str (), msv);
          assert (!msv.empty ());
          for (unsigned int j = 0; j < msv.size (); ++j){
            assert (isdigit (msv [j]));
          }
          _syncControl.push_back (static_cast <unsigned int> (atoi (msv.c_str ())));
        }

        getConfigParameter (config, "physics_wait_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semPhysicsWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "physics_post_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semPhysicsPostIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "collision_wait_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semCollisionWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "collision_post_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semCollisionPostIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "intersection_wait_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semIntersectionWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "intersection_post_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semIntersectionPostIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "graphics_wait_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semGraphicsWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "graphics_post_index", mStr);
        assert (!mStr.empty ());
        for (unsigned int i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semGraphicsPostIndex = atoi (mStr.c_str ());
      }

      /*************************** INITIALIZE OPENGL PARAMETERS ***************************/

      // initialize octant specific GL buffer ID's
      _glIndexBufferId.reserve (_submesh.size ());
      _glTexCoordBufferId.reserve (_submesh.size ());
      _gl2DTextureId.reserve (_submesh.size ());
      _glRenderVertexArrayId.reserve (2*_submesh.size ());

      GLuint tmpu = 0;
      _glIndexBufferId.resize (_submesh.size (), tmpu);
      _glTexCoordBufferId.resize (_submesh.size (), tmpu);
      _gl2DTextureId.resize (_submesh.size (), tmpu);
      _glRenderVertexArrayId.resize (2*_submesh.size (), tmpu);

      // get names of the shading programs
      getConfigParameter (config, "normal_shader", _glProgramName [0]);
      assert (!_glProgramName [0].empty ());

      getConfigParameter (config, "color_shader1", _glProgramName [1]);
      assert (!_glProgramName [1].empty ());

      getConfigParameter (config, "color_shader2", _glProgramName [2]);
      assert (!_glProgramName [2].empty ());

      for (int i = 0; i < 3; ++i){
        _glProgram [i] = 0;
      }

      string texStr;
      getConfigParameter (config, "texture", texStr);

      if (!texStr.empty ()) {

        // initialize texture-related variables
        _glTextureFlag = true;

        Texture3D tex3d;
        int status = 0;

        // load and read 3D texture file
        {
          string texInfoFile;
          getConfigParameter (config, "textureinfo", texInfoFile);
          assert (!texInfoFile.empty ());

          FILE *txfp = fopen (texInfoFile.c_str (), "r");
          assert (txfp);
          status = fscanf (txfp, "%u %u %u", &(tex3d._dimension [0]), &(tex3d._dimension [1]), &(tex3d._dimension [2]));
          assert (status != 0);
          assert (tex3d._dimension [0] > 0 && tex3d._dimension [1] > 0 && tex3d._dimension [2] > 0);

#ifdef SF_DOUBLE_PRECISION
          status = fscanf (txfp, "%lf %lf %lf\n", &(tex3d._aspectRatio [0]), &(tex3d._aspectRatio [1]), &(tex3d._aspectRatio [2]));
#else
          status = fscanf (txfp, "%f %f %f\n", &(tex3d._aspectRatio [0]), &(tex3d._aspectRatio [1]), &(tex3d._aspectRatio [2]));
#endif
          assert (status != 0);
          assert (tex3d._aspectRatio [0] > 0. && tex3d._aspectRatio [1] > 0. && tex3d._aspectRatio [2] > 0.);
          fclose (txfp);

          unsigned int size = 4*tex3d._dimension [0]*tex3d._dimension [1]*tex3d._dimension [2];
          tex3d._rgba.resize (size);

          txfp = fopen (texStr.c_str (), "rb");
          assert (txfp);
          status = fread (&(tex3d._rgba [0]), sizeof (unsigned char), size, txfp);
          assert (status != 0);
          texStr.clear ();
        }

        // initialize 3D vertex texture coordinates (is properly populated in initTextureAtlas (..) function)
        _texCoords3D.resize (_vertices [0].size (), vec3 (2., 2., 2.));

        // initialize non-texture related buffers
        initGLBufferObjects ();

        // initialize texture-related objects
        string atlasShader;
        getConfigParameter (config, "atlas_shader", atlasShader);
        assert (!atlasShader.empty ());
        unsigned int scale = 0;
        {
          string scaleStr;
          getConfigParameter (config, "atlas_scale", scaleStr);
          assert (!scaleStr.empty ());
          for (unsigned int i = 0; i < scaleStr.length (); ++i){
            if (!isdigit (scaleStr [i])){
              PRINT ("fatal error: atlas scale %s in %s not a number", scaleStr.c_str (), config.c_str ());
              exit (EXIT_FAILURE);
            }
          }
          scale = static_cast <unsigned int> (atoi (scaleStr.c_str ()));
        }
        initGLTextureObjects (scale, atlasShader, tex3d);
      }
      else {
        string cStr;
        getConfigParameter (config, "color", cStr);

        if (!cStr.empty ()){

          unsigned int first = cStr.find_first_of (" ", 0);
          unsigned int last = cStr.find_last_of (" ", cStr.size () - 1);

          string red, green, blue;
          red.append (cStr, 0, first);
          green.append (cStr, first + 1, last - first - 1);
          blue.append (cStr, last + 1, cStr.size () - last);

          assert (!red.empty () && !green.empty ()  && !blue.empty ());
          for (unsigned int i = 0; i < red.size (); ++i){
            assert (isdigit (red [i]) || red [i] == '.');
          }
          for (unsigned int i = 0; i < green.size (); ++i){
            assert (isdigit (green [i]) || green [i] == '.');
          }
          for (unsigned int i = 0; i < blue.size (); ++i){
            assert (isdigit (blue [i]) || blue [i] == '.');
          }

          _glColor [0] = static_cast <real> (atof (red.c_str ()));
          _glColor [1] = static_cast <real> (atof (green.c_str ()));
          _glColor [2] = static_cast <real> (atof (blue.c_str ()));
        }
        else {
          for (int i = 0; i < 3; ++i){
            _glColor [i] = .5;
          }
        }

        initGLBufferObjects ();
      }

      // initialize GPU programs and update view volume
      initGPUPrograms ();

      GL_Window *disp = driver._display.get ();
      if (_glNumLights){
        _glLightDir1 = &(disp->_lightDir1 [0]);
        _glLightAmb1 = &(disp->_lightAmb1 [0]);
        _glLightDiff1 = &(disp->_lightDiff1 [0]);
        _glLightSpec1 = disp->_lightSpec1;
        _glLightExp1 = disp->_lightExp1;
      }
      if (_glNumLights > 1){
        _glLightDir2 = &(disp->_lightDir2 [0]);
        _glLightAmb2 = &(disp->_lightAmb2 [0]);
        _glLightDiff2 = &(disp->_lightDiff2 [0]);
        _glLightSpec2 = disp->_lightSpec2;
        _glLightExp2 = disp->_lightExp2;
      }
      for (int i = 0; i < 3; ++i){
        if (disp->_bbox._v [0]._v [i] > _bbox._v [0]._v [i]){
          disp->_bbox._v [0]._v [i] = _bbox._v [0]._v [i];
        }
      }
      for (int i = 0; i < 3; ++i){
        if (disp->_bbox._v [1]._v [i] < _bbox._v [1]._v [i]){
          disp->_bbox._v [1]._v [i] = _bbox._v [1]._v [i];
        }
      }

      /*************************** INITIALIZE RESOURCE FUNCTION POINTERS ***************************/
      if (_glTextureFlag){
        draw = &texturedDraw;
      } else {
        draw = &plainDraw;
      }
      reprogram = &reloadPrograms;

      /*************************** INITIALIZE GLX RELATED PARAMETERS ***************************
      _glContext = glXGetCurrentContext ();
      _glDisplay = glXGetCurrentDisplay ();
      _glDrawable = glXGetCurrentDrawable ();

      /*************************** INITIALIZE CUDA RELATED PARAMETERS ***************************
      getConfigParameter (config, "compute_program", _cudaPTXFileName);
      assert (!_cudaPTXFileName.empty ());

      getConfigParameter (config, "compute_function", _cudaKernelFuncName);
      assert (!_cudaKernelFuncName.empty ());
      */

		} // end - Mesh::Mesh (const string &config, Driver &driver)

    // destructor method
    Mesh::~Mesh ()
    {
      /*
      // clean up
      cuGraphicsUnregisterResource (_cudaVertexBufferId [0]);
      cuGraphicsUnregisterResource (_cudaVertexBufferId [1]);

      glXMakeCurrent (_glDisplay, _glDrawable, _cudaGLContext);
      glXDestroyContext (_glDisplay, _cudaGLContext);

      glXDestroyContext (_glDisplay, _glContext);
      cuCtxDestroy (_cudaContext);
      */
    }

    // mesh's run method
    void
    Mesh::run ()
    {
      {
        unsigned int ctr = 0, nedges = 0;
        Submesh *sm;
        for (unsigned int i = 0; i < _submesh.size (); ++i){
          sm = _submesh [i].get ();
          nedges += sm->_edges.size ();
          for (unsigned int j = 0; j < sm->_edges.size (); ++j){
            ctr += sm->_edges [j]._numOwners;
          }
        }
        PRINT ("Total edges: %u Avg owners: %g\n", nedges, static_cast <float> (ctr)/ static_cast <float> (nedges));

        ctr = 0;
        for (unsigned int i = 0; i < _vertexInfo.size (); ++i){
          for (unsigned int j = 0; j < _vertexInfo [i]._numSubmeshes; ++j){
            ctr += (_vertexInfo [i]._owners [j] [1] - 1);
          }
        }
        PRINT ("Vertex incidence: %g\n", static_cast <float> (ctr)/ static_cast <float> (_vertexInfo.size ()));
      }
      /*
      size_t numBytes;
#ifndef NDEBUG
      cudaError_t error;
      GLenum gerror;
#endif

      // make a new context from existing resources
      XVisualInfo *visualInfo = glXChooseVisual (_glDisplay, 0, GLX_ATTRIBUTE_LIST);
#ifndef NDEBUG
      checkGLError (gerror);
#endif

      // get the context
      _cudaGLContext = glXCreateContext (_glDisplay, visualInfo, _glContext, GL_TRUE);
#ifndef NDEBUG
      checkGLError (gerror);
#endif
      glXMakeCurrent (_glDisplay, _glDrawable, _cudaGLContext);
#ifndef NDEBUG
      checkGLError (gerror);
#endif

      // initialize CUDA
      initCUDA ();

      // variables for CUDA kernel
      int threadsPerBlock = 32;
      unsigned int problemsize = _vertices [0].size ();
      void *kernelArgs1 [] = { &_devVertexBufferPtr [0], &_devVertexBufferPtr [1], &problemsize};
      void *kernelArgs2 [] = { &_devVertexBufferPtr [1], &_devVertexBufferPtr [0], &problemsize};
      */

      // go to computation loop
      while (true) {

        // wait for physics end to get control
        _syncControl [_semPhysicsWaitIndex].wait ();

        /*
      // push cuda context into stack
      cuCtxPushCurrent (_cudaContext);
#ifndef NDEBUG
      checkCUDAError (error);
#endif

      // map GL buffers
      cuGraphicsMapResources (2, _cudaVertexBufferId, 0);
#ifndef NDEBUG
      checkCUDAError (error);
#endif

      // retrieve device pointers from mapped buffers
      cuGraphicsResourceGetMappedPointer (&_devVertexBufferPtr [0], &numBytes, _cudaVertexBufferId [0]);
#ifndef NDEBUG
      checkCUDAError (error);
#endif

      cuGraphicsResourceGetMappedPointer (&_devVertexBufferPtr [1], &numBytes, _cudaVertexBufferId [1]);
      #ifndef NDEBUG
      checkCUDAError (error);
      #endif

      // launch CUDA kernel
      if (_glBufferFlag){
        //cuLaunchKernel (_cudaKernelFunc, problemsize/ threadsPerBlock, 1, 1, threadsPerBlock, 1, 1, 0, 0, kernelArgs1, 0);
      } else {
        //cuLaunchKernel (_cudaKernelFunc, problemsize/ threadsPerBlock, 1, 1, threadsPerBlock, 1, 1, 0, 0, kernelArgs2, 0);
      }

      // unmap GL buffers so that OpenGL can use them
      cuGraphicsUnmapResources (2, _cudaVertexBufferId, 0);
      #ifndef NDEBUG
      checkCUDAError (error);
      #endif

      // pop context
      cuCtxPopCurrent (&_cudaContext);
        */

        // toggle swap buffer flag
        _glBufferFlag = !_glBufferFlag; // should be the last line in this segment

        // release data
        _syncControl [_semPhysicsPostIndex].post ();
      }
    }

    // method to re-adjust vertices that collide with blades
    void
    Mesh::adjustVertices (vector <vec> &curr, vector <vec> &prev, vector <unsigned int> &inds, vector <vec> &normals1, vector <vec> &normals2)
    {
      // gather the colliding vertices from every partition
      Submesh *sm;
      for (unsigned int i = 0; i < _submesh.size (); ++i){
        sm = _submesh [i].get ();
        for (unsigned int j = 0; j < sm->_partitions.size (); ++j){
          if (!sm->_partitions [j]._collidingVertices.empty ()){
            while (!sm->_partitions [j]._collidingVertices.empty ()){
              _collidingVertices.push_front (sm->_partitions [j]._collidingVertices.front ());
              sm->_partitions [j]._collidingVertices.pop_front ();
            }
          }
        }
      }
      if (_collidingVertices.empty ()){
        return;
      }

      _collidingVertices.sort ();
      _collidingVertices.unique ();

      vec ed;
      unsigned int ind, cInd;
      bool finishedFlag = false, notOkFlag = false, surfaceFlag = false, conditionFlag = false;

      while (!_collidingVertices.empty ()){

        finishedFlag = false;
        conditionFlag = false;
        ind = _collidingVertices.front ();

        // get the first cell in the list and determine if it belongs to object surface
        sm = _submesh [_vertexInfo [ind]._owners [0][0]].get ();
        for (unsigned i = 0; i < 4; ++i){
          if (sm->_cells [_vertexInfo [ind]._owners [0][2]]._index [i] == ind){
            surfaceFlag = sm->_cells [_vertexInfo [ind]._owners [0][2]].testExternalVertexFlag (i);
            break;
          }
        }
        if (!surfaceFlag){
          conditionFlag = true;
        }

        for (unsigned int i = 0; i < _vertexInfo [ind]._numSubmeshes; ++i){

          sm = _submesh [_vertexInfo [ind]._owners [i][0]].get ();
          for (unsigned int j = 2; j < _vertexInfo [ind]._owners [i][1]; ++j){

            cInd = _vertexInfo [ind]._owners [i][j];

            for (unsigned int k = 0; k < 4; ++k){

              conditionFlag |= sm->_cells [cInd].testExternalVertexFlag (k);

              if (conditionFlag && sm->_cells [cInd]._index [k] != ind){

                notOkFlag = false;
                ed = _curr->at (sm->_cells [cInd]._index [k]) - _curr->at (ind);

                for (unsigned int l = 0; l < normals1.size (); ++l){
                  if (ABS(ed.dot (normals1 [l])) > 1. - EPSILON || ABS(ed.dot (normals2 [l])) > 1. - EPSILON){
                    notOkFlag = true;
                    break;
                  }
                }
                if (!notOkFlag){
                  ed *= .2;
                  _curr->at (ind) += ed;
                  finishedFlag = true;
                  break;
                }
              } // end - if (conditionFlag && sm->_cells [_vertexInfo [ind]._owners [i][j]]._index [k] != ind)
            } // end - for (unsigned int k = 0; k < 4; ++k)
            if (finishedFlag){
              break;
            }
          } // end - for (unsigned int j = 2; j < _vertexInfo [ind]._owners [i][1]; ++j)
          if (finishedFlag){
            break;
          }
        } // end - for (unsigned int i = 0; i < _vertexInfo [ind]._numSubmeshes; ++i)

        _collidingVertices.pop_front ();
      }
    }

/*
    // function to initialize CUDA-related parameters
    void
    Mesh::initCUDA ()
    {
      // initialize CUDA driver API
      cuInit (0);

      // check if any device is CUDA capable
      int numCudaDevices = 0;
      cuDeviceGetCount (&numCudaDevices);
      if (!numCudaDevices){
        PRINT ("Error: Could not find CUDA devices\n");
      }
      cudaError_t error;

      // get handle for device 0
      CUdevice cudaDevice;
      cuDeviceGet (&cudaDevice, 0);
      checkCUDAError (error);

      // create context
      cuGLCtxCreate (&_cudaContext, CU_CTX_BLOCKING_SYNC, cudaDevice);
      checkCUDAError (error);

      //cuCtxPushCurrent (_cudaContext);

      // create and load module
      CUmodule cudaModule;
      string cudaPTXFileContents;
      readCudaPTXFile (_cudaPTXFileName, cudaPTXFileContents);
      cuModuleLoadDataEx (&cudaModule, cudaPTXFileContents.c_str (), 0, 0, 0);
      checkCUDAError (error);

      // load module function
      cuModuleGetFunction (&_cudaKernelFunc, cudaModule, _cudaKernelFuncName.c_str ());
      checkCUDAError (error);

      // register GL buffers to be used
      cuGraphicsGLRegisterBuffer (&_cudaVertexBufferId [0], _glVertexBufferId [0], CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE);
      checkCUDAError (error);

      cuGraphicsGLRegisterBuffer (&_cudaVertexBufferId [1], _glVertexBufferId [1], CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE);
      checkCUDAError (error);

      // release current context so that it can be used by while loop
      //cuCtxPopCurrent (&_cudaContext);
    }
*/
		// private method to initialize non-texture related OpenGL buffer objects
		bool
		Mesh::initGLBufferObjects ()
		{
		  GLenum error;

		  glGenBuffers (2, _glVertexBufferId);
      checkGLError (error);

		  glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [0]);
      checkGLError (error);
		  glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*_vertices [0].size (), &(_vertices [0][0]), GL_DYNAMIC_DRAW);
      checkGLError (error);

		  glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [1]);
      checkGLError (error);
		  glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*_vertices [1].size (), &(_vertices [1][0]), GL_DYNAMIC_DRAW);
      checkGLError (error);
		  glBindBuffer (GL_ARRAY_BUFFER, 0);

      glGenBuffers (_glIndexBufferId.size (), &(_glIndexBufferId [0]));
      checkGLError (error);
		  for (unsigned int i = 0; i < _glIndexBufferId.size (); ++i){
		    glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glIndexBufferId [i]);
		    checkGLError (error);
		    glBufferData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int)*_numFaces [i], &(_faceIndices [i][0]), GL_DYNAMIC_DRAW);
		    checkGLError (error);
		  }
		  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);

		  /*************************** INITIALIZE NORMAL CALCULATION BUFFERS ***************************/
      unsigned int width = static_cast <unsigned int> (trunc (ceil (sqrt (static_cast <double> (_numSurfaceVertices)))));
      unsigned int height = static_cast <unsigned int> (trunc (floor (sqrt (static_cast <double> (_numSurfaceVertices)))));

      unsigned int pow2 = 1;
      while (width > pow2){
        pow2 *= 2;
      }
      width = pow2;
      pow2 = 1;
      while (height > pow2){
        pow2 *= 2;
      }
      height = pow2;
      while (width*height > _numSurfaceVertices){
        height /= 2;
      }
      if (width*height < _numSurfaceVertices){
        height *= 2;
      }
      while (width*height > _numSurfaceVertices){
        width /= 2;
      }
      if (width*height < _numSurfaceVertices){
        width *= 2;
      }

      _glNormalFramebufferDimensions [0] = static_cast <GLuint> (width);
      _glNormalFramebufferDimensions [1] = static_cast <GLuint> (height);

      glGenTextures (1, &_glNormalTextureId);
      checkGLError (error);
      glBindTexture (GL_TEXTURE_2D, _glNormalTextureId);
      checkGLError (error);

      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA32F, _glNormalFramebufferDimensions [0], _glNormalFramebufferDimensions [1], 0, GL_RGBA, GL_FLOAT, 0);
      checkGLError (error);

      glGenFramebuffers (1, &_glNormalFramebufferId);
      checkGLError (error);
      glBindFramebuffer (GL_FRAMEBUFFER, _glNormalFramebufferId);
      checkGLError (error);
      glFramebufferTexture2D (GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, _glNormalTextureId, 0);
      checkGLError (error);
      glBindFramebuffer (GL_FRAMEBUFFER, 0);

      // initialize texture coordinates for every surface vertex into the normal texture
      vector <vec2> texcoords;
      texcoords.resize (_numSurfaceVertices, vec2::ZERO);

      unsigned int xcount = 0, ycount = 0;
      for (unsigned int i = 0; i < _numSurfaceVertices; ++i){
        texcoords [i] = vec2 (static_cast <real> (xcount)/ static_cast <real> (width), static_cast <real> (ycount)/ static_cast <real> (height));
        ++xcount;
        if (xcount >= width){
          xcount = 0;
          ++ycount;
        }
      }

      // advance each coordinate by half pixel to align with GL texture
      real xstep = 1./ static_cast <real> (2*_glNormalFramebufferDimensions [0]);
      real ystep = 1./ static_cast <real> (2*_glNormalFramebufferDimensions [1]);
      for (unsigned int i = 0; i < _numSurfaceVertices; ++i){
        texcoords [i]._v [0] += xstep;
        texcoords [i]._v [1] += ystep;
      }

      // generate the GL buffers
      glGenBuffers (1, &_glNormalTexCoordBufferId);
      checkGLError (error);
      glBindBuffer (GL_ARRAY_BUFFER, _glNormalTexCoordBufferId);
      checkGLError (error);
      glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*texcoords.size (), &(texcoords [0]), GL_STATIC_DRAW);
      checkGLError (error);
      glBindBuffer (GL_ARRAY_BUFFER, 0);

      // add faces to the GL normal index buffer
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        _glNumFaces += _faceIndices [i].size ();
      }

      vector <unsigned int> indices;
      indices.reserve (_glNumFaces);
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        for (unsigned int j = 0; j < _faceIndices [i].size (); ++j){
          indices.push_back (_faceIndices [i] [j]);
        }
      }
      glGenBuffers (1, &_glNormalIndexBufferId);
      checkGLError (error);
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glNormalIndexBufferId);
      checkGLError (error);
      glBufferData (GL_ELEMENT_ARRAY_BUFFER, _glNumFaces * sizeof (unsigned int), &(indices [0]), GL_STATIC_DRAW);
      checkGLError (error);
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);

      return true;
		}

    // method to initialize all GPU programs
    bool
    Mesh::initGPUPrograms ()
    {
      /*************************** INITIALIZE NORMAL RENDERING PROGRAM ***************************/
      if (!initGPUProgram (true, *_glslPrefixString, _glProgramName [0], _glProgram [0])){
        PRINT ("error: could not initialize %s\n", _glProgramName [0].c_str ());
        return false;
      }

      GLenum error;

      // initialize vertex array for normal rendering program
      glUseProgram (_glProgram [0]);
      checkGLError (error);

      GLint vertLocation = glGetAttribLocation (_glProgram [0], "vertex");
      assert (vertLocation > -1);
      GLint texCoordLocation = glGetAttribLocation (_glProgram [0], "inTexCoord");
      assert (texCoordLocation > -1);

      glBindFragDataLocation (_glProgram [0], 0, "fragColor");
      checkGLError (error);

      glGenVertexArrays (2, _glNormalVertexArrayId);
      checkGLError (error);

      for (int i = 0; i < 2; ++i){

        glBindVertexArray (_glNormalVertexArrayId [i]);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [i]);
        checkGLError (error);
        glVertexAttribPointer (vertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (vertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glNormalTexCoordBufferId);
        checkGLError (error);
        glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (texCoordLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      glUseProgram (0);

      /*************************** INITIALIZE EXTERNAL SURFACE RENDERING PROGRAM ***************************/
      if (!initGPUProgram (false, *_glslPrefixString, _glProgramName [1], _glProgram [1])){
        PRINT ("error: could not initialize %s\n", _glProgramName [1].c_str ());
        return false;
      }

      // initialize vertex array for extermal surface rendering program
      glUseProgram (_glProgram [1]);
      checkGLError (error);

      vertLocation = glGetAttribLocation (_glProgram [1], "vertex");
      assert (vertLocation > -1);
      texCoordLocation = glGetAttribLocation (_glProgram [1], "normalTexCoord");
      assert (texCoordLocation > -1);
      _glModelviewMatrixLocation1 = glGetUniformLocation (_glProgram [1], "modelview");
      assert (_glModelviewMatrixLocation1 > -1);
      _glProjectionMatrixLocation1 = glGetUniformLocation (_glProgram [1], "projection");
      assert (_glProjectionMatrixLocation1 > -1);
      _glNormalTextureLocation1 = glGetUniformLocation (_glProgram [1], "normalTexture");
      assert (_glNormalTextureLocation1 > -1);

      if (_glEnvTextureId){
        _glEnvTextureLocation = glGetUniformLocation (_glProgram [1], "envTexture");
        assert (_glEnvTextureLocation > -1);
      }

      if (_glNumLights){
        _glLightDirLocation1 = glGetUniformLocation (_glProgram [1], "lightDir1");
        assert (_glLightDirLocation1 > -1);
        _glLightAmbLocation1 = glGetUniformLocation (_glProgram [1], "lightAmbient1");
        assert (_glLightAmbLocation1 > -1);
        _glLightDiffLocation1 = glGetUniformLocation (_glProgram [1], "lightDiffuse1");
        assert (_glLightDiffLocation1 > -1);
        _glLightSpecLocation1 = glGetUniformLocation (_glProgram [1], "lightSpecular1");
        assert (_glLightSpecLocation1 > -1);
        _glLightExpLocation1 = glGetUniformLocation (_glProgram [1], "lightExp1");
        assert (_glLightExpLocation1 > -1);
      }

      if (_glNumLights > 1){
        _glLightDirLocation2 = glGetUniformLocation (_glProgram [1], "lightDir2");
        assert (_glLightDirLocation2 > -1);
        _glLightAmbLocation2 = glGetUniformLocation (_glProgram [1], "lightAmbient2");
        assert (_glLightAmbLocation2 > -1);
        _glLightDiffLocation2 = glGetUniformLocation (_glProgram [1], "lightDiffuse2");
        assert (_glLightDiffLocation2 > -1);
        _glLightSpecLocation2 = glGetUniformLocation (_glProgram [1], "lightSpecular2");
        assert (_glLightSpecLocation2 > -1);
        _glLightExpLocation2 = glGetUniformLocation (_glProgram [1], "lightExp2");
        assert (_glLightExpLocation2 > -1);
      }

      glBindFragDataLocation (_glProgram [1], 0, "fragColor");
      checkGLError (error);

      glGenVertexArrays (2*_submesh.size (), &(_glRenderVertexArrayId [0]));
      checkGLError (error);

      for (unsigned int i = 0; i < _faceIndices.size (); ++i){

        glBindVertexArray (_glRenderVertexArrayId [i]);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [0]);
        checkGLError (error);
        glVertexAttribPointer (vertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (vertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glNormalTexCoordBufferId);
        checkGLError (error);
        glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (texCoordLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      for (unsigned int i = _faceIndices.size (); i < 2*_faceIndices.size (); ++i){

        glBindVertexArray (_glRenderVertexArrayId [i]);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [1]);
        checkGLError (error);
        glVertexAttribPointer (vertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (vertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glNormalTexCoordBufferId);
        checkGLError (error);
        glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (texCoordLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      if (_glTextureFlag){

        GLint ctexCoordLocation = glGetAttribLocation (_glProgram [1], "inColorTexCoord");
        assert (ctexCoordLocation > -1);
        _glColorTextureLocation1 = glGetUniformLocation (_glProgram [1], "colorTexture");
        assert (_glColorTextureLocation1 > -1);
        _glTexCoordTextureLocation1 = glGetUniformLocation (_glProgram [1], "texCoordTexture");
        assert (_glTexCoordTextureLocation1 > -1);

        for (unsigned int i = 0; i < _faceIndices.size (); ++i){

          glBindVertexArray (_glRenderVertexArrayId [i]);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
          checkGLError (error);
          glVertexAttribPointer (ctexCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (ctexCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, 0);

          glBindVertexArray (_glRenderVertexArrayId [i + _faceIndices.size ()]);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
          checkGLError (error);
          glVertexAttribPointer (ctexCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (ctexCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, 0);
          glBindVertexArray (0);
        }

      } else {
        _glColorLocation1 = glGetUniformLocation (_glProgram [1], "color");
        assert (_glColorLocation1 > -1);
      }
      glUseProgram (0);

      /*************************** INITIALIZE CUT SURFACE RENDERING PROGRAM ****************************/
      if (!initGPUProgram (true, *_glslPrefixString, _glProgramName [2], _glProgram [2])){
        return false;
      }
      glUseProgram (_glProgram [2]);
      checkGLError (error);

      vertLocation = glGetAttribLocation (_glProgram [2], "vertex");
      assert (vertLocation > -1);
      checkGLError (error);
      _glModelviewMatrixLocation2 = glGetUniformLocation (_glProgram [2], "modelview");
      assert (_glModelviewMatrixLocation2 > -1);
      checkGLError (error);
      _glProjectionMatrixLocation2 = glGetUniformLocation (_glProgram [2], "projection");
      assert (_glProjectionMatrixLocation2 > -1);
      checkGLError (error);

      glBindFragDataLocation (_glProgram [2], 0, "fragColor");
      checkGLError (error);

      Submesh *sm;
      for (unsigned int i = 0; i < _submesh.size (); ++i){
        sm = _submesh [i].get ();

        glGenVertexArrays (1, &(sm->_glInRenderVertexArrayId));
        checkGLError (error);
        glBindVertexArray (sm->_glInRenderVertexArrayId);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, sm->_glInVertexBufferId);
        checkGLError (error);
        glVertexAttribPointer (vertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (vertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      if (_glTextureFlag){

        _glColorTextureLocation2 = glGetUniformLocation (_glProgram [2], "colorTexture");
        assert (_glColorTextureLocation2 > -1);
        _glTexCoordTextureLocation2 = glGetUniformLocation (_glProgram [2], "texCoordTexture");
        assert (_glTexCoordTextureLocation2 > -1);

        GLint texCoordCoordLocation = glGetAttribLocation (_glProgram [2], "inTexCoordCoord");
        assert (texCoordCoordLocation > -1);
        GLint surfaceFlagLocation = glGetAttribLocation (_glProgram [2], "inSurfaceFlag");
        assert (surfaceFlagLocation > -1);
        texCoordLocation = glGetAttribLocation (_glProgram [2], "inTexCoord");
        assert (texCoordLocation > -1);

        for (unsigned int i = 0; i < _submesh.size (); ++i){
          sm = _submesh [i].get ();

          glBindVertexArray (sm->_glInRenderVertexArrayId);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, sm->_glIn2DTexCoordBufferId);
          checkGLError (error);
          glVertexAttribPointer (texCoordCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (texCoordCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, sm->_glInVertexStatusBufferId);
          checkGLError (error);
          glVertexAttribPointer (surfaceFlagLocation, 1, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (surfaceFlagLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, sm->_glIn3DTexCoordBufferId);
          checkGLError (error);
          glVertexAttribPointer (texCoordLocation, 3, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (texCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, 0);
          glBindVertexArray (0);
        }
      } else {
        _glColorLocation2 = glGetUniformLocation (_glProgram [2], "color");
      }

      glUseProgram (0);

      return true;
    }

    // private method to initialize texture-related objects
    bool
    Mesh::initGLTextureObjects (unsigned int atlasScaleFactor, const string &atlasShader, const Texture3D &texture)
    {
      // store original number of faces and get three outer rings for each submesh triangles
      {
        vector <vector <unsigned int> > extraFaces;
        extraFaces.reserve (_faceIndices.size ());
        extraFaces.resize (_faceIndices.size ());

        for (unsigned int i = 0; i < extraFaces.size (); ++i){
          getFaceRings (i, _faceIndices, extraFaces [i]);
        }
        for (unsigned int i = 0; i < extraFaces.size (); ++i){
          for (unsigned int j = 0; j < extraFaces [i].size (); ++j){
            _faceIndices [i].push_back (extraFaces [i][j]);
          }
          extraFaces [i].clear ();
        }

        for (unsigned int i = 0; i < extraFaces.size (); ++i){
          getFaceRings (i, _faceIndices, extraFaces [i]);
        }
        for (unsigned int i = 0; i < extraFaces.size (); ++i){
          for (unsigned int j = 0; j < extraFaces [i].size (); ++j){
            _faceIndices [i].push_back (extraFaces [i][j]);
          }
          extraFaces [i].clear ();
        }
      }

      // initialize 2D texture coordinates for each submesh using Tutte's parametric coordinates method
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        calculateParametricCoordinates (_numSurfaceVertices, _vertices [0], _faceIndices [i], _submesh [i].get ()->_meshSurfaceVertexTexCoords);
      }

      // determine scale factors to be used for rasterizing charts
      vector <real> area2d (_faceIndices.size ());
      {
        // calculate the resolution of surface pixels
        int first, second;
        unsigned int numPixels = 0;
        unsigned int offset1 = 4*texture._dimension[0]*texture._dimension[1];
        unsigned int offset2 = 4*texture._dimension[0];
        for (unsigned int i = 0; i < texture._dimension [2]; ++i){
          for (unsigned int j = 0; j < texture._dimension [1]; ++j){
            first = -1;
            second = -1;
            for (unsigned int k = 0; k < texture._dimension [0]; ++k){
              if (texture._rgba [offset1*i + offset2*j + 4*k + 3] > .5){
                first = static_cast <int> (k);
                break;
              }
            }
            if (first >= 0){
              for (int k = static_cast <int> (texture._dimension [0]) - 1; k >=0; --k){
                if (texture._rgba [offset1*i + offset2*j + 4*k + 3] > .5){
                  second = k;
                  break;
                }
              }
              numPixels += first < second ? 2 : 1;
            }
          }
        }

        // calculate scaling factors
        vec te, e1, e2;
        vec2 e3, e4;
        real pixelArea = 0.;
        vector <vec2> *texCoordptr = NULL;
        vector <real>area3d (_faceIndices.size ());

        for (unsigned int i = 0; i < _faceIndices.size (); ++i){
          area2d [i] = 0.;
          area3d [i] = 0.;
          texCoordptr = &(_submesh [i].get ()->_meshSurfaceVertexTexCoords);

          for (unsigned int j = 0; j < _faceIndices [i].size (); j += 3){
            e1 = _vertices [0][_faceIndices [i][j + 1]] - _vertices [0][_faceIndices [i][j]];
            e2 = _vertices [0][_faceIndices [i][j + 2]] - _vertices [0][_faceIndices [i][j]];
            e1.fast_cross (te, e2);
            area3d [i] += te.length ();

            e3 = (*texCoordptr) [_faceIndices[i][j + 1]] - (*texCoordptr) [_faceIndices[i][j]];
            e4 = (*texCoordptr) [_faceIndices[i][j + 2]] - (*texCoordptr) [_faceIndices[i][j]];
            e1 = vec (e3, 0.
#ifdef SF_VECTOR4_ENABLED
                      , 1.
#endif
                      );
            e2 = vec (e4, 0.
#ifdef SF_VECTOR4_ENABLED
                      , 1.
#endif
                      );
            e1.fast_cross (te, e2);
            area2d [i] += te.length ();
          }
          pixelArea += area3d [i];
        }
        pixelArea /= static_cast <real> (numPixels);

        // after this step area2d contains the scale factors for each octant
        for (unsigned int i = 0; i < _faceIndices.size (); ++i){
          area2d [i] = static_cast <real> (sqrt (area3d [i]/ (area2d [i] * pixelArea)));
        }
      }

      // rasterize each flattened 2D submesh and generate rectangular charts
      rasterizeCharts (atlasScaleFactor, atlasShader, texture, area2d);

      // update the 3D texture coordinates for internal vertices
      real min [3] = {2., 2., 2.};
      real max [3] = {-1., -1., -1.};
      for (unsigned int i = 0; i < _texCoords3D.size (); ++i){
        if (_texCoords3D [i]._v [0] < 1.5){
          for (int j = 0; j < 3; ++j){
            if (min [j] > _texCoords3D [i]._v [j]){
              min [j] = _texCoords3D [i]._v [j];
            }
            else if (max [j] < _texCoords3D [i]._v [j]){
              max [j] = _texCoords3D [i]._v [j];
            }
          }
        }
      }
      for (int i = 0; i < 3; ++i){
        max [i] -= min [i];
      }
      real bmin [3] = {_bbox._v [0]._v[0], _bbox._v [0]._v[1], _bbox._v [0]._v[2]};
      real bmax [3] = {_bbox._v [1]._v[0], _bbox._v [1]._v[1], _bbox._v [1]._v[2]};
      for (int j = 0; j < 3; ++j){
        bmax [j] -= bmin [j];
      }
      for (int j = 0; j < 3; ++j){
        bmax [j] = 1./ bmax [j];
      }

      for (unsigned int i = 0; i < _texCoords3D.size (); ++i){
        if (_texCoords3D [i]._v [0] > 1.5){
          for (int j = 0; j < 3; ++j){
            _texCoords3D [i]._v [j] = bmax [j] * (_vertices [0][i]._v [j] - bmin [j]);
            _texCoords3D [i]._v [j] *= max [j];
            _texCoords3D [i]._v [j] += min [j];
          }
        }
      }

      // resize the indices (done to undo the additions made to the indices at the very start of this method)
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        _faceIndices [i].resize (_numFaces [i]);
      }

      return true;
    }

    // private method to rasterize flattened submeshes
    void
    Mesh::rasterizeCharts (unsigned int atlasScale, const string &atlasShader, const Texture3D &texture, const vector <real> &scales)
    {
      // scale vertices to account for aspect ratios of the texture dataset
      vector <vec> normalizedVerts (_numSurfaceVertices, vec::ZERO);
      scaleVertices (&(texture._aspectRatio [0]), _vertices [0], _bbox, normalizedVerts);

      // locally calculate normals
      vector <vec> normals(_numSurfaceVertices, vec::ZERO);
      calculateVertexNormals (normalizedVerts, _faceIndices, normals);

      // scale normals so that they go from [-1:1] to [0, 1]
      for (unsigned int i = 0; i < normals.size (); ++i){
        normals [i] *= .5;
        normals [i] += .5;
      }

      // load GL shader file
      GLuint program = 0;
      initGPUProgram (false, *_glslPrefixString, atlasShader, program);
      assert (program);

      // initialize GL texture atlas and texture coordinate objects
      GLenum error;

      glGenTextures (_faceIndices.size (), &(_gl2DTextureId [0]));
      checkGLError (error);
      glGenBuffers (_faceIndices.size (), &(_glTexCoordBufferId [0]));
      checkGLError (error);

      // variables to be used inside the loop
      vec2 tmpCoord;
      Submesh *sptr = NULL;
      vector <vec2> *texCoordptr = NULL;
      int index, fIndex, dim, intCoord [2];
      real mag, scale, offset, alpha, alphasum, delta [2];
      vector <bool> changeFlag (_numSurfaceVertices);

      // generate rasterized charts for every octant
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){

        // calculate the size of each chart in pixel-dimensions
        index = static_cast <int> (ceil (static_cast <double> (scales [i])));
        dim = 1;
        while (dim < index){
          dim *= 2;
        }
        dim *= static_cast <int> (atlasScale);

        // offset the texture coordinates by size of half a pixel
        for (unsigned int j = 0; j < _numSurfaceVertices; ++j){
          changeFlag [j] = false;
        }
        offset = 1./ static_cast <real> (2*dim);
        scale = 1. - 2. * offset;

        sptr = _submesh [i].get ();
        texCoordptr = &(sptr->_meshSurfaceVertexTexCoords);
        for (unsigned int j = 0; j < _faceIndices [i].size (); ++j){
          if (!changeFlag [_faceIndices [i][j]]) {
            changeFlag [_faceIndices [i][j]] = true;
            for (int k = 0; k < 2; ++k){
              (*texCoordptr) [_faceIndices [i][j]]._v [k] *= scale;
              (*texCoordptr) [_faceIndices [i][j]]._v [k] += offset;
            }
          }
        }

        // generate texture atlases for coordinate and normal data
        vector <GLfloat> rgbaData, coData, noData;
        coData.resize (4*dim*dim, 0.);
        initTextureAtlas (program, dim, normalizedVerts, sptr->_meshSurfaceVertexTexCoords, _faceIndices [i], coData);

        noData.resize (4*dim*dim, 0.);
        initTextureAtlas (program, dim, normals, sptr->_meshSurfaceVertexTexCoords, _faceIndices [i], noData);

        // normalize incoming normals
        for (int j = 0; j < 4*dim*dim; j += 4){
          if (noData [j + 3] > .5) {
            mag = 0.;
            for (int k = 0; k < 3; ++k){
              mag += noData [j + k] * noData [j + k];
            }
            mag = 1./ static_cast <real> (sqrt (mag));
            for (int k = 0; k < 3; ++k){
              noData [j + k] *= mag;
            }
            // rescale normals from [0:1] range to [-1:1]
            for (int k = 0; k < 3; ++k){
              noData [j + k] *= 2.;
              noData [j + k] -= 1.;
            }
          }
        }

        rgbaData.resize (4*dim*dim, 0.);
        raytraceThroughVolumef (dim, coData, noData, texture, rgbaData);

        coData.clear ();
        noData.clear ();

        // calculate 3D texture coordinates for every vertex from the data obtained above
        scale = 1./ scale;
        for (unsigned int j = 0; j < _faceIndices [i].size () - _numFaces [i]; ++j) {
          fIndex = _faceIndices [i][j];
          if (_texCoords3D [fIndex]._v [0] > 1.5){ // do this calculation once for each vertex

            // get texture coordinate needed to access texture atlas
            tmpCoord = (*texCoordptr) [j];
            tmpCoord *= scale;
            tmpCoord -= offset;
            tmpCoord *= static_cast <real> (dim - 1);

            // generate rounded integer coordinate
            for (int k = 0; k < 2; ++k){
              intCoord [k] = static_cast <int> (floor (tmpCoord._v [k]));
              delta [k] = tmpCoord._v [k] - intCoord [k];
            }

            // do bilinear interpolation to get final 3D texture coordinate
            index = 4* (dim*intCoord [1] + intCoord [0]);
            alphasum = alpha = (1. - delta [0]) * (1. - delta [1]) * rgbaData [index + 3];
            for (int k = 0; k < 3; ++k){
              _texCoords3D [fIndex]._v [k] = alpha * rgbaData [index + k];
            }

            if (intCoord [0] < dim - 1){
              alpha = delta [0] * (1. - delta [1]) * rgbaData [index + 4 + 3];
              for (int k = 0; k < 3; ++k){
                _texCoords3D [fIndex]._v [k] += alpha * rgbaData [index + 4 + k];
              }
              alphasum += alpha;
            }

            if (intCoord [1] < dim - 1){
              alpha = (1. - delta [0]) * delta [1] * rgbaData [index + 4*dim + 3];
              for (int k = 0; k < 3; ++k){
                _texCoords3D [fIndex]._v [k] += alpha * rgbaData [index + 4*dim + k];
              }
              alphasum += alpha;
            }

            if (intCoord [0] < dim - 1 && intCoord [1] < dim - 1){
              alpha = delta [0] * delta [1] * rgbaData [index + 4*dim + 4 + 3];
              for (int k = 0; k < 3; ++k){
                _texCoords3D [fIndex]._v [k] += alpha * rgbaData [index + 4*dim + 4 + k];
              }
              alphasum += alpha;
            }

            if (alphasum > EPSILON){
              alphasum = 1./ alphasum;
              for (int k = 0; k < 3; ++k){
                _texCoords3D [fIndex]._v [k] *= alphasum;
              }
            }
          } // end - if (_texCoords3D [_faceIndices [i][j]]._v [0] > 1.5)
        } // end - for (unsigned int j = 0; j < _faceIndices [i].size () - _numFaces [i]; ++j)

        // generate 2D texture atlas containing texture coordinates
        glBindTexture (GL_TEXTURE_2D, _gl2DTextureId [i]);
        checkGLError (error);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, dim, dim, 0, GL_RGBA, GL_FLOAT, &(rgbaData [0]));
        checkGLError (error);

        rgbaData.clear ();

        // generate 2D texture coordinate objects
        glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*_numSurfaceVertices, &(texCoordptr->at (0)), GL_STATIC_DRAW);
        checkGLError (error);

      } // end - for (unsigned int i = 0; i < _faceIndices.size (); ++i)

      // cleanup
      glBindTexture (GL_TEXTURE_2D, 0);
      glBindBuffer (GL_ARRAY_BUFFER, 0);

      // populate 3D texture
      glGenTextures (1, &_gl3DTextureId);
      checkGLError (error);
      glBindTexture (GL_TEXTURE_3D, _gl3DTextureId);
      checkGLError (error);
      glTexParameteri (GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri (GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri (GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri (GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri (GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
      glTexImage3D (GL_TEXTURE_3D, 0, GL_RGBA, texture._dimension [0], texture._dimension [1], texture._dimension [2], 0, GL_RGBA, GL_UNSIGNED_BYTE, &(texture._rgba [0]));
      checkGLError (error);
      glBindTexture (GL_TEXTURE_3D, 0);
    }

	} // end - namespace XFE
} // end - namespace SF
