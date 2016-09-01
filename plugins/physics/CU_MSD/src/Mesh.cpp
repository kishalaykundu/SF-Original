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
 * The mesh class for the CU_MSD library. This is derived class from
 * Resource.
 */

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "Preprocess.h"

#define SF_MSD_CUDA_THREADS_PER_BLOCK 256

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
#include "Mesh.h"

extern "C" {
#include <cuda_runtime_api.h>
}

namespace SF {

	namespace MSD {

	  static int GLX_ATTRIBUTE_LIST [] = {GLX_RGBA, None};

	  // static function to reload GPU programs
	  static void reloadPrograms (Resource & r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);
	    mptr->_glReprogramFlag = true;
	  }

    // inline function to draw normals
    static void drawNormals (Mesh *mptr)
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

		// function to calculate vertex forces using Hooke's Law
		static void calcForces (Mesh *mptr)
		{
#ifndef NDEBUG
      GLenum error;
#endif

      glUseProgram (mptr->_glProgram [2]);
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
      glBindFramebuffer (GL_FRAMEBUFFER, mptr->_glForceFrameBufferId);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glPushAttrib (GL_VIEWPORT_BIT);
      glViewport (0, 0, mptr->_glForceFrameBufferDimensions [0], mptr->_glForceFrameBufferDimensions [1]);
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
        glBindVertexArray (mptr->_glForceVertexArrayId [0]);
      } else {
        glBindVertexArray (mptr->_glForceVertexArrayId [1]);
      }
#ifndef NDEBUG
      checkGLError (error);
#endif

        glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glSpringIndexBufferId);
        glDrawElements (GL_LINES, mptr->_numSprings, GL_UNSIGNED_INT, 0);
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

      /*************************** CALCULATION OF NORMALS HAPPENS HERE ****************************/
      drawNormals (mptr);

      /************************ RENDERING OF EXTERNAL SURFACE HAPPENS HERE ************************/
#ifndef NDEBUG
      GLenum error;
#endif

      glUseProgram (mptr->_glProgram [1]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation, 1, false, mptr->_glProjection);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniform3f (mptr->_glColorLocation, mptr->_glColor [0], mptr->_glColor [1], mptr->_glColor [2]);
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
      glUniform1i (mptr->_glNormalTextureLocation, 0);
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

      /*************************** CALCULATION OF FORCES HAPPENS HERE *****************************/
			calcForces (mptr);

      mptr->_syncControl [mptr->_semGraphicsPostIndex].post ();
	  }

    // drawing function textured datasets
	  static void texturedDraw (Resource &r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);

	    mptr->_syncControl [mptr->_semGraphicsWaitIndex].wait ();

	    // reload program if needed
	    if (mptr->_glReprogramFlag) {
        mptr->initGPUPrograms ();
        mptr->_glReprogramFlag = false;
	    }

      /**************************** CALCULATION OF NORMAL HAPPENS HERE ****************************/
      drawNormals (mptr);

      /************************ RENDERING OF EXTERNAL SURFACE HAPPENS HERE ************************/
#ifndef NDEBUG
      GLenum error;
#endif
      glUseProgram (mptr->_glProgram [1]);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glUniformMatrix4fv (mptr->_glModelviewMatrixLocation, 1, false, mptr->_glModelview);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniformMatrix4fv (mptr->_glProjectionMatrixLocation, 1, false, mptr->_glProjection);
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
      glBindTexture (GL_TEXTURE_2D, mptr->_glNormalTextureId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glUniform1i (mptr->_glNormalTextureLocation, 0);
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
        glActiveTexture (GL_TEXTURE2);
        glBindTexture (GL_TEXTURE_2D, mptr->_glTextureId [i]);
#ifndef NDEBUG
      checkGLError (error);
#endif
        glUniform1i (mptr->_glColorTextureLocation, 2);
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

      /*************************** CALCULATION OF FORCES HAPPENS HERE *****************************/
			calcForces (mptr);

	    mptr->_syncControl [mptr->_semGraphicsPostIndex].post ();
	  }

		// protected constructor and assignment functions
		Mesh::Mesh () { }
		Mesh::Mesh (const Mesh &m) { }
		Mesh& Mesh::operator = (const Mesh &m) { return *this; }

		// only legitimate constructor
		Mesh::Mesh (const string &config, Driver &driver)
		: _semPhysicsWaitIndex (-1), _semPhysicsPostIndex (-1),
		  _semGraphicsWaitIndex (-1), _semGraphicsPostIndex (-1),
		  _numSurfaceVertices (0), _numTotalVertices (0), _curr (&(_vertices [0])), _prev (&(_vertices [1])), _mass (NULL), _numSprings (0),
		  _past (boost::posix_time::microsec_clock::universal_time ()), _present (boost::posix_time::microsec_clock::universal_time ()),
		  _glBufferFlag (false), _glTextureFlag (false), _glReprogramFlag (false),
		  _glNormalFramebufferId (0), _glNormalTexCoordBufferId (0), _glNormalTextureId (0), _glNumFaces (0), _glNormalIndexBufferId (0),
		  _glSpringIndexBufferId (0), _glForceFrameBufferId (0), _glForceTexCoordBufferId (0), _glForceTextureId (0),
			_glEnvTextureId (driver._display.get ()->_glEnvTextureId),
		  _glModelviewMatrixLocation (-1), _glProjectionMatrixLocation (-1), _glNormalTextureLocation (-1), _glColorTextureLocation (-1), _glColorLocation (-1),
		  _glLightDirLocation1 (-1), _glLightAmbLocation1 (-1), _glLightDiffLocation1 (-1), _glLightSpecLocation1 (-1), _glLightExpLocation1 (-1),
		  _glLightDirLocation2 (-1), _glLightAmbLocation2 (-1), _glLightDiffLocation2 (-1), _glLightSpecLocation2 (-1), _glLightExpLocation2 (-1),
		  _glEnvTextureLocation (-1),
		  _glNumLights (driver._display.get ()->_numLights),
		  _glLightDir1 (NULL), _glLightAmb1 (NULL), _glLightDiff1 (NULL), _glLightSpec1 (0.), _glLightExp1 (0.),
		  _glLightDir2 (NULL), _glLightAmb2 (NULL), _glLightDiff2 (NULL), _glLightSpec2 (0.), _glLightExp2 (0.),
		  _glModelview (&(driver._display.get ()->_modelview [0])), _glProjection (&(driver._display.get ()->_projection [0])),
		  _glslPrefixString (&(driver._display.get ()->_glslPrefixString))
		{
			assert (!config.empty ());

			_owner = boost::shared_ptr <string> (new string ("CudaMsd"));

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

        unsigned int numPartitions = 1;
				for (unsigned int i = 0; i < depth; ++i){
					numPartitions *= 8;
				}

				// generate file prefix
				string prefix (folder);
				prefix.append (name);

				/*************************** READ NODE (VERTEX) FILE ***************************/
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
				_numTotalVertices = _vertices [0].size ();

				/*************************** READ VERTEX MASS RECIPROCAL FILE ***************************/
				file = prefix;
				file.append (".lm");

				fp = fopen (file.c_str (), "r");
				assert (fp);

				status = fscanf (fp, "%d\n", &(tmpd [0]));
				assert (status != 0);
				if (tmpd [0] <= 0){
          PRINT ("fatal error: invalid number of vertex masses \'%d\' in %s\n", tmpd [0], file.c_str ());
          exit (EXIT_FAILURE);
				}
				{
				  unsigned int nmass = static_cast <unsigned int> (tmpd [0]);
				  _mass = new real [nmass];

				  for (unsigned int i = 0; i < nmass; ++i){
#ifdef SF_DOUBLE_PRECISION
            status = fscanf (fp, "%lf\n", &(_mass [i]));
#else
            status = fscanf (fp, "%f\n", &(_mass [i]));
#endif
            assert (status != 0);
				  }
				}


				fclose (fp);

				/*************************** READ EDGE ELEMENT FILES ***************************/
				{
					file = prefix;
					file.append (".edge");

					fp = fopen (file.c_str (), "r");
					assert (fp);

					int status = fscanf (fp, "%d\n", &(tmpd [0]));
					assert (status != 0);
					if (tmpd [0] <= 0){
		        PRINT ("fatal error: invalid number of vertices \'%d\' in %s\n", tmpd [0], file.c_str ());
		        exit (EXIT_FAILURE);
					}

#ifndef NDEBUG
            int vertSize = static_cast <int> (_vertices [0].size ());
#endif

					_numSprings = static_cast <unsigned int> (tmpd [0]);

					vector <unsigned int> springIndices (_numSprings * 2);
					for (unsigned int i = 0; i < _numSprings; ++i){
						status = fscanf (fp, "%d %d\n", &(tmpd [0]), &(tmpd [1]));
						assert (status != 0);
						assert (tmpd [0] >= 0 && tmpd [0] < vertSize);
						assert (tmpd [1] >= 0 && tmpd [1] < vertSize);
						springIndices [2*i] = tmpd [0];
						springIndices [2*i + 1] = tmpd [1];
					}

					fclose (fp);
				}

				/*************************** READ TRIANGULAR ELEMENT FILES ***************************/
				_numFaces.resize (numPartitions, 0);
				_faceIndices.resize (numPartitions, vector <unsigned int> ());

				for (unsigned int i = 0; i < numPartitions; ++i){

					// generate file name
					file = prefix;
					file.append (".");

					char indexStr [4];
					sprintf (indexStr, "%u", i);
					file.append (indexStr);
					file.append (".tri");

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
						_numFaces [i] *= 3;
					}
					fclose (fp);
				} // end - for (unsigned int i = 0; i < numSubmeshes; ++i)

				++_numSurfaceVertices; // this is done because before this _numSurfaceVertices contains biggest index of triangles
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
      GLuint tmpu = 0;
      _glIndexBufferId.resize (_faceIndices.size (), tmpu);
      _glTexCoordBufferId.resize (_faceIndices.size (), tmpu);
      _glTextureId.resize (_faceIndices.size (), tmpu);
      _glRenderVertexArrayId.resize (_faceIndices.size (), tmpu);

      // get names of the shading programs
      getConfigParameter (config, "normal_shader", _glProgramName [0]);
      assert (!_glProgramName [0].empty ());

      getConfigParameter (config, "color_shader", _glProgramName [1]);
      assert (!_glProgramName [1].empty ());

      getConfigParameter (config, "force_shader", _glProgramName [2]);
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

      /**************************** INITIALIZE GLX RELATED PARAMETERS *****************************/
      _glContext = glXGetCurrentContext ();
      _glDisplay = glXGetCurrentDisplay ();
      _glDrawable = glXGetCurrentDrawable ();

      /**************************** INITIALIZE CUDA RELATED PARAMETERS ****************************/
      getConfigParameter (config, "compute_program", _cuPTXFileName);
      assert (!_cuPTXFileName.empty ());

      getConfigParameter (config, "compute_function0", _cuKernelFuncName [0]);
      assert (!_cuKernelFuncName [0].empty ());

      getConfigParameter (config, "compute_function1", _cuKernelFuncName [1]);
      assert (!_cuKernelFuncName [1].empty ());

		} // end - Mesh::Mesh (const string &config, Driver &driver)

		// class destructor
		Mesh::~Mesh ()
		{
		}

		// CUDA-enabled run method (run within individual thread)
		void
		Mesh::run ()
		{
		  // initialize all pertinent CUDA resources
		  {
        CUresult status;

        XVisualInfo *visualInfo = glXChooseVisual (_glDisplay, 0, GLX_ATTRIBUTE_LIST);
        GLenum error;
        checkGLError (error);

        _cuglContext = glXCreateContext (_glDisplay, visualInfo, _glContext, GL_TRUE);
        checkGLError (error);
        if (!glXMakeCurrent (_glDisplay, _glDrawable, _cuglContext)){
          checkGLError (error);
          PRINT ("glXMakeCurrent error\n");
        }

        status = cuInit (0);
        printCUResult (status);

        // check if any device is CUDA capable
        int numCudaDevices = 0;
        status = cuDeviceGetCount (&numCudaDevices);
        printCUResult (status);
        if (!numCudaDevices){
          PRINT ("Error: Could not find CUDA devices\n");
        }

        // get handle for device 0
        CUdevice cudaDevice;
        status = cuDeviceGet (&cudaDevice, 0);
        printCUResult (status);

        // calculate CUDA block, grid and threads-per-block sizes
        CUdevprop devProps;
        status = cuDeviceGetProperties (&devProps, cudaDevice);
        printCUResult (status);

        // create context
        status = cuGLCtxCreate (&_cuContext, CU_CTX_SCHED_BLOCKING_SYNC, cudaDevice);
        printCUResult (status);
      PRINT ("HERE\n");

        status = cuCtxPushCurrent (_cuContext);
        printCUResult (status);
      PRINT ("HERE\n");

        status = cuCtxSynchronize ();
        printCUResult (status);
      PRINT ("HERE\n");

        // load CUDA module and get required function handles
        CUmodule cuMod;
        status = cuModuleLoad (&cuMod, _cuPTXFileName.c_str ());
        printCUResult (status);
        status = cuModuleGetFunction (&_cuKernelFunc [0], cuMod, _cuKernelFuncName [0].c_str ());
        printCUResult (status);
        status = cuModuleGetFunction (&_cuKernelFunc [1], cuMod, _cuKernelFuncName [1].c_str ());
        printCUResult (status);

        // register the vertex buffers and force-framebuffer
        cuGraphicsGLRegisterBuffer (&_cuVertexBufferId [0], _glVertexBufferId [0], CU_GRAPHICS_REGISTER_FLAGS_NONE);
        printCUResult (status);
        cuGraphicsGLRegisterBuffer (&_cuVertexBufferId [1], _glVertexBufferId [1], CU_GRAPHICS_REGISTER_FLAGS_NONE);
        printCUResult (status);
        cuGraphicsGLRegisterBuffer (&_cuForceTexCoordBufferId, _glForceTexCoordBufferId, CU_GRAPHICS_REGISTER_FLAGS_NONE);
        printCUResult (status);
        cuGraphicsGLRegisterImage (&_cuForceTextureId, _glForceTextureId, GL_TEXTURE_2D, CU_GRAPHICS_REGISTER_FLAGS_READ_ONLY);
        printCUResult (status);

        // get a handle for the texture in the CUDA module and sync it with the texture reference on host
        status = cuModuleGetTexRef (&_devTexRef, cuMod, "forceTexture");
        printCUResult (status);
        status = cuTexRefSetArray (_devTexRef, _devForceTexturePtr, CU_TRSA_OVERRIDE_FORMAT);
        printCUResult (status);
        status = cuTexRefSetAddressMode (_devTexRef, 0, CU_TR_ADDRESS_MODE_WRAP);
        printCUResult (status);
        status = cuTexRefSetAddressMode (_devTexRef, 1, CU_TR_ADDRESS_MODE_WRAP);
        printCUResult (status);
        status = cuTexRefSetFilterMode (_devTexRef, CU_TR_FILTER_MODE_LINEAR);
        printCUResult (status);
        status = cuTexRefSetFlags (_devTexRef, CU_TRSF_NORMALIZED_COORDINATES);
        printCUResult (status);
        status = cuTexRefSetFormat (_devTexRef, CU_AD_FORMAT_FLOAT, 4);
        printCUResult (status);

        status = cuParamSetTexRef (_cuKernelFunc [0], CU_PARAM_TR_DEFAULT, _devTexRef);
        printCUResult (status);
        status = cuParamSetTexRef (_cuKernelFunc [1], CU_PARAM_TR_DEFAULT, _devTexRef);
        printCUResult (status);

        // initialize the mass-reciprocal device pointer (exclusively used by CUDA)
        status = cuMemAlloc (&_devMassArrayPtr, _numTotalVertices * sizeof (real));
        printCUResult (status);
        status = cuMemcpyHtoD (_devMassArrayPtr, _mass, _numTotalVertices * sizeof (real));
        printCUResult (status);
        delete [] _mass;

        // release context
        cuCtxPopCurrent (&_cuContext);
		  }

      PRINT ("HERE\n");
      // set kernel arguments used by CUDA function (CUDA 4.x notation)
      size_t argBufferSize = 0;
      char argBuffer1 [256];
      char argBuffer2 [256];

      *((CUdeviceptr *) &argBuffer1 [argBufferSize]) = _devVertexBufferPtr [0];
      *((CUdeviceptr *) &argBuffer2 [argBufferSize]) = _devVertexBufferPtr [1];
      argBufferSize += sizeof (_devVertexBufferPtr [0]);
      *((CUdeviceptr *) &argBuffer1 [argBufferSize]) = _devVertexBufferPtr [1];
      *((CUdeviceptr *) &argBuffer2 [argBufferSize]) = _devVertexBufferPtr [0];
      argBufferSize += sizeof (_devVertexBufferPtr [0]);
      *((CUdeviceptr *) &argBuffer1 [argBufferSize]) = _devForceTexCoordBufferPtr;
      *((CUdeviceptr *) &argBuffer2 [argBufferSize]) = _devForceTexCoordBufferPtr;
      argBufferSize += sizeof (_devForceTexCoordBufferPtr);
      *((CUdeviceptr *) &argBuffer1 [argBufferSize]) = _devMassArrayPtr;
      *((CUdeviceptr *) &argBuffer2 [argBufferSize]) = _devMassArrayPtr;
      argBufferSize += sizeof (_devMassArrayPtr);

      void *_cuConfig1 [] = {CU_LAUNCH_PARAM_BUFFER_POINTER, argBuffer1, CU_LAUNCH_PARAM_BUFFER_SIZE, &argBufferSize, CU_LAUNCH_PARAM_END};
      void *_cuConfig2 [] = {CU_LAUNCH_PARAM_BUFFER_POINTER, argBuffer2, CU_LAUNCH_PARAM_BUFFER_SIZE, &argBufferSize, CU_LAUNCH_PARAM_END};

      size_t numBytes = 0;
      size_t problemSize = (_numTotalVertices + SF_MSD_CUDA_THREADS_PER_BLOCK - 1)/ SF_MSD_CUDA_THREADS_PER_BLOCK;

      CUresult status;
      cudaError_t error;
      unsigned int numIters = 0;

      while (true){

        // wait for physics end to get control
        _syncControl [_semPhysicsWaitIndex].wait ();

        // calculate time-differences
        _past = _present;
        _present = boost::posix_time::microsec_clock::universal_time ();
        _deltaT0 = _deltaT1;
        _deltaT1 = _present - _past;

        // push CUDA context
        status = cuCtxPushCurrent (_cuContext);
#ifndef NDEBUG
        printCUResult (status);
#endif

        // map GL vertex buffers and texture
        status = cuGraphicsMapResources (2, _cuVertexBufferId, 0);
#ifndef NDEBUG
        printCUResult (status);
#endif
        status = cuGraphicsMapResources (1, &_cuForceTextureId, 0);
#ifndef NDEBUG
        printCUResult (status);
#endif

        // retrieve device pointers to access the mapped graphics resources
        status = cuGraphicsResourceGetMappedPointer (&(_devVertexBufferPtr [0]), &numBytes, _cuVertexBufferId [0]);
#ifndef NDEBUG
        printCUResult (status);
#endif
        status = cuGraphicsResourceGetMappedPointer (&(_devVertexBufferPtr [1]), &numBytes, _cuVertexBufferId [1]);
#ifndef NDEBUG
        printCUResult (status);
#endif
        status = cuGraphicsResourceGetMappedPointer (&_devForceTexCoordBufferPtr, &numBytes, _cuForceTexCoordBufferId);
#ifndef NDEBUG
        printCUResult (status);
#endif
        status = cuGraphicsSubResourceGetMappedArray (&_devForceTexturePtr, _cuForceTextureId, 0, 0);
#ifndef NDEBUG
        printCUResult (status);
#endif

        // copy constants to the device
        double deltaT0 = static_cast <double> (_deltaT0.total_nanoseconds ())*1.e-9;
        double deltaT1 = static_cast <double> (_deltaT1.total_nanoseconds ())*1.e-9;

        real factor0 = deltaT1/ deltaT0;
        real factor1 = deltaT1 * deltaT1;

        error = cudaMemcpyToSymbol ("factor0", &factor0, sizeof (real), 0, cudaMemcpyHostToDevice);
#ifndef NDEBUG
        checkCUDAError (error);
#endif
        error = cudaMemcpyToSymbol ("factor1", &factor1, sizeof (real), 0, cudaMemcpyHostToDevice);
#ifndef NDEBUG
        checkCUDAError (error);
#endif
        error = cudaMemcpyToSymbol ("N", &_numTotalVertices, sizeof (unsigned int), 0, cudaMemcpyHostToDevice);
#ifndef NDEBUG
        checkCUDAError (error);
#endif

        // launch kernel
        if (numIters > 1){
          if (_glBufferFlag){
            status = cuLaunchKernel (_cuKernelFunc [1], problemSize, 1, 1, SF_MSD_CUDA_THREADS_PER_BLOCK, 1, 1, 0, 0, NULL, _cuConfig1);
#ifndef NDEBUG
            printCUResult (status);
#endif
          } else {
            status = cuLaunchKernel (_cuKernelFunc [1], problemSize, 1, 1, SF_MSD_CUDA_THREADS_PER_BLOCK, 1, 1, 0, 0, NULL, _cuConfig2);
#ifndef NDEBUG
            printCUResult (status);
#endif
          }
        } else {
          if (_glBufferFlag){
            status = cuLaunchKernel (_cuKernelFunc [0], problemSize, 1, 1, SF_MSD_CUDA_THREADS_PER_BLOCK, 1, 1, 0, 0, NULL, _cuConfig1);
#ifndef NDEBUG
            printCUResult (status);
#endif
          } else {
            status = cuLaunchKernel (_cuKernelFunc [0], problemSize, 1, 1, SF_MSD_CUDA_THREADS_PER_BLOCK, 1, 1, 0, 0, NULL, _cuConfig2);
#ifndef NDEBUG
            printCUResult (status);
#endif
          }
          ++numIters;
        }

        // wait for everything to synchronize at this point
        status = cuCtxSynchronize ();
#ifndef NDEBUG
        printCUResult (status);
#endif

        // unmap GL vertex buffers and texture
        status = cuGraphicsUnmapResources (2, _cuVertexBufferId, 0);
#ifndef NDEBUG
        printCUResult (status);
#endif
        status = cuGraphicsUnmapResources (1, &_cuForceTextureId, 0);
#ifndef NDEBUG
        printCUResult (status);
#endif


        // pop CUDA context
        cuCtxPopCurrent (&_cuContext);

        // swap buffers
        vector <vec> *tmp = _curr;
        _curr = _prev;
        _prev = tmp;

        // toggle swap buffer flag
        _glBufferFlag = !_glBufferFlag; // should be the last line in this segment

        // release data
        _syncControl [_semPhysicsPostIndex].post ();
      }
		}

    // method to do cleanup CUDA related resources (run from within threads)
    void
    Mesh::cleanup ()
    {
    }

		// method to initialize all the GPU programs
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
      _glModelviewMatrixLocation = glGetUniformLocation (_glProgram [1], "modelview");
      assert (_glModelviewMatrixLocation > -1);
      _glProjectionMatrixLocation = glGetUniformLocation (_glProgram [1], "projection");
      assert (_glProjectionMatrixLocation > -1);
      _glNormalTextureLocation = glGetUniformLocation (_glProgram [1], "normalTexture");
      assert (_glNormalTextureLocation > -1);

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

      glGenVertexArrays (2*_faceIndices.size (), &(_glRenderVertexArrayId [0]));
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

        texCoordLocation = glGetAttribLocation (_glProgram [1], "inColorTexCoord");
        assert (texCoordLocation > -1);
        _glColorTextureLocation = glGetUniformLocation (_glProgram [1], "colorTexture");
        assert (_glColorTextureLocation > -1);

        for (unsigned int i = 0; i < _faceIndices.size (); ++i){

          glBindVertexArray (_glRenderVertexArrayId [i]);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
          checkGLError (error);
          glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (texCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, 0);

          glBindVertexArray (_glRenderVertexArrayId [i + _faceIndices.size ()]);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
          checkGLError (error);
          glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
          checkGLError (error);
          glEnableVertexAttribArray (texCoordLocation);
          checkGLError (error);

          glBindBuffer (GL_ARRAY_BUFFER, 0);
          glBindVertexArray (0);
        }

      } else {
        _glColorLocation = glGetUniformLocation (_glProgram [1], "color");
        assert (_glColorLocation > -1);
      }
      glUseProgram (0);

      /*************************** INITIALIZE FORCE CALCULATION PROGRAM ***************************/
      if (!initGPUProgram (true, *_glslPrefixString, _glProgramName [2], _glProgram [2])){
        PRINT ("error: could not initialize %s\n", _glProgramName [2].c_str ());
        return false;
      }

      // initialize vertex array for force calculation program
      glUseProgram (_glProgram [2]);
      checkGLError (error);

      vertLocation = glGetAttribLocation (_glProgram [2], "vertex");
      assert (vertLocation > -1);
      GLint restVertLocation = glGetAttribLocation (_glProgram [2], "restVertex");
      assert (vertLocation > -1);
      texCoordLocation = glGetAttribLocation (_glProgram [2], "inTexCoord");
      assert (texCoordLocation > -1);

      glBindFragDataLocation (_glProgram [2], 0, "fragColor");
      checkGLError (error);

      glGenVertexArrays (2, _glForceVertexArrayId);
      checkGLError (error);

      for (int i = 0; i < 2; ++i){

        glBindVertexArray (_glForceVertexArrayId [i]);
        checkGLError (error);
        glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [i]);
        checkGLError (error);
        glVertexAttribPointer (vertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (vertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glRestVertexBufferId);
        checkGLError (error);
        glVertexAttribPointer (restVertLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (restVertLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, _glForceTexCoordBufferId);
        checkGLError (error);
        glVertexAttribPointer (texCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
        checkGLError (error);
        glEnableVertexAttribArray (texCoordLocation);
        checkGLError (error);

        glBindBuffer (GL_ARRAY_BUFFER, 0);
        glBindVertexArray (0);
      }

      glUseProgram (0);

      return true;
		}

		// private method to initialize the force-calculation related GL buffer objects
		bool
		Mesh::initGLForceBufferObjects (const vector <unsigned int> &springs)
		{
      unsigned int width = static_cast <unsigned int> (trunc (ceil (sqrt (static_cast <double> (_numSprings)))));
      unsigned int height = static_cast <unsigned int> (trunc (floor (sqrt (static_cast <double> (_numSprings)))));

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
      while (width*height > _numSprings){
        height /= 2;
      }
      if (width*height < _numSprings){
        height *= 2;
      }
      while (width*height > _numSprings){
        width /= 2;
      }
      if (width*height < _numSprings){
        width *= 2;
      }

      _glForceFrameBufferDimensions [0] = static_cast <GLuint> (width);
      _glForceFrameBufferDimensions [1] = static_cast <GLuint> (height);

      GLenum error;

      glGenTextures (1, &_glForceTextureId);
      checkGLError (error);
      glBindTexture (GL_TEXTURE_2D, _glForceTextureId);
      checkGLError (error);

      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA32F, _glForceFrameBufferDimensions [0], _glForceFrameBufferDimensions [1], 0, GL_RGBA, GL_FLOAT, 0);
      checkGLError (error);

      glGenFramebuffers (1, &_glForceFrameBufferId);
      checkGLError (error);
      glBindFramebuffer (GL_FRAMEBUFFER, _glForceFrameBufferId);
      checkGLError (error);
      glFramebufferTexture2D (GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, _glForceTextureId, 0);
      checkGLError (error);
      glBindFramebuffer (GL_FRAMEBUFFER, 0);

      // initialize texture coordinates for every vertex into the force texture
      vector <vec2> texcoords;
      texcoords.resize (_numSprings, vec2::ZERO);

      unsigned int xcount = 0, ycount = 0;
      for (unsigned int i = 0; i < _numSprings; ++i){
        texcoords [i] = vec2 (static_cast <real> (xcount)/ static_cast <real> (width), static_cast <real> (ycount)/ static_cast <real> (height));
        ++xcount;
        if (xcount >= width){
          xcount = 0;
          ++ycount;
        }
      }

      // advance each coordinate by half pixel to align with GL texture
      real xstep = 1./ static_cast <real> (2*_glForceFrameBufferDimensions [0]);
      real ystep = 1./ static_cast <real> (2*_glForceFrameBufferDimensions [1]);
      for (unsigned int i = 0; i < _numSprings; ++i){
        texcoords [i]._v [0] += xstep;
        texcoords [i]._v [1] += ystep;
      }

      // generate the GL buffers
      glGenBuffers (1, &_glForceTexCoordBufferId);
      checkGLError (error);
      glBindBuffer (GL_ARRAY_BUFFER, _glForceTexCoordBufferId);
      checkGLError (error);
      glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*texcoords.size (), &(texcoords [0]), GL_STATIC_DRAW);
      checkGLError (error);
      glBindBuffer (GL_ARRAY_BUFFER, 0);

      // add springs to the GL force index buffer
      glGenBuffers (1, &_glSpringIndexBufferId);
      checkGLError (error);
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glSpringIndexBufferId);
      checkGLError (error);
      glBufferData (GL_ELEMENT_ARRAY_BUFFER, 2 * _numSprings * sizeof (unsigned int), &(springs [0]), GL_STATIC_DRAW);
      checkGLError (error);
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);

      _numSprings *= 2;

      return true;
		}

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

		  glGenBuffers (1, &_glRestVertexBufferId);
		  checkGLError (error);
		  glBindBuffer (GL_ARRAY_BUFFER, _glRestVertexBufferId);
		  checkGLError (error);
		  glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*_vertices [0].size (), &(_vertices [0][0]), GL_STATIC_DRAW);
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

    // private method to initialize texture-related objects
    bool
    Mesh::initGLTextureObjects (unsigned int atlasScaleFactor, const string &atlasShader, const Texture3D &texture)
    {
      // store original number of faces and get three outer rings for each submesh triangles
      {
        vector <vector <unsigned int> > extraFaces;
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
      vector <vector <vec2> > texCoords (_faceIndices.size (), vector <vec2> (_numSurfaceVertices));
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        calculateParametricCoordinates (_numSurfaceVertices, _vertices [0], _faceIndices [i], texCoords [i]);
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
          texCoordptr = &(texCoords [i]);

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
      rasterizeCharts (atlasScaleFactor, atlasShader, texture, area2d, texCoords);

      // resize the indices (done to undo the additions made to the indices at the very start of this method)
      for (unsigned int i = 0; i < _faceIndices.size (); ++i){
        _faceIndices [i].resize (_numFaces [i]);
      }

      return true;
    }

    // private method to rasterize flattened submeshes
    void
    Mesh::rasterizeCharts (unsigned int atlasScale, const string &atlasShader, const Texture3D &texture, const vector <real> &scales, vector <vector <vec2> > &texCoords)
    {
      // scale vertices to account for aspect ratios of the texture dataset
      //vector <vec> normalizedVerts (_numSurfaceVertices, vec::ZERO);
      //scaleVertices (&(texture._aspectRatio [0]), _vertices [0], _bbox, normalizedVerts);

      vector <vec> normalizedVerts (_vertices [0]);
      normalizedVerts.resize (_numSurfaceVertices);

      // locally calculate normals
      vector <vec> normals (_numSurfaceVertices, vec::ZERO);
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

      glGenTextures (_faceIndices.size (), &(_glTextureId [0]));
      checkGLError (error);
      glGenBuffers (_faceIndices.size (), &(_glTexCoordBufferId [0]));
      checkGLError (error);

      // variables to be used inside the loop
      vec2 tmpCoord;
      vector <vec2> *texCoordptr = NULL;
      int index, dim;
      real mag, scale, offset;
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

        texCoordptr = &(texCoords [i]);
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
        vector <GLubyte> rgbaData;
        vector <GLfloat> coData, noData;
        coData.resize (4*dim*dim, 0.);
        initTextureAtlas (program, dim, normalizedVerts, texCoords [i], _faceIndices [i], coData);

        noData.resize (4*dim*dim, 0.);
        initTextureAtlas (program, dim, normals, texCoords [i], _faceIndices [i], noData);

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

        rgbaData.resize (4*dim*dim, 0);
        raytraceThroughVolumeb (dim, coData, noData, texture, rgbaData);

        coData.clear ();
        noData.clear ();

        // generate 2D texture atlas containing texture coordinates
        glBindTexture (GL_TEXTURE_2D, _glTextureId [i]);
        checkGLError (error);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, dim, dim, 0, GL_RGBA, GL_UNSIGNED_BYTE, &(rgbaData [0]));
        checkGLError (error);

        rgbaData.clear ();

        // generate 2D texture coordinate objects
        glBindBuffer (GL_ARRAY_BUFFER, _glTexCoordBufferId [i]);
        checkGLError (error);
        glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*_numSurfaceVertices, &(texCoords [0]), GL_STATIC_DRAW);
        checkGLError (error);

      } // end - for (unsigned int i = 0; i < _faceIndices.size (); ++i)

      // cleanup
      glBindTexture (GL_TEXTURE_2D, 0);
      glBindBuffer (GL_ARRAY_BUFFER, 0);
    }

	} // end - namespace MSD
} // end - namespace SF

