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
 * The mesh class for the Rigid library.
 */

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glx.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#endif
}

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

#include "mat4x4.h"
#include "aabb.h"
#include "GL/common.h"

#include "ThreadControl.h"
#include "Driver.h"
#include "Display.h"

#include "Common.h"
#include "Mesh.h"

namespace SF {

  namespace RM {

    // static function to transform the mesh
    static void toggleTransformFlag (Resource &r)
    {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);
	    mptr->_transformFlag = !mptr->_transformFlag;
    }

	  // static function to reload GPU programs
	  static void reloadPrograms (Resource &r)
	  {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);
	    mptr->_glReprogramFlag = true;
	  }

    // static function to read OFF style file
    static bool readOFFMeshFile (const string &file, vector <vec> &verts, vector <unsigned int> &indices)
    {
      assert (!file.empty ());

      FILE *fp = fopen (file.c_str (), "r");
      assert (fp);

      char header [32];
      int status = fscanf (fp, "%s\n", header);
      assert (status != 0);

      int nverts = 0, nfaces = 0, nedges = 0;
      status = fscanf (fp, "%d %d %d\n", &nverts, &nfaces, &nedges);
      assert (status != 0);
      assert (nverts > 0 && (nfaces > 0 || nedges > 0));

      // read vertex data
      real tmpr [SF_VECTOR_SIZE] = {0., 0., 0.
#ifdef SF_VECTOR4_ENABLED
        , 1.
#endif
      };
      verts.reserve (static_cast <size_t> (nverts));
      for (int i = 0; i < nverts; ++i){
#ifdef SF_DOUBLE_PRECISION
        status = fscanf (fp, "%lf %lf %lf\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#else
        status = fscanf (fp, "%f %f %f\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#endif
        assert (status != 0);
        verts.push_back (vec (tmpr));
      }

      int dummy;
      int tmpd [3];
      if (nfaces){
        for (int i = 0; i < nfaces; ++i){
          status = fscanf (fp, "%d %d %d %d\n", &dummy, &(tmpd [0]), &(tmpd [1]), &(tmpd [2]));
          assert (status != 0);
          assert (tmpd [0] >= 0 && tmpd [0] < nverts);
          assert (tmpd [1] >= 0 && tmpd [1] < nverts);
          assert (tmpd [2] >= 0 && tmpd [2] < nverts);
          for (int j = 0; j < 3; ++j){
            indices.push_back ( static_cast <unsigned int> (tmpd [j]));
          }
        }
      } else if (nedges){
        for (int i = 0; i < nedges; ++i){
          status = fscanf (fp, "%d %d %d\n", &dummy, &(tmpd [0]), &(tmpd [1]));
          assert (status != 0);
          assert (tmpd [0] >= 0 && tmpd [0] < nverts);
          assert (tmpd [1] >= 0 && tmpd [1] < nverts);
          for (int j = 0; j < 2; ++j){
            indices.push_back ( static_cast <unsigned int> (tmpd [j]));
          }
        }
      }
      fclose (fp);
      return true;
    }

    // static function to draw mesh
    static void plainDraw (Resource &r)
    {
	    Mesh* mptr = dynamic_cast <Mesh *> (&r);

	    mptr->_syncControl [mptr->_semGraphicsWaitIndex].wait ();
#ifndef NDEBUG
      GLenum error;
#endif

#ifndef NDEBUG
        checkGLError (error);
#endif
	    // reload program if needed
	    if (mptr->_glReprogramFlag) {
        mptr->initGPUPrograms ();
        mptr->_glReprogramFlag = false;
	    }

      /************************ RENDERING OF NORMALS HAPPENS HERE ************************/

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
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glDrawElements (GL_TRIANGLES, mptr->_numFaces [0], GL_UNSIGNED_INT, 0);
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

      /************************ RENDERING OF EXTERNAL SURFACE HAPPENS HERE ************************/
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

      if (mptr->_glBufferFlag){
        glBindVertexArray (mptr->_glRenderVertexArrayId [1]);
      } else {
        glBindVertexArray (mptr->_glRenderVertexArrayId [0]);
      }
#ifndef NDEBUG
      checkGLError (error);
#endif

      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, mptr->_glIndexBufferId);
#ifndef NDEBUG
      checkGLError (error);
#endif
      glDrawElements (GL_TRIANGLES, mptr->_numFaces [0], GL_UNSIGNED_INT, 0);
#ifndef NDEBUG
      checkGLError (error);
#endif

      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
      glBindVertexArray (0);
      glBindTexture (GL_TEXTURE_2D, 0);
      glBindTexture (GL_TEXTURE_CUBE_MAP, 0);
      glUseProgram (0);

	    mptr->_syncControl [mptr->_semGraphicsPostIndex].post ();
    }

		// protected constructor and assignment functions
		Mesh::Mesh () { }
		Mesh::Mesh (const Mesh &m) { }
		Mesh& Mesh::operator = (const Mesh &m) { return *this; }

		// only legitimate constructor
		Mesh::Mesh (const string &config, Driver &driver)
		: _transformFlag (false),
		_semPhysicsWaitIndex (-1), _semPhysicsPostIndex (-1),
		_semCollisionWaitIndex (-1), _semCollisionPostIndex (-1),
    _semIntersectionWaitIndex (-1), _semIntersectionPostIndex (-1),
    _semGraphicsWaitIndex (-1), _semGraphicsPostIndex (-1),
    _numSurfaceVertices (0), _curr (&(_vertices [0])), _prev (&(_vertices [1])),
    _bladeCurr (NULL), _bladePrev (NULL), _bladeIndices (NULL),
    _glBufferFlag (false), _glReprogramFlag (false),
    _glNormalFramebufferId (0), _glNormalTexCoordBufferId (0), _glNormalTextureId (0),
    _glEnvTextureId (driver._display.get ()->_glEnvTextureId), _glIndexBufferId (0), _glEnvTextureLocation (-1),
    _glModelviewMatrixLocation (-1), _glProjectionMatrixLocation (-1), _glNormalTextureLocation (-1), _glColorLocation (-1),
    _glLightDirLocation1 (-1), _glLightAmbLocation1 (-1), _glLightDiffLocation1 (-1), _glLightSpecLocation1 (-1), _glLightExpLocation1 (-1),
    _glLightDirLocation2 (-1), _glLightAmbLocation2 (-1), _glLightDiffLocation2 (-1), _glLightSpecLocation2 (-1), _glLightExpLocation2 (-1),
    _glNumLights (driver._display.get ()->_numLights),
    _glLightDir1 (NULL), _glLightAmb1 (NULL), _glLightDiff1 (NULL), _glLightSpec1 (0.), _glLightExp1 (0.),
    _glLightDir2 (NULL), _glLightAmb2 (NULL), _glLightDiff2 (NULL), _glLightSpec2 (0.), _glLightExp2 (0.),
    _glModelview (&(driver._display.get ()->_modelview [0])), _glProjection (&(driver._display.get ()->_projection [0])),
    _glslPrefixString (&(driver._display.get ()->_glslPrefixString))
		{
			assert (!config.empty ());

			_owner = boost::shared_ptr <string> (new string ("Rigid"));

      {
        string name;
				getConfigParameter (config, "name", name);
				assert (!name.empty ());
				_name = boost::shared_ptr <string> (new string (name));

				string file;
				getConfigParameter (config, "data_file", file);
				assert (!file.empty ());

				// initialize faceIndices structure
				_faceIndices.reserve (1);
				_faceIndices.push_back (vector <unsigned int> ());
				readOFFMeshFile (file, _vertices [0], _faceIndices [0]);

				_numSurfaceVertices = _vertices [0].size ();

				_numFaces.reserve (1);
				_numFaces.push_back (_faceIndices [0].size ());

				vec3 min (_vertices [0][0]._v [0], _vertices [0][0]._v [1], _vertices [0][0]._v [2]);
				vec3 max (min);

				for (size_t i = 1; i < _vertices [0].size (); ++i){
          for (int j = 0; j < 3; ++j){
            if (min._v [j] > _vertices [0][i]._v [j]){
              min._v [j] = _vertices [0][i]._v [j];
            } else if (max._v [j] < _vertices [0][i]._v [j]){
              max._v [j] = _vertices [0][i]._v [j];
            }
          }
				}
				_bbox = aabb (min, max);
      }
      _vertices [1] = _vertices [0];

      // add blade-related parameters if body is of
      {
        string type;
        getConfigParameter (config, "type", type);
        if (!type.compare ("cut")){

          string bladeFile;
          getConfigParameter (config, "cut_data", bladeFile);
          assert (!bladeFile.empty ());

          _bladeCurr = new vector <vec> ();
          _bladeIndices = new vector <unsigned int> ();

          readOFFMeshFile (bladeFile, *_bladeCurr, *_bladeIndices);

          _bladePrev = new vector <vec> (*_bladeCurr);
        }
      }

      // add any additional displacement vector to node
      {
        mat4x4 m (1., 0., 0., 0.,
                  0., 0., 1., 0.,
                  0., -1., 0., 0.,
                  0., 0., 0., 1.);
        for (unsigned int i = 0; i < _numSurfaceVertices; ++i){
          _vertices [0] [i] = m * _vertices [0] [i];
        }
        for (unsigned int i = 0; i < _numSurfaceVertices; ++i){
          _vertices [1] [i] = m * _vertices [1] [i];
        }
        for (unsigned int i = 0; i < _bladeCurr->size (); ++i){
          _bladeCurr->at (i) = m * _bladeCurr->at (i);
        }
        for (unsigned int i = 0; i < _bladePrev->size (); ++i){
          _bladePrev->at (i) = m * _bladePrev->at (i);
        }
        string dispStr;
        getConfigParameter (config, "displacement_vector", dispStr);
        if (!dispStr.empty ()){
          size_t first = dispStr.find_first_of (" ", 0);
          size_t last = dispStr.find_last_of (" ", dispStr.size () - 1);

          string x, y, z;
          x.append (dispStr, 0, first);
          y.append (dispStr, first + 1, last - first - 1);
          z.append (dispStr, last + 1, dispStr.size () - last);

          assert (!x.empty () && !y.empty ()  && !z.empty ());
          for (size_t i = 0; i < x.size (); ++i){
            assert (isdigit (x [i]) || x [i] == '.' || x [i] == '-');
          }
          for (size_t i = 0; i < y.size (); ++i){
            assert (isdigit (y [i]) || y [i] == '.' || y [i] == '-');
          }
          for (size_t i = 0; i < z.size (); ++i){
            assert (isdigit (z [i]) || z [i] == '.' || z [i] == '-');
          }

          real disp [3] = {static_cast <real> (atof (x.c_str ())), static_cast <real> (atof (y.c_str ())), static_cast <real> (atof (z.c_str ()))};

          for (size_t i = 0; i < _numSurfaceVertices; ++i){
            for (int j = 0; j < 3; ++j){
              _vertices [0][i]._v [j] += disp [j];
            }
          }

          for (size_t i = 0; i < _numSurfaceVertices; ++i){
            for (int j = 0; j < 3; ++j){
              _vertices [1][i]._v [j] += disp [j];
            }
          }

          for (size_t i = 0; i < _bladeCurr->size (); ++i){
            for (int j = 0; j < 3; ++j){
              _bladeCurr->at (i)._v [j] += disp [j];
            }
          }

          for (size_t i = 0; i < _bladePrev->size (); ++i){
            for (int j = 0; j < 3; ++j){
              _bladePrev->at (i)._v [j] += disp [j];
            }
          }
        }
      }

      /*************************** INITIALIZE THREAD CONTROL PARAMETERS ***************************/
      {
        string mStr;
        getConfigParameter (config, "num_mutexes", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }

        int numMutex = atoi (mStr.c_str ());
        string msp ("mutex_startval"), msv;
        for (int i = 0; i < numMutex; ++i){
          mStr = msp;
          stringstream ss;
          ss << i + 1;
          mStr.append (ss.str ());

          getConfigParameter (config, mStr.c_str (), msv);
          assert (!msv.empty ());
          for (size_t j = 0; j < msv.size (); ++j){
            assert (isdigit (msv [j]));
          }
          _syncControl.push_back (static_cast <unsigned int> (atoi (msv.c_str ())));
        }

        getConfigParameter (config, "physics_wait_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semPhysicsWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "physics_post_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semPhysicsPostIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "collision_wait_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semCollisionWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "collision_post_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semCollisionPostIndex = atoi (mStr.c_str ());

        if (_bladeCurr){
          getConfigParameter (config, "intersection_wait_index", mStr);
          assert (!mStr.empty ());
          for (size_t i = 0; i < mStr.size (); ++i){
            assert (isdigit (mStr [i]));
          }
          _semIntersectionWaitIndex = atoi (mStr.c_str ());

          getConfigParameter (config, "intersection_post_index", mStr);
          assert (!mStr.empty ());
          for (size_t i = 0; i < mStr.size (); ++i){
            assert (isdigit (mStr [i]));
          }
          _semIntersectionPostIndex = atoi (mStr.c_str ());
        }

        getConfigParameter (config, "graphics_wait_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semGraphicsWaitIndex = atoi (mStr.c_str ());

        getConfigParameter (config, "graphics_post_index", mStr);
        assert (!mStr.empty ());
        for (size_t i = 0; i < mStr.size (); ++i){
          assert (isdigit (mStr [i]));
        }
        _semGraphicsPostIndex = atoi (mStr.c_str ());
      }

      /*************************** INITIALIZE OPENGL PARAMETERS ***************************/
      string colorStr;
      getConfigParameter (config, "color", colorStr);
      if (!colorStr.empty ()){
        size_t first = colorStr.find_first_of (" ", 0);
        size_t last = colorStr.find_last_of (" ", colorStr.size () - 1);

        string red, green, blue;
        red.append (colorStr, 0, first);
        green.append (colorStr, first + 1, last - first - 1);
        blue.append (colorStr, last + 1, colorStr.size () - last);

        assert (!red.empty () && !green.empty ()  && !blue.empty ());
        for (size_t i = 0; i < red.size (); ++i){
          assert (isdigit (red [i]) || red [i] == '.');
        }
        for (size_t i = 0; i < green.size (); ++i){
          assert (isdigit (green [i]) || green [i] == '.');
        }
        for (size_t i = 0; i < blue.size (); ++i){
          assert (isdigit (blue [i]) || blue [i] == '.');
        }

        _glColor [0] = static_cast <real> (atof (red.c_str ()));
        _glColor [1] = static_cast <real> (atof (green.c_str ()));
        _glColor [2] = static_cast <real> (atof (blue.c_str ()));

      } else {
        for (int i = 0; i < 3; ++i){
          _glColor [i] = .5;
        }
      }

      for (int i = 0; i < 2; ++i){
        _glNormalFramebufferDimensions [i] = 0;
        _glNormalVertexArrayId [i] = 0;
        _glVertexBufferId [i] = 0;
        _glRenderVertexArrayId [i] = 0;
      }

      initGLBufferObjects ();

      // get names of the shading programs
      getConfigParameter (config, "normal_shader", _glProgramName [0]);
      assert (!_glProgramName [0].empty ());

      getConfigParameter (config, "color_shader", _glProgramName [1]);
      assert (!_glProgramName [1].empty ());

      for (int i = 0; i < 2; ++i){
        _glProgram [i] = 0;
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
      draw = &plainDraw;
      reprogram = &reloadPrograms;
      transform = &toggleTransformFlag;
		}

    // mesh's run method
    void
    Mesh::run ()
    {
      vector <vec> *tmp;

      while (true) {

        // wait for physics end to get control
        _syncControl [_semPhysicsWaitIndex].wait ();

        // swap vertex buffers
        tmp = _curr;
        _curr = _prev;
        _prev = tmp;

        if (_bladeCurr){
          tmp = _bladeCurr;
          _bladeCurr = _bladePrev;
          _bladePrev = tmp;
        }

        // do useful stuff
        if (_transformFlag){
          move ();
        }

        _glBufferFlag = !_glBufferFlag; // should be the last line in this segment

        // release data
        _syncControl [_semPhysicsPostIndex].post ();
      }
    }

    void
    Mesh::move ()
    {
      vec moveVec (-.02, 0., 0.);
      for (unsigned int i = 0; i < _curr->size (); ++i){
        _curr->at (i) = _prev->at (i) + moveVec;
      }
      if (_bladeCurr){
        for (unsigned int i = 0; i < _bladeCurr->size (); ++i){
          _bladeCurr->at (i) = _bladePrev->at (i) + moveVec;
        }
      }
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

      /*************************** INITIALIZE SURFACE RENDERING PROGRAM ***************************/
      if (!initGPUProgram (false, *_glslPrefixString, _glProgramName [1], _glProgram [1])){
        PRINT ("error: could not initialize %s\n", _glProgramName [1].c_str ());
        return false;
      }

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
      _glColorLocation = glGetUniformLocation (_glProgram [1], "color");

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

      glGenVertexArrays (2, _glRenderVertexArrayId);
      checkGLError (error);

      for (size_t i = 0; i < 2; ++i){

        glBindVertexArray (_glRenderVertexArrayId [i]);
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

      return true;
    }

    // private method to initialize texture-related objects
    bool
    Mesh::initGLBufferObjects ()
    {
      GLenum error;

      glGenBuffers (2, _glVertexBufferId);
      checkGLError (error);

      glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [0]);
      checkGLError (error);
      glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*_numSurfaceVertices, &(_vertices [0][0]), GL_DYNAMIC_DRAW);
      checkGLError (error);

		  glBindBuffer (GL_ARRAY_BUFFER, _glVertexBufferId [1]);
      checkGLError (error);
		  glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*_vertices [1].size (), &(_vertices [1][0]), GL_DYNAMIC_DRAW);
      checkGLError (error);
      glBindBuffer (GL_ARRAY_BUFFER, 0);

      glGenBuffers (1, &(_glIndexBufferId));
      checkGLError (error);
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, _glIndexBufferId);
      checkGLError (error);
      glBufferData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int)*_numFaces [0], &(_faceIndices [0][0]), GL_DYNAMIC_DRAW);
      checkGLError (error);
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
      texcoords.reserve (_numSurfaceVertices);
      texcoords.resize (_numSurfaceVertices, vec2::ZERO);

      unsigned int xcount = 0, ycount = 0;
      for (size_t i = 0; i < _numSurfaceVertices; ++i){
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
      for (size_t i = 0; i < _numSurfaceVertices; ++i){
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

      return true;
    }
  }
}
