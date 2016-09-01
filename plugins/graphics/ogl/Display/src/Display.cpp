/**
 * @file GL_Window.cpp
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
 * The GL_Window class for Simulate framework. This contains all the OpenGL
 * and GLUT specific routines used for rendering.
 */

#include <cmath>

extern "C" {
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libpng12/png.h>
#include <libpng12/pngconf.h>
}

#include <vector>
#include <string>
#include <sstream>

#include "Preprocess.h"
#include "vec3.h"
#include "aabb.h"

#include "Driver.h"
#include "Resource.h"
#include "GL/common.h"
#include "Display.h"
#include "GL_Display.h"

using namespace std;

namespace SF {

	static const GLfloat OGL_HITHER = 1.;
	static const GLfloat OGL_YON = 1500.;
	static const GLfloat OGL_FOV = 45.;

	// local variables (to speed up computations)
	static bool _wireFlag = false;
	static bool _fpsShowFlag = false;
	static bool _helpShowFlag = false;

#ifdef SF_DOUBLE_PRECISION
	GLdouble _result [16];
#else
	GLfloat _result [16];
#endif

	static real _deltaX, _deltaY;
	static real _mag, _cosVal, _sinVal;
	static real _tr00, _tr01, _tr02, _tr03, _tr11, _tr12, _tr13, _tr22, _tr23;

	// copy of the display class displayed
	static GL_Window *_disp = NULL;

  // function to change endian-ness of float value
	inline GLfloat changeEndian (GLfloat &num)
	{
		uint32_t mem = *(uint32_t *) (&num);
		uint32_t result = (mem >> 24) | ((mem >> 8) & 0x0000ff00) | ((mem << 8) & 0x00ff0000) | (mem << 24);
		return *(GLfloat *) (&result);
	}

	// static function to read ppm-format cube-map image files
	static int readCubeMapFile (const bool swapFlag, const string &imgFile, vector <GLfloat> &rgb)
	{
		assert (rgb.empty ());

		FILE *fp = fopen (imgFile.c_str (), "rb");
		if (!fp){
			PRINT ("GL error: could not open %s\n", imgFile.c_str ());
			exit (EXIT_FAILURE);
		}

		char dummy [8];
		int d1 = 0, d2 = 0;
		float d3 = 0.;
		int status = fscanf (fp, "%s\n", &(dummy [0]));
		assert (status);
		status = fscanf (fp, "%d\n", &d1);
		assert (status);
		status = fscanf (fp, "%d\n", &d2);
		assert (status);
		status = fscanf (fp, "%f\n", &d3);
		assert (status);
		assert (status);
		assert (d1 == d2);

		rgb.push_back (0.);
		rgb.resize (3*d1*d1);
		unsigned int nbytes = fread (&(rgb [0]), sizeof (GLfloat), 3*d1*d1, fp);
		if (static_cast <int> (nbytes) != 3*d1*d1){
      PRINT ("Warning: Number of bytes read (%u) not same as dimensions specified (%d) in %s\n", nbytes, 3*d1*d1, imgFile.c_str ());
		}

		fclose (fp);

		if (swapFlag){

			GLfloat val;
			for (unsigned int i = 0; i < rgb.size (); ++i){
				val = changeEndian (rgb.at (i));
				rgb.at (i) = val;
			}
		}

		return d1;
	}

	// function to read framebuffer and output to sample.png file on Desktop
	static void printScreen ()
	{
		// wait for all previous GL commands to be finished
		glFinish ();

		// read framebuffer
		vector <GLubyte> rgb (3 * _disp->_windowWidth * _disp->_windowHeight);
		assert (rgb.size () == static_cast <size_t> (3 * _disp->_windowWidth * _disp->_windowHeight));

		GLenum error;

		glReadBuffer (0);
		checkGLError (error);

		glReadPixels (0, 0, _disp->_windowWidth, _disp->_windowHeight, GL_RGB, GL_UNSIGNED_BYTE, &(rgb [0]));
		checkGLError (error);

		glFinish ();

		// write out to PNG file
		FILE *fp = fopen ("~/Desktop/sample.png", "wb");
		assert (fp);

		// initialize png-specific functions
		png_structp png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if (!png_ptr){
			PRINT ("GL error: could not allocate png_ptr\n");
			fclose (fp);
			return;
		}

		png_infop info_ptr = png_create_info_struct (png_ptr);
		if (!info_ptr){
			PRINT ("GL error: could not allocate png info structure\n");
			png_destroy_write_struct (&png_ptr, static_cast <png_infopp> (NULL));
			fclose (fp);
			return;
		}

		if (setjmp (png_jmpbuf (png_ptr))){
			PRINT ("GL error: could not write png header\n");
			png_destroy_write_struct (&png_ptr, &info_ptr);
			fclose (fp);
			return;
		}

		png_init_io (png_ptr, fp);

		png_set_filter (png_ptr, 0, PNG_FILTER_NONE);
		png_bytep *row_headers = new png_bytep [_disp->_windowHeight];

		png_set_IHDR (png_ptr, info_ptr, _disp->_windowWidth, _disp->_windowHeight, 8,
				PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

		for (int i = 0; i < _disp->_windowHeight; ++i){
			row_headers [_disp->_windowHeight - i - 1] = (png_bytep)(&(rgb [3 * _disp->_windowWidth * i]));
		}

		png_write_info (png_ptr, info_ptr);
		png_write_image (png_ptr, row_headers);

		delete [] row_headers;

		png_write_end (png_ptr, NULL);

		// clean up
		png_destroy_write_struct (&png_ptr, &info_ptr);

		fclose (fp);
	}

	// input config-file parsing function
	static bool getProperty (const string &cfgfile, const char *nodeName, const char *property, string &result)
	{
		if (!result.empty ()){
			result.clear ();
		}

		xmlDocPtr doc = xmlReadFile (cfgfile.c_str (), NULL, 0);
		if (!doc){
			PRINT ("error opening file %s\n", cfgfile.c_str ());
			return false;
		}

		xmlNodePtr node = xmlDocGetRootElement (doc);
		if (strcmp (reinterpret_cast <const char *> (node->name), "SFGLConfig")){
			PRINT ("error: root element in %s not of \'SFGLConfig\' type", node->name);
			return false;
		}

		node = node->children;
		node = node->next;

		while (node){
			if (!strcmp (reinterpret_cast <const char *> (node->name), nodeName )){
				char *rstr = reinterpret_cast <char *> (xmlGetProp (node, reinterpret_cast <const xmlChar *> (property)));
				result = string (rstr);
				free (rstr);
				rstr = NULL;
				break;
			}
			node = node->next;
			node = node->next;
		}

		// clean up
		xmlFreeDoc (doc);
		xmlCleanupParser ();

		if (result.empty ()){
			return false;
		}

		return true;
	}

	// function for keyboard binding
	static void keys (unsigned char k, int a, int b)
	{
		switch (k){

			// reload programs for objects by numbers
			case '1':
				if (!_disp->_drawables.empty ()){
					_disp->_drawables [0]->reprogram (* (_disp->_drawables [0].get ()));
				}
				break;
			case '2':
				if (_disp->_drawables.size () > 1){
					_disp->_drawables [1]->reprogram (* (_disp->_drawables [1].get ()));
				}
				break;
			case '3':
				if (_disp->_drawables.size () > 2){
					_disp->_drawables [2]->reprogram (* (_disp->_drawables [2].get ()));
				}
				break;
			case '4':
				if (_disp->_drawables.size () > 3){
					_disp->_drawables [3]->reprogram (* (_disp->_drawables [3].get ()));
				}
				break;
			case '5':
				if (_disp->_drawables.size () > 4){
					_disp->_drawables [4]->reprogram (* (_disp->_drawables [4].get ()));
				}
				break;

			// transformation flag
			case 'm': case 'M':
				if (_disp->_moveToggleCounter < _disp->_moveables.size ()){
					_disp->_moveables [_disp->_moveToggleCounter]->transform (* (_disp->_moveables [_disp->_moveToggleCounter].get ()));
				}
				break;

			// toggle transformation flags
			case 't': case 'T':
				if (++_disp->_moveToggleCounter >= _disp->_moveables.size ()){
					_disp->_moveToggleCounter = 0;
				}
				break;

			// wire-frame rendering toggle key
			case 'w': case 'W':
				if (!_wireFlag){
					_wireFlag = true;
					glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
				}
				else {
					_wireFlag = false;
					glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
				}
				break;

			// print screen key
			case 'p': case 'P':
				printScreen ();
				break;

			// display help
			case 'h': case 'H':
				_helpShowFlag = !_helpShowFlag;
				break;

			// display frames per second
			case 'f': case 'F':
				_fpsShowFlag = !_fpsShowFlag;
				break;

			// quit
			case 'q': case 'Q': case 27:
        _disp->_parent->cleanup ();
				exit (EXIT_SUCCESS);
		}
	}

	// display function
	static void display ()
	{
		// clear from previous frame
		glClearColor (_disp->_background [0], _disp->_background [1], _disp->_background [2], 0.);
		glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (_fpsShowFlag){
			frameStart ();
		}

		for (unsigned int i = 0; i < _disp->_drawables.size (); ++i){
			_disp->_drawables [i]->draw (* (_disp->_drawables [i].get ()));
		}

		if (_fpsShowFlag){
			frameEnd (GLUT_BITMAP_HELVETICA_12, 1., 0., 0., .89, .035);
		}
		if (_helpShowFlag){
			displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, .95, "KEY BINDINGS:");
			displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, .91, "f: GL_Window FramesPerSec");
			displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, .87, "w: GL_Window wireframe");
			displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, .83, "p: Print window (to sample.png)");
			displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, .79, "Shader reload keys:");
			float pos = .79;
			string _dispStr;
			for (unsigned int i = 0; i < _disp->_drawables.size (); ++i){

			  pos -= .04;

				stringstream ss;
				ss << i + 1;
				_dispStr = string (ss.str ());
				_dispStr.append (": ");
				_dispStr.append (* (_disp->_drawables [i].get ()->_name.get ()));

				displayOnScreen (GLUT_BITMAP_HELVETICA_18, .25, 0., 0., .02, pos, _dispStr.c_str ());
			}
		}

    // push frame-buffers to display queue: last line of function
    glutSwapBuffers ();
	}

	// idle function
	static void idle ()
	{
		display ();
	}

	// window resize function
	static void resize (int w, int h)
	{
		_disp->_windowWidth = w;
		_disp->_windowHeight = h;

		glViewport (0, 0, w, h);
		_disp->updateProjection ();

		char wStr [16];
		sprintf (wStr, "%d x %d", w, h);

		string title ("RunSim GL_Window System - ");
		title.append (wStr);

		glutSetWindowTitle (title.c_str ());

		display ();
	}

	// mouse function
	static void mouse (int b, int s, int x, int y)
	{
		_disp->_mouseButton = b;
		_disp->_mouseX = x;
		_disp->_mouseY = y;
	}

	// mouse motion function
	static void motion (int x, int y)
	{
		// record change
		_deltaX = x - _disp->_mouseX;
		_deltaY = y - _disp->_mouseY;

		// update position
		_disp->_mouseX = x;
		_disp->_mouseY = y;


		if (_disp->_mouseButton == GLUT_LEFT_BUTTON){ // rotate

			_deltaX *= 0.01; // to damp motion
			_deltaY *= 0.01; // to damp motion

			_mag = static_cast <real> (sqrt (_deltaX*_deltaX + _deltaY*_deltaY));
			_deltaX /= _mag;
			_deltaY /= _mag;
			_cosVal = static_cast <real> (cos (180.* _mag/ _disp->_windowWidth));
			_sinVal = static_cast <real> (sqrt (1. - _cosVal*_cosVal));

			_tr00 = _deltaY*_deltaY*(1. - _cosVal) + _cosVal;
			_tr01 = _deltaY*_deltaX*(1. - _cosVal);
			_tr02 = _deltaX*_sinVal;
			_tr11 = _deltaX*_deltaX*(1. - _cosVal) + _cosVal;
			_tr12 = -_deltaY*_sinVal;
			_tr22 = _cosVal;

			_result [0] = _tr00 * _disp->_modelview [0] + _tr01 * _disp->_modelview [1] + _tr02 * _disp->_modelview [2];
			_result [1] = _tr01 * _disp->_modelview [0] + _tr11 * _disp->_modelview [1] + _tr12 * _disp->_modelview [2];
			_result [2] = -_tr02 * _disp->_modelview [0] - _tr12 * _disp->_modelview [1] + _tr22 * _disp->_modelview [2];
			_result [3] = _disp->_modelview [3];

			_result [4] = _tr00 * _disp->_modelview [4] + _tr01 * _disp->_modelview [5] + _tr02 * _disp->_modelview [6];
			_result [5] = _tr01 * _disp->_modelview [4] + _tr11 * _disp->_modelview [5] + _tr12 * _disp->_modelview [6];
			_result [6] = -_tr02 * _disp->_modelview [4] - _tr12 * _disp->_modelview [5] + _tr22 * _disp->_modelview [6];
			_result [7] = _disp->_modelview [7];

			_result [8] = _tr00 * _disp->_modelview [8] + _tr01 * _disp->_modelview [9] + _tr02 * _disp->_modelview [10];
			_result [9] = _tr01 * _disp->_modelview [8] + _tr11 * _disp->_modelview [9] + _tr12 * _disp->_modelview [10];
			_result [10] = -_tr02 * _disp->_modelview [8] - _tr12 * _disp->_modelview [9] + _tr22 * _disp->_modelview [10];
			_result [11] = _disp->_modelview [11];

			_result [12] = _tr00 * _disp->_modelview [12] + _tr01 * _disp->_modelview [13] + _tr02 * _disp->_modelview [14];
			_result [13] = _tr01 * _disp->_modelview [12] + _tr11 * _disp->_modelview [13] + _tr12 * _disp->_modelview [14];
			_result [14] = -_tr02 * _disp->_modelview [12] - _tr12 * _disp->_modelview [13] + _tr22 * _disp->_modelview [14];
			_result [15] = _disp->_modelview [15];

			memcpy (_disp->_modelview, _result, 16*sizeof (real));
		}
		else if (_disp->_mouseButton == GLUT_RIGHT_BUTTON){ // zoom

			_tr23 = -_disp->_cameraScales [2] * _deltaY;

			_result [0] = _disp->_modelview [2] + _tr23 * _disp->_modelview [3];
			_result [1] = _disp->_modelview [6] + _tr23 * _disp->_modelview [7];
			_result [2] = _disp->_modelview [10] + _tr23 * _disp->_modelview [11];
			_result [3] = _disp->_modelview [14] + _tr23 * _disp->_modelview [15];

			for (int i = 0; i < 4; ++i){
				_disp->_modelview [4*i + 2] = _result [i];
			}
			for (int i = 0; i < 16; ++i){
				_disp->_modelview [i] /= _disp->_modelview [15];
			}
		}
		else if (_disp->_mouseButton == GLUT_MIDDLE_BUTTON){ // pan

			_tr03 = -_disp->_cameraScales [0] * _deltaX;
			_tr13 = -_disp->_cameraScales [1] * _deltaY;

			_result [0] = _disp->_modelview [0] + _tr03 * _disp->_modelview [3];
			_result [1] = _disp->_modelview [4] + _tr03 * _disp->_modelview [7];
			_result [2] = _disp->_modelview [8] + _tr03 * _disp->_modelview [11];
			_result [3] = _disp->_modelview [12] + _tr03 * _disp->_modelview [15];
			_result [4] = _disp->_modelview [1] + _tr13 * _disp->_modelview [3];
			_result [5] = _disp->_modelview [5] + _tr13 * _disp->_modelview [7];
			_result [6] = _disp->_modelview [9] + _tr13 * _disp->_modelview [11];
			_result [7] = _disp->_modelview [13] + _tr13 * _disp->_modelview [15];

			for (int i = 0; i < 4; ++i){
				_disp->_modelview [4*i] = _result [i];
			}
			for (int i = 0; i < 4; ++i){
				_disp->_modelview [4*i + 1] = _result [4 + i];
			}
			for (int i = 0; i < 16; ++i){
				_disp->_modelview [i] /= _disp->_modelview [15];
			}
		}
	}

	// overloaded constructor
	GL_Window::GL_Window (int &argc, char **argv, const string &config)
	: _parent (NULL), _moveToggleCounter (0),
	  _windowX (0), _windowY (0), _windowWidth (0), _windowHeight (0),
	  _mouseX (0), _mouseY (0), _mouseButton (0),
	  _numLights (0), _lightSpec1 (0.), _lightExp1 (0.), _lightSpec2 (0.), _lightExp2 (0.)
	{
		for (int i = 0; i < 3; ++i){
			_background [i] = 1.;
			_cameraPosition [i] = 0.;
			_cameraScales [i] = 1.;

			_lightDir1 [i] = 0.;
			_lightAmb1 [i] = 0.;
			_lightDiff1 [i] = 0.;

			_lightDir2 [i] = 0.;
			_lightAmb2 [i] = 0.;
			_lightDiff2 [i] = 0.;
		}
		for(int i = 0; i < 16; ++i){
			_projection [i] = 0.;
		}
		for(int i = 0; i < 16; ++i){
			_modelview [i] = 0.;
		}
		_projection [0] = _projection [5] = _projection [10] = _projection [15] = 1.;
		_modelview [0] = _modelview [5] = _modelview [10] = _modelview [15] = 1.;

		// get parameters from input configuration file
		string inStr;
		if (getProperty (config, "dimensions", "width", inStr)){
			_windowWidth = atoi (inStr.c_str ());
		}
		assert (_windowWidth);

		if (getProperty (config, "dimensions", "height", inStr)){
			_windowHeight = atoi (inStr.c_str ());
		}
		assert (_windowHeight);

		if (getProperty (config, "background", "color", inStr)){
			unsigned int first = inStr.find_first_of (" ", 0);
			unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

			string red, green, blue;
			red.append (inStr, 0, first);
			green.append (inStr, first + 1, last - first - 1);
			blue.append (inStr, last + 1, inStr.size () - last);

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

			_background [0] = static_cast <real> (atof (red.c_str ()));
			_background [1] = static_cast <real> (atof (green.c_str ()));
			_background [2] = static_cast <real> (atof (blue.c_str ()));
		}

		if (getProperty (config, "camera", "position", inStr)){
			unsigned int first = inStr.find_first_of (" ", 0);
			unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

			string x, y, z;
			x.append (inStr, 0, first);
			y.append (inStr, first + 1, last - first - 1);
			z.append (inStr, last + 1, inStr.size () - last);

      assert (!x.empty () && !y.empty ()  && !z.empty ());
      for (unsigned int i = 0; i < x.size (); ++i){
        assert (isdigit (x [i]) || x [i] == '.');
      }
      for (unsigned int i = 0; i < y.size (); ++i){
        assert (isdigit (y [i]) || y [i] == '.');
      }
      for (unsigned int i = 0; i < z.size (); ++i){
        assert (isdigit (z [i]) || z [i] == '.');
      }

			_cameraPosition [0] = static_cast <real> (atof (x.c_str ()));
			_cameraPosition [1] = static_cast <real> (atof (y.c_str ()));
			_cameraPosition [2] = static_cast <real> (atof (z.c_str ()));
		}

    if (getProperty (config, "light0", "specular", inStr)){

      ++_numLights;

      for (unsigned int i = 0; i < inStr.length (); ++i){
        assert (isdigit (inStr [i]) || inStr [i] == '.');
      }
      _lightSpec1 = static_cast <real> (atof (inStr.c_str ()));

      if (getProperty (config, "light0", "exp", inStr)){
        for (unsigned int i = 0; i < inStr.length (); ++i){
          assert (isdigit (inStr [i]) || inStr [i] == '.');
        }
        _lightExp1 = static_cast <real> (atof (inStr.c_str ()));
      }

      if (getProperty (config, "light0", "direction", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.' || x [i] == '-');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.' || y [i] == '-');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.' || z [i] == '-');
        }

        _lightDir1 [0] = static_cast <real> (atof (x.c_str ()));
        _lightDir1 [1] = static_cast <real> (atof (y.c_str ()));
        _lightDir1 [2] = static_cast <real> (atof (z.c_str ()));
      }

      if (getProperty (config, "light0", "ambient", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.');
        }

        _lightAmb1 [0] = static_cast <real> (atof (x.c_str ()));
        _lightAmb1 [1] = static_cast <real> (atof (y.c_str ()));
        _lightAmb1 [2] = static_cast <real> (atof (z.c_str ()));
      }

      if (getProperty (config, "light0", "diffuse", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.');
        }

        _lightDiff1 [0] = static_cast <real> (atof (x.c_str ()));
        _lightDiff1 [1] = static_cast <real> (atof (y.c_str ()));
        _lightDiff1 [2] = static_cast <real> (atof (z.c_str ()));
      }
    }

    if (getProperty (config, "light1", "specular", inStr)){

      ++_numLights;

      for (unsigned int i = 0; i < inStr.length (); ++i){
        assert (isdigit (inStr [i]) || inStr [i] == '.');
      }
      _lightSpec2 = static_cast <real> (atof (inStr.c_str ()));

      if (getProperty (config, "light1", "exp", inStr)){
        for (unsigned int i = 0; i < inStr.length (); ++i){
          assert (isdigit (inStr [i]) || inStr [i] == '.');
        }
        _lightExp2 = static_cast <real> (atof (inStr.c_str ()));
      }

      if (getProperty (config, "light1", "direction", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.' || x [i] == '-');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.' || y [i] == '-');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.' || z [i] == '-');
        }

        _lightDir2 [0] = static_cast <real> (atof (x.c_str ()));
        _lightDir2 [1] = static_cast <real> (atof (y.c_str ()));
        _lightDir2 [2] = static_cast <real> (atof (z.c_str ()));
      }

      if (getProperty (config, "light1", "ambient", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.');
        }

        _lightAmb2 [0] = static_cast <real> (atof (x.c_str ()));
        _lightAmb2 [1] = static_cast <real> (atof (y.c_str ()));
        _lightAmb1 [2] = static_cast <real> (atof (z.c_str ()));
      }

      if (getProperty (config, "light1", "diffuse", inStr)){
        unsigned int first = inStr.find_first_of (" ", 0);
        unsigned int last = inStr.find_last_of (" ", inStr.size () - 1);

        string x, y, z;
        x.append (inStr, 0, first);
        y.append (inStr, first + 1, last - first - 1);
        z.append (inStr, last + 1, inStr.size () - last);

        assert (!x.empty () && !y.empty ()  && !z.empty ());
        for (unsigned int i = 0; i < x.size (); ++i){
          assert (isdigit (x [i]) || x [i] == '.');
        }
        for (unsigned int i = 0; i < y.size (); ++i){
          assert (isdigit (y [i]) || y [i] == '.');
        }
        for (unsigned int i = 0; i < z.size (); ++i){
          assert (isdigit (z [i]) || z [i] == '.');
        }

        _lightDiff2 [0] = static_cast <real> (atof (x.c_str ()));
        _lightDiff2 [1] = static_cast <real> (atof (y.c_str ()));
        _lightDiff2 [2] = static_cast <real> (atof (z.c_str ()));
      }
    }

		// initialize rendering context
    GLenum error;

		glutInit (&argc, argv);
		checkGLError (error);
		glutInitWindowPosition (_windowX, _windowY);
		checkGLError (error);
		glutInitWindowSize (_windowWidth, _windowHeight);
		checkGLError (error);
		glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
		checkGLError (error);
		glutCreateWindow ("Simulate OpenGL GL_Window");
		checkGLError (error);

		// enable GL states
		glShadeModel (GL_SMOOTH);
		glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

		glEnable (GL_DEPTH_TEST);
		glDepthFunc (GL_LESS);

		glEnable (GL_NORMALIZE);
		glEnable (GL_POLYGON_SMOOTH);

		glEnable (GL_CULL_FACE);
		glCullFace (GL_BACK);

		// store the OpenGL version string
		{
			const GLubyte *version = glGetString (GL_VERSION);
			PRINT ("GL version: %s\n", version);

			char vch [4];
			vch [0] = static_cast <char> (version [0]);
			vch [1] = static_cast <char> (version [2]);
			vch [2] = '0';
			vch [3] ='\0';

			int versionNum = atoi (vch);
			if (versionNum < 330){
        PRINT ("Warning: Current GL version too old to support Ashwini");
			}

			_glslPrefixString = string ("#version ");
			_glslPrefixString.append (vch);
			_glslPrefixString.append ("\n");
		}

		// optionally create environment cube-map
		if (getProperty (config, "environmentmap", "home", inStr)){

			glGenTextures (1, &_glEnvTextureId);
			checkGLError (error);

			glBindTexture (GL_TEXTURE_CUBE_MAP, _glEnvTextureId);
			checkGLError (error);

			glTexParameterf (GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameterf (GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameterf (GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameterf (GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP);
			glTexParameterf (GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP);

			bool endianSwap = false;
			string tmpStr;
			if (getProperty (config, "environmentmap", "endian_swap", tmpStr)){
				if (!tmpStr.compare ("yes")){
					endianSwap = true;
				}
			}

			// read +X envmap file
			{
				getProperty (config, "environmentmap", "posx", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}
			// read -X envmap file
			{
				getProperty (config, "environmentmap", "negx", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}

			// read +Y envmap file
			{
				getProperty (config, "environmentmap", "posy", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}
			// read -Y envmap file
			{
				getProperty (config, "environmentmap", "negy", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}

			// read +Z envmap file
			{
				getProperty (config, "environmentmap", "posz", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}
			// read -Z envmap file
			{
				getProperty (config, "environmentmap", "negz", tmpStr);
				string filename (inStr);
				filename.append (tmpStr);
				vector <GLfloat> rgb;
				int dim = readCubeMapFile (endianSwap, filename, rgb);
				glTexImage2D (GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGB32F, dim, dim, 0, GL_RGB, GL_FLOAT, &(rgb [0]));
				checkGLError (error);
			}
		} // end - if (getProperty (config, "environmentmap", "home", inStr))

		// GLUT function binding
		glutDisplayFunc( display );
		checkGLError (error);
		glutIdleFunc( idle );
		checkGLError (error);
		glutReshapeFunc( resize );
		checkGLError (error);
		glutKeyboardFunc( keys );
		checkGLError (error);
		glutMouseFunc( mouse );
		checkGLError (error);
		glutMotionFunc( motion );
		checkGLError (error);
	}

	// destructor
	GL_Window::~GL_Window ()
	{
    PRINT ("here\n");
	}

	// run function
	void
	GL_Window::run ()
	{
		updateModelview ();
		updateProjection ();

		_disp = const_cast <SF::GL_Window *> (this);

		glutMainLoop ();
	}

	// method to add drawable resources
	void
	GL_Window::addDrawables (const boost::shared_ptr <Resource> &r)
	{
		_drawables.push_back (r);
	}

	// method to add moveable resources
	void
	GL_Window::addMoveables (const boost::shared_ptr <Resource> &r)
	{
		_moveables.push_back (r);
	}

	// method to update projection matrix
	void
	GL_Window::updateProjection ()
	{
		real top = OGL_HITHER * static_cast <real> (tan (OGL_FOV * M_PI/ 360.));
		real bottom = -top;
		real aspect = static_cast <real> (_windowWidth)/ _windowHeight;
		real left = bottom * aspect;
		real right = top * aspect;

		_projection [0] = 2 * OGL_HITHER / (right - left);
		_projection [2] =  (right + left)/ (right - left);
		_projection [5] = 2 * OGL_HITHER / (top - bottom);
		_projection [6] =  (top + bottom)/ (top - bottom);
		_projection [10] = - (OGL_YON + OGL_HITHER)/ (OGL_YON - OGL_HITHER);
		_projection [11] = -2 * OGL_YON * OGL_HITHER/ (OGL_YON - OGL_HITHER);
		_projection [14] = -1.0;
	}

	// method to update modelview matrix
	void
	GL_Window::updateModelview ()
	{
		// calculate the look-at position
		vec3 at ((_bbox._v [0]._v [0] + _bbox._v [1]._v [0])* 0.5, (_bbox._v [0]._v [1] + _bbox._v [1]._v [1])* 0.5, _bbox._v [0]._v [2]);

		// get the maximum XY extent of the view volume
		real maxExtent = 0.0;
		if (_bbox._v [1]._v [0] - _bbox._v [0]._v [0] > _bbox._v [1]._v [1] - _bbox._v [0]._v [1]){
			maxExtent = _bbox._v [1]._v [0] - _bbox._v [0]._v [0];
		} else {
			maxExtent = _bbox._v [1]._v [1] - _bbox._v [0]._v [1];
		}
		maxExtent *= 0.5; // we use half of the max-extent for our calculations

		// set camera position such that it forms a 60deg angle with
		vec3 camera (at);
		real distance = maxExtent/ static_cast <real> (tan (M_PI * OGL_FOV / 360.));
		camera._v [2] += distance;

		// form view, up and right vectors
		vec3 view (camera - at);
		view.normalize ();
		vec3 right = vec3 (0., 1., 0.).ncross (view);
		vec3 up = view.cross (right);

		// assign model view matrix factors
		_modelview [0] = right._v [0];
		_modelview [4] = right._v [1];
		_modelview [8] = right._v [2];

		_modelview [1] = up._v [0];
		_modelview [5] = up._v [1];
		_modelview [9] = up._v [2];

		_modelview [2] = view._v [0];
		_modelview [6] = view._v [1];
		_modelview [10] = view._v [2];

		_modelview [12] = -right.dot (camera);
		_modelview [13] = -up.dot (camera);
		_modelview [14] = -view.dot (camera);

		_cameraScales [0] = 0.01*_modelview [12];
		_cameraScales [1] = 0.01*_modelview [13];
		_cameraScales [2] = 0.01*_modelview [14];
	}
}
