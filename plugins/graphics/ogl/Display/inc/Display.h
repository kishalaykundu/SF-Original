/**
 * @file Display.h
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

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
}

using namespace std;
using namespace boost;

namespace SF {

	class aabb;
	class Driver;
	class Resource;

	class GL_Window {

	public:

    // pointer to the parent driver
    Driver *_parent;

		// resources participating in display
		vector <boost::shared_ptr <Resource> > _drawables;

		// resources participating in keyboard-driven transformations
		unsigned int _moveToggleCounter;
		vector <boost::shared_ptr <Resource> > _moveables;

		// view volume
		aabb _bbox;

		// background color
		real _background [3];

		// column-major transform matrices
		real _projection [16], _modelview [16];

		// OGL-window attributes
		int _windowX, _windowY;
		int _windowWidth, _windowHeight;

		// mouse attributes
		int _mouseX, _mouseY, _mouseButton;

		// camera attributes
		real _cameraPosition[3];
		real _cameraScales[3];

    // global lighting info (max 2)
    int _numLights;

    // light one information
    real _lightDir1 [3];
    real _lightAmb1 [3];
    real _lightDiff1 [3];
    real _lightSpec1;
    real _lightExp1;

    // light two information
    real _lightDir2 [3];
    real _lightAmb2 [3];
    real _lightDiff2 [3];
    real _lightSpec2;
    real _lightExp2;

		// shader prefix string
		string _glslPrefixString;

		// global environment map
		GLuint _glEnvTextureId;

	public:
		GL_Window (int &argc, char **argv, const string &config);
		~GL_Window ();

	private:
		GL_Window () { }
		GL_Window (const GL_Window &d) { }
		GL_Window& operator = (const GL_Window &d) { return *this; }

	public:
		void run ();
		void addDrawables (const boost::shared_ptr <Resource> &r);
		void addMoveables (const boost::shared_ptr <Resource> &r);

		void updateProjection ();
		void updateModelview ();
	};
}
