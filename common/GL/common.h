/**
 * @file common.h
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
 * Common functions used by OpenGL based functions.
 */
#pragma once

#include <cstdio>

#include <string>

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
}

using namespace std;

namespace SF {

  // function to check OpenGL error and print string format
	inline void checkGLErrorPrivate (GLenum &error, const string &file, unsigned int line)
	{
		if ((error = glGetError ()) != GL_NO_ERROR){
			fprintf ( stdout, "%s[%u]:\tGL error: %s\n", basename (file.c_str ()), line, gluErrorString (error));
		}
	}

	#define checkGLError(error) (checkGLErrorPrivate (error, __FILE__, __LINE__))

  // function to initialize a GL shader
  bool initGLShader (const string &header, const string &file, GLuint &shaderId, const GLenum shaderType);

  // function to initialize GPU program
  bool initGPUProgram (const bool geometryFlag, const string &header, const string &file, GLuint &programId);
}
