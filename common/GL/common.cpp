/**
 * @file common.cpp
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

#include <cassert>

#include <string>
#include <iostream>
#include <fstream>

#include "Preprocess.h"

extern "C" {
#if defined( __APPLE__ ) || defined( MACOSX )
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#include <OpenGL/glx.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glx.h>
#endif
}

#include "common.h"

namespace SF {

  // function to initialize a GL shader
  bool initGLShader (const string &header, const string &file, GLuint &shaderId, const GLenum shaderType)
  {
    GLenum error;

    shaderId = glCreateShader (shaderType);
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    ifstream is;
    is.open (file.c_str (), ios::binary);
    assert (!is.fail ());
    is.seekg (0, ios::end);
    size_t length = is.tellg ();
    is.seekg (0, ios::beg);
    char *src = new char [length + 1];
    is.read (src, length);
    is.close ();
    src [length] = '\0';

    string code (header);
    code.append (src);
    delete [] src;

    const char *source = code.c_str ();

    glShaderSource (shaderId, 1, &source, NULL);
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    glCompileShader (shaderId);
    checkGLError (error);

    GLint compileStatus;
    glGetShaderiv (shaderId, GL_COMPILE_STATUS, &compileStatus);
    if (error != GL_NO_ERROR || compileStatus == GL_FALSE){

      char log [1024];
      GLsizei loglength;
      glGetShaderInfoLog (shaderId, 1024, &loglength, log);
      PRINT ("Shader Info Log [%s]\n%s\n", file.c_str (), log);
      return false;
    }
    return true;
  }

  // function to initialize GPU program
  bool initGPUProgram (const bool geometryFlag, const string &header, const string &file, GLuint &programId)
  {
    GLenum error;

    if (programId){
      glDeleteProgram (programId);
      checkGLError (error);
      if (error != GL_NO_ERROR){
        return false;
      }
    }

    programId = glCreateProgram ();
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    GLuint shaderId;

    // attach vertex shader
    string shader (file + ".vs");
    if ( !initGLShader (header, shader, shaderId, GL_VERTEX_SHADER)){
      PRINT ("error while initializing %s\n", shader.c_str ());
      return false;
    }
    glAttachShader (programId, shaderId);
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    // optionally attach geometry shader
    if (geometryFlag){
      shader = file;
      shader.append (".gs");

      if ( !initGLShader (header, shader, shaderId, GL_GEOMETRY_SHADER)){
        PRINT ("error while initializing %s\n", shader.c_str ());
        return false;
      }
      glAttachShader (programId, shaderId);
      checkGLError (error);
      if (error != GL_NO_ERROR){
        return false;
      }
    }

    // attach fragment shader
    shader = file;
    shader.append (".fs");

    if (!initGLShader (header, shader, shaderId, GL_FRAGMENT_SHADER)){
      PRINT ("error while initializing %s\n", shader.c_str ());
      return false;
    }
    glAttachShader (programId, shaderId);
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    // link program object
    glLinkProgram (programId);
    checkGLError (error);
    if (error != GL_NO_ERROR){
      return false;
    }

    return true;
  }
}
