/**
 * @file Preprocess.h
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
 * Definition for Simulation Framework prerequisites
 */

#pragma once

extern "C" {
#include "config.h"
}

// linux-specific preprocessor command to include GLExt
#if defined( __linux__ )
#define GL_GLEXT_PROTOTYPES 1
#define GLX_GLXEXT_PROTOTYPES 1
#endif

#ifndef SF_NO_PRINT
#include <cstdio>
#include <cstring>
// print debug message, prefixed by the file and line number
#define PRINT(format, ...) (fprintf (stdout, "%s[%d]:\t", basename( __FILE__ ), __LINE__), fprintf (stdout, format, ##__VA_ARGS__))
#else
// no debug messages
#define PRINT(format, ...) ((void) 0)
#endif

namespace SF {

	// set precision
#	ifndef SF_DOUBLE_PRECISION
	typedef float real;
#	else
	typedef double real;
#	endif

	// user defined constants
#	ifndef SF_DOUBLE_PRECISION
#		define EPSILON 1e-6
# else
#		define EPSILON 1e-9
# endif

	// user-defined macros
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef SIGN
#define SIGN(x) ((x) < 0 ? -1 : 1)
#endif

	// user-defined macros
#ifndef MAX
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#endif

#define EXPORT __attribute__((visibility("default")))
#define HIDDEN __attribute__((visibility("hidden")))
}
