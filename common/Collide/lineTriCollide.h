/**
 * @file lineTriCollide.h
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
 * Functions for line segment-triangle collision detections. It also
 * returns the local co-ordinate of the point in the line-segment.
 */

#pragma once

#include "Preprocess.h"
#include "mat3x3.h"

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

namespace SF {

  // list of functions to detect line triangle intersection
  bool pointInTriangle (vec &p, vec &t1, vec &t2, vec &t3, vec &normal, bool planeTestFlag = true);

  bool lineLineCollide (vec &l11, vec &l12, vec &l21, vec &l22);

  bool lineTriCollide (real &eu, vec &l1, vec &l2, vec &t1, vec &t2, vec &t3, vec &normal);
}
