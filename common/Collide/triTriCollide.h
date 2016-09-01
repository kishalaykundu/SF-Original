/**
 * @file triTriCollide.h
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
 * Functions for fast triangle-triangle collision detections.
 */

#pragma once

#include "Preprocess.h"

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

namespace SF {

  /**
    * Fast Triangle-Triangle Collision detection technique:
    * 1. Trivial rejection test:
    * 		i.  See if all vertices of triangle 2 lie on the same side of triangle 1. If true, return "No Collision".
    * 		ii. See if all vertices of triangle 1 lie on the same side of triangle 2. If true, return "No Collision".
    * 2. Compute the direction of intersection line (cross product between 2 triangle normals).
    * 3. Compute the largest component (axis) of the intersection line direction vector.
    * 4. Project the 3 components of each triangle to plane perpendicular to the largest component.
    */
  bool triTriCollide (vec &n1, vec &v0, vec &v1, vec &v2, vec &n2, vec &u0, vec &u1, vec &u2, vec &e1);
}
