/**
 * @file Intersect.h
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
 * Functions for fast tetrahedron-blade intersection response. Part of
 * the CU_XFEM library.
 */

#pragma once

#include <vector>

#include "Preprocess.h"

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

#include "Cut.h"
#include "Cell.h"
#include "Collide/triTriCollide.h"

using namespace std;

namespace SF {
  namespace XFE {

    inline void SORT (real &a, real &b, int &index)
    {
      if (a > b){
        real c = a;
        a = b;
        b = c;
        index = 1;
      }
    }

    inline void intersect (vec &v0, vec &v1, vec &v2, real vv0, real vv1, real vv2, real d0, real d1, real d2, real &isect0, real &isect1, vec &isectpt0, vec &isectpt1)
    {
      real tmp = d0/ (d0 - d1);
      isect0 = vv0 + (vv1 - vv0) * tmp;
      isectpt0 = v0 + (v1 - v0) * tmp;

      tmp = d0/ (d0 - d2);
      isect1 = vv0 + (vv2 - vv0) * tmp;
      isectpt1 = v0 + (v2 - v0) * tmp;
    }

    static bool computeIntersection (vec &vert0, vec &vert1, vec &vert2, real vv0, real vv1, real vv2, real d0, real d1, real d2,
                                     real d0d1, real d0d2, real &isect0, real &isect1, vec &isectpt0, vec &isectpt1)
    {
      if (d0d1 > 0.){
        intersect (vert2, vert0, vert1, vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpt0, isectpt1);
      }
      else if (d0d2 > 0.){
        intersect (vert1, vert0, vert2, vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpt0, isectpt1);
      }
      else if (d1*d2 > 0. || d0 != 0.){
        intersect (vert0, vert1, vert2, vv0, vv1, vv2, d0, d1, d2, isect0, isect1, isectpt0, isectpt1);
      }
      else if (d1 != 0.){
        intersect (vert1, vert0, vert2, vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpt0, isectpt1);
      }
      else if (d2 != 0.){
        intersect (vert2, vert0, vert1, vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpt0, isectpt1);
      }
      else {
        return true;
      }
      return false;
    }

    /** Fast Triangle-Triangle Collision detection technique:
    * 1. Trivial rejection test:
    * 		i.  See if all vertices of triangle 2 lie on the same side of triangle 1. If true, return "No Collision".
    * 		ii. See if all vertices of triangle 1 lie on the same side of triangle 2. If true, return "No Collision".
    * 2. Compute the direction of intersection line (cross product between 2 triangle normals).
    * 3. Compute the largest component (axis) of the intersection line direction vector.
    * 4. Project the 3 components of each triangle to plane perpendicular to the largest component.
    */
    static bool
    triTriIntersect (vec &n1, vec &v0, vec &v1, vec &v2, vec &n2, vec &u0, vec &u1, vec &u2, vec &e1, vec &l1, vec &l2)
    {
      /****************** STEP 1  (i) ******************/
      real d1 = - (n1.dot (v0));

      real du0 = n1.dot (u0) + d1;
      if (ABS (du0) < EPSILON) {
        du0 = 0.;
      }
      real du1 = n1.dot (u1) + d1;
      if (ABS (du1) < EPSILON) {
        du1 = 0.;
      }
      real du2 = n1.dot (u2) + d1;
      if (ABS (du2) < EPSILON) {
        du2 = 0.;
      }

      real du0du1 = du0 * du1;
      real du0du2 = du0 * du2;

      // same non-zero sign on all of them - no intersection
      if (du0du1 > 0. && du0du2 > 0.) {
        return false;
      }

      /****************** STEP 1  (ii) ******************/
      real d2 = - (n2.dot (u0));

      real dv0 = n2.dot (v0) + d2;
      if (ABS (dv0) < EPSILON) {
        dv0 = 0.;
      }
      real dv1 = n2.dot (v1) + d2;
      if (ABS (dv1) < EPSILON) {
        dv1 = 0.;
      }
      real dv2 = n2.dot (v2) + d2;
      if (ABS (dv2) < EPSILON) {
        dv2 = 0.;
      }

      real dv0dv1 = dv0 * dv1;
      real dv0dv2 = dv0 * dv2;

      // same non-zero sign on all of them - no intersection
      if (dv0dv1 > 0. && dv0dv2 > 0.) {
        return false;
      }

      /****************** STEP 2 ******************/
      n1.fast_cross (e1, n2);

      /****************** STEP 3 ******************/
      unsigned int index = 0;
      real max = ABS (e1._v [0]);
      real bb = ABS (e1._v [1]);
      real cc = ABS (e1._v [2]);
      if (bb > max) {
        max = bb;
        index = 1;
      }
      if (cc > max) {
        index = 2;
      }

      /****************** STEP 4 ******************/
      real vp0 = v0._v [index];
      real vp1 = v1._v [index];
      real vp2 = v2._v [index];

      real up0 = u0._v [index];
      real up1 = u1._v [index];
      real up2 = u2._v [index];

      // compute interval for triangle 1
      real isect1 [2];
      vec isectptA1, isectptA2;
      if (computeIntersection (v0, v1, v2, vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, isect1 [0], isect1 [1], isectptA1, isectptA2)){
        return false;
      }

      // compute interval for triangle 2
      real isect2 [2];
      vec isectptB1, isectptB2;
      computeIntersection (u0, u1, u2, up0, up1, up2, du0, du1, du2, du0du1, du0du2, isect2 [0], isect2 [1], isectptB1, isectptB2);

      int smallest1 = 0;
      SORT (isect1 [0], isect1 [1], smallest1);

      int smallest2 = 0;
      SORT (isect2 [0], isect2 [1], smallest2);

      if (isect1 [1] < isect2 [0] || isect2 [1] < isect1 [0]){
        return false;
      }

      // set proper intersection points
      if (isect2 [0] < isect1 [0]){
        l1 = smallest1 > 0 ? isectptA2 : isectptA1;
        if (isect2 [1] < isect1 [1]){
          l2 = smallest2 > 0 ? isectptB1 : isectptB2;
        } else {
          l2 = smallest1 > 0 ? isectptA1 : isectptA2;
        }
      } else {
        l1 = smallest2 > 0 ? isectptB2 : isectptB1;
        if (isect2 [1] > isect1 [1]){
          l2 = smallest1 > 0 ? isectptA1 : isectptA2;
        } else {
          l2 = smallest2 > 0 ? isectptB1 : isectptB2;
        }
      }

      return true;
    }

  }
}
