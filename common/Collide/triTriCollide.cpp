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

#include "triTriCollide.h"

namespace SF {

  // set of functions to check intersection between two triangles (courtesy of Tomas Akenine-Moller)
  inline bool pointInTriTest (vec& v0, vec& u0, vec& u1, vec& u2, unsigned int i0, unsigned int i1)
	{
		real a = u1._v [i1] - u0._v [i1];
		real b = -(u1._v [i0] - u0._v [i0]);
		real c = -a * u0._v [i0] - b * u0._v [i1];
		real dd0 = a * v0._v [i0] + b * v0._v [i1] + c;

		a = u2._v [i1] - u1._v [i1];
		b = -(u2._v [i0] - u1._v [i0]);
		c = -a * u1._v [i0] - b * u1._v [i1];
		real dd1 = a * v0._v [i0] + b * v0._v [i1] + c;

		a = u0._v [i1] - u2._v [i1];
		b = -(u0._v [i0] - u2._v [i0]);
		c = -a * u2._v [i0] - b * u2._v [i1];
		real dd2 = a * v0._v [i0] + b * v0._v [i1] + c;

		if (dd0*dd1 > 0. && dd0*dd2 > 0.){
			return true;
		}
		return false;
	}

	inline bool edgeEdgeTest (real Ax, real Ay, vec& v0, vec& u0, vec& u1, unsigned int i0, unsigned int i1)
	{
		real Bx = u0._v [i0] - u1._v [i0];
		real By = u0._v [i1] - u1._v [i1];
		real Cx = v0._v [i0] - u0._v [i0];
		real Cy = v0._v [i1] - u0._v [i1];
		real f = Ay * Bx - Ax * By;
		real d = By * Cx - Bx * Cy;

		if ((f > 0. && d >= 0. && d <= f) || (f < 0. && d <= 0. && d >= f)) {
			real e = Ax * Cy - Ay * Cx;
			if (f > 0.){
				if (e >= 0. && e <= f) {
					return true;
				}
			} else {
				if (e <= 0. && e >= f) {
					return true;
				}
			}
		}

		return false;
	}

	inline bool edgeTriEdgeTest (vec& v0, vec& v1, vec& u0, vec& u1, vec& u2, unsigned int i0, unsigned int i1)
	{
		real Ax = v1._v [i0] - v0._v [i0];
		real Ay = v1._v [i1] - v0._v [i1];
		if (edgeEdgeTest (Ax, Ay, v0, u0, u1, i0, i1)) {
			return true;
		}
		if (edgeEdgeTest (Ax, Ay, v0, u1, u2, i0, i1)) {
			return true;
		}
		if (edgeEdgeTest(Ax, Ay, v0, u2, u0, i0, i1)) {
			return true;
		}
		return false;
	}

	inline bool coplanarTriTri (vec &e1, vec& n1, vec& v0, vec& v1, vec& v2, vec& u0, vec& u1, vec& u2)
	{
	  unsigned int i0, i1;

		// project normal to axis-aligned plane that maximizes area of triangle1 and compute indices
		e1._v [0] = ABS (n1._v [0]);
		e1._v [1] = ABS (n1._v [1]);
		e1._v [2] = ABS (n1._v [2]);
		if (e1._v [0] > e1._v [1]) {
			if (e1._v [0] > e1._v [2]) { // x-component greatest
				i0 = 1;
				i1 = 2;
			} else { // z-component greatest
				i0 = 0;
				i1 = 1;
			}
		}
		else {
			if (e1._v [2] > e1._v [1]) { // z-component greatest
				i0 = 0;
				i1 = 1;
			} else { // y-component greatest
				i0 = 0;
				i1 = 2;
			}
		}
		// test edges of triangle 1 against edges of triangle 2
		if (edgeTriEdgeTest (v0, v1, u0, u1, u2, i0, i1)) {
			return true;
		}
		if (edgeTriEdgeTest (v1, v2, u0, u1, u2, i0, i1)) {
			return true;
		}
		if (edgeTriEdgeTest (v2, v0, u0, u1, u2, i0, i1)) {
			return true;
		}

		// test if triangle 1 is totally contained in triangle 2
		if (pointInTriTest (v0, u0, u1, u2, i0, i1)) {
			return true;
		}
		if (pointInTriTest (v0, u0, u1, u2, i0, i1)){
			return true;
		}

		return false;
	}

	inline bool computeInterval (real vv0, real vv1, real vv2, real d0, real d1, real d2, real d0d1, real d0d2, real& a, real& b, real& c, real& x0, real& x1)
	{
		if (d0d1 > 0.) { // d0d2 <= 0, i.e d0, d1 are on the same side, d2 on the other or on the plane
			a = vv2;
			b = (vv0 - vv2) * d2;
			c = (vv1 - vv2) * d2;
			x0 = d2 - d0;
			x1 = d2 - d1;
			return true;
		}
		else if (d0d2 > 0.) { // d0d1 <= 0
			a = vv1;
			b = (vv0 - vv1) * d1;
			c = (vv2 - vv1) * d1;
			x0 = d1 - d0;
			x1 = d1 - d2;
			return true;
		}
		else if (d1*d2 > 0. || d0 != 0.) {
			a = vv0;
			b = (vv1 - vv0) * d0;
			c = (vv2 - vv0) * d0;
			x0 = d0 - d1;
			x1 = d0 - d2;
			return true;
		}
		else if (d1 != 0.) {
			a = vv1;
			b = (vv0 - vv1) * d1;
			c = (vv2 - vv1) * d1;
			x0 = d1 - d0;
			x1 = d1 - d2;
			return true;
		}
		else if (d2 != 0.) {
			a = vv2;
			b = (vv0 - vv2) * d2;
			c = (vv1 - vv2) * d2;
			x0 = d2 - d0;
			x1 = d2 - d1;
			return true;
		}
		return false;
	}

  /**
    * Fast Triangle-Triangle Collision detection technique:
    * 1. Trivial rejection test:
    * 		i.  See if all vertices of triangle 2 lie on the same side of triangle 1. If true, return "No Collision".
    * 		ii. See if all vertices of triangle 1 lie on the same side of triangle 2. If true, return "No Collision".
    * 2. Compute the direction of intersection line (cross product between 2 triangle normals).
    * 3. Compute the largest component (axis) of the intersection line direction vector.
    * 4. Project the 3 components of each triangle to plane perpendicular to the largest component.
    */
  bool triTriCollide (vec &n1, vec &v0, vec &v1, vec &v2, vec &n2, vec &u0, vec &u1, vec &u2, vec &e1)
  {
		/****************** STEP 1  (i) ******************/
		real d1 = -n1.dot  (v0);

		real du0 = n1.dot  (u0) + d1;
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
		real d2 = -n2.dot (u0);

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
			max = cc;
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
		real a, b, c, x0, x1;
		if (!computeInterval (vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, a, b, c, x0, x1)) {
			return coplanarTriTri (e1, n1, v0, v1, v2, u0, u1, u2);
		}

		// compute interval for triangle 2
		real d, e, f, y0, y1;
		if (!computeInterval (up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f, y0, y1)) {
			return coplanarTriTri (e1, n1, v0, v1, v2, u0, u1, u2);
		}

		real xx = x0 * x1;
		real yy = y0 * y1;
		real xxyy = xx * yy;

		real tempr = a * xxyy;

		real isect1 [2];
		isect1 [0] = tempr + b * x1 * yy;
		isect1 [1] = tempr + c * x0 * yy;

		tempr = d * xxyy;

		real isect2 [2];
		isect2 [0] = tempr + e * xx * y1;
		isect2 [1] = tempr + f * xx * y0;

		if (isect1 [0] > isect1 [1]) {
			tempr = isect1 [0];
			isect1 [0] = isect1 [1];
			isect1 [1] = tempr;
		}
		if (isect2 [0] > isect2 [1]) {
			tempr = isect2 [0];
			isect2 [0] = isect2 [1];
			isect2 [1] = tempr;
		}

		if (isect1 [1] < isect2 [0] || isect2 [1] < isect1 [0]) {
			return false;
		}

    return true;
  }

}
