/**
 * @file aabb.h
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
 * The axis-aligned bounding box class for Simulate Framework
 */

#pragma once

#include "Preprocess.h"
#include "vec3.h"
#include "vec4.h"

namespace SF {

	static inline bool axis_test (const int ind1, const int ind2, const real a, const real b, const real fa, const real fb,
			const vec &halflength, const vec &v0, const vec &v1, real& min, real& max)
	{
		real p0 = a*v0._v[1] - b*v0._v[2];
		real p1 = a*v1._v[1] - b*v1._v[2];
		if (p0 < p1 ){
			min = p0;
			max = p1;
		} else {
			min = p1;
			max = p0;
		}
		real rad = fa*halflength._v[ind1] + fb*halflength._v[ind2];
		return !(min > rad || max < -rad);
	}

	static inline void min_max (real x0, real x1, real x2, real& min, real& max)
	{
		min = max = x0;
		if (x1 < min){
			min = x1;
		}
		if (x1 > max){
			max = x1;
		}
		if (x2 < min){
			min = x2;
		}
		if (x2 > max){
			max = x2;
		}
	}

	static inline bool overlap (const vec &halflength, const vec &normal, const vec &vert)
	{
		real val;
		vec vmin, vmax;

		for (int i = 0; i < 3; ++i){
			val = normal._v[i];
			vmin._v[i] = val > 0. ? - (halflength._v[i] + vert._v[i]) : halflength._v[i] - vert._v[i];
			vmax._v[i] = val < 0. ? - (halflength._v[i] + vert._v[i]) : halflength._v[i] - vert._v[i];
		}
		return (vmin.dot (normal) <= EPSILON) & (vmax.dot (normal) >= EPSILON);
	}

	class aabb {

	public:
		vec _v[2];
		vec _center;
		vec _halflength;

	public:
		// default constructor
		inline aabb ()
		{
			_v[0] = vec::ZERO;
			_v[1] = vec::ZERO;
			_center = vec::ZERO;
			_halflength = vec::ZERO;
		}
		// copy constructor
		inline aabb (const aabb &bv)
		{
			memcpy (_v, bv._v, 2*SF_VECTOR_SIZE*sizeof (real));
			memcpy (_center._v, bv._center._v, SF_VECTOR_SIZE*sizeof (real));
			memcpy (_halflength._v, bv._halflength._v, SF_VECTOR_SIZE*sizeof (real));
		}
		// overloaded constructor
		inline aabb (const vec &v1, const vec &v2)
		{
			memcpy (_v[0]._v, v1._v, SF_VECTOR_SIZE*sizeof (real));
			memcpy (_v[1]._v, v2._v, SF_VECTOR_SIZE*sizeof (real));
			update ();
		}
		// overloaded constructor
		inline aabb (const real* v1, const real* v2)
		{
			memcpy (_v[0]._v, v1, SF_VECTOR_SIZE*sizeof (real));
			memcpy (_v[1]._v, v2, SF_VECTOR_SIZE*sizeof (real));
			update ();
		}
		// destructor
		inline ~aabb () { }

		// assignment operator
		inline aabb &operator = (const aabb &bv)
		{
			memcpy (_v, bv._v, 2*SF_VECTOR_SIZE*sizeof (real));
			memcpy (_center._v, bv._center._v, SF_VECTOR_SIZE*sizeof (real));
			memcpy (_halflength._v, bv._halflength._v, SF_VECTOR_SIZE*sizeof (real));
			return *this;
		}

		// corner accessor
		inline vec operator [] (const int i) const
		{
			assert (i >= 0 && i < 8);
			switch (i){
				case 0:
					return vec (_v[0]);
				case 1:
					return vec (_v[1]._v[0], _v[0]._v[1], _v[0]._v[2]);
				case 2:
					return vec (_v[0]._v[0], _v[1]._v[1], _v[0]._v[2]);
				case 3:
					return vec (_v[1]._v[0], _v[1]._v[1], _v[0]._v[2]);
				case 4:
					return vec (_v[0]._v[0], _v[0]._v[1], _v[1]._v[2]);
				case 5:
					return vec (_v[1]._v[0], _v[0]._v[1], _v[1]._v[2]);
				case 6:
					return vec (_v[0]._v[0], _v[1]._v[1], _v[1]._v[2]);
				case 7:
					return vec (_v[1]);
			}
		}

		// update method to update center and halflength
		inline void update ()
		{
			_center = _v[0] + _v[1];
			_center *= 0.5;
			_halflength = _center - _v[0];
		}

		// bounding box inside another bounding box test
		inline bool inside (const aabb &bv) const
		{
			return (bv._v[0]._v[0] >= _v[0]._v[0]) & (bv._v[1]._v[0] <= _v[1]._v[0]) & (bv._v[0]._v[1] >= _v[0]._v[1]) &
					(bv._v[1]._v[1] <= _v[1]._v[1]) & (bv._v[0]._v[2] >= _v[0]._v[2]) & (bv._v[1]._v[2] <= _v[1]._v[2]);
		}

		// bounding box collision test
		inline bool collide (const aabb &bv) const
		{
			if (inside (bv) || bv.inside (*this)){
				return true;
			}
			if (bv._v[0]._v[0] > _v[1]._v[0]){
				return false;
			}
			if (bv._v[0]._v[1] > _v[1]._v[1]){
				return false;
			}
			if (bv._v[0]._v[2] > _v[1]._v[2]){
				return false;
			}
			if (bv._v[1]._v[0] < _v[0]._v[0]){
				return false;
			}
			if (bv._v[1]._v[1] < _v[0]._v[1]){
				return false;
			}
			if (bv._v[1]._v[2] < _v[0]._v[2]){
				return false;
			}
			return true;
		}

		// vertex-bounding box collision test
		inline bool collide (const vec3 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (vec._v[i] - _v[0]._v[i] < -EPSILON){
					return false;
				}
			}
			for (int i = 0; i < 3; ++i){
				if (vec._v[i] - _v[1]._v[i] > EPSILON){
					return false;
				}
			}
			return true;
		}

		// vertex-bounding box collision test
		inline bool collide (const vec4 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (vec._v[i] - _v[0]._v[i] < -EPSILON){
					return false;
				}
			}
			for (int i = 0; i < 3; ++i){
				if (vec._v[i] - _v[1]._v[i] > EPSILON){
					return false;
				}
			}
			return true;
		}

		/*
		 * Collision method for triangle-AABB (courtesy of Tomas Akenine Moller)
		 * Uses separating axis theorem to test overlap between triangle and box
		 * Need to test for overlap in these directions:
		 * 1. the {x,y,z}-directions (actually, since we use the AABB of the triangle
		 *    we do not even need to test these)
		 * 2. normal of the triangle
		 * 3. crossproduct(edge from tri, {x,y,z}-direction)
		 * This gives 3x3=9 more tests
		 */
		inline bool collide (const vec &vec0, const vec &vec1, const vec &vec2) const
		{
			// move everything so that the boxcenter is in (0,0,0)
			vec v0 (vec0 - _center);
			vec v1 (vec1 - _center);
			vec v2 (vec2 - _center);

			// compute triangle edges
			vec e0 (v1 - v0); /* tri edge 0 */
			vec e1 (v2 - v1); /* tri edge 1 */
			vec e2 (v0 - v2); /* tri edge 2 */

			/* Bullet 3: test for the 9 cases first (this was faster) */
			real fex = ABS (e0._v[0]);
			real fey = ABS (e0._v[1]);
			real fez = ABS (e0._v[2]);

			real min = 0, max = 0;

			if (!axis_test (1, 2, e0._v[2], e0._v[1], fez, fey, _halflength, v0, v2, min, max)){
				return false;
			}
			if (!axis_test (0, 2, e0._v[2], e0._v[0], fez, fex, _halflength, v0, v2, min, max)){
				return false;
			}
			if (!axis_test (0, 1, e0._v[1], e0._v[0], fey, fex, _halflength, v1, v2, min, max)){
				return false;
			}

			fex = ABS (e1[0]);
			fey = ABS (e1[1]);
			fez = ABS (e1[2]);

			if (!axis_test (1, 2, e1._v[2], e1._v[1], fez, fey, _halflength, v0, v2, min, max)){
				return false;
			}
			if (!axis_test (0, 2, e1._v[2], e1._v[0], fez, fex, _halflength, v0, v2, min, max)){
				return false;
			}
			if (!axis_test (0, 1, e1._v[1], e1._v[0], fey, fex, _halflength, v0, v1, min, max)){
				return false;
			}

			fex = ABS (e2[0]);
			fey = ABS (e2[1]);
			fez = ABS (e2[2]);

			if (!axis_test (1, 2, e2._v[2], e2._v[1], fez, fey, _halflength, v0, v1, min, max)){
				return false;
			}
			if (!axis_test (0, 2, e2._v[2], e2._v[0], fez, fex, _halflength, v0, v1, min, max)){
				return false;
			}
			if (!axis_test (0, 1, e2._v[1], e2._v[0], fey, fex, _halflength, v1, v2, min, max)){
				return false;
			}

			/* Bullet 1:
			 * first test overlap in the {x,y,z}-directions
			 * find min, max of the triangle each direction, and test for overlap in
			 * that direction -- this is equivalent to testing a minimal AABB around
			 * the triangle against the AABB
			 */

			// test in X-direction
			min_max (v0._v[0], v1._v[0], v2._v[0], min, max);
			if (min > _halflength._v[0] || max < -_halflength._v[0]){
				return false;
			}

			// test in Y-direction
			min_max (v0._v[1], v1._v[1], v2._v[1], min, max);
			if (min > _halflength._v[1] || max < -_halflength._v[1]){
				return false;
			}

			// test in Z-direction
			min_max (v0._v[2], v1._v[2], v2._v[2], min, max);
			if (min > _halflength._v[2] || max < -_halflength._v[2]){
				return false;
			}

			/* Bullet 2:
			 * test if the box intersects the plane of the triangle
			 * compute plane equation of triangle: normal*x + d = 0
			 */
			vec normal;
			e0.fast_cross (normal, e1);
			return overlap (_halflength, normal, v0);
		}
	};
}
