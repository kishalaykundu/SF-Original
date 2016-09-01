/**
 * @file vec3.h
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
 * The 3-element vector class for Simulate Framework
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#include "Preprocess.h"
#include "vec2.h"

namespace SF {

	class vec3 {

	public:
		real _v[3];

	// methods described later
	public:
		static const vec3 ZERO;
		static const vec3 UNIT;
		static const vec3 UNIT_X;
		static const vec3 UNIT_Y;
		static const vec3 UNIT_Z;
		static const vec3 NEG_UNIT_X;
		static const vec3 NEG_UNIT_Y;
		static const vec3 NEG_UNIT_Z;

	public:
		// default contructor
		inline vec3 ()
		{
			for (int i = 0; i < 3; ++i){
				_v[i] = 0.;
			}
		}
		inline ~vec3 (){ }
		// overloaded constructor
		inline vec3 (const real* vec)
		{
			memcpy (_v, vec, 3 * sizeof (real));
		}
		// overloaded constructor
		inline vec3 (const vec3 &vec)
		{
			memcpy (_v, vec._v, 3 * sizeof (real));
		}
		// overloaded constructor
		inline vec3 (const real x, const real y, const real z)
		{
			_v[0] = x;
			_v[1] = y;
			_v[2] = z;
		}
		inline vec3 (const vec2 &vec, const real z)
		{
			memcpy (_v, vec._v, 2 * sizeof (real));
			_v[2] = z;
		}

		// assignment operator
		inline vec3 &operator = (const vec3 &vec)
		{
			memcpy (_v, vec._v, 3 * sizeof (real));
			return *this;
		}
		// assignment to val operator
		inline vec3 &operator = (const real val)
		{
			_v[0] = _v[1] = _v[2] = val;
			return *this;
		}

		// element mutator
		inline real &operator [] (const int i)
		{
			assert (i >= 0 && i < 3);
			return _v[i];
		}

		// largest coordinate
    inline real &largest_coord ()
    {
    	real a0 = ABS (_v[0]);
    	real a1 = ABS (_v[1]);
    	real a2 = ABS (_v[2]);
    	if (a0 > a1){
    		if (a0 > a2){
    			return _v[0];
    		} else {
    			return _v[2];
    		}
    	} else if (a1 > a2){
    		return _v[1];
    	} else {
    		return _v[2];
    	}
    }

		// smallest coordinate
    inline real &smallest_coord ()
    {
    	real a0 = fabs (_v[0]);
    	real a1 = fabs (_v[1]);
    	real a2 = fabs (_v[2]);
    	if (a0 < a1){
    		if (a0 < a2){
    			return _v[0];
    		} else {
    			return _v[2];
    		}
    	} else if (a1 < a2){
    		return _v[1];
    	} else {
    		return _v[2];
    	}
    }

		// unary negation operator
		inline vec3 operator - () const
		{
			return vec3 (-_v[0], -_v[1], -_v[2]);
		}

		// binary equality operator
		inline bool operator == (const vec3 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (ABS (_v[i] - vec._v[i]) >= EPSILON){
					return false;
				}
			}
			return true;
		}
		// binary inequality operator
		inline bool operator != (const vec3 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (ABS (_v[i] - vec._v[i]) >= EPSILON){
					return true;
				}
			}
			return false;
		}

		// non-mutating addition operator
		inline vec3 operator + (const vec3 &vec) const
		{
			return vec3 (_v[0] + vec._v[0], _v[1] + vec._v[1], _v[2] + vec._v[2]);
		}
		// non-mutating subtraction operator
		inline vec3 operator - (const vec3 &vec) const
		{
			return vec3 (_v[0] - vec._v[0], _v[1] - vec._v[1], _v[2] - vec._v[2]);
		}
		// non-mutating per-element multiplication operator
		inline vec3 operator * (const vec3 &vec) const
		{
			return vec3 (_v[0] * vec._v[0], _v[1] * vec._v[1], _v[2] * vec._v[2]);
		}
 		// scaling operator
		inline vec3 operator * (const real val) const
		{
			return vec3 (_v[0] * val, _v[1] * val, _v[2] * val);
		}
 		// scaling operator
		inline vec3 operator / (const real val) const
		{
			assert (ABS (val) >= EPSILON);
			real reciprocal = 1/ val;
			return vec3 (_v[0] * reciprocal, _v[1] * reciprocal, _v[2] * reciprocal);
		}
		// non-mutating per-element division operator
		inline vec3 operator / (const vec3 &vec) const
		{
			for (unsigned int i = 0; i < 3; ++i){
				assert (ABS (vec._v[i]) >= EPSILON);
			}
			return vec3 (_v[0]/ vec._v[0], _v[1]/ vec._v[1], _v[2]/ vec._v[2]);
		}

		// mutating addition operator
		inline vec3 &operator += (const vec3 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] += vec._v[i];
			}
			return *this;
		}
		// mutating addition operator
		inline vec3 &operator += (const real val)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] += val;
			}
			return *this;
		}
		// mutating subtraction operator
		inline vec3 &operator -= (const vec3 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] -= vec._v[i];
			}
			return *this;
		}
		// mutating subtraction operator
		inline vec3 &operator -= (const real val)
			{
				for (int i = 0; i < 3; ++i){
					_v[i] -= val;
				}
				return *this;
			}
		// mutating multiplication operator
		inline vec3 &operator *= (const vec3 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] *= vec._v[i];
			}
			return *this;
		}
		// mutating scaling operator
		inline vec3 &operator *= (const real val)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] *= val;
			}
			return *this;
		}
		// mutating division operator
		inline vec3 &operator /= (const vec3 &vec)
		{
			for (int i = 0; i < 3; ++i){
				assert (ABS (vec._v[i]) >= EPSILON);
			}
			for (int i = 0; i < 3; ++i){
				_v[i] /= vec._v[i];
			}
			return *this;
		}
		// mutating scaling operator
		inline vec3 &operator /= (const real val)
		{
			assert (ABS (val) >= EPSILON);
			real reciprocal = 1/ val;
			for (int i = 0; i < 3; ++i){
				_v[i] *=  reciprocal;
			}
			return *this;
		}

		// reflection around XY
		inline void reflect_xy (){ _v[2] = -_v[2]; }
		// reflection around YZ
		inline void reflect_zx (){ _v[1] = -_v[1]; }
		// reflection around YZ
		inline void reflect_yz (){ _v[0] = -_v[0]; }

		// dot product
		inline real dot (const vec3 &vec)
		{
			return _v[0] * vec._v[0] + _v[1] * vec._v[1] + _v[2] * vec._v[2];
		}
		// angle
		inline real angle (const vec3 &vec)
		{
			real mag = this->length ()* vec.length ();
			assert ( mag >= EPSILON);
			return this->dot (vec)/ mag;
		}

		// cross product
		inline vec3 cross (const vec3 &vec)
		{
			return vec3 (_v[1]*vec._v[2] - _v[2]*vec._v[1], _v[2]*vec._v[0] - _v[0]*vec._v[2], _v[0]*vec._v[1] - _v[1]*vec._v[0]);
		}
		// fast cross product (no Vector constructor needed)
		inline void fast_cross (vec3 &prod, const vec3 &vec)
		{
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];
		}

		// normalized cross product
		inline vec3 ncross (const vec3 &vec)
		{
			vec3 prod;
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];

			real mag = prod._v[0] * prod._v[0] + prod._v[1] * prod._v[1] + prod._v[2] * prod._v[2];
			assert (mag >= EPSILON);

			real invmag = 1/ static_cast< real > (sqrt (mag));
			for (unsigned int i = 0; i < 3; ++i){
				prod._v[i] *= invmag;
			}
			return prod;
		}
		// fast normalized cross product (no Vector constructor needed)
		inline void fast_ncross (vec3 &prod, const vec3 &vec)
		{
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];
			real mag = prod._v[0] * prod._v[0] + prod._v[1] * prod._v[1] + prod._v[2] * prod._v[2];

			assert (mag >= EPSILON);
			real invmag = 1./ static_cast< real > (sqrt (mag));
			for (unsigned int i = 0; i < 3; ++i){
				prod._v[i] *= invmag;
			}
		}

		// length
		inline real length () const
		{
			return static_cast< real > (sqrt (_v[0]*_v[0] +  _v[1]*_v[1] + _v[2]*_v[2]));
		}
		// squared length
		inline real square_length () const
		{
			return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2];
		}

		// normalization method
		inline void normalize ()
		{
			real mag = length ();
			assert (mag >= EPSILON);
			mag = 1./ mag;
			for (int i = 0; i < 3; ++i){
				_v[i] *= mag;
			}
		}

		// distance from another vector
 		inline real distance (const vec3 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			real t3 = vec._v[2] - _v[2];
			return static_cast< real > (sqrt (t1*t1 + t2*t2 + t3*t3));
		}
		// squared distance from another vector
		inline real square_dist (const vec3 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			real t3 = vec._v[2] - _v[2];
			return t1*t1 + t2*t2 + t3*t3;
		}

	};
}

#ifdef SF_VECTOR3_ENABLED
#ifndef vec
	const unsigned int SF_VECTOR_SIZE = 3;
	typedef class SF::vec3 vec;
#endif
#endif
