/**
 * @file vec4.h
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
 * The 4-element vector class for Simulate Framework
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#include "Preprocess.h"
#include "vec2.h"
#include "vec3.h"

namespace SF {

	class vec4 {

	public:
		real _v[4];

	// methods described later
	public:
		static const vec4 ZERO;
		static const vec4 UNIT;
		static const vec4 UNIT_X;
		static const vec4 UNIT_Y;
		static const vec4 UNIT_Z;
		static const vec4 NEG_UNIT_X;
		static const vec4 NEG_UNIT_Y;
		static const vec4 NEG_UNIT_Z;

	public:
		// default constructor/ destructor
		inline vec4 ()
		{
			for (int i = 0; i < 3; ++i){
				_v[i] = 0.;
			}
			_v[3] = 1.;
		}
		inline ~vec4 (){ }

		// overloaded constructors
		inline vec4 (const real* vec)
		{
			memcpy (_v, vec, 4 * sizeof (real));
		}
		inline vec4 (const vec4 &vec)
		{
			memcpy (_v, vec._v, 4 * sizeof (real));
		}
		inline vec4 (const real x, const real y, const real z, const real w)
		{
			_v[0] = x;
			_v[1] = y;
			_v[2] = z;
			_v[3] = w;
		}
		inline vec4 (const real x, const real y, const real z)
		{
			_v[0] = x;
			_v[1] = y;
			_v[2] = z;
			_v[3] = 1.;
		}
		inline vec4 (const vec3 &vec, const real w)
		{
			memcpy (_v, vec._v, 3 * sizeof (real));
			_v[3] = w;
		}
		inline vec4 (const vec2 &vec, const real z, const real w)
		{
			memcpy (_v, vec._v, 2 * sizeof (real));
			_v[2] = z;
			_v[3] = w;
		}
		inline vec4 (const vec3 &vec)
		{
			memcpy (_v, vec._v, 3 * sizeof (real));
			_v[3] = 1.;
		}

		// assignment operators
		inline vec4 &operator = (const vec4 &vec)
		{
			memcpy (_v, vec._v, 4 * sizeof (real));
			return *this;
		}

		inline vec4 &operator = (const real val)
		{
			_v[0] = _v[1] = _v[2] = val; _v[3] = 1.0;
			return *this;
		}

		// element mutator
		inline real &operator [] (const int i)
		{
			assert (i >= 0 && i < 4);
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
    	real a0 = ABS (_v[0]);
    	real a1 = ABS (_v[1]);
    	real a2 = ABS (_v[2]);
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
		inline vec4 operator - () const
		{
			return vec4 (-_v[0], -_v[1], -_v[2], _v[3]);
		}

		// binary comparison operators
		inline bool operator == (const vec4 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (ABS (_v[i] - vec._v[i]) >= EPSILON){
					return false;
				}
			}
			return true;
		}
		inline bool operator != (const vec4 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				if (ABS (_v[i] - vec._v[i]) >= EPSILON){
					return true;
				}
			}
			return false;
		}

		// non-mutating mathematical operators
		inline vec4 operator + (const vec4 &vec) const
		{
			return vec4 (_v[0] + vec._v[0], _v[1] + vec._v[1], _v[2] + vec._v[2], 1.0);
		}
		inline vec4 operator - (const vec4 &vec) const
		{
			return vec4 (_v[0] - vec._v[0], _v[1] - vec._v[1], _v[2] - vec._v[2], 1.0);
		}
		inline vec4 operator * (const vec4 &vec) const
		{
			return vec4 (_v[0] * vec._v[0], _v[1] * vec._v[1], _v[2] * vec._v[2], 1.0);
		}
		inline vec4 operator * (const real val) const
		{
			return vec4 (_v[0] * val, _v[1] * val, _v[2] * val, 1.0);
		}
		inline vec4 operator / (const real val) const
		{
			assert (ABS (val) >= EPSILON);
			real reciprocal = 1/ val;
			return vec4 (_v[0] * reciprocal, _v[1] * reciprocal, _v[2] * reciprocal, 1.0);
		}
		inline vec4 operator / (const vec4 &vec) const
		{
			for (int i = 0; i < 3; ++i){
				assert (ABS (vec._v[i]) >= EPSILON);
			}
			return vec4 (_v[0]/ vec._v[0], _v[1]/ vec._v[1], _v[2]/ vec._v[2], 1.0);
		}

		// mutating mathematical operators
		inline vec4 &operator += (const vec4 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] += vec._v[i];
			}
			return *this;
		}
		inline vec4 &operator += (const real val)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] += val;
			}
			return *this;
		}
		inline vec4 &operator -= (const vec4 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] -= vec._v[i];
			}
			return *this;
		}
		inline vec4 &operator -= (const real val)
			{
				for (int i = 0; i < 3; ++i){
					_v[i] -= val;
				}
				return *this;
			}
		inline vec4 &operator *= (const vec4 &vec)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] *= vec._v[i];
			}
			return *this;
		}
		inline vec4 &operator *= (const real val)
		{
			for (int i = 0; i < 3; ++i){
				_v[i] *= val;
			}
			return *this;
		}
		inline vec4 &operator /= (const vec4 &vec)
		{
			for (int i = 0; i < 3; ++i){
				assert (ABS (vec._v[i]) >= EPSILON);
			}
			for (unsigned int i = 0; i < 3; ++i){
				_v[i] /= vec._v[i];
			}
			return *this;
		}
		inline vec4 &operator /= (const real val)
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
		inline real dot (const vec4 &vec)
		{
			return _v[0] * vec._v[0] + _v[1] * vec._v[1] + _v[2] * vec._v[2];
		}
		// angle
		inline real angle (const vec4 &vec)
		{
			real mag = this->length ()* vec.length ();
			assert ( mag >= EPSILON);
			return this->dot (vec)/ mag;
		}

		// cross product
		inline vec4 cross (const vec4 &vec)
		{
			return vec4 (_v[1]*vec._v[2] - _v[2]*vec._v[1], _v[2]*vec._v[0] - _v[0]*vec._v[2], _v[0]*vec._v[1] - _v[1]*vec._v[0], 1.);
		}
		// fast cross product (no constructor needed)
		inline void fast_cross (vec4 &prod, const vec4 &vec)
		{
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];
		}

		// normalized cross product
		inline vec4 ncross (const vec4 &vec)
		{
			vec4 prod;
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];
			real mag = prod._v[0] * prod._v[0] + prod._v[1] * prod._v[1] + prod._v[2] * prod._v[2];

			assert (mag >= EPSILON);
			real invmag = 1./ sqrt (mag);
			for (int i = 0; i < 3; ++i){
				prod._v[i] *= invmag;
			}
			return prod;
		}
		// fast normalized cross product (no constructor needed)
		inline void fast_ncross (vec4 &prod, const vec4 &vec)
		{
			prod._v[0] = _v[1]*vec._v[2] - _v[2]*vec._v[1];
			prod._v[1] = _v[2]*vec._v[0] - _v[0]*vec._v[2];
			prod._v[2] = _v[0]*vec._v[1] - _v[1]*vec._v[0];
			real mag = prod._v[0] * prod._v[0] + prod._v[1] * prod._v[1] + prod._v[2] * prod._v[2];
      mag = static_cast< real > (sqrt (mag));
			assert (mag > EPSILON);
			mag = 1./ mag;
			for (int i = 0; i < 3; ++i){
				prod._v[i] *= mag;
			}
		}

		// length
		inline real length () const
		{
			return static_cast <real> (sqrt (_v[0]*_v[0] +  _v[1]*_v[1] + _v[2]*_v[2]));
		}
		// squared length
		inline real square_length () const
		{
			return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2];
		}

		// distance from another vector
		inline real distance (const vec4 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			real t3 = vec._v[2] - _v[2];
			return static_cast< real > (sqrt (t1*t1 + t2*t2 + t3*t3));
		}
		// squared distance from another vector
		inline real square_dist (const vec4 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			real t3 = vec._v[2] - _v[2];
			return t1*t1 + t2*t2 + t3*t3;
		}

		// normalize
		inline void normalize ()
		{
		  real mag = static_cast <real> (sqrt (_v[0]*_v[0] +  _v[1]*_v[1] + _v[2]*_v[2]));
			assert (mag >= EPSILON);
			mag = 1./ mag;
			for (int i = 0; i < 3; ++i){
				_v[i] *= mag;
			}
		}
	};
}

#ifdef SF_VECTOR4_ENABLED
#ifndef vec
	const unsigned int SF_VECTOR_SIZE = 4;
	typedef class SF::vec4 vec;
#endif
#endif
