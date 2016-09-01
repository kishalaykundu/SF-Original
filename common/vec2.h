/**
 * @file vec2.h
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
 * The 2-element vector class for Simulate Framework
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#include "Preprocess.h"

namespace SF {

	class vec2 {

	public:
		real _v[2];

	// constant vectors
	public:
		static const vec2 ZERO;
		static const vec2 UNIT;
		static const vec2 UNIT_X;
		static const vec2 UNIT_Y;
		static const vec2 NEG_UNIT_X;
		static const vec2 NEG_UNIT_Y;

	public:
		// default contructor
		inline vec2 ()
		{
			for (int i = 0; i < 2; ++i){
				_v[i] = 0.;
			}
		}
		inline ~vec2 (){ }
		// overloaded constructor
		inline vec2 (const real *vec){ memcpy (_v, vec, 2*sizeof (real)); }
		// overloaded constructor
		inline vec2 (const vec2 &vec){ memcpy (_v, vec._v, 2*sizeof (real)); }
		// overloaded constructor
		inline vec2 (const real x, const real y)
		{
			_v[0] = x;
			_v[1] = y;
		}

		// assignment operator
		inline vec2 &operator = (const vec2 &vec)
		{
			memcpy (_v, vec._v, 2*sizeof (real));
			return *this;
		}
		// assignment to val operator
		inline vec2 &operator = (const real val)
		{
			_v[0] = val;
			_v[1] = val;
			return *this;
		}

		// element mutator
		inline real &operator [] (const int i)
		{
			assert (i >= 0 && i < 2);
			return _v[i];
		}

		// largest coordinate
    inline real &largest_coord ()
    {
    	real a0 = ABS (_v[0]);
    	real a1 = ABS (_v[1]);
			return a0 > a1? _v[0] : _v[1];
    }

		// smallest coordinate
    inline real &smallest_coord ()
    {
    	real a0 = ABS (_v[0]);
    	real a1 = ABS (_v[1]);
			return a0 < a1? _v[0] : _v[1];
    }

		// unary negation operator
		inline vec2 operator - () const
		{
			return vec2 (-_v[0], -_v[1]);
		}

		// binary equality operator
		inline bool operator == (const vec2 &vec) const
		{
			for (int i = 0; i < 2; ++i){
				if (ABS (_v[i] - vec._v[i]) > EPSILON){
					return false;
				}
			}
			return true;
		}
		// binary inequality operator
		inline bool operator != (const vec2 &vec) const
		{
			for (int i = 0; i < 2; ++i){
				if (ABS (_v[i] - vec._v[i]) > EPSILON){
					return true;
				}
			}
			return false;
		}

		// non-mutating addition operator
		inline vec2 operator + (const vec2 &vec) const
		{
			return vec2 (_v[0] + vec._v[0], _v[1] + vec._v[1]);
		}
		// non-mutating subtraction operator
		inline vec2 operator - (const vec2 &vec) const
		{
			return vec2 (_v[0] - vec._v[0], _v[1] - vec._v[1]);
		}
		// non-mutating per-element multiplication operator
		inline vec2 operator * (const vec2 &vec) const
		{
			return vec2 (_v[0] * vec._v[0], _v[1] * vec._v[1]);
		}
 		// scaling operator
		inline vec2 operator * (const real val) const
		{
			return vec2 (_v[0] * val, _v[1] * val);
		}
 		// scaling operator
		inline vec2 operator / (const real val) const
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1/ val;
			return vec2 (_v[0] * reciprocal, _v[1] * reciprocal);
		}
		// non-mutating per-element division operator
		inline vec2 operator / (const vec2 &vec) const
		{
			for (int i = 0; i < 2; ++i){
				assert (ABS (vec._v[i]) > EPSILON);
			}
			return vec2 (_v[0]/ vec._v[0], _v[1]/ vec._v[1]);
		}

		// mutating addition operator
		inline vec2 &operator += (const vec2 &vec)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] += vec._v[i];
			}
			return *this;
		}
		// mutating addition operator
		inline vec2 &operator += (const real val)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] += val;
			}
			return *this;
		}
		// mutating subtraction operator
		inline vec2 &operator -= (const vec2 &vec)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] -= vec._v[i];
			}
			return *this;
		}
		// mutating subtraction operator
		inline vec2 &operator -= (const real val)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] -= val;
			}
			return *this;
		}
		// mutating multiplication operator
		inline vec2 &operator *= (const vec2 &vec)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] *= vec._v[i];
			}
			return *this;
		}
		// mutating scaling operator
		inline vec2 &operator *= (const real val)
		{
			for (int i = 0; i < 2; ++i){
				_v[i] *= val;
			}
			return *this;
		}
		// mutating division operator
		inline vec2 &operator /= (const vec2 &vec)
		{
			for (int i = 0; i < 2; ++i){
				assert (ABS (vec._v[i]) > EPSILON);
			}
			for (int i = 0; i < 2; ++i){
				_v[i] /= vec._v[i];
			}
			return *this;
		}
		// mutating scaling operator
		inline vec2 &operator /= (const real val)
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1/ val;
			for (int i = 0; i < 2; ++i){
				_v[i] *=  reciprocal;
			}
			return *this;
		}

		// reflection around X
		inline void reflect_x (){ _v[1] = -_v[1]; }
		// reflection around Y
		inline void reflect_y (){ _v[0] = -_v[0]; }

		// dot product
		inline real dot (const vec2 &vec)
		{
			return _v[0] * vec._v[0] + _v[1] * vec._v[1];
		}
		// angle
		inline real angle (const vec2 &vec)
		{
			real mag = this->length ()* vec.length ();
			assert ( mag > EPSILON);
			return this->dot (vec)/ mag;
		}

		// length
		inline real length () const
		{
			return static_cast< real > (sqrt (_v[0]*_v[0] +  _v[1]*_v[1]));
		}
		// squared length
		inline real square_length () const
		{
			return _v[0]*_v[0] + _v[1]*_v[1];
		}

		// distance from another vector
		inline real distance (const vec2 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			return static_cast< real > (sqrt (t1*t1 + t2*t2));
		}
		// squared distance from another vector
		inline real square_dist (const vec2 &vec) const
		{
			real t1 = vec._v[0] - _v[0];
			real t2 = vec._v[1] - _v[1];
			return t1*t1 + t2*t2;
		}

		// normalize
		inline void normalize ()
		{
			assert (length () > EPSILON);
			real invmag = 1/ length ();
			for (int i = 0; i < 2; ++i){
				_v[i] *= invmag;
			}
		}

	};
}
