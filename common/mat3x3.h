/**
 * @file mat3x3.h
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
 * The 3 by 3 matrix class for Simulate Framework
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#include "Preprocess.h"
#include "vec3.h"

namespace SF {

	class mat3x3 {

	public:
		real _m[9];

		// methods described later
	public:
		static const mat3x3 ZERO;
		static const mat3x3 ONE;
		static const mat3x3 IDENTITY;

	public:
		inline mat3x3 ()
		{
			for (int i = 0; i < 9; ++i){
				_m[i] = 0.;
			}
		}
		// overloaded constructor
		inline mat3x3 (const real* &arr)
		{
			memcpy (_m, arr, 9*sizeof (real));
		}
		// overloaded constructor
		inline mat3x3 (const mat3x3 &mat)
		{
			memcpy (_m, mat._m, 9*sizeof (real));
		}
		// overloaded constructor
		inline mat3x3 (const real a00, const real a01, const real a02,
				const real a10, const real a11, const real a12,
				const real a20, const real a21, const real a22)
		{
			_m[0] = a00; _m[1] = a01; _m[2] = a02;
			_m[3] = a10; _m[4] = a11; _m[5] = a12;
			_m[6] = a20; _m[7] = a21; _m[8] = a22;
		}

		// assignment operator
		inline mat3x3 &operator = (const mat3x3 &mat)
		{
			memcpy (_m, mat._m, 9*sizeof (real));
			return *this;
		}

		// element access and mutation operator
		inline real &operator () (const int i, const int j)
		{
			assert (i >= 0 && i < 3);
			assert (j >= 0 && j < 3);
			return _m[3*i + j];
		}

		// non-mutating addition operator
		inline mat3x3 operator + (const mat3x3 &mat) const
		{
			return mat3x3 (_m[0]+mat._m[0], _m[1]+mat._m[1], _m[2]+mat._m[2],
												_m[3]+mat._m[3], _m[4]+mat._m[4], _m[5]+mat._m[5],
												_m[6]+mat._m[6], _m[7]+mat._m[7], _m[8]+mat._m[8]);
		}
		// non-mutating subtraction operator
		inline mat3x3 operator - (const mat3x3 &mat) const
		{
			return mat3x3 (_m[0]-mat._m[0], _m[1]-mat._m[1], _m[2]-mat._m[2],
												_m[3]-mat._m[3], _m[4]-mat._m[4], _m[5]-mat._m[5],
												_m[6]-mat._m[6], _m[7]-mat._m[7], _m[8]-mat._m[8]);
		}
		// non-mutating per-element multiplication operator
		inline mat3x3 operator * (const mat3x3 &mat) const
		{
			return mat3x3 (_m[0]*mat._m[0] + _m[1]*mat._m[3] + _m[2]*mat._m[6],
												_m[0]*mat._m[1] + _m[1]*mat._m[4] + _m[2]*mat._m[7],
												_m[0]*mat._m[2] + _m[1]*mat._m[5] + _m[2]*mat._m[8],
												_m[3]*mat._m[0] + _m[4]*mat._m[3] + _m[5]*mat._m[6],
												_m[3]*mat._m[1] + _m[4]*mat._m[4] + _m[5]*mat._m[7],
												_m[3]*mat._m[2] + _m[4]*mat._m[5] + _m[5]*mat._m[8],
												_m[6]*mat._m[0] + _m[7]*mat._m[3] + _m[8]*mat._m[6],
												_m[6]*mat._m[1] + _m[7]*mat._m[4] + _m[8]*mat._m[7],
												_m[6]*mat._m[2] + _m[7]*mat._m[5] + _m[8]*mat._m[8]);
		}
		// vector multiplication
		inline vec3 operator * (vec3 &v)
		{
			return vec3 (v._v[0]*_m[0] + v._v[1]*_m[1] + v._v[2]*_m[2],
											v._v[0]*_m[3] + v._v[1]*_m[4] + v._v[2]*_m[5],
											v._v[0]*_m[6] + v._v[1]*_m[7] + v._v[2]*_m[8]);
		}
 		// scaling operator
		inline mat3x3 operator * (const real val) const
		{
			return mat3x3 (_m[0]*val, _m[1]*val, _m[2]*val,
												_m[3]*val, _m[4]*val, _m[5]*val,
												_m[6]*val, _m[7]*val, _m[8]*val);
		}
 		// scaling operator
		inline mat3x3 operator / (const real val) const
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1/ val;
			return mat3x3 (_m[0]*reciprocal, _m[1]*reciprocal, _m[2]*reciprocal,
												_m[3]*reciprocal, _m[4]*reciprocal, _m[5]*reciprocal,
												_m[6]*reciprocal, _m[7]*reciprocal, _m[8]*reciprocal);
		}

		// mutating addition operator
		inline mat3x3 &operator += (const mat3x3 &mat)
		{
			for (int i = 0; i < 9; ++i){
				_m[i] += mat._m[i];
			}
			return *this;
		}
		// mutating addition operator
		inline mat3x3 &operator += (const real val)
		{
			for (int i = 0; i < 9; ++i){
				_m[i] += val;
			}
			return *this;
		}
		// mutating subtraction operator
		inline mat3x3 &operator -= (const mat3x3 &mat)
		{
			for (int i = 0; i < 9; ++i){
				_m[i] -= mat._m[i];
			}
			return *this;
		}
		// mutating subtraction operator
		inline mat3x3 &operator -= (const real val)
		{
			for (int i = 0; i < 9; ++i){
				_m[i] -= val;
			}
			return *this;
		}
		// mutating multiplication operator
		inline mat3x3 &operator *= (const mat3x3 &mat)
		{
			real tmp[9];
			tmp[0] = _m[0]*mat._m[0] + _m[1]*mat._m[3] + _m[2]*mat._m[6];
			tmp[1] = _m[0]*mat._m[1] + _m[1]*mat._m[4] + _m[2]*mat._m[7];
			tmp[2] = _m[0]*mat._m[2] + _m[1]*mat._m[5] + _m[2]*mat._m[8];
			tmp[3] = _m[3]*mat._m[0] + _m[4]*mat._m[3] + _m[5]*mat._m[6];
			tmp[4] = _m[3]*mat._m[1] + _m[4]*mat._m[4] + _m[5]*mat._m[7];
			tmp[5] = _m[3]*mat._m[2] + _m[4]*mat._m[5] + _m[5]*mat._m[8];
			tmp[6] = _m[6]*mat._m[0] + _m[7]*mat._m[3] + _m[8]*mat._m[6];
			tmp[7] = _m[6]*mat._m[1] + _m[7]*mat._m[4] + _m[8]*mat._m[7];
			tmp[8] = _m[6]*mat._m[2] + _m[7]*mat._m[5] + _m[8]*mat._m[8];
			memcpy (_m, tmp, 9*sizeof (real));
			return *this;
		}
		// mutating scaling operator
		inline mat3x3 &operator *= (const real val)
		{
			for (int i = 0; i < 9; ++i){
				_m[i] *= val;
			}
			return *this;
		}
		// mutating scaling operator
		inline mat3x3 &operator /= (const real val)
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1/ val;
			for (int i = 0; i < 9; ++i){
				_m[i] *=  reciprocal;
			}
			return *this;
		}

		// transpose operator
		inline void transpose ()
		{
			real tmp;
			tmp = _m[1]; _m[1] = _m[3]; _m[3] = tmp;
			tmp = _m[2]; _m[2] = _m[6]; _m[6] = tmp;
			tmp = _m[5]; _m[5] = _m[7]; _m[7] = tmp;
		}

		// determinant calculator
		inline real determinant ()
		{
			return _m[0]* (_m[4]*_m[8] - _m[5]*_m[7]) + _m[1]* (_m[5]*_m[6] - _m[3]*_m[8]) + _m[2]* (_m[3]*_m[7] - _m[4]*_m[6]);
		}

		// inverse operator
		inline void invert ()
		{
			real tmp[9];
			tmp[0] = _m[4]*_m[8] - _m[5]*_m[7];
			tmp[1] = _m[2]*_m[7] - _m[1]*_m[8];
			tmp[2] = _m[1]*_m[5] - _m[2]*_m[4];
			tmp[3] = _m[5]*_m[6] - _m[3]*_m[8];
			tmp[4] = _m[0]*_m[8] - _m[2]*_m[6];
			tmp[5] = _m[2]*_m[3] - _m[0]*_m[5];
			tmp[6] = _m[3]*_m[7] - _m[4]*_m[6];
			tmp[7] = _m[1]*_m[6] - _m[0]*_m[7];
			tmp[8] = _m[0]*_m[4] - _m[1]*_m[3];

			real det = determinant ();
			assert (ABS (det) > EPSILON);
			det = 1./ det;
			for (int i = 0; i < 9; ++i){
				_m[i] = det*tmp[i];
			}
		}
	};
}
