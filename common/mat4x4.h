/**
 * @file mat4x4.h
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
 * The 4 by 4 matrix class for Simulate Framework
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#include "Preprocess.h"
#include "vec4.h"

namespace SF {

	class mat4x4 {

	public:
		real _m[16];

		// methods described later
	public:
		static const mat4x4 ZERO;
		static const mat4x4 ONE;
		static const mat4x4 IDENTITY;

	public:
		inline mat4x4()
		{
			for (int i = 0; i < 15; ++i){
				_m[i] = 0.;
			}
			_m[15] = 1.;
		}
		// overloaded constructor
		inline mat4x4 (const real* &arr)
		{
			memcpy (_m, arr, 16*sizeof (real));
		}
		// overloaded constructor
		inline mat4x4 (const mat4x4 &mat)
		{
			memcpy (_m, mat._m, 16*sizeof (real));
		}
		// overloaded constructor
		inline mat4x4 (const real a00, const real a01, const real a02, const real a03,
				const real a10, const real a11, const real a12, const real a13,
				const real a20, const real a21, const real a22, const real a23,
				const real a30, const real a31, const real a32, const real a33)
		{
			_m[0] = a00; _m[1] = a01; _m[2] = a02; _m[3] = a03;
			_m[4] = a10; _m[5] = a11; _m[6] = a12; _m[7] = a13;
			_m[8] = a20; _m[9] = a21; _m[10] = a22; _m[11] = a23;
			_m[12] = a30; _m[13] = a31; _m[14] = a32; _m[15] = a33;
		}
		// overloaded constructor
		inline mat4x4 (const real a00, const real a01, const real a02,
				const real a10, const real a11, const real a12,
				const real a20, const real a21, const real a22)
		{
			_m[0] = a00; _m[1] = a01; _m[2] = a02; _m[3] = 0.;
			_m[4] = a10; _m[5] = a11; _m[6] = a12; _m[7] = 0.;
			_m[8] = a20; _m[9] = a21; _m[10] = a22; _m[11] = 0.;
			_m[12] = 0.0; _m[13] = 0.0; _m[14] = 0.0; _m[15] = 1.;
		}

		// assignment operator
		inline mat4x4 &operator = (const mat4x4 &mat)
		{
			memcpy (_m, mat._m, 16*sizeof (real));
			return *this;
		}

		// element access and mutation operator
		inline real &operator() (const int i, const int j)
		{
			assert (i >= 0 && i < 4);
			assert (j >= 0 && j < 4);
			return _m[4*i + j];
		}

		// non-mutating addition operator
		inline mat4x4 operator + (const mat4x4 &mat) const
		{
			return mat4x4(
					_m[0]+mat._m[0], _m[1]+mat._m[1], _m[2]+mat._m[2], _m[3]+mat._m[3],
					_m[4]+mat._m[4], _m[5]+mat._m[5], _m[6]+mat._m[6], _m[7]+mat._m[7],
					_m[8]+mat._m[8], _m[9]+mat._m[9], _m[10]+mat._m[10], _m[11]+mat._m[11],
					_m[12]+mat._m[12], _m[13]+mat._m[13], _m[14]+mat._m[14], _m[15]+mat._m[15]);
		}
		// non-mutating subtraction operator
		inline mat4x4 operator - (const mat4x4 &mat) const
		{
			return mat4x4(
					_m[0]-mat._m[0], _m[1]-mat._m[1], _m[2]-mat._m[2], _m[3]-mat._m[3],
					_m[4]-mat._m[4], _m[5]-mat._m[5], _m[6]-mat._m[6], _m[7]-mat._m[7],
					_m[8]-mat._m[8], _m[9]-mat._m[9], _m[10]-mat._m[10], _m[11]-mat._m[11],
					_m[12]-mat._m[12], _m[13]-mat._m[13], _m[14]-mat._m[14], _m[15]-mat._m[15]);
		}
		// non-mutating per-element multiplication operator
		inline mat4x4 operator * (const mat4x4 &mat) const
		{
			return mat4x4(
					_m[0]*mat._m[0] + _m[1]*mat._m[4] + _m[2]*mat._m[8] + _m[3]*mat._m[12],
					_m[0]*mat._m[1] + _m[1]*mat._m[5] + _m[2]*mat._m[9] + _m[3]*mat._m[13],
					_m[0]*mat._m[2] + _m[1]*mat._m[6] + _m[2]*mat._m[10] + _m[3]*mat._m[14],
					_m[0]*mat._m[3] + _m[1]*mat._m[7] + _m[2]*mat._m[11] + _m[3]*mat._m[15],
					_m[4]*mat._m[0] + _m[5]*mat._m[4] + _m[6]*mat._m[8] + _m[7]*mat._m[12],
					_m[4]*mat._m[1] + _m[5]*mat._m[5] + _m[6]*mat._m[9] + _m[7]*mat._m[13],
					_m[4]*mat._m[2] + _m[5]*mat._m[6] + _m[6]*mat._m[10] + _m[7]*mat._m[14],
					_m[4]*mat._m[3] + _m[5]*mat._m[7] + _m[6]*mat._m[11] + _m[7]*mat._m[15],
					_m[8]*mat._m[0] + _m[9]*mat._m[4] + _m[10]*mat._m[8] + _m[11]*mat._m[12],
					_m[8]*mat._m[1] + _m[9]*mat._m[5] + _m[10]*mat._m[9] + _m[11]*mat._m[13],
					_m[8]*mat._m[2] + _m[9]*mat._m[6] + _m[10]*mat._m[10] + _m[11]*mat._m[14],
					_m[8]*mat._m[3] + _m[9]*mat._m[7] + _m[10]*mat._m[11] + _m[11]*mat._m[15],
					_m[12]*mat._m[0] + _m[13]*mat._m[4] + _m[14]*mat._m[8] + _m[15]*mat._m[12],
					_m[12]*mat._m[1] + _m[13]*mat._m[5] + _m[14]*mat._m[9] + _m[15]*mat._m[13],
					_m[12]*mat._m[2] + _m[13]*mat._m[6] + _m[14]*mat._m[10] + _m[15]*mat._m[14],
					_m[12]*mat._m[3] + _m[13]*mat._m[7] + _m[14]*mat._m[11] + _m[15]*mat._m[15]);
		}
		// vector multiplication
		inline vec4 operator * (vec4 &v)
		{
			real w = v._v[0]*_m[12] + v._v[1]*_m[13] + v._v[2]*_m[14] + v._v[3]*_m[15];
			assert (ABS (w) > EPSILON);
			w = 1./w;
			return vec4(
					 (v._v[0]*_m[0] + v._v[1]*_m[1] + v._v[2]*_m[2] + v._v[3]*_m[3]) * w,
					 (v._v[0]*_m[4] + v._v[1]*_m[5] + v._v[2]*_m[6] + v._v[3]*_m[7]) * w,
					 (v._v[0]*_m[8] + v._v[1]*_m[9] + v._v[2]*_m[10] + v._v[3]*_m[11]) * w,
					1.);
		}
 		// scaling operator
		inline mat4x4 operator * (const real val) const
		{
			return mat4x4(
					_m[0]*val, _m[1]*val, _m[2]*val, _m[3]*val,
					_m[4]*val, _m[5]*val, _m[6]*val, _m[7]*val,
					_m[8]*val, _m[9]*val, _m[10]*val, _m[11]*val,
					_m[12]*val, _m[13]*val, _m[14]*val, _m[15]*val);
		}
 		// scaling operator
		inline mat4x4 operator / (const real val) const
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1./ val;
			return mat4x4(
					_m[0]*reciprocal, _m[1]*reciprocal, _m[2]*reciprocal, _m[3]*reciprocal,
					_m[4]*reciprocal, _m[5]*reciprocal, _m[6]*reciprocal, _m[7]*reciprocal,
					_m[8]*reciprocal, _m[9]*reciprocal, _m[10]*reciprocal, _m[11]*reciprocal,
					_m[12]*reciprocal, _m[13]*reciprocal, _m[14]*reciprocal, _m[15]*reciprocal);
		}

		// mutating addition operator
		inline mat4x4 &operator += (const mat4x4 &mat)
		{
			for (int i = 0; i < 16; ++i){
				_m[i] += mat._m[i];
			}
			return *this;
		}
		// mutating addition operator
		inline mat4x4 &operator += (const real val)
		{
			for (int i = 0; i < 16; ++i){
				_m[i] += val;
			}
			return *this;
		}
		// mutating subtraction operator
		inline mat4x4 &operator -= (const mat4x4 &mat)
		{
			for (int i = 0; i < 16; ++i){
				_m[i] -= mat._m[i];
			}
			return *this;
		}
		// mutating subtraction operator
		inline mat4x4 &operator -= (const real val)
		{
			for (int i = 0; i < 16; ++i){
				_m[i] -= val;
			}
			return *this;
		}
		// mutating multiplication operator
		inline mat4x4 &operator *= (const mat4x4 &mat)
		{
			real tmp[16];
			tmp[0] = _m[0]*mat._m[0] + _m[1]*mat._m[4] + _m[2]*mat._m[8] + _m[3]*mat._m[12];
			tmp[1] = _m[0]*mat._m[1] + _m[1]*mat._m[5] + _m[2]*mat._m[9] + _m[3]*mat._m[13];
			tmp[2] = _m[0]*mat._m[2] + _m[1]*mat._m[6] + _m[2]*mat._m[10] + _m[3]*mat._m[14];
			tmp[3] = _m[0]*mat._m[3] + _m[1]*mat._m[7] + _m[2]*mat._m[11] + _m[3]*mat._m[15];
			tmp[4] = _m[4]*mat._m[0] + _m[5]*mat._m[4] + _m[6]*mat._m[8] + _m[7]*mat._m[12];
			tmp[5] = _m[4]*mat._m[1] + _m[5]*mat._m[5] + _m[6]*mat._m[9] + _m[7]*mat._m[13];
			tmp[6] = _m[4]*mat._m[2] + _m[5]*mat._m[6] + _m[6]*mat._m[10] + _m[7]*mat._m[14];
			tmp[7] = _m[4]*mat._m[3] + _m[5]*mat._m[7] + _m[6]*mat._m[11] + _m[7]*mat._m[15];
			tmp[8] = _m[8]*mat._m[0] + _m[9]*mat._m[4] + _m[10]*mat._m[8] + _m[11]*mat._m[12];
			tmp[9] = _m[8]*mat._m[1] + _m[9]*mat._m[5] + _m[10]*mat._m[9] + _m[11]*mat._m[13];
			tmp[10] = _m[8]*mat._m[2] + _m[9]*mat._m[6] + _m[10]*mat._m[10] + _m[11]*mat._m[14];
			tmp[11] = _m[8]*mat._m[3] + _m[9]*mat._m[7] + _m[10]*mat._m[11] + _m[11]*mat._m[15];
			tmp[12] = _m[12]*mat._m[0] + _m[13]*mat._m[4] + _m[14]*mat._m[8] + _m[15]*mat._m[12];
			tmp[13] = _m[12]*mat._m[1] + _m[13]*mat._m[5] + _m[14]*mat._m[9] + _m[15]*mat._m[13];
			tmp[14] = _m[12]*mat._m[2] + _m[13]*mat._m[6] + _m[14]*mat._m[10] + _m[15]*mat._m[14];
			tmp[15] = _m[12]*mat._m[3] + _m[13]*mat._m[7] + _m[14]*mat._m[11] + _m[15]*mat._m[15];
			memcpy (_m, tmp, 16*sizeof (real));
			return *this;
		}
		// mutating scaling operator
		inline mat4x4 &operator *= (const real val)
		{
			for (int i = 0; i < 16; ++i){
				_m[i] *= val;
			}
			return *this;
		}
		// mutating scaling operator
		inline mat4x4 &operator /= (const real val)
		{
			assert (ABS (val) > EPSILON);
			real reciprocal = 1./ val;
			for (int i = 0; i < 16; ++i){
				_m[i] *=  reciprocal;
			}
			return *this;
		}

		// transpose operator
		inline void transpose ()
		{
			real tmp = _m[1]; _m[1] = _m[4]; _m[4] = tmp;
			tmp = _m[2]; _m[2] = _m[8]; _m[8] = tmp;
			tmp = _m[3]; _m[3] = _m[12]; _m[12] = tmp;
			tmp = _m[6]; _m[6] = _m[9]; _m[9] = tmp;
			tmp = _m[7]; _m[7] = _m[13]; _m[13] = tmp;
			tmp = _m[11]; _m[11] = _m[14]; _m[14] = tmp;
		}

		// determinant
		inline real determinant ()
		{
			return
			_m[3] * _m[6] * _m[9] * _m[12] - _m[2] * _m[7] * _m[9] * _m[12] -
			_m[3] * _m[5] * _m[10] * _m[12]+_m[1] * _m[7] * _m[10] * _m[12] +
			_m[2] * _m[5] * _m[11] * _m[12]-_m[1] * _m[6] * _m[11] * _m[12] -
			_m[3] * _m[6] * _m[8] * _m[13]+_m[2] * _m[7] * _m[8] * _m[13] +
			_m[3] * _m[4] * _m[10] * _m[13]-_m[0] * _m[7] * _m[10] * _m[13] -
			_m[2] * _m[4] * _m[11] * _m[13]+_m[0] * _m[6] * _m[11] * _m[13] +
			_m[3] * _m[5] * _m[8] * _m[14]-_m[1] * _m[7] * _m[8] * _m[14] -
			_m[3] * _m[4] * _m[9] * _m[14]+_m[0] * _m[7] * _m[9] * _m[14] +
			_m[1] * _m[4] * _m[11] * _m[14]-_m[0] * _m[5] * _m[11] * _m[14] -
			_m[2] * _m[5] * _m[8] * _m[15]+_m[1] * _m[6] * _m[8] * _m[15] +
			_m[2] * _m[4] * _m[9] * _m[15]-_m[0] * _m[6] * _m[9] * _m[15] -
			_m[1] * _m[4] * _m[10] * _m[15]+_m[0] * _m[5] * _m[10] * _m[15];
		}

		// inverse operator
		inline void invert()
		{
			real tmp[16];
			tmp[0] = _m[6]*_m[11]*_m[13] - _m[7]*_m[10]*_m[13] + _m[7]*_m[9]*_m[14] - _m[5]*_m[11]*_m[14] - _m[6]*_m[9]*_m[15] + _m[5]*_m[10]*_m[15];
			tmp[1] = _m[3]*_m[10]*_m[13] - _m[2]*_m[11]*_m[13] - _m[3]*_m[9]*_m[14] + _m[1]*_m[11]*_m[14] + _m[2]*_m[9]*_m[15] - _m[1]*_m[10]*_m[15];
			tmp[2] = _m[2]*_m[7]*_m[13] - _m[3]*_m[6]*_m[13] + _m[3]*_m[5]*_m[14] - _m[1]*_m[7]*_m[14] - _m[2]*_m[5]*_m[15] + _m[1]*_m[6]*_m[15];
			tmp[3] = _m[3]*_m[6]*_m[9] - _m[2]*_m[7]*_m[9] - _m[3]*_m[5]*_m[10] + _m[1]*_m[7]*_m[10] + _m[2]*_m[5]*_m[11] - _m[1]*_m[6]*_m[11];
			tmp[4] = _m[7]*_m[10]*_m[12] - _m[6]*_m[11]*_m[12] - _m[7]*_m[8]*_m[14] + _m[4]*_m[11]*_m[14] + _m[6]*_m[8]*_m[15] - _m[4]*_m[10]*_m[15];
			tmp[5] = _m[2]*_m[11]*_m[12] - _m[3]*_m[10]*_m[12] + _m[3]*_m[8]*_m[14] - _m[0]*_m[11]*_m[14] - _m[2]*_m[8]*_m[15] + _m[0]*_m[10]*_m[15];
			tmp[6] = _m[3]*_m[6]*_m[12] - _m[2]*_m[7]*_m[12] - _m[3]*_m[4]*_m[14] + _m[0]*_m[7]*_m[14] + _m[2]*_m[4]*_m[15] - _m[0]*_m[6]*_m[15];
			tmp[7] = _m[2]*_m[7]*_m[8] - _m[3]*_m[6]*_m[8] + _m[3]*_m[4]*_m[10] - _m[0]*_m[7]*_m[10] - _m[2]*_m[4]*_m[11] + _m[0]*_m[6]*_m[11];
			tmp[8] = _m[5]*_m[11]*_m[12] - _m[7]*_m[9]*_m[12] + _m[7]*_m[8]*_m[13] - _m[4]*_m[11]*_m[13] - _m[5]*_m[8]*_m[15] + _m[4]*_m[9]*_m[15];
			tmp[9] = _m[3]*_m[9]*_m[12] - _m[1]*_m[11]*_m[12] - _m[3]*_m[8]*_m[13] + _m[0]*_m[11]*_m[13] + _m[1]*_m[8]*_m[15] - _m[0]*_m[9]*_m[15];
			tmp[10] = _m[1]*_m[7]*_m[12] - _m[3]*_m[5]*_m[12] + _m[3]*_m[4]*_m[13] - _m[0]*_m[7]*_m[13] - _m[1]*_m[4]*_m[15] + _m[0]*_m[5]*_m[15];
			tmp[11] = _m[3]*_m[5]*_m[8] - _m[1]*_m[7]*_m[8] - _m[3]*_m[4]*_m[9] + _m[0]*_m[7]*_m[9] + _m[1]*_m[4]*_m[11] - _m[0]*_m[5]*_m[11];
			tmp[12] = _m[6]*_m[9]*_m[12] - _m[5]*_m[10]*_m[12] - _m[6]*_m[8]*_m[13] + _m[4]*_m[10]*_m[13] + _m[5]*_m[8]*_m[14] - _m[4]*_m[9]*_m[14];
			tmp[13] = _m[1]*_m[10]*_m[12] - _m[2]*_m[9]*_m[12] + _m[2]*_m[8]*_m[13] - _m[0]*_m[10]*_m[13] - _m[1]*_m[8]*_m[14] + _m[0]*_m[9]*_m[14];
			tmp[14] = _m[2]*_m[5]*_m[12] - _m[1]*_m[6]*_m[12] - _m[2]*_m[4]*_m[13] + _m[0]*_m[6]*_m[13] + _m[1]*_m[4]*_m[14] - _m[0]*_m[5]*_m[14];
			tmp[15] = _m[1]*_m[6]*_m[8] - _m[2]*_m[5]*_m[8] + _m[2]*_m[4]*_m[9] - _m[0]*_m[6]*_m[9] - _m[1]*_m[4]*_m[10] + _m[0]*_m[5]*_m[10];

			real det = determinant();
			assert (ABS (det) > EPSILON);
			det = 1./ det;
			for (int i = 0; i < 16; ++i){
				_m[i] = det*tmp[i];
			}
		}
	};
}
