/**
 * @file mat4x4.cpp
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

#include "mat4x4.h"

namespace SF {

	const mat4x4 mat4x4::ZERO (0., 0., 0., 0.,
			0., 0., 0., 0.,
			0., 0., 0., 0.,
			0., 0., 0., 0.);
	const mat4x4 mat4x4::ONE (1., 1., 1., 1.,
			1., 1., 1., 1.,
			1., 1., 1., 1.,
			1., 1., 1., 1.);
	const mat4x4 mat4x4::IDENTITY (1., 0., 0., 0.,
			0., 1., 0., 0.,
			0., 0., 1., 0.,
			0., 0., 0., 1.);

}
