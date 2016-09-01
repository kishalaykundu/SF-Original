/**
 * @file mat3x3.cpp
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

#include "mat3x3.h"

namespace SF {

	const mat3x3 mat3x3::ZERO (0., 0., 0., 0., 0., 0., 0., 0., 0.);
	const mat3x3 mat3x3::ONE (1., 1., 1., 1., 1., 1., 1., 1., 1.);
	const mat3x3 mat3x3::IDENTITY (1., 0., 0., 0., 1., 0., 0., 0., 1.);

}
