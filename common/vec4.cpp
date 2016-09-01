/**
 * @file vec4.cpp
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

#include "vec4.h"

namespace SF {

	const vec4 vec4::ZERO (0., 0., 0., 1.);
	const vec4 vec4::UNIT (1., 1., 1., 1.);

	const vec4 vec4::UNIT_X (1., 0., 0., 1.);
	const vec4 vec4::UNIT_Y (0., 1., 0., 1.);
	const vec4 vec4::UNIT_Z (0., 0., 1., 1.);

	const vec4 vec4::NEG_UNIT_X (-1., 0., 0., 1.);
	const vec4 vec4::NEG_UNIT_Y (0., -1., 0., 1.);
	const vec4 vec4::NEG_UNIT_Z (0., 0., -1., 1.);

}
