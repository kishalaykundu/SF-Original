/**
 * @file vec2.cpp
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

#include "vec2.h"

namespace SF {
	const vec2 vec2::ZERO (0., 0.);
	const vec2 vec2::UNIT (1., 1.);

	const vec2 vec2::UNIT_X (1., 0.);
	const vec2 vec2::UNIT_Y (0., 1.);

	const vec2 vec2::NEG_UNIT_X (-1., 0.);
	const vec2 vec2::NEG_UNIT_Y (0., -1.);

}
