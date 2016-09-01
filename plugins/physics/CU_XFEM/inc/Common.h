/**
 * @file Common.h
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
 * Common functions for the CU_XFEM library.
 */

#pragma once

#include <climits>
#include <string>

using namespace std;

namespace SF {
	namespace XFE {

		typedef struct fBits {
			bool _cbit;
			size_t _cfrom;
			size_t _cto;

			inline fBits ()
			: _cbit (false), _cfrom (UINT_MAX), _cto (0)
			{ }

			inline void reset ()
			{
				_cbit = false;
				_cfrom = UINT_MAX;
				_cto = 0;
			}

		} FaceChangeStruct;

		bool getConfigParameter (const string &cfile, const char *param, string &result);
	}
}
