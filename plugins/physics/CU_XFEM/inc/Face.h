/**
 * @file Face.h
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
 * The face info structure class for the CU_XFEM library.
 */

#pragma once

#include <climits>

 namespace SF {
  namespace XFE {
    class Face {

      public:
        unsigned int _owner; // index of owning cell
        unsigned char _index; // face-index of owning cell (0-012, 1-023, 2-031, 3-032)

      public:
        inline Face ()
        : _owner (UINT_MAX), _index (0xFF) { }

        ~Face () { }

        inline Face (unsigned int o, unsigned char i)
        : _owner (o), _index (i) { }

        inline Face (const Face &f)
        : _owner (f._owner), _index (f._index) { }

        inline Face & operator = (const Face &f)
        {
          memcpy (&_owner, &(f._owner), sizeof (unsigned int));
          memcpy (&_index, &(f._index), sizeof (unsigned char));
          return *this;
        }
    };
  }
 }
