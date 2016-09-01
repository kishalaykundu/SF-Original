/**
 * @file Vertex.h
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
 * The vertex topology class for the CU_XFEM library.
 */

#pragma once

#include <cstdlib>

#include "Preprocess.h"

namespace SF {
  namespace XFE {

    class Vertex {

      public:
        bool _testFlag;
        unsigned int _numSubmeshes;
        /**
        * Each array contains three distinct numbers:
        * 1. 1st element: Submesh index
        * 2. 2nd element: index of the last element
        * 3. Array of owning cells from the submesh
        */
        unsigned int **_owners;

        Vertex () : _testFlag (false), _numSubmeshes (0), _owners (NULL) { }
        ~Vertex ()
        {
          if (_numSubmeshes){
            for (unsigned int i = 0; i < _numSubmeshes; ++i){
              delete [] _owners [i];
            }
            free (_owners);
          }
        }

        inline void
        allocateSubmeshSpace (unsigned int index, unsigned int size)
        {
          if (!_numSubmeshes){
            _owners = (unsigned int **)malloc (sizeof (unsigned int *));
          } else {
            _owners = (unsigned int **)realloc (_owners, (_numSubmeshes + 1 )* sizeof (unsigned int *));
          }
          _owners [_numSubmeshes] = new unsigned int [size + 2];
          _owners [_numSubmeshes] [0] = index;
          _owners [_numSubmeshes] [1] = 2;

          ++_numSubmeshes;
        }

        inline void
        addOwner (unsigned int index, unsigned int owner)
        {
          for (unsigned int i = 0; i < _numSubmeshes; ++i){
            if (_owners [i][0] == index){
              _owners [i][_owners [i][1]] = owner;
              ++_owners [i][1];
              break;
            }
          }
        }

        inline void
        setCollisionFlag ()
        {
          _testFlag = true;
        }
        inline bool
        testCollisionFlag ()
        {
          return _testFlag;
        }
        inline void
        reset ()
        {
          _testFlag = false;
        }
    };
  }
}
