/**
 * @file Edge.h
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
 * The edge topology class for the CU_XFEM library.
 */

#pragma once

#include <cstdlib>

#include "Preprocess.h"

namespace SF {
  namespace XFE {

    class Edge {

      public:
        real _u;
        bool _testFlag;
        unsigned int _firstVertex;
        unsigned int _numOwners;
        unsigned int *_owner;

        Edge ()
        : _u (0), _testFlag (false), _firstVertex (0), _numOwners (0), _owner (NULL) { }

        ~Edge () { delete [] _owner; }

        // overloaded constructor
        Edge (unsigned int v, unsigned int n, unsigned int *o)
        : _u (0), _testFlag (false), _firstVertex (v), _numOwners (n)
        {
          _owner = new unsigned int [_numOwners];
          memcpy (_owner, o, _numOwners * sizeof (unsigned int));
        }

        // copy constructor
        Edge (const Edge &e)
        : _u (e._u), _testFlag (e._testFlag), _firstVertex (e._firstVertex), _numOwners (e._numOwners)
        {
          if (_numOwners){
            _owner = new unsigned int [_numOwners];
            memcpy (_owner, e._owner, _numOwners * sizeof (unsigned int));
          }
        }

        // assignment operator
        Edge& operator = (const Edge &e)
        {
          _u = e._u;
          _testFlag = e._testFlag;
          _firstVertex = e._firstVertex;
          if (_numOwners){
            delete [] _owner;
          }
          _numOwners = e._numOwners;
          if (_numOwners){
            _owner = new unsigned int [_numOwners];
            memcpy (_owner, e._owner, _numOwners * sizeof (unsigned int));
          }
          return *this;
        }

        // test if collision has been recorded
        inline void
        setCollisionFlag ()
        {
          _testFlag = true;
        }
        inline bool
        testCollisionFlag () const
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
