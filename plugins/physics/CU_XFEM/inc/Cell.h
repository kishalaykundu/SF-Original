/**
 * @file Cell.h
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
 * The cell topology class for the CU_XFEM library.
 */

#pragma once

#include <climits>

 namespace SF {
  namespace XFE {
    const unsigned int CELL_BIT_ARRAY_SIZE = 3;

    class Cell {

      public:
        /**
        * State bits:
        * 0-3   : External Vertex Flags
        * 4-7   : Inside Surface Vertex Flags
        * 8-11  : External Face Flags
        * 12-15 : Internal Face Flags
        * 16-19 : Vertex Collision Flags
        * 20    : Flag to signify that vertex has collided
        * 21    : Flag to signify that edge has collided
        * 22    : Cell Examined Flag
        * 23    : Cell finalized flag - completely severred
        */
        unsigned char _bitset [CELL_BIT_ARRAY_SIZE];

        // index into the cut vector
        int _cutIndex;

        // indices to vertices
        unsigned int _index [4];
        /**
        * Storage convention:
        * 0 - face 0,1,2
        * 1 - face 0,2,3
        * 2 - face 0,3,1
        * 3 - face 1,3,2
        */
        int _neighbor [4];  // contains cell indices of neighbors

        /**
        * Storage convention:
        * 0 - edge 0,1
        * 1 - edge 0,2
        * 2 - edge 0,3
        * 3 - edge 1,2
        * 4 - edge 1,3
        * 5 - edge 2,3
        */
        unsigned int _edgeIndex [6];

      public:
        inline Cell ()
        : _cutIndex (-1)
        {
          for (unsigned int i = 0; i < CELL_BIT_ARRAY_SIZE; ++i){
            _bitset [i] = 0x00;
          }
          for (unsigned int i = 0; i < 4; ++i){
            _index [i] = UINT_MAX;
            _neighbor [i] = -1;
          }
          for (unsigned int i = 0; i < 6; ++i){
            _edgeIndex [i] = UINT_MAX;
          }
        }

        ~Cell () { }

        // overloaded constructor
        inline Cell (const unsigned int ind [])
        : _cutIndex (-1)
        {
          for (unsigned int i = 0; i < CELL_BIT_ARRAY_SIZE; ++i){
            _bitset [i] = 0x00;
          }
          memcpy (_index, ind, 4*sizeof (unsigned int));
        }

        // copy constructor
        inline Cell (const Cell &c)
        : _cutIndex (c._cutIndex)
        {
          memcpy (_bitset, c._bitset, CELL_BIT_ARRAY_SIZE*sizeof (unsigned char));
          memcpy (_index, c._index, 4*sizeof (unsigned int));
          memcpy (_neighbor, c._neighbor, 4*sizeof (int));
          memcpy (_edgeIndex, c._edgeIndex, 6*sizeof (unsigned int));
        }
        // assignment operator
        inline Cell& operator = (const Cell &c)
        {
          memcpy (_bitset, c._bitset, CELL_BIT_ARRAY_SIZE*sizeof (unsigned char));
          _cutIndex = c._cutIndex;
          memcpy (_index, c._index, 4*sizeof (unsigned int));
          memcpy (_neighbor, c._neighbor, 4*sizeof (int));
          memcpy (_edgeIndex, c._edgeIndex, 6*sizeof (unsigned int));
          return *this;
        }

        // method to add neighbors
        inline void addNeighbors (const int ind [])
        {
          memcpy (_neighbor, ind, 4*sizeof (int));
        }

        // external vertex flag functions
        inline void setExternalVertexFlag (unsigned int index)
        {
          assert (index < 4);
          _bitset [0] |= 0x01 << index;
        }
        inline unsigned int numExternalVertexBits () const
        {
          return (1 && (_bitset [0] & 0x01)) + (1 && (_bitset [0] & 0x02)) + (1 && (_bitset [0] & 0x04)) + (1 && (_bitset [0] & 0x08));
        }
        inline bool testExternalVertexFlag (unsigned int index)
        {
          assert (index < 4);
          return _bitset [0] & (0x01 << index);
        }

        // internal vertex flag functions
        inline void setInternalVertexFlag (unsigned int index)
        {
          assert (index < 4);
          _bitset [0] |= 0x10 << index;
        }
        inline unsigned int numInternalVertexBits () const
        {
          return (1 && (_bitset [0] & 0x10)) + (1 && (_bitset [0] & 0x20)) + (1 && (_bitset [0] & 0x40)) + (1 && (_bitset [0] & 0x80));
        }

        // external face flag functions
        inline void setExternalFaceFlag (unsigned int index)
        {
          assert (index < 4);
          _bitset [1] |= 0x01 << index;
        }
        inline bool testExternalFaceFlag (unsigned int index) const
        {
          assert (index < 4);
          return _bitset [1] & (0x01 << index);
        }
        inline bool testAnyExternalFaceFlag () const
        {
          return (_bitset [1] & 0x01) || (_bitset [1] & 0x02) || (_bitset [1] & 0x04) || (_bitset [1] & 0x08);
        }
        inline unsigned int numExternalFaceBits () const
        {
          return (1 && (_bitset [1] & 0x01)) + (1 && (_bitset [1] & 0x02)) + (1 && (_bitset [1] & 0x04)) + (1 && (_bitset [1] & 0x08));
        }

        // internal face flag functions
        inline void setInternalFaceFlag (unsigned int index)
        {
          assert (index < 4);
          _bitset [1] |= 0x10 << index;
        }
        inline unsigned int numInternalFaceBits () const
        {
          return (1 && (_bitset [1] & 0x10)) + (1 && (_bitset [1] & 0x20)) + (1 && (_bitset [1] & 0x40)) + (1 && (_bitset [1] & 0x80));
        }

        // vertex collision flag functions
        inline void setVertexCollisionFlag (unsigned int index)
        {
          assert (index < 4);
          _bitset [2] |= 0x01 << index;
          _bitset [2] |= 0x10;
        }
        inline bool testVertexCollisionFlag (unsigned int index) const
        {
          assert (index < 4);
          return _bitset [2] & (0x01 << index);
        }
        inline unsigned int numVertexCollisionBits () const
        {
          return (1 && (_bitset [2] & 0x01)) + (1 && (_bitset [2] & 0x02)) + (1 && (_bitset [2] & 0x04)) + (1 && (_bitset [2] & 0x08));
        }

        inline void setThisVertexCollisionFlag (unsigned int index)
        {
          for (unsigned int i = 0; i < 4; ++i){
            if (_index [i] == index){
              _bitset [2] |= 0x01 << i;
              break;
            }
          }
          _bitset [2] |= 0x10;
        }

        // edge collision flag
        inline void setEdgeCollisionFlag ()
        {
          _bitset [2] |= 0x10;
        }

        // test if any collision has occured
        inline bool testAnyCollisionFlag () const
        {
          return (_bitset [2] & 0x10) || (_bitset [2] & 0x20);
        }

        // cell exam flag functions
        inline void setCellExamFlag ()
        {
          _bitset [2] |= 0x40;
        }
        inline bool testCellExamFlag () const
        {
          return (_bitset [2] & 0x80) || (_bitset [2] & 0x40);
        }

        // set flag to certify that cell has been finalized
        inline void finalize ()
        {
          _bitset[2] |= 0x80;
        }
        inline bool testCellFinalizeFlag () const
        {
          return _bitset [2] & 0x80;
        }

        // reset collision related flags
        inline void reset ()
        {
          _bitset [2] &= 0x80;
        }
    };
  }
 }
