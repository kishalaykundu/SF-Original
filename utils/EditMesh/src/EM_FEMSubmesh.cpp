/**
 * @file EM_FEMSubmesh.cpp
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
 * Finite Element sub-mesh for Edit Mesh application. Derived from Mesh
 */
#include <cstdlib>
#include <climits>
#include <vector>
#include <list>

#include "Preprocess.h"
#include "crc32.h"
#include "aabb.h"
#include "EM_common.h"
#include "EM_FEMSubmesh.h"

using namespace std;

namespace SF {

  static bool isthere;
  static bool indices [4];
  static bool defaults [4] = {false, false, false, false};

  // static function to test if triangle belongs to cell
  static inline bool
  faceBelongsToCell (int *f, int *c, unsigned int &ind)
  {
    memcpy (indices, defaults, 4*sizeof (bool));
    for (unsigned int i = 0; i < 3; ++i){
      isthere = false;
      for (unsigned int j = 0; j < 4; ++j){
        if (f [i] == c [j]){
          indices [j] = true;
          isthere = true;
        }
      }
      if (!isthere){
        return false;
      }
    }
    ind = 0;
    if (!indices [0]){
      ind = 3;
    } else if (!indices [1]){
      ind = 1;
    }
    else if (!indices [2]){
      ind = 2;
    }
    return true;
  }

	// default constructor
	FEMSubmesh::FEMSubmesh () { }

	// overloaded constructor
	FEMSubmesh::FEMSubmesh (const vec &min, const vec &max)
	: _bbox (min, max)
	{ }

	// destructor
	FEMSubmesh::~FEMSubmesh () { }

	// method to generate topological info
	void
	FEMSubmesh::generateTopology (const vector <int> &faces)
	{
		// generate cell topology
		generateCellTopology ();

		// generate edge topology
		generateEdgeTopology ();

		// separate faces into surface and internal types
		{
			vector <bool> flags;
			flags.reserve (_efaces.size ()/ 3);
			flags.push_back (false);
			flags.resize (_efaces.size ()/ 3, false);

			for (size_t i = 0; i < _efaces.size (); i +=3){
				for (size_t j = 0; j < faces.size (); j += 3){
					if (sameTriangle (&(_efaces[i]), &(faces[j]))){
						flags.at (i/ 3) = true;
					}
				}
			}

			vector <int> tvecs;
			for (size_t i = 0; i < _efaces.size (); i += 3){
				if (!flags.at (i/ 3)){
					for (int j = 0; j < 3; ++j){
						_ifaces.push_back( _efaces.at (i + j));
					}
				} else {
					for (int j = 0; j < 3; ++j){
						tvecs.push_back( _efaces.at (i + j));
					}
				}
			}
			_efaces.swap (tvecs);
		}

		// initialize face topology
		generateFaceTopology (_eftop, _efaces);
		generateFaceTopology (_iftop, _ifaces);

		// initialize face owner structure
		generateFaceOwnerInfo ();
	}

	// protected method to generate face owner info
	void
	FEMSubmesh::generateFaceOwnerInfo ()
	{
	  _efown.resize (_eftop.size ());

    unsigned int index;
	  for (unsigned int i = 0; i < _efown.size (); ++i){
      for (unsigned int j = 0; j < _cells.size ()/ 4; ++j){
        if (faceBelongsToCell (&(_efaces [3*i]), &(_cells [4*j]), index)){
          _efown [i]._owner = j;
          _efown [i]._index = index;
          break;
        }
      }
	  }

	  _ifown.resize (_iftop.size ());
	  for (unsigned int i = 0; i < _ifown.size (); ++i){
      for (unsigned int j = 0; j < _cells.size ()/ 4; ++j){
        if (faceBelongsToCell (&(_ifaces [3*i]), &(_cells [4*j]), index)){
          _ifown [i]._owner = j;
          _ifown [i]._index = index;
          break;
        }
      }
	  }
	}

	// protected method to generate edge-related topological info
	void
	FEMSubmesh::generateEdgeTopology ()
	{
		vector <list <Edge> > edges;
		list <Edge> *eptr;
		list <Edge>::iterator iter;

		edges.reserve (USHRT_MAX + 1);
		edges.push_back (list <Edge> ());
		edges.resize (USHRT_MAX + 1, list <Edge> ());

		bool inflag = true;
		uint32_t index;
		int inds[2];
		char tmpstr[16], finalstr[32];

		for (size_t i = 0; i < _cells.size (); i += 4){
			for (int j = 0; j < 6; ++j){

				switch (j){
					case 0:
						inds[0] = _cells.at (i);
						inds[1] = _cells.at (i + 1);
						break;
					case 1:
						inds[0] = _cells.at (i);
						inds[1] = _cells.at (i + 2);
						break;
					case 2:
						inds[0] = _cells.at (i);
						inds[1] = _cells.at (i + 3);
						break;
					case 3:
						inds[0] = _cells.at (i + 1);
						inds[1] = _cells.at (i + 2);
						break;
					case 4:
						inds[0] = _cells.at (i + 1);
						inds[1] = _cells.at (i + 3);
						break;
					case 5:
						inds[0] = _cells.at (i + 2);
						inds[1] = _cells.at (i + 3);
						break;
				}

				// sort the indices
				if (inds[0] > inds[1]){
					int tmpi = inds[0];
					inds[0] = inds[1];
					inds[1] = tmpi;
				}

				// form string by concatenating three indices
				sprintf (finalstr, "%d", inds[0]);
				sprintf (tmpstr, "%d", inds[1]);
				assert (strlen (finalstr) + strlen (tmpstr) < 32);
				strcat (finalstr, tmpstr);

				// get crc32 hash-code and mask it to below 65536
				index = crc32 (finalstr, strlen (finalstr));
				index &= 0x0000ffff;

				inflag = false;
				eptr = &(edges.at (index));
				iter = eptr->begin ();

				while (iter != eptr->end ()){
					if (*iter == &(inds[0])){
						inflag = true;
						iter->add (i/ 4);
						break;
					}
					++iter;
				}
				if (!inflag){
					eptr->push_back( Edge (i/ 4, &(inds[0])));
				}
			}
		}
		// copy edges to edge-vector
		for (int i = 0; i <= USHRT_MAX; ++i){

			eptr = &(edges.at (i));

			while (!eptr->empty ()){
				_edges.push_back (eptr->front ());
				eptr->pop_front ();
			}
		}
	}

	// protected method to generate cell-related topological info
	void
	FEMSubmesh::generateCellTopology ()
	{
		_ctop.reserve (_cells.size ()/ 4);
		_ctop.push_back (Cell ());
		_ctop.resize (_cells.size ()/ 4, Cell ());

		vector < list < Trig > > faces;
		list < Trig >* fptr;
		list < Trig >::iterator iter;
		faces.reserve (USHRT_MAX + 1);
		faces.push_back (list < Trig > ());
		faces.resize (USHRT_MAX + 1, list < Trig > ());

		bool inflag = false;
		uint32_t index;
		int sorted_inds[3], orig_inds[3];
		char tmpstr1[16], tmpstr2[16], finalstr[48];

		for (size_t i = 0; i < _ctop.size (); ++i){
			for (int j = 0; j < 4; ++j){

				// make copy of 3 face indices for each cell
				switch (j){
					case 0: // 0-1-2
						orig_inds[0] = sorted_inds[0] = _cells.at (4*i);
						orig_inds[1] = sorted_inds[1] = _cells.at (4*i + 1);
						orig_inds[2] = sorted_inds[2] = _cells.at (4*i + 2);
						break;
					case 1: // 0-2-3
						orig_inds[0] = sorted_inds[0] = _cells.at (4*i);
						orig_inds[1] = sorted_inds[1] = _cells.at (4*i + 2);
						orig_inds[2] = sorted_inds[2] = _cells.at (4*i + 3);
						break;
					case 2: // 0-3-1
						orig_inds[0] = sorted_inds[0] = _cells.at (4*i);
						orig_inds[1] = sorted_inds[1] = _cells.at (4*i + 3);
						orig_inds[2] = sorted_inds[2] = _cells.at (4*i + 1);
						break;
					case 3: // 1-3-2
						orig_inds[0] = sorted_inds[0] = _cells.at (4*i + 1);
						orig_inds[1] = sorted_inds[1] = _cells.at (4*i + 3);
						orig_inds[2] = sorted_inds[2] = _cells.at (4*i + 2);
				}

				qsort (&(sorted_inds[0]), 3, sizeof (int), compare);

				// form string by concatenating three indices
				sprintf (finalstr, "%d", sorted_inds[0]);
				sprintf (tmpstr1, "%d", sorted_inds[1]);
				sprintf (tmpstr2, "%d", sorted_inds[2]);
				assert (strlen (finalstr) + strlen (tmpstr1) + strlen (tmpstr2) < 48);
				strcat (finalstr, tmpstr1);
				strcat (finalstr, tmpstr2);

				// get crc32 hash-code and mask it to below 65536
				index = crc32 (finalstr, strlen (finalstr));
				index &= 0x0000ffff;

				inflag = false;
				fptr = &(faces.at (index));
				iter = fptr->begin ();

				// check for existing entry
				while (iter != fptr->end ()){
					if (*iter == &(orig_inds[0])){
						inflag = true;
						_ctop.at (i)._neighbors[j] = iter->_owner;
						_ctop.at (iter->_owner)._neighbors[iter->_face] = i;

						fptr->erase (iter);
						break;
					}
					++iter;
				}

				// new entry
				if (!inflag){
					fptr->push_back (Trig (i, j, &(orig_inds[0])));
				}
			} // end - for (int j = 0; j < 4; ++j)
		} // end - for (size_t i = 0; i < cells.size (); ++i)

		// faces now contains list of surface triangles
		for (size_t i = 0; i < faces.size (); ++i){

			fptr = &(faces.at (i));

			while (!fptr->empty ()){
				for (int j = 0; j < 3; ++j){
					_efaces.push_back (fptr->front ()._indices[j]);
				}
				fptr->pop_front ();
			}
		}
	}

}
