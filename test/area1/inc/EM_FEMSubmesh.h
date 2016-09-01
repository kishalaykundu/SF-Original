/**
 * @file EM_FEMSubmesh.h
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

#pragma once

#include <cstdlib>
#include <vector>

#include "Preprocess.h"
#include "vec4.h"
#include "EM_common.h"

using namespace std;

namespace SF {

	typedef struct edge_t {
		int _indices[2];
		int _nOwners;
		int *_owners;

		inline edge_t ()
		: _nOwners (0), _owners (NULL)
		{
			for (int i = 0; i < 2; ++i){
				_indices[i] = -1;
			}
		}

		inline edge_t (const int o, const int *i)
		: _nOwners (1)
		{
			_owners = (int *) malloc (sizeof (int));
			_owners[0] = o;
			memcpy (_indices, i, 2*sizeof (int));
		}

		inline edge_t (const edge_t &e)
		:_nOwners (e._nOwners)
		{
			_owners = (int *) malloc (_nOwners*sizeof (int));
			memcpy (_indices, e._indices, 2*sizeof (int));
			memcpy (_owners, e._owners, _nOwners*sizeof (int));
		}

		inline bool operator == (const int *f) const
		{
			bool flags[2] = {false, false};
			for(int i = 0; i < 2; ++i){
				for (int j = 0; j < 2; ++j){
					flags[i] |= (_indices[i] == f[j]);
				}
			}
			return flags[0] & flags[1];
		}

		inline void add (const int o)
		{
			if (_nOwners){
				_owners = (int *) realloc (_owners, (_nOwners + 1)* sizeof (int));
			} else {
				_owners = (int *) malloc (_nOwners*sizeof (int));
			}
			_owners[_nOwners] = o;
			++_nOwners;
		}

		~edge_t ()
		{
			free (_owners); _owners = NULL;
		}
	} Edge;

	class aabb;

	class FEMSubmesh {

	public:
		aabb _bbox;

		// edges
		vector <Edge> _edges;

		// internal faces
		vector <int> _ifaces;
		vector <Face> _iftop;

		// external faces
		vector <int> _efaces;
		vector <Face> _eftop;

		// cells
		vector <int> _cells;
		vector <Cell> _ctop;

	public:
		FEMSubmesh ();
		~FEMSubmesh ();

		// overloaded constructor
		FEMSubmesh (const vec &min, const vec &max);

		// method to generate topological info
		void generateTopology (const vector <int> &faces);

	protected:
		void generateEdgeTopology ();
		void generateCellTopology ();
	};
}
