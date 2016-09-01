/**
 * @file EM_FEMMesh.h
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
 * Finite Element mesh for Edit Mesh application. Derived from Mesh
 */

#pragma once

#include <cstdlib>
#include "EM_Mesh.h"

namespace SF {

  typedef struct oinfo_t {
    unsigned int _submesh;
    unsigned int _cellIndex;
  } OwnerInfo;

  typedef struct vinfo_t {
    unsigned int _nOwners;
    OwnerInfo *_owner;

    vinfo_t ()
    : _nOwners (0), _owner (NULL) { }

    ~vinfo_t ()
    {
      if (_nOwners){
        free (_owner);
      }
    }

    void addOwner (unsigned int submesh, unsigned int cellIndex)
    {
      if (_nOwners){
        _owner = (OwnerInfo *) realloc (_owner, (_nOwners + 1)* sizeof (OwnerInfo));
      } else {
        _owner = (OwnerInfo *) malloc (sizeof (OwnerInfo));
      }
      _owner [_nOwners]._submesh = submesh;
      _owner [_nOwners]._cellIndex = cellIndex;
      ++_nOwners;
    }
  } Vertex;

	class FEMSubmesh;

	class FEMMesh : public Mesh {

	protected:
    vector <Vertex> _vertInfo;
		vector <FEMSubmesh> _submesh;

	public:
		FEMMesh ();
		~FEMMesh ();

		void process (const int);

	protected:
		void shuffleVertices ();
		void writeElementsToFiles (const string &folder, const string &prefix) const;
	};
}
