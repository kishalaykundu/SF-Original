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

#include "EM_Mesh.h"

namespace SF {

	class FEMSubmesh;

	class FEMMesh : public Mesh {

	protected:
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
