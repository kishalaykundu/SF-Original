/**
 * @file EM_MSDMesh.h
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
 * Mass spring mesh for Edit Mesh application. Derived from Mesh
 */

#pragma once

#include <vector>

#include "EM_Mesh.h"

using namespace std;

namespace SF {

	class aabb;

	class MSDMesh : public Mesh {

	private:
		vector <aabb> _bbox;

		vector <vector <int> > _trigs;
		vector <vector <Face> > _ftop;

		vector <int> _edges;

	public:
		MSDMesh ();
		~MSDMesh ();

		void process (const int);

	protected:
		void shuffleVertices (const vector <aabb> &bvs);
		void generateEdgeList ();
		void writeElementsToFiles (const string &folder, const string &prefix) const;
	};
}
