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
 * Mass spring mesh for Edit Mesh application
 */

#pragma once

#include <vector>
#include <string>

#include "vec3.h"

using namespace std;

namespace SF {

	class Mesh {

	public:
		vector< vec3 > _vertices;
		vector< int > _cells;
		vector <int> _faces;

	public:
		Mesh () { }
		virtual ~Mesh () { }

		// method to process
		virtual void process (const int) = 0;

		// method to output to files (to be customized by derived classes)
		void writeToFiles (const string &folder, const string &prefix) const
		{
			writeNodesToFile (folder, prefix);
			writeElementsToFiles (folder, prefix);
		}

	protected:
		// method to write geometric elements to files
		virtual void writeElementsToFiles (const string &folder, const string &prefix) const = 0;

		// method to write vertices out .node file
		void writeNodesToFile (const string &folder, const string &prefix) const
		{
			assert (!folder.empty ());
			assert (!prefix.empty ());
			assert (!_vertices.empty ());

			string nodeFile (folder + prefix);
			nodeFile.append (".node");

			FILE* fp = fopen (nodeFile.c_str (), "w");
			assert (fp);

			fprintf (fp, "%lu\n", _vertices.size ());

			for (size_t i = 0; i < _vertices.size (); ++i){
				fprintf (fp, "%g %g %g\n", _vertices.at (i)._v[0], _vertices.at (i)._v[1], _vertices.at (i)._v[2]);
			}

			fclose (fp);
		}
	};
}
