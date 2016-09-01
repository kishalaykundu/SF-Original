/**
 * @file EM_FEMMesh.cpp
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
#include <climits>

#include <sstream>

#include "Preprocess.h"
#include "vec3.h"
#include "vec4.h"
#include "aabb.h"
#include "EM_FEMSubmesh.h"
#include "EM_FEMMesh.h"

namespace SF {

	// function to decide which submesh a cell belongs to
	// race to 2. (default: submesh the 1st vertex belongs to)
	static inline size_t getCellSubmeshIndex (const vector <FEMSubmesh> &submesh, const vec3 &v1, const vec3 &v2, const vec3 &v3, const vec3 &v4)
	{
		size_t result = submesh.size ();
		for (size_t i = 0; i < submesh.size (); ++i){
			if (submesh.at (i)._bbox.collide (v1)){
				result = i;
				break;
			}
		}

		int counter = 0;
		aabb *bptr;
		for (size_t i = 0; i < submesh.size (); ++i){

			counter = 0;
			bptr = const_cast <aabb*> (&(submesh.at (i)._bbox));
			if (bptr->collide (v1)){
				++counter;
			}
			if (bptr->collide (v2)){
				++counter;
				if (counter > 1){
					result = i;
					break;
				}
			}
			if (bptr->collide (v3)){
				++counter;
				if (counter > 1){
					result = i;
					break;
				}
			}
			if (bptr->collide (v4)){
				++counter;
				if (counter > 1){
					result = i;
					break;
				}
			}
		}
		return result;
	}

	// default constructor
	FEMMesh::FEMMesh () { }

	// destructor
	FEMMesh::~FEMMesh () { }

	// method to process data to break them to a FEM format
	void
	FEMMesh::process (const int depth)
	{
		{
			// calculate bounding box
			vec min (_vertices.at (0));
			vec max (min);

			for (size_t i = 1; i < _vertices.size (); ++i){
				for (int j = 0; j < 3; ++j){
					min._v[j] = min._v[j] > _vertices.at (i)._v[j] ? _vertices.at (i)._v[j] : min._v[j];
				}
				for (int j = 0; j < 3; ++j){
					max._v[j] = max._v[j] < _vertices.at (i)._v[j] ? _vertices.at (i)._v[j] : max._v[j];
				}
			}
			aabb bv (min, max);

			// calculate number of submeshes
			size_t numSubmeshes = 1;
			for (int i = 1; i <= depth; ++i){
				numSubmeshes *= 8;
			}
			_submesh.reserve (numSubmeshes);

			// initialize submeshes
			vec step (max - min);
			size_t factor = 1;
			for (int i = 1; i <= depth; ++i){
				factor *= 2;
			}
			step *= 1./ static_cast <real> (factor);

			vec bmin, bmax;
			for (int i = 0; i <= depth; ++i){
				for (int j = 0; j <= depth; ++j){
					for (int k = 0; k <= depth; ++k){
						bmin = min;
						bmin._v[2] += i * step._v[2];
						bmin._v[1] += j * step._v[1];
						bmin._v[0] += k * step._v[0];
						bmax = bmin + step;
						_submesh.push_back (FEMSubmesh (bmin, bmax));
					}
				}
			}

			shuffleVertices ();

			size_t submeshIndex;
			for (size_t i = 0; i < _cells.size (); i += 4){
				submeshIndex = getCellSubmeshIndex (_submesh,
						_vertices.at (_cells.at (i)), _vertices.at (_cells.at (i+1)), _vertices.at (_cells.at (i+2)), _vertices.at (_cells.at (i+3)));
				assert (submeshIndex < numSubmeshes);
				for (int j = 0; j < 4; ++j){
					_submesh.at (submeshIndex)._cells.push_back (_cells.at (i+j));
				}
			}
			_cells.clear ();
		}

		for (size_t i = 0; i < _submesh.size (); ++i){
			_submesh.at (i).generateTopology (_faces);
		}
	}

	// protected method to shuffle cell vertices
	void
	FEMMesh::shuffleVertices ()
	{
		vector <bool> sflag;
		sflag.reserve (_vertices.size ());
		sflag.push_back (false);
		sflag.resize (_vertices.size (), false);

		for (size_t i = 0; i < _faces.size (); ++i){
			sflag.at (_faces.at (i)) = true;
		}

		// contains  and current counters offsets for each bv
		vector <int> sOffset;
		sOffset.reserve (_submesh.size ());
		sOffset.push_back (0);
		sOffset.resize (_submesh.size (), 0);

		vector <int> iOffset;
		iOffset.reserve (_submesh.size ());
		iOffset.push_back (0);
		iOffset.resize (_submesh.size (), 0);

		for (size_t i = 0; i < _vertices.size (); ++i){
			for (size_t j = 0; j < _submesh.size (); ++j){

				if (_submesh.at (j)._bbox.collide (_vertices.at (i))){
					if (sflag.at (i)){
						++ (sOffset.at (j));
					} else {
						++ (iOffset.at (j));
					}
					break;
				}
			}
		}

		// do cumulative indices
		if (_submesh.size () > 1){
			int final = sOffset.at (sOffset.size () - 1);

			for (size_t i = sOffset.size () - 1; i > 0; --i){
				sOffset.at (i) = sOffset.at (i - 1);
			}
			sOffset.at (0) = 0;
			for (size_t i = 1; i < sOffset.size (); ++i){
				sOffset.at (i) += sOffset.at (i - 1);
			}

			for (size_t i = iOffset.size () - 1; i > 0; --i){
				iOffset.at (i) = iOffset.at (i - 1);
			}
			iOffset.at (0) = 0;
			for (size_t i = 1; i < iOffset.size (); ++i){
				iOffset.at (i) += iOffset.at (i - 1);
			}

			for (size_t i = 0; i < iOffset.size (); ++i){
				iOffset.at (i) += final;
			}
		}
		int addendum = sOffset.at (sOffset.size () - 1);
		for (size_t i = 0; i < iOffset.size (); ++i){
			iOffset.at (i) += addendum;
		}

		vector <int> sCounter;
		sCounter.reserve (_submesh.size ());
		sCounter.push_back (0);
		sCounter.resize (_submesh.size (), 0);

		vector <int> iCounter;
		iCounter.reserve (_submesh.size ());
		iCounter.push_back (0);
		iCounter.resize (_submesh.size (), 0);

		vector <size_t> newIndices;
		newIndices.reserve (_vertices.size ());
		newIndices.push_back (UINT_MAX);
		newIndices.resize (_vertices.size (), UINT_MAX);

		for (size_t i = 0; i < _vertices.size (); ++i){
			for (size_t j = 0; j < _submesh.size (); ++j){

				if (_submesh.at (j)._bbox.collide (_vertices.at (i))){
					if (sflag.at (i)){
						newIndices.at (i) = sOffset.at (j) + sCounter.at (j);
						++ (sCounter.at (j));
					} else {
						newIndices.at (i) = iOffset.at (j) + iCounter.at (j);
						++ (iCounter.at (j));
					}
					break;
				}
			}
		}

		// shuffle actual vertices
		{
			vector <vec3> verts;
			verts.reserve (_vertices.size ());
			verts.push_back (vec3::ZERO);
			verts.resize (_vertices.size(), vec3::ZERO);

			for (size_t i = 0; i < newIndices.size (); ++i){
				verts.at (newIndices.at (i)) = _vertices.at (i);
			}
			_vertices.swap (verts);
		}

		for (size_t i = 0; i < _cells.size (); ++i){
			_cells.at (i) = newIndices.at (_cells.at (i));
		}

		for (size_t i = 0; i < _faces.size (); ++i){
			_faces.at (i) = newIndices.at (_faces.at (i));
		}
	}

	// method to write elements to files
	void
	FEMMesh::writeElementsToFiles (const string &folder, const string &prefix) const
	{
		string fname;
		FILE* fp = NULL;

		// write cells to files
		{
			vector <int> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".tet.ele");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <int>* > (&(_submesh. at(i)._cells));
				fprintf (fp, "%lu\n", cptr->size ()/ 4);

				for (size_t j = 0; j < cptr->size (); j += 4){
					fprintf (fp, "%d %d %d %d\n", cptr->at (j), cptr->at (j+1), cptr->at (j+2), cptr->at (j+3));
				}

				fclose (fp);
			}
		}

		// write cell topology to files
		{
			vector <Cell> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".tet.top");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <Cell>* > (&(_submesh. at(i)._ctop));
				fprintf (fp, "%lu\n", cptr->size ());

				for (size_t j = 0; j < cptr->size (); ++j){
					fprintf (fp, "%d %d %d %d\n", cptr->at (j)._neighbors[0], cptr->at (j)._neighbors[1], cptr->at (j)._neighbors[2], cptr->at (j)._neighbors[3]);
				}

				fclose (fp);
			}
		}

		// write external triangle elements to files
		{
			vector <int> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".trio.ele");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <int>* > (&(_submesh. at(i)._efaces));
				fprintf (fp, "%lu\n", cptr->size ()/ 3);

				for (size_t j = 0; j < cptr->size (); j += 3){
					fprintf (fp, "%d %d %d\n", cptr->at (j), cptr->at (j+1), cptr->at (j+2));
				}

				fclose (fp);
			}
		}

		// write external triangle topology to files
		{
			vector <Face> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".trio.top");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <Face>* > (&(_submesh. at(i)._eftop));
				fprintf (fp, "%lu\n", cptr->size ());

				for (size_t j = 0; j < cptr->size (); ++j){
					fprintf (fp, "%d %d %d\n", cptr->at (j)._neighbors[0], cptr->at (j)._neighbors[1], cptr->at (j)._neighbors[2]);
				}

				fclose (fp);
			}
		}

		// write internal triangle elements to files
		{
			vector <int> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".trii.ele");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <int>* > (&(_submesh. at(i)._ifaces));
				fprintf (fp, "%lu\n", cptr->size ()/ 3);

				for (size_t j = 0; j < cptr->size (); j += 3){
					fprintf (fp, "%d %d %d\n", cptr->at (j), cptr->at (j+1), cptr->at (j+2));
				}

				fclose (fp);
			}
		}

		// write internal triangle topology to files
		{
			vector <Face> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".trii.top");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <Face>* > (&(_submesh. at(i)._iftop));
				fprintf (fp, "%lu\n", cptr->size ());

				for (size_t j = 0; j < cptr->size (); ++j){
					fprintf (fp, "%d %d %d\n", cptr->at (j)._neighbors[0], cptr->at (j)._neighbors[1], cptr->at (j)._neighbors[2]);
				}

				fclose (fp);
			}
		}

		// write edge elements and topology to files
		{
			vector <Edge> *cptr;
			for (size_t i = 0; i < _submesh.size (); ++i){

				fname = folder + prefix;
				fname.append (".");

				stringstream ss;
				ss << i;
				fname.append (ss.str ());
				fname.append (".edge.ele");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				cptr = const_cast <vector <Edge>* > (&(_submesh. at(i)._edges));
				fprintf (fp, "%lu\n", cptr->size ());

				for (size_t j = 0; j < cptr->size (); ++j){
					fprintf (fp, "%d %d\n", cptr->at (j)._indices[0], cptr->at (j)._indices[1]);
				}

				fclose (fp);

				fname = folder + prefix;
				fname.append (".");
				fname.append (ss.str ());
				fname.append (".edge.top");

				fp = fopen (fname.c_str (), "w");
				assert (fp);

				fprintf (fp, "%lu\n", cptr->size ());

				for (size_t j = 0; j < cptr->size (); ++j){
					fprintf (fp, "%d", cptr->at (j)._nOwners);
					for (int k = 0; k < cptr->at (j)._nOwners; ++k){
						fprintf (fp, " %d", cptr->at (j)._owners[k]);
					}
					fprintf (fp, "\n");
				}

				fclose (fp);
			}
		}

	}
}
