/**
 * @file EM_common.cpp
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
 * Common definitions for the Edit Mesh application
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
extern "C"{
#include <stdint.h>
}

#include <vector>
#include <list>
#include <string>

#include "crc32.h"
#include "vec3.h"
#include "EM_common.h"

using namespace std;

namespace SF {

 	static inline bool sameCell(int *t1, int *t2)
 	{
 		bool flags[4] = {false, false, false, false};
 		for(int i = 0; i < 4; ++i){
 			for (int j = 0; j < 4; ++j){
 				flags[i] |= (t1[i] == t2[j]);
 			}
 		}
 		return flags[0] & flags[1] & flags[2] & flags[3];
 	}

 	static inline bool sameOrder (const int *t1, const int *t2)
 	{
 		int index = -1;
 		for (int i = 0; i < 3; ++i){
 			if (t1[0] == t2[i]){
 				index = i;
 			}
 		}
 		index = index > 1 ? 0 : index + 1;
 		return t1[1] == t2[index];
 	}

	// function to read mesh file
 	static void detectAnomalies (const vector <int> &indices)
 	{
 		size_t size = indices.size () / 4;

 		// detect degenerate cells
 		for (size_t i = 0; i < size; ++i){
 			i <<= 2;
 			assert (indices.at (i) != indices.at (i + 1));
 			assert (indices.at (i) != indices.at (i + 2));
 			assert (indices.at (i) != indices.at (i + 3));
 			assert (indices.at (i + 1) != indices.at (i + 2));
 			assert (indices.at (i + 1) != indices.at (i + 3));
 			assert (indices.at (i + 2) != indices.at (i + 3));
 			i >>= 2;
 		}
 		// detect duplicate cells
 		vector < vector < Cell > > chash;
 		vector <Cell>* hptr;

 		chash.reserve( USHRT_MAX + 1);
 		chash.push_back (vector <Cell> ());
 		chash.resize (USHRT_MAX + 1, vector <Cell> ());

 		int tmpi[4];
 		size_t index;
		char fstr[64], tstr1[16], tstr2[16], tstr3[16];

 		for (size_t i = 0; i < size; ++i){
 			i <<= 2;
 			for (int j = 0; j < 4; ++j){
 				tmpi[j] = indices.at (i + j);
 			}
 			i >>= 2;
 			qsort (&( tmpi[0]), 4, sizeof (int), compare);

			sprintf (fstr, "%d", tmpi[0]);
			sprintf (tstr1, "%d", tmpi[1]);
			strcat (fstr, tstr1);
			sprintf (tstr2, "%d", tmpi[2]);
			strcat (fstr, tstr2);
			sprintf (tstr3, "%d", tmpi[3]);
			strcat (fstr, tstr3);

 			index = crc32 (fstr, strlen (fstr));
 			index &= 0x0000ffff;

 			hptr = &( chash.at (index));
 			for (size_t j = 0; j < hptr->size (); ++j){
 				assert (!sameCell (tmpi, hptr->at (j)._neighbors));
 			}
 			hptr->push_back( Cell (tmpi));
 		}
 	}
	void readMesh (const string &folder, const string &prefix, vector < vec3 > &verts, vector < int > &indices)
	{
		assert (!folder.empty ());
		assert (!prefix.empty ());

		string filename (folder + prefix);
		filename.append (".tet");

		FILE* fp = fopen (filename.c_str (), "r");
		if (!fp){
      fprintf (stderr, "error: could not open %s\n", filename.c_str ());
      exit (EXIT_FAILURE);
		}

		int numVerts = 0, numIndices = 0;
		int status = fscanf (fp, "%d %d\n", &numVerts, &numIndices);
		assert (status);
		assert (numVerts);
		assert (numIndices);

		verts.reserve (numVerts);
		indices.reserve (4*numIndices);

		float tmpf[3];
		for (int i = 0; i < numVerts; ++i){
			status = fscanf (fp, "%f %f %f\n", &(tmpf[0]), &(tmpf[1]), &(tmpf[2]));
			assert (status);
			verts.push_back (vec3 (tmpf));
		}

		int tmpi[4], mini = numVerts;
		for (int i = 0; i < numIndices; ++i){
			status = fscanf (fp, "%d %d %d %d\n", &(tmpi[0]), &(tmpi[1]), &(tmpi[2]), &(tmpi[3]));
			assert (status);
			for (int j = 0; j < 4; ++j){
				mini = mini < tmpi[j] ? mini: tmpi[j];
				indices.push_back (tmpi[j]);
			}
		}
		fclose (fp);
		assert (mini >= 0);

		if (mini){
#ifndef NDEBUG
		  int vsize = static_cast <int> (verts.size ());
#endif
			for (size_t i = 0; i < indices.size (); ++i){
				indices.at (i) -= mini;
				assert (indices.at (i) < vsize);
			}
		}

		detectAnomalies (indices);
	}

	// function to scale vertices given an extent file
	void processVertices (const string &extentfile, const float aspect_ratio[3], vector < vec3 > &vertices, int startcode, int &startvertex)
	{
		// get the extents of the mesh data set
		vec3 min (vertices.at (0));
		vec3 max (min);

		int minindex[3] = {0, 0, 0}, maxindex[3] = {0, 0, 0};

		for (size_t i = 1; i < vertices.size (); ++i){
			for (int j = 0; j < 3; ++j){
				if (min._v[j] > vertices.at (i)._v[j]){
					min._v[j] = vertices.at (i)._v[j];
					minindex[j] = static_cast< int > (i);
				}
				else if (max._v[j] < vertices.at (i)._v[j]){
					max._v[j] = vertices.at (i)._v[j];
					maxindex[j] = static_cast< int >(i);
				}
			}
		}
		startvertex = startcode < 3 ? minindex[startcode] : maxindex[startcode - 3];

		if (extentfile.empty ()){
			return;
		}

		// get the extents of the file
		FILE* fp = fopen (extentfile.c_str (), "r");
		assert(fp);

		vec3 from, to;
		int status = fscanf (fp, "%f %f %f\n", &(from._v[0]), &(from._v[1]), &(from._v[2]));
		assert (status);
		status = fscanf (fp, "%f %f %f\n", &(to._v[0]), &(to._v[1]), &(to._v[2]));
		assert (status);

		fclose (fp);

		// change the extent properties to account for aspect ratio
		to -= from;
		for (int i = 0; i < 3; ++i){
			from._v[i] *= aspect_ratio[i];
		}
		for (int i = 0; i < 3; ++i){
			to._v[i] *= aspect_ratio[i];
		}

		// scale to between 0. and 1.
		max -= min;
		for (int i = 0; i < 3; ++i){
			max._v[i] = 1./max._v[i];
		}
		for (size_t i = 0; i < vertices.size (); ++i){
			vertices.at (i) -= min;
			vertices.at (i) *= max;
			vertices.at (i) *= to;
			vertices.at (i) += from;
		}
	}

	// function to get index starting cell from which topological algorithms begin
	int getStartingCell (const int startvert, const vector < int > &indices)
	{
		int result = -1;
		size_t index;

		for (size_t i = 0; i < indices.size ()/4; ++i){
			index = 4*i;
			if (indices.at (index) == startvert || indices.at (index + 1) == startvert ||
					indices.at (index + 2) == startvert || indices.at (index + 3) == startvert){
				result = static_cast< int > (i);
				break;
			}
		}
		assert (result >= 0);

		return result;
	}

	// function to generate face-topological info for triangular mesh
	void generateFaceTopology(vector <Face> &top, const vector <int> &faces)
	{
		top.reserve (faces.size ()/ 3);
		top.push_back (Face ());
		top.resize (faces.size ()/ 3, Face ());

		vector <list <Line> > edges;
		list <Line> *eptr;
		list <Line>::iterator iter;

		edges.reserve (USHRT_MAX + 1);
		edges.push_back (list <Line> ());
		edges.resize (USHRT_MAX + 1, list <Line> ());

		bool inflag = true;
		uint32_t index;
		int tmpi, inds[2];
		char tmpstr[16], finalstr[32];

		for (size_t i = 0; i < top.size (); ++i){
			for (int j = 0; j < 3; ++j){

				if (j > 1){
					inds[0] = faces.at (3*i + 2);
					inds[1] = faces.at (3*i);
				} else {
					inds[0] = faces.at (3*i + j);
					inds[1] = faces.at (3*i + j + 1);
				}

				// sort the indices
				if (inds[0] > inds[1]){
					tmpi = inds[0];
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

				// check for existing entry
				while (iter != eptr->end ()){
					if (*iter == &(inds[0])){
						inflag = true;
						top.at (i)._neighbors[j] = iter->_owner;
						top.at (iter->_owner)._neighbors[iter->_face] = i;

						eptr->erase (iter);
						break;
					}
					++iter;
				}

				// new entry
				if (!inflag){
					eptr->push_back (Line (i, j, &(inds[0])));
				}

			} // end - for (int j = 0; j < 3; ++j)
		} // end - for (size_t i = 0; i < top.size (); ++i)
	}

	// function to check ordering consistency for mesh cells
	static void checkOrder (const int start, const vector < Cell > &cells, const vector < int > &indices, vector < bool > &testflags, vector < bool > &flipflags)
	{
		list < int > ilist;
		ilist.push_back (start);

		while (!ilist.empty ()) {

			// get current index and pop the front
			size_t index = static_cast <size_t> (ilist.front ());
			ilist.pop_front ();

			// examine test-flags for all legit neighbors
			bool flag = true;
			for (int i = 0; i < 4; ++i){
				if (cells.at (index)._neighbors[i] >= 0){
					flag &= testflags.at (cells.at (index)._neighbors[i]);
				}
			}

			if (!flag){

			}
			size_t nindex;
			int self, t1[3], t2[3];

			// check handed-ness of each unchecked neighbor
			for (int i = 0; i < 4; ++i){
				if (cells.at (index)._neighbors[i] >= 0 && !testflags.at (static_cast < size_t > (cells.at (index)._neighbors[i]))){

					ilist.push_back (cells.at (index)._neighbors[i]);

					self = -1;
					nindex = static_cast <size_t> (cells.at (index)._neighbors[i]);
					for (int j = 0; j < 4; ++j){
						if (cells.at (nindex)._neighbors[j] == static_cast <int> (index)){
							self = j;
						}
					}
					assert (self >= 0);

					index <<= 2;
					assert (index + 3 < indices.size ());
					switch (i){
						case 0:
							t1[0] = indices.at (index);
							t1[1] = indices.at (index + 1);
							t1[2] = indices.at (index + 2);
							break;
						case 1:
							t1[0] = indices.at (index);
							t1[1] = indices.at (index + 2);
							t1[2] = indices.at (index + 3);
							break;
						case 2:
							t1[0] = indices.at (index);
							t1[1] = indices.at (index + 3);
							t1[2] = indices.at (index + 1);
							break;
						case 3:
							t1[0] = indices.at (index + 1);
							t1[1] = indices.at (index + 3);
							t1[2] = indices.at (index + 2);
					}
					index >>= 2;

					nindex <<= 2;
					assert (nindex + 3 < indices.size ());
					switch (self){
						case 0:
							t2[0] = indices.at (nindex);
							t2[1] = indices.at (nindex + 1);
							t2[2] = indices.at (nindex + 2);
							break;
						case 1:
							t2[0] = indices.at (nindex);
							t2[1] = indices.at (nindex + 2);
							t2[2] = indices.at (nindex + 3);
							break;
						case 2:
							t2[0] = indices.at (nindex);
							t2[1] = indices.at (nindex + 3);
							t2[2] = indices.at (nindex + 1);
							break;
						case 3:
							t2[0] = indices.at (nindex + 1);
							t2[1] = indices.at (nindex + 3);
							t2[2] = indices.at (nindex + 2);
					}
					nindex >>= 2;

					if (sameOrder (&(t1[0]), &(t2[0]))){
						flipflags.at (nindex) = flipflags.at (index) == false ? true : false;
					}
					else {
						flipflags.at (nindex) = flipflags.at (index) == true ? true : false;
					}
					testflags.at (static_cast < size_t > (cells.at (index)._neighbors[i])) = true;
				}
			}
		}
	}
	void orderCells (bool rflag, int start, vector < int > &indices, vector <int> &trigs)
	{
		vector < Cell > cells;
		cells.reserve (indices. size()/ 4);
		cells.push_back (Cell ());
		cells.resize (indices. size()/ 4, Cell ());

		// generate neighborhood information for cells
		{
			vector < list < Trig > > faces;
			list < Trig >* fptr;
			list < Trig >::iterator iter;
			faces.reserve (USHRT_MAX + 1);
			faces.push_back (list < Trig > ());
			faces.resize (USHRT_MAX + 1, list < Trig > ());

			bool inflag = false;
			uint32_t index;
			int sorted_inds[3];
			char tmpstr1[16], tmpstr2[16], finalstr[48];

			for (size_t i = 0; i < cells.size (); ++i){
				for (int j = 0; j < 4; ++j){

					// make copy of 3 face indices for each cell
					switch (j){
						case 0: // 0-1-2
							sorted_inds[0] = indices.at (4*i);
							sorted_inds[1] = indices.at (4*i + 1);
							sorted_inds[2] = indices.at (4*i + 2);
							break;
						case 1: // 0-2-3
							sorted_inds[0] = indices.at (4*i);
							sorted_inds[1] = indices.at (4*i + 2);
							sorted_inds[2] = indices.at (4*i + 3);
							break;
						case 2: // 0-3-1
							sorted_inds[0] = indices.at (4*i);
							sorted_inds[1] = indices.at (4*i + 3);
							sorted_inds[2] = indices.at (4*i + 1);
							break;
						case 3: // 1-3-2
							sorted_inds[0] = indices.at (4*i + 1);
							sorted_inds[1] = indices.at (4*i + 3);
							sorted_inds[2] = indices.at (4*i + 2);
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
						if (*iter == &(sorted_inds[0])){
							inflag = true;
							cells.at (i)._neighbors[j] = iter->_owner;
							cells.at (iter->_owner)._neighbors[iter->_face] = i;

							fptr->erase (iter);
							break;
						}
						++iter;
					}

					// new entry
					if (!inflag){
						fptr->push_back (Trig (i, j, &(sorted_inds[0])));
					}
				} // end - for (int j = 0; j < 4; ++j)
			} // end - for (size_t i = 0; i < cells.size (); ++i)


			// faces now contains list of surface triangles
			for (size_t i = 0; i < faces.size (); ++i){

				fptr = &(faces.at (i));

				while (!fptr->empty ()){
					for (int j = 0; j < 3; ++j){
						trigs.push_back (fptr->front ()._indices[j]);
					}
					fptr->pop_front ();
				}
			}

		} // end - generate neighborhood information for cells

		// initialize vector flags
		vector < bool > testflags, flipflags;

		testflags.reserve (cells.size ());
		testflags.push_back (false);
		testflags.resize (cells.size (), false);

		flipflags.reserve (cells.size ());
		flipflags.push_back (false);
		flipflags.resize (cells.size (), false);

		testflags.at (static_cast <size_t> (start)) = true;
		checkOrder (start, cells, indices, testflags, flipflags);

		cells.clear ();
		testflags.clear ();

		int counter = 0, tmpi;
		for (size_t i = 0; i < flipflags.size (); ++i){
			if (flipflags.at (i)){
				++counter;
				tmpi = indices.at (4*i);
				indices.at (4*i) = indices.at (4*i + 1);
				indices.at (4*i + 1) = tmpi;
			}
		}

    if (rflag){
      ++counter;
    }

		fprintf (stdout, "%lu cells (%d flipped)\n%lu surface triangles\n", indices.size ()/ 4, counter, trigs.size ()/ 3);
	}
}
