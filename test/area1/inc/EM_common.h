/**
 * @file EM_common.h
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

#pragma once

#include <vector>
#include <string>

using namespace std;

namespace SF {

	class vec3;

	typedef struct line_t {
		int _owner;
		int _face;
		int _indices[2];

		inline line_t ()
		: _owner (-1), _face (-1)
		{
			for (int i = 0; i < 2; ++i){
				_indices[i] = -1;
			}
		}

		inline line_t (const int o, const int f, const int *i)
		: _owner (o), _face (f)
		{
			memcpy (_indices, i, 2*sizeof (int));
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
	} Line;

	typedef struct trig_t {
		int _owner;
		int _face;
		int _indices[3];

		inline trig_t ()
		: _owner (-1), _face (-1)
		{
			for (int i = 0; i < 3; ++i){
				_indices[i] = -1;
			}
		}

		inline trig_t (const int o, const int f, const int *i)
		: _owner (o), _face (f)
		{
			memcpy (_indices, i, 3*sizeof (int));
		}

		inline bool operator == (const int *f) const
		{
			bool flags[3] = {false, false, false};
			for(int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					flags[i] |= (_indices[i] == f[j]);
				}
			}
			return flags[0] & flags[1] & flags[2];
		}
	} Trig;

	typedef struct face_t {
		int _neighbors[3];

		inline face_t ()
		{
			for (int i = 0; i < 3; ++i){
				_neighbors[i] = -1;
			}
		}
		inline face_t (const face_t &c)
		{
			memcpy (_neighbors, c._neighbors, 3*sizeof (int));
		}

		inline face_t ( const int *n)
		{
			memcpy (_neighbors, n, 3*sizeof (int));
		}
	} Face;

	typedef struct cell_t {
		int _neighbors[4];

		inline cell_t ()
		{
			for (int i = 0; i < 4; ++i){
				_neighbors[i] = -1;
			}
		}
		inline cell_t (const cell_t &c)
		{
			memcpy (_neighbors, c._neighbors, 4*sizeof (int));
		}

		inline cell_t ( const int *n)
		{
			memcpy (_neighbors, n, 4*sizeof (int));
		}
	} Cell;

	// comparison function used by qsort
 	static int compare (const void *a, const void *b)
	{
		return (*(int *) a - *(int *) b);
	}

 	inline bool sameTriangle(const int *t1, const int *t2)
 	{
 		bool flags[3] = {false, false, false};
 		for(int i = 0; i < 3; ++i){
 			for (int j = 0; j < 3; ++j){
 				flags[i] |= (t1[i] == t2[j]);
 			}
 		}
 		return flags[0] & flags[1] & flags[2];
 	}

	// function to load a mesh from a file
	void readMesh (const string &folder, const string &prefix, vector <vec3> &verts, vector <int> &indices);

	// function to scale vertices based on an extent file
	void processVertices (const string &efile, const float asp[3], vector < vec3 > &verts, int startcode, int &startvert);

	// function to get index starting cell from which topological algorithms begin
	int getStartingCell (const int startvert, const vector < int > &indices);

	// function to generate face-topological info for triangular mesh
	void generateFaceTopology(vector <Face> &top, const vector <int> &faces);

	// function to examine ordering consistency in cells
	void orderCells (const int start, vector < int > &indices, vector <int> &faces);
}
