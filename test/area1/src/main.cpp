/**
 * @file main.cpp
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
 */

#include <cstdlib>
#include <cstring>
extern "C" {
#include <sys/stat.h>
}

#include <iostream>
#include <sstream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "Preprocess.h"
#include "vec3.h"
#include "EM_common.h"
#include "EM_Mesh.h"
#include "EM_FEMMesh.h"
#include "EM_MSDMesh.h"

using namespace std;
using namespace boost;
using namespace SF;

static inline void makeOutputFolder (const int formatcode, const int max_depth, string &folder)
{
	switch (formatcode){
		case 0:
			folder.append ("fem/");
			break;
		case 1:
			folder.append ("msd/");
			break;
	}

	struct stat st;
	if (stat (folder.c_str (), &st) != 0){
		mkdir (folder.c_str (), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	}

	stringstream ss;
	ss << max_depth;
	folder.append (ss.str ());
	folder.append ("/");

	if (stat (folder.c_str (), &st) != 0){
		mkdir (folder.c_str (), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	}
}

static void display_usage ()
{
	cerr << "usage: edit-mesh <arg> <option>" << endl;
	cerr << "List of Arguments and Options: <arg> <list-of-opts>" << endl;

	cerr << "Mandatory argument:" << endl;
	cerr << "\t-i[--input-dir] <folder_name>\tfolder name" << endl;
	cerr << "\t-p[--file-prefix] <file_prefix>\tfile prefix: of tetrahedron mesh file" << endl;
	cerr << "\t-f[--format] <format>\t\toutput format (\"fem\" or \"msd\")" << endl;
	cerr << "\t-d[--depth] <depth>\t\tdepth of recursion for mesh sub-division" << endl;
	cerr << "\t\t\t\t\tThe total number of sub-divisions is 8^depth" << endl;

	cerr << "Optional arguments:" << endl;
	cerr << "\t-e[--ext-file] <f> <x> <y> <z>\t<f>: file with mesh extents, <x> <y> <z>: aspect ratio (Default: none)" << endl;
	cerr << "\t-r[--reverse]\t\t\treverse flag: reverses orientation of starting tetrahedron" << endl;
	cerr << "\t-xyz[--start-axis] <opt>\t<opt> valid inputs - \"x\", \"y\", \"z\", \"X\", \"Y\" or \"Z\" (Default: x)" << endl;
	cerr << "\t\t\t\t\tStart axis signifies tetrahedron with min/max vertex in specified" << endl;
	cerr << "\t\t\t\t\taxis. This is used for vertex ordering. (x - Min, X - Max, ...)" << endl;
}

int
main (int argc, char* argv [])
{
	// sanity check
	if (argc < 2 || !strcmp (argv[1], "-h") || !strcmp (argv[1], "--help")){
		display_usage ();
		exit (EXIT_FAILURE);
	}

	// input variables
	int max_depth = -1;
	string folder, prefix;

	shared_ptr< Mesh > mesh;

	{
		// variables not needed by topological algorithms
		int startcell = -1;
		bool reverseflag = false;
		int formatcode = -1;
		int startcode = 0;
		string extentfile;
		float aspect_ratio[3] = {1., 1., 1.};

		// read command line arguments
		int index = 1;
		while (index < argc){
			if (!strcmp (argv[index], "-i") || !strcmp (argv[index], "--input-dir")){
				folder = string (argv[++index]);

				// conditionally append "/"
				if (folder .at(folder.size () - 1) != '/'){
					folder.append ("/");
				}
			}
			else if (!strcmp (argv[index], "-p") || !strcmp (argv[index], "--file-prefix")){
				prefix = string (argv[++index]);
			}
			else if (!strcmp (argv[index], "-f") || !strcmp (argv[index], "--format")){
				++index;
				if ( !strcmp (argv[index], "fem")){
					formatcode = 0;
					mesh = shared_ptr< Mesh > (new FEMMesh ());
				}
				else if (!strcmp (argv[index], "msd")){
					formatcode = 1;
					mesh = shared_ptr< Mesh > (new MSDMesh ());
				}
				else {
					cerr << "error: could not recognize format: " << argv[index] << "...aborting" << endl << endl;
					display_usage ();
					exit (EXIT_FAILURE);
				}
			}
			else if (!strcmp (argv[index], "-d") || !strcmp (argv[index], "--max-depth")){
				++index;
				for (int i = 0; i < strlen (argv[index]); ++i){
					assert (isdigit (argv[index][i]));
				}
				max_depth = atoi (argv[index]);
			}
			else if (!strcmp (argv[index], "-e") || !strcmp (argv[index], "--ext-file")){
				extentfile = string (argv[++index]);

				for (int i = 0; i < 3; ++i){
					++index;
					for (int j = 0; j < strlen (argv[index]); ++j){
						assert (isdigit (argv[index][j]) || argv[index][j] == '.');
					}
					aspect_ratio[i] = atof (argv[index]);
				}
			}
			else if (!strcmp (argv[index], "-r") || !strcmp (argv[index], "--reverse")){
				reverseflag = true;
			}
			else if (!strcmp (argv[index], "-xyz") || !strcmp (argv[index], "--start-axis")){
				++index;
				if (!strcmp (argv[index], "x")){
					startcode = 0;
				}
				else if (!strcmp (argv[index], "y")){
					startcode = 1;
				}
				else if (!strcmp (argv[index], "z")){
					startcode = 2;
				}
				else if (!strcmp (argv[index], "X")){
					startcode = 3;
				}
				else if (!strcmp (argv[index], "Y")){
					startcode = 4;
				}
				else if (!strcmp (argv[index], "Z")){
					startcode = 5;
				}
				else {
					cerr << "error: could not recognize start axis: " << argv[index] << "...aborting" << endl << endl;
					display_usage ();
					exit (EXIT_FAILURE);
				}
			}
			++index;
		}

		// sanity check
		assert (!folder.empty ());
		assert (!prefix.empty ());
		assert (max_depth >= 0);

		// read in file
		Mesh* mptr = mesh.get ();
		readMesh (folder, prefix, mptr->_vertices, mptr->_cells);

		// get index of starting vertex and optionally scale vertices to the extents given by the extent file
		int startvertex = -1;
		processVertices (extentfile, aspect_ratio, mptr->_vertices, startcode, startvertex);
		assert (startvertex >= 0);

		// get index of the starting cell
		startcell = getStartingCell (startvertex, mptr->_cells);

		// if reverseflag is ON, change orientation of starting cell
		if (reverseflag){
			int ind = 4*startcell;
			int tmp = mptr->_cells.at (ind);
			mptr->_cells.at (ind) = mptr->_cells.at (ind + 1);
			mptr->_cells.at(ind + 1) = tmp;
		}

		// formulate the output folder-name
		makeOutputFolder (formatcode, max_depth, folder);

		orderCells (startcell, mptr->_cells, mptr->_faces);
	}

	mesh.get ()->process (max_depth);
	mesh.get ()->writeToFiles (folder, prefix);
}
