/**
 * @file EM_MSDMesh.cpp
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

#include <cstring>
#include <climits>

#include <sstream>
#include <vector>
#include <list>

#include "crc32.h"
#include "aabb.h"
#include "EM_common.h"
#include "EM_MSDMesh.h"

using namespace std;

namespace SF {

	typedef struct edge_t {
		int _indices[2];

		inline edge_t ()
		{
			for (int i = 0; i < 2; ++i){
				_indices[i] = -1;
			}
		}

		inline edge_t (const int *i)
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
	} Edge;

	// static function to detect which submesh a triangle belongs to
	static size_t getFaceSubmeshIndex (const vector <aabb> &bvs, const vec3 &v1, const vec3 &v2, const vec3 &v3)
	{
		size_t result = bvs.size ();
		for (size_t i = 0; i < bvs.size (); ++i){
			if (bvs.at (i).collide (v1)){
				result = i;
				break;
			}
		}
		int counter = 0;
		for (size_t i = 0; i < bvs.size (); ++i){

			counter = 0;
			if (bvs. at(i).collide (v1)){
				++counter;
			}

			if (bvs. at(i).collide (v2)){
				++counter;
				if (counter > 1){
					result = i;
					break;
				}
			}

			if (bvs. at(i).collide (v3)){
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
	MSDMesh::MSDMesh () { }

	// destructor
	MSDMesh::~MSDMesh () { }

	// method to process data to break them to an MSD format
	void
	MSDMesh::process (const int depth)
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

    // expand the bounding box
    for (int j = 0; j < 3; ++j){
      min._v [j] -= 1.0;
    }
    for (int j = 0; j < 3; ++j){
      max._v [j] += 1.0;
    }
		aabb bv (min, max);

		// calculate number of submeshes
		size_t numSubmeshes = 1;
		for (int i = 1; i <= depth; ++i){
			numSubmeshes *= 8;
		}

		vector <aabb> bboxes;
		bboxes.reserve (numSubmeshes);

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
					bboxes.push_back (aabb (bmin, bmax));
				}
			}
		}

		shuffleVertices (bboxes);

		generateEdgeList ();
		fprintf (stdout, "%lu edges\n", _edges.size ()/ 2);

		calcMassReciprocal ();

		// generate topology information
		_trigs.reserve (numSubmeshes);
		_trigs.push_back (vector <int> ());
		_trigs.resize (numSubmeshes, vector <int> ());

		_ftop.reserve (numSubmeshes);
		_ftop.push_back (vector <Face> ());
		_ftop.resize (numSubmeshes, vector <Face> ());

		size_t index;
		for (size_t i = 0; i < _faces.size (); i += 3){
			index = getFaceSubmeshIndex (bboxes, _vertices.at (_faces.at (i)), _vertices.at (_faces.at (i+1)), _vertices.at (_faces.at (i+2)));
			for (int j = 0; j < 3; ++j){
				_trigs. at(index).push_back (_faces.at (i + j));
			}
		}
		_faces.clear ();

		for (size_t i = 0; i < numSubmeshes; ++i){
			generateFaceTopology (_ftop.at (i), _trigs.at (i));
		}
	}

	// protected method to shuffle vertices such that surface vertices are in front and
	// vertices in a submesh are congruent
	void
	MSDMesh::shuffleVertices (const vector <aabb> &bvs)
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
		sOffset.reserve (bvs.size ());
		sOffset.push_back (0);
		sOffset.resize (bvs.size (), 0);

		vector <int> iOffset;
		iOffset.reserve (bvs.size ());
		iOffset.push_back (0);
		iOffset.resize (bvs.size (), 0);

		for (size_t i = 0; i < _vertices.size (); ++i){
			for (size_t j = 0; j < bvs.size (); ++j){

				if (bvs.at (j).collide (_vertices.at (i))){
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
		if (bvs.size () > 1){
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
		sCounter.reserve (bvs.size ());
		sCounter.push_back (0);
		sCounter.resize (bvs.size (), 0);

		vector <int> iCounter;
		iCounter.reserve (bvs.size ());
		iCounter.push_back (0);
		iCounter.resize (bvs.size (), 0);

		vector <size_t> newIndices;
		newIndices.reserve (_vertices.size ());
		newIndices.push_back (UINT_MAX);
		newIndices.resize (_vertices.size (), UINT_MAX);

		for (size_t i = 0; i < _vertices.size (); ++i){
			for (size_t j = 0; j < bvs.size (); ++j){

				if (bvs.at (j).collide (_vertices.at (i))){
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

	// protected method to generate list of unique edges
	void
	MSDMesh::generateEdgeList ()
	{
		vector <list <Edge> > edges;
		list <Edge> *eptr;
		list <Edge>::iterator iter;

		edges.reserve (USHRT_MAX + 1);
		edges.push_back (list <Edge> ());
		edges.resize (USHRT_MAX + 1, list <Edge> ());

		bool inflag = true;
		uint32_t index;
		int tmpi, inds[2];
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
					tmpi = inds[0];
					inds[0] = inds[1];
					inds[1] = tmpi;
				}

				// form string by concatenating three indices
				sprintf (finalstr, "%x", inds[0]);
				sprintf (tmpstr, "%x", inds[1]);
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
						break;
					}
					++iter;
				}
				if (!inflag){
					eptr->push_back( Edge (&(inds[0])));
				}
			}
		}
		// copy edges to edge-vector
		for (int i = 0; i <= USHRT_MAX; ++i){

			eptr = &(edges.at (i));

			while (!eptr->empty ()){
				for (int k = 0; k < 2; ++k){
					_edges.push_back (eptr->front ()._indices[k]);
				}
				eptr->pop_front ();
			}
		}
	}

  // protected method to calculate the reciprocal mass for every vertex
  void
  MSDMesh::calcMassReciprocal ()
  {
    _mass.resize (_vertices.size (), 0.);

    real volume;
    vec3 a, b, c;
    for (size_t i = 0; i < _cells.size (); i += 4){
      a = _vertices [_cells [i + 1]] - _vertices [_cells [i]];
      b = _vertices [_cells [i + 2]] - _vertices [_cells [i]];
      c = _vertices [_cells [i + 3]] - _vertices [_cells [i]];

      volume = a.dot (b.cross (c))/ 24.;
      volume = volume < 0. ? -volume: volume;

      for (unsigned int j = 0; j < 4; ++j){
        _mass [_cells [i + j]] += volume;
      }
    }

    for (size_t i = 0; i < _mass.size (); ++i){
      _mass [i] = 1./ _mass [i];
    }
  }

	// method to write elements to files
	void
	MSDMesh::writeElementsToFiles (const string &folder, const string &prefix) const
	{
		// write out edge file
		string fname (folder + prefix);
		fname.append (".edge");

		FILE* fp = fopen (fname.c_str (), "w");
		assert (fp);

		fprintf (fp, "%lu\n", _edges.size ()/ 2);
		for (size_t i = 0; i < _edges.size (); i += 2){
			fprintf (fp, "%d %d\n", _edges.at (i), _edges.at (i + 1));
		}

		fclose (fp);

	  // write mass reciprocal to file
	  fname = folder + prefix;
	  fname.append (".lm");

	  fp = fopen (fname.c_str (), "w");
	  assert (fp);

	  fprintf (fp, "%lu\n", _mass.size ());
	  for (size_t i = 0; i < _mass.size (); ++i){
      fprintf (fp, "%f\n", _mass.at (i));
	  }

	  fclose (fp);

		// write out triangle files
		vector <int> *tptr;
		for (size_t i = 0; i < _trigs.size (); ++i){
			fname = folder + prefix;
			fname.append (".");

			stringstream ss;
			ss << i;
			fname.append (ss.str ());
			fname.append (".tri");

			fp = fopen (fname.c_str (), "w");
			assert (fp);

			tptr = const_cast <vector <int>* > (&(_trigs.at (i)));
			fprintf (fp, "%lu\n", tptr->size ()/ 3);
			for (size_t j = 0; j < tptr->size (); j += 3){
				fprintf (fp, "%d %d %d\n", tptr->at (j), tptr->at (j + 1), tptr->at (j + 2));
			}

			fclose (fp);
		}
	}
}
