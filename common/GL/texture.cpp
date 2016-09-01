/**
 * @file texture.cpp
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
 * Common texture related functions used by SF modules
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>

extern "C" {
#include <libpng12/png.h>
#include <libpng12/pngconf.h>
}

#include <vector>
#include <list>

#include <Eigen/Core>
#include <Eigen/LU>

#include "crc32.h"
#include "vec3.h"
#include "mat3x3.h"

#include "common.h"
#include "texture.h"

using namespace std;
using namespace Eigen;

namespace SF {

  static const real RAY_SCALE = .3;
  static const real ALPHA_THRESHOLD = .9;
  static const real ALPHA_DISTANCE = .05;
  static const real SCALE_CONSTANT = 1./ 255.;

  // write unsigned byte RGB data to png file
  void writeRGBToPng (const char *prefix, int index, int dim, const vector <GLubyte> &rgb)
  {
    string filename (prefix);
    {
      char indStr [32];
      sprintf (indStr, "%u.png", index);
      filename.append (indStr);
    }

    FILE *fp = fopen (filename.c_str (), "wb");
    assert (fp);

		// initialize png-specific functions
		png_structp png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if (!png_ptr){
			PRINT ("GL error: could not allocate png_ptr\n");
			fclose (fp);
			return;
		}

		png_infop info_ptr = png_create_info_struct (png_ptr);
		if (!info_ptr){
			PRINT ("GL error: could not allocate png info structure\n");
			png_destroy_write_struct (&png_ptr, static_cast <png_infopp> (NULL));
			fclose (fp);
			return;
		}

		if (setjmp (png_jmpbuf (png_ptr))){
			PRINT ("GL error: could not write png header\n");
			png_destroy_write_struct (&png_ptr, &info_ptr);
			fclose (fp);
			return;
		}

		png_init_io (png_ptr, fp);

		png_set_filter (png_ptr, 0, PNG_FILTER_NONE);
		png_bytep *row_headers = new png_bytep [dim];

		png_set_IHDR (png_ptr, info_ptr, dim, dim, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

		for (int i = 0; i < dim; ++i){
			row_headers [dim - i - 1] = (png_bytep) (&(rgb [3*dim*i]));
		}

		png_write_info (png_ptr, info_ptr);
		png_write_image (png_ptr, row_headers);

		delete [] row_headers;

		png_write_end (png_ptr, NULL);

		// clean up
		png_destroy_write_struct (&png_ptr, &info_ptr);

		fclose (fp);
  }
  // write unsigned RGBA data to png file (drop alpha channel)
  void writeRGBAToPng (const char *prefix, int index, int dim, const vector <GLubyte> &rgba)
  {
    vector <GLubyte> rgb (3*dim*dim);
    for (int i = 0; i < dim*dim; ++i){
      for (int j = 0; j < 3; ++j){
        rgb [3*i + j] = rgba [4*i + j];
      }
    }
    writeRGBToPng (prefix, index, dim, rgb);
  }
  // write float RGBA data to png file (normalize data to between 0 and 255)
  void writeFloatToPng (const char *prefix, int index, int dim, const vector <GLfloat> &rgba)
  {
    GLfloat min [3], max [3];
    bool first = true;
    for (int i = 0; i < 4*dim*dim; i += 4){
      if (rgba [i + 3] > .5){
        if (!first) {
          for (int j = 0; j < 3; ++j){
            if (min [j] > rgba [i + j]){
              min [j] = rgba [i + j];
            } else if (max [j] < rgba [i + j]){
              max [j] = rgba [i + j];
            }
          }
        } else {
          memcpy (min, &(rgba [i]), 3*sizeof (GLfloat));
          memcpy (max, min, 3*sizeof (GLfloat));
          first = false;
        }
      }
    }
    for (int i = 0; i < 3; ++i){
      max [i] -= min [i];
    }
    for (int i = 0; i < 3; ++i){
      max [i] = 255./ max [i]; // scales values to between 0 and 255
    }

    vector <GLubyte> rgbu;
    rgbu.resize (3*dim*dim, 0x00);
    for (int i = 0; i < dim*dim; ++i){
      if (rgba [4*i + 3] > .5){
        for (int j = 0; j < 3; ++j){
          rgbu [3*i + j] = static_cast <unsigned char> (floor (max [j]* (rgba [4*i + j] - min [j])));
        }
      }
    }
    writeRGBToPng (prefix, index, dim, rgbu);
  }

  // function to get topology information for a triangular mesh
  void initTopologyInfo (const vector <unsigned int> &faces, vector <FaceEdge> &edges, vector <FaceNeighbor> &neighbors)
  {
    // initialize neighbors vector
    neighbors.reserve (faces.size ()/3);
    neighbors.resize (faces.size ()/3, FaceNeighbor ());

    // make hash (vector) tables of size 2^16
    vector <list <FaceEdge> > elists; // edges
    elists.resize (USHRT_MAX + 1, list <FaceEdge> ());

    vector <list <unsigned int> > olists; // edge owners
    olists.resize (USHRT_MAX + 1, list <unsigned int> ());

    vector <list <unsigned int> > ilists; // edge owner index
    ilists.resize (USHRT_MAX + 1, list <unsigned int> ());

    bool inflag = false;
    uint32_t index = 0;
    unsigned int inds [2];
    char tmpStr [9], finalStr [32];

    list <FaceEdge>::iterator eiter;
    list <unsigned int>::iterator oiter;
    list <unsigned int>::iterator iiter;

    for (unsigned int i = 0; i < faces.size ()/ 3; ++i){
      for (unsigned int j = 0; j < 3; ++j){

        // make a copy of two edges and sort them
        if (j < 2){
          inds [0] = faces [3*i + j];
          inds [1] = faces [3*i + j + 1];
        } else {
          inds [0] = faces [3*i + 2];
          inds [1] = faces [3*i];
        }

        // sort indices - to get consistent edge ordering
        if (inds [0] > inds [1]){
          unsigned int tmpi = inds [0];
          inds [0] = inds [1];
          inds [1] = tmpi;
        }

        // form a string by concatenating the two indices
        sprintf (tmpStr, "%u", inds [0]);
        assert (strlen (tmpStr) < 9);
        sprintf (finalStr, "%u", inds [1]);
        assert (strlen (finalStr) + strlen (tmpStr) < 32);
        strcat (finalStr, tmpStr);

        // get crc32 hash code and mask it to below 65536
        index = crc32 (finalStr, strlen (finalStr));
        index &= 0x0000ffff;

        eiter = elists [index].begin ();
        oiter = olists [index].begin ();
        iiter = ilists [index].begin ();

        // look for existence of edge in hash-table
        inflag = false;
        while (eiter != elists [index].end ()){
          if (eiter->_v [0] != inds [0]){
            ++eiter;
            ++oiter;
            ++iiter;
          } else if (eiter->_v [1] != inds [1]){
            ++eiter;
            ++oiter;
            ++iiter;
          } else {

            neighbors [i]._v [j] = static_cast <int> (*oiter);
            neighbors [*oiter]._v [*iiter] = static_cast <int> (i);

            edges.push_back (*eiter);

            elists [index].erase (eiter);
            olists [index].erase (oiter);
            ilists [index].erase (iiter);

            inflag = true;
            break;
          }
        }

        if (!inflag) {
          elists [index].push_back (FaceEdge (inds [0], inds [1]));
          olists [index].push_back (i);
          ilists [index].push_back (j);
        }
      } // end - for (int j = 0; j < 3; ++j)
    } // end - for (unsigned int i = 0; i < faces.size ()/3; ++i)

    for (unsigned int i = 0; i <= USHRT_MAX; ++i){
      while (!elists [i].empty ()){
        edges.push_back (elists [i].front ());
        elists [i].pop_front ();
      }
    }
  }

  // function to get one-ring of faces around face-mesh
  static inline bool edgeBelongsToFace (const unsigned int *edge, const unsigned int *face1, const unsigned int *face2)
  {
    bool trueflag = false;
    for (int i = 0; i < 3; ++i){
      if (edge [0] == face1 [i]){
        trueflag = true;
      }
    }
    if (!trueflag){
      return false;
    }

    trueflag = false;
    for (int i = 0; i < 3; ++i){
      if (edge [1] == face1 [i]){
        trueflag = true;
      }
    }
    if (!trueflag){
      return false;
    }

    // at this point we know edge belongs to face1 and need to check if face1 and face2 are same
    return (face1[0] != face2[0]) || (face1[1] != face2[1]) || (face1[2] != face2[2]);
  }
  void getFaceRings (unsigned int index, const vector <vector <unsigned int> > &faceIndices, vector <unsigned int> &faces)
  {
    assert (index < faceIndices.size ());

    vector <FaceEdge> edges;
    vector <FaceNeighbor> neighbors;

    initTopologyInfo (faceIndices [index], edges, neighbors);
    edges.clear ();

    unsigned int uind [2];
    bool yesflag = false, noflag = true;
    vector <unsigned int> faceVectorIndex;
    vector <unsigned int> faceIndex;

    for (unsigned int i = 0; i < neighbors.size (); ++i){
      for (unsigned int j = 0; j < 3; ++j){

        // only examine faces with no neighbors (border triangles)
        if (neighbors [i]._v [j] < 0){

          // gather edge info
          if (j > 1){
            uind [0] = faceIndices [index][3*i + 2];
            uind [1] = faceIndices [index][3*i];
          } else {
            uind [0] = faceIndices [index][3*i + j];
            uind [1] = faceIndices [index][3*i + j + 1];
          }

          yesflag = false;

          // examine neighboring sets for triangle neighbor
          for (unsigned int k = 0; k < index; ++k){

            for (unsigned int l = 0; l < faceIndices [k].size ()/ 3; ++l){
              if (edgeBelongsToFace (&(uind [0]), &(faceIndices [k][3*l]), &(faceIndices [index][3*i]))){

                // if not already in index, do so
                noflag = true;
                for (unsigned int a = 0; a < faceVectorIndex.size (); ++a){
                  if (k == faceVectorIndex [a] && l == faceIndex [a]){
                    noflag = false;
                    break;
                  }
                }
                if (noflag){
                  faceVectorIndex.push_back (k);
                  faceIndex.push_back (l);
                  for (unsigned int a = 0; a < 3; ++a){
                    faces.push_back (faceIndices[k][3*l + a]);
                  }
                }
                yesflag = true;
                break;
              }
            }
            if (yesflag){
              break;
            }
          }
          if (!yesflag){
            for (unsigned int k = index + 1; k < faceIndices.size (); ++k){
              for (unsigned int l = 0; l < faceIndices [k].size ()/ 3; ++l){
                if (edgeBelongsToFace (&(uind [0]), &(faceIndices [k][3*l]), &(faceIndices [index][3*i]))){

                  // if not already in index, do so
                  noflag = true;
                  for (unsigned int a = 0; a < faceVectorIndex.size (); ++a){
                    if (k == faceVectorIndex [a] && l == faceIndex [a]){
                      noflag = false;
                      break;
                    }
                  }
                  if (noflag){
                    faceVectorIndex.push_back (k);
                    faceIndex.push_back (l);
                    for (unsigned int a = 0; a < 3; ++a){
                      faces.push_back (faceIndices[k][3*l + a]);
                    }
                  }
                  yesflag = true;
                  break;
                }
              }
              if (yesflag){
                break;
              }
            }
          }
          assert (yesflag);
        } // end - if (neighbors [i]._v [j] < 0)
      } // end - for (int j = 0; j < 3; ++j)
    } // end - for (unsigned int i = 0; i < neighbors.size (); ++i)
  }

  // function to calculate parametric coordinates using Tutte's method (used to calculate texture coordinates)
  static inline void getVertexSubset (const vector <vec> &inVerts, const vector <unsigned int> &inFaces,
                                      vector <vec3> &outVerts, vector <unsigned int> &outFaces, vector <unsigned int> &uniqueVertIndices)
  {
    // get unique list of indices (also tells the size of outverts)
    list <unsigned int> tmpi;
    for (unsigned int i = 0; i < inFaces.size (); ++i){
      tmpi.push_back (inFaces [i]);
    }
    tmpi.sort ();
    tmpi.unique ();

    uniqueVertIndices.reserve (tmpi.size ());
    while (!tmpi.empty ()){
      uniqueVertIndices.push_back (tmpi.front ());
      tmpi.pop_front ();
    }

    // gather the relevant vertices in new array
    outVerts.reserve (uniqueVertIndices.size ());
    for (unsigned int i = 0; i < uniqueVertIndices.size (); ++i){
      outVerts.push_back (vec3 (inVerts [uniqueVertIndices [i]]._v [0], inVerts [uniqueVertIndices [i]]._v [1], inVerts [uniqueVertIndices [i]]._v [2]));
    }

    // output the new indices
#ifndef NDEBUG
    bool checkFlag = false;
#endif

    for (unsigned int i = 0; i < inFaces.size (); ++i){
#ifndef NDEBUG
      checkFlag = false;
#endif
      for (unsigned int j = 0; j < uniqueVertIndices.size (); ++j){
        if (inFaces [i] == uniqueVertIndices [j]){
          outFaces.push_back (j);
#ifndef NDEBUG
          checkFlag = true;
#endif
          break;
        }
      }
#ifndef NDEBUG
      assert (checkFlag);
#endif
    }
  }
  void calculateParametricCoordinates (unsigned int numSurfaceVerts, const vector <vec> &vertices, const vector <unsigned int> &indices, vector <vec2> &texCoords)
  {
    if (indices.empty ()){
      return;
    }

    // get a local subset of vertices used by current set of indices (also changes indices)
    vector <vec3> verts;
    vector <unsigned int> faces;
    vector <unsigned int> uniqueVertIndices;
    getVertexSubset (vertices, indices, verts, faces, uniqueVertIndices);

    // get neighboring information for the faces
    vector <FaceEdge> edges;
    vector <FaceNeighbor> neighbors;
    initTopologyInfo (faces, edges, neighbors);

    // transform 3D vertices to 2D form and identify border vertices
    vector <vec2> tmpTexCoords;
    vector <bool> borderFlag;
    unsigned int numBorderVerts = 0;
    {
      // calculate area-weighted average normal
      vec3 normal (vec3::ZERO), tmpv, e1, e2;
      for (unsigned int i = 0; i < faces.size (); i += 3){
        e1 = verts [faces [i + 1]] - verts [faces [i]];
        e2 = verts [faces [i + 2]] - verts [faces [i]];
        e1.fast_cross (tmpv, e2);
        normal += tmpv;
      }
      normal.normalize ();

      // transform vertices to normal plane
      mat3x3 transform (mat3x3::ZERO);
      transform (0, 0) = normal._v [1] * normal._v [1] + normal._v [2] * normal._v [2];
      transform (0, 1) = -normal._v [0] * normal._v [1];
      transform (0, 2) = -normal._v [0] * normal._v [2];
      transform (1, 0) = transform (0, 1);
      transform (1, 1) = normal._v [0] * normal._v [0] + normal._v [2] * normal._v [2];
      transform (1, 2) = -normal._v [1] * normal._v [2];
      transform (2, 0) = transform (0, 2);
      transform (2, 1) = transform (1, 2);
      transform (2, 2) = normal._v [0] * normal._v [0] + normal._v [1] * normal._v [1];

      for (unsigned int i = 0; i < verts.size (); ++i){
        verts [i] = transform * verts [i];
      }

      // rotate transformed vertices to cardinal plane closest to normal plane
      unsigned int maxi = 0;
      real maxval = ABS (normal._v [0]);
      if (ABS (normal._v [1]) > maxval){
        maxi = 1;
        maxval = ABS (normal._v [1]);
      }
      if (ABS (normal._v [2]) > maxval){
        maxi = 2;
      }
      vec3 axis (vec3::ZERO);
      if (normal._v [maxi] < 0.){
        axis._v [maxi] = -1.;
      } else {
        axis._v [maxi] = 1.;
      }

      vec3 rotationAxis;
      normal.fast_cross (rotationAxis, axis);

      real sinval = rotationAxis.length ();
      real cosval = static_cast <real> (sqrt (1. - sinval * sinval));

      rotationAxis.normalize ();

      real u2 = rotationAxis._v [0]*rotationAxis._v [0];
      real v2 = rotationAxis._v [1]*rotationAxis._v [1];
      real w2 = rotationAxis._v [2]*rotationAxis._v [2];
      real uvcos = rotationAxis._v [0]*rotationAxis._v [1]*(1. - cosval);
      real uwcos = rotationAxis._v [0]*rotationAxis._v [2]*(1. - cosval);
      real vwcos = rotationAxis._v [1]*rotationAxis._v [2]*(1. - cosval);
      real usin = rotationAxis._v [0]*sinval;
      real vsin = rotationAxis._v [1]*sinval;
      real wsin = rotationAxis._v [2]*sinval;
      transform = mat3x3 (u2 + (1. - u2)*cosval, uvcos - wsin, uwcos + vsin,
                          uvcos + wsin, v2 + (1. - v2)*cosval, vwcos - usin,
                          uwcos - vsin, vwcos + usin, w2 + (1. - w2)*cosval);

      for (unsigned int i = 0; i < verts.size (); ++i){
        verts [i] = transform * verts [i];
      }

      // initialize temporary texture coordinates
      tmpTexCoords.reserve (verts.size ());
      switch (maxi){
        case 0:
          for (unsigned int i = 0; i < verts.size (); ++i){
            tmpTexCoords.push_back (vec2 (verts [i]._v [1], verts [i]._v [2]));
          }
          break;
        case 1:
          for (unsigned int i = 0; i < verts.size (); ++i){
            tmpTexCoords.push_back (vec2 (verts [i]._v [0], verts [i]._v [2]));
          }
          break;
        case 2:
          for (unsigned int i = 0; i < verts.size (); ++i){
            tmpTexCoords.push_back (vec2 (verts [i]._v [0], verts [i]._v [1]));
          }
      }

      // tag border vertices as such and make a count
      borderFlag.resize (verts.size (), false);
      for (unsigned int i = 0; i < faces.size ()/ 3; ++i){
        for (unsigned int j = 0; j < 3; ++j){
          if (neighbors [i]._v [j] < 0){
            switch (j){
              case 0:
                borderFlag [faces [3*i]] = true;
                borderFlag [faces [3*i + 1]] = true;
                break;
              case 1:
                borderFlag [faces [3*i + 1]] = true;
                borderFlag [faces [3*i + 2]] = true;
                break;
              case 2:
                borderFlag [faces [3*i]] = true;
                borderFlag [faces [3*i + 2]] = true;
            }
          }
        }
      }
      for (unsigned int i = 0; i < borderFlag.size (); ++i){
        if (borderFlag [i]){
          ++numBorderVerts;
        }
      }

      // scale border 2D texture coordinates to 0-1 range
      real min [2] = {tmpTexCoords [0]._v [0], tmpTexCoords [0]._v [1]};
      real max [2] = {min [0], min [1]};
      for (unsigned int i = 1; i < tmpTexCoords.size (); ++i){
        for (unsigned int j = 0; j < 2; ++j){
          if (min [j] > tmpTexCoords [i]._v [j]){
            min [j] = tmpTexCoords [i]._v [j];
          } else if (max [j] < tmpTexCoords [i]._v [j]){
            max [j] = tmpTexCoords [i]._v [j];
          }
        }
      }
      for (int j = 0; j < 2; ++j){
        max [j] -= min [j];
      }
      for (int j = 0; j < 2; ++j){
        max [j] = 1./ max [j];
      }
      for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
        if (borderFlag [i]){
          for (int j = 0; j < 2; ++j){
            tmpTexCoords [i]._v [j] -= min [j];
            tmpTexCoords [i]._v [j] *= max [j];
          }
        }
      }
    }

    // run convex hull on border vertices and read in the output
    vector <vec2> convexCoords;
    {
      // write border vertices to temp file
      FILE *fp = fopen ("./.tmpQHullInput", "w");
      assert (fp);
      fprintf (fp, "2\n%u\n", numBorderVerts);
      for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
        if (borderFlag [i]){
          fprintf (fp, "%g %g\n", tmpTexCoords [i]._v [0], tmpTexCoords [i]._v [1]);
        }
      }
      fclose (fp);

      // run qhull on the file and print output to temp file
      int status = system ("qconvex Qc p < ./.tmpQHullInput > ./.tmpQHullOutput");
      assert (!status);

      // read temp file outputted by qhull
      unsigned int numVerts = 0;
      fp = fopen ("./.tmpQHullOutput", "r");
      assert (fp);
      status = fscanf (fp, "%u\n", &numVerts); // this is the cardinality of each point (2 in our case)
      assert (status != 0);
      status = fscanf (fp, "%u\n", &numVerts);
      assert (status != 0);
      assert (numVerts);

      convexCoords.reserve (numVerts);

      real tmpr [2];
      for (unsigned int i = 0; i < numVerts; ++i){
#ifdef SF_DOUBLE_PRECISION
        status = fscanf (fp, "%lf %lf\n", &(tmpr [0]), &(tmpr [1]));
#else
        status = fscanf (fp, "%f %f\n", &(tmpr [0]), &(tmpr [1]));
#endif
        assert (status != 0);
        convexCoords.push_back (vec2 (tmpr));
      }
      fclose (fp);

      // remove temp files
      status = system ("rm -f ./.tmpQHull*");
      assert (!status);
    }

    // update boundary flags after convex hull operation
    {
      bool present = true;
      for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
        if (borderFlag [i]){
          present = false;
          for (unsigned int j = 0; j < convexCoords.size (); ++j){
            if (tmpTexCoords [i] == convexCoords [j]){
              present = true;
              break;
            }
          }
          borderFlag [i] = present;
        }
      }
    }

    // calculate degree of incidence for all vertices
    vector <real> lambda;
    lambda.resize (verts.size (), 0.);
    for (unsigned int i = 0; i < edges.size (); ++i){
      lambda  [edges [i]._v [0]] += 1.;
      lambda  [edges [i]._v [1]] += 1.;
    }
    for (unsigned int i = 0; i < lambda.size (); ++i){
      assert (lambda [i] > 0.);
      lambda [i] = 1./ lambda [i];
    }

    // calculate new indices separating border and internal vertices
    vector <unsigned int> newIndices (verts.size ());
    {
      unsigned int counter1 = 0, counter2 = 0;
      for (unsigned int i = 0; i < borderFlag.size (); ++i){
        if (borderFlag [i]){
          newIndices [i] = counter1;
          ++counter1;
        } else {
          newIndices [i] = counter2;
          ++counter2;
        }
      }
    }

    // initialize eigen matrix and vector
    unsigned int numInsideVerts = tmpTexCoords.size () - convexCoords.size ();
#ifdef SF_DOUBLE_PRECISION
    MatrixXd A (MatrixXd::Identity (numInsideVerts, numInsideVerts));
    VectorXd B (VectorXd::Zero (numInsideVerts));
    VectorXd C (VectorXd::Zero (numInsideVerts));
#else
    MatrixXf A (MatrixXf::Identity (numInsideVerts, numInsideVerts));
    VectorXf B (VectorXf::Zero (numInsideVerts));
    VectorXf C (VectorXf::Zero (numInsideVerts));
#endif

    unsigned int index1, index2;
    for (unsigned int i = 0; i < edges.size (); ++i){
      index1 = edges [i]._v [0];
      index2 = edges [i]._v [1];
      assert (index1 < verts.size ());
      assert (index2 < verts.size ());

      if (!borderFlag [index1] && !borderFlag [index2]){
        assert (newIndices [index1] < numInsideVerts);
        assert (newIndices [index2] < numInsideVerts);
        A (newIndices [index1], newIndices [index2]) = -lambda [index1];
        A (newIndices [index2], newIndices [index1]) = -lambda [index2];
      }
      else if (!borderFlag [index1]){
        assert (newIndices [index1] < numInsideVerts);
        assert (newIndices [index2] < convexCoords.size ());
        B [newIndices [index1]] += convexCoords [newIndices [index2]]._v [0] * lambda [index1];
        C [newIndices [index1]] += convexCoords [newIndices [index2]]._v [1] * lambda [index1];
      }
      else if (!borderFlag [index2]){
        assert (newIndices [index2] < numInsideVerts);
        assert (newIndices [index1] < convexCoords.size ());
        B [newIndices [index2]] += convexCoords [newIndices [index1]]._v [0] * lambda [index2];
        C [newIndices [index2]] += convexCoords [newIndices [index1]]._v [1] * lambda [index2];
      }
    }

    // do LU decomposition to obtain jittered internal vertices
#ifdef SF_DOUBLE_PRECISION
    FullPivLU <MatrixXd> luOfA (A);
    VectorXd x1;
#else
    FullPivLU <MatrixXf> luOfA (A);
    VectorXf x1;
#endif
    x1 = luOfA.solve (B);

#ifdef SF_DOUBLE_PRECISION
    VectorXd x2;
#else
    VectorXf x2;
#endif
    x2 = luOfA.solve (C);

    // copy to temporary 2D texture coordinate buffer
    index1 = 0;
    for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
      if (!borderFlag [i]){
        tmpTexCoords [i]._v [0] = x1 [index1];
        tmpTexCoords [i]._v [1] = x2 [index1];
        ++index1;
      }
    }

    // scale it again to between 0 and 1
    real min [2] = {tmpTexCoords [0]._v [0], tmpTexCoords [0]._v [1]};
    real max [2] = {min [0], min [1]};
    for (unsigned int i = 1; i < tmpTexCoords.size (); ++i){
      for (int j = 0; j < 2; ++j){
        if (min [j] > tmpTexCoords [i]._v [j]){
          min [j] = tmpTexCoords [i]._v [j];
        } else if (max [j] < tmpTexCoords [i]._v [j]){
          max [j] = tmpTexCoords [i]._v [j];
        }
      }
    }
    for (int j = 0; j < 2; ++j){
      max [j] -= min [j];
    }
    for (int j = 0; j < 2; ++j){
      max [j] = 1./ max [j];
    }
    for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
      for (int j = 0; j < 2; ++j){
        tmpTexCoords [i]._v [j] -= min [j];
        tmpTexCoords [i]._v [j] *= max [j];
      }
    }

    // copy the resulting values to destination
    for (unsigned int i = 0; i < texCoords.size (); ++i){
      for (int j = 0; j < 2; ++j){
        texCoords [i]._v [j] = 0.;
      }
    }
    for (unsigned int i = 0; i < tmpTexCoords.size (); ++i){
      assert (uniqueVertIndices [i] < texCoords.size ());
      texCoords [uniqueVertIndices [i]] = tmpTexCoords [i];
    }
  }

  // function to scale vertices with proper aspect ratio to between 0 and 1
  void scaleVertices (const real *aspect, const vector <vec> &src, const aabb &bv, vector <vec> &dest)
  {
    //real scale [3] = {1./ aspect [0], 1./ aspect [1], 1./ aspect [2]};
    for (unsigned int i = 0; i < dest.size (); ++i){
      for (int j = 0; j < 3; ++j){
        dest [i]._v [j] = /*scale [j] */(src [i]._v [j] - bv._v [0]._v [j]);
      }
    }
    real min [3] = {dest [0]._v [0], dest [0]._v [1], dest [0]._v [2]};
    real max [3] = {min [0], min [1], min [2]};
    for (unsigned int i = 1; i < dest.size (); ++i){
      for (int j = 0; j < 3; ++j){
        if (min [j] > dest [i]._v [j]){
          min [j] = dest [i]._v [j];
        } else if (max [j] < dest [i]._v [j]){
          max [j] = dest [i]._v [j];
        }
      }
    }
    for (int i = 0; i < 3; ++i){
      max [i] -= min [i];
    }
    for (int i = 0; i < 3; ++i){
      max [i] = 1./ max [i];
    }
    for (unsigned int i = 0; i < dest.size (); ++i){
      for (int j = 0; j < 3; ++j){
        dest [i]._v [j] -= min [j];
        dest [i]._v [j] *= max [j];
      }
    }
  }

  // function to calculate area-weighted vertex normals for a mesh
  void calculateVertexNormals (const vector <vec> &verts, const vector <vector <unsigned int> > &faces, vector <vec> &normals)
  {
    vec e1, e2, normal;
    for (unsigned int i = 0; i < faces.size (); ++i){
      for (unsigned int j = 0; j < faces [i].size (); j += 3){
        e1 = verts [faces [i][j + 1]] - verts [faces [i][j]];
        e2 = verts [faces [i][j + 2]] - verts [faces [i][j]];
        e1.fast_cross (normal, e2);

        for (int k = 0; k < 3; ++k){
          normals [faces [i][j + k]] += normal;
        }
      }
    }
    for (unsigned int i = 0; i < normals.size (); ++i){
      normals [i].normalize ();
    }
  }

  // function to calculate texture atlases
  void initTextureAtlas (GLuint program, int dim, const vector <vec> &verts, const vector <vec2> &texCoords,
                         const vector <unsigned int> &faces, vector <GLfloat> &rgbaData)
  {
    GLenum error;

    glUseProgram (program);
    checkGLError (error);

    // get shader variable locations
    GLint vertLocation = glGetAttribLocation (program, "vertex");
    assert (vertLocation > -1);
    GLint texCoordLocation = glGetAttribLocation (program, "texCoord");
    assert (texCoordLocation > -1);
    glBindFragDataLocation (program, 0, "fragColor");
    checkGLError (error);

    glDisable (GL_CULL_FACE);

    glClampColor (GL_CLAMP_VERTEX_COLOR, GL_FALSE);
    glClampColor (GL_CLAMP_READ_COLOR, GL_FALSE);
    glClampColor (GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);

    // initialize framebuffer texture
    GLuint textureId;
    glGenTextures (1, &textureId);
    checkGLError (error);
    glBindTexture (GL_TEXTURE_2D, textureId);
    checkGLError (error);
    glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA32F, dim, dim, 0, GL_RGBA, GL_FLOAT, 0);
    checkGLError (error);

    // initialize framebuffer and bind texture
    GLuint fboId;
    glGenFramebuffers (1, &fboId);
    checkGLError (error);
    glBindFramebuffer (GL_FRAMEBUFFER, fboId);
    checkGLError (error);
    glFramebufferTexture2D (GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureId, 0);

    // initialize index buffer
    GLuint indexId;
    glGenBuffers (1, &indexId);
    checkGLError (error);
    glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, indexId);
    checkGLError (error);
    glBufferData (GL_ELEMENT_ARRAY_BUFFER, sizeof (unsigned int)*faces.size (), &(faces [0]), GL_STATIC_DRAW);
    checkGLError (error);

    // initialize render vertex array
    GLuint arrayId;
    glGenVertexArrays (1, &arrayId);
    checkGLError (error);
    glBindVertexArray (arrayId);
    checkGLError (error);

    GLuint vertexId;
    glGenBuffers (1, &vertexId);
    checkGLError (error);
    glBindBuffer (GL_ARRAY_BUFFER, vertexId);
    checkGLError (error);
    glBufferData (GL_ARRAY_BUFFER, 2*sizeof (real)*texCoords.size (), &(texCoords [0]), GL_STATIC_DRAW);
    checkGLError (error);
    glVertexAttribPointer (vertLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray (vertLocation);

    GLuint texCoordId;
    glGenBuffers (1, &texCoordId);
    checkGLError (error);
    glBindBuffer (GL_ARRAY_BUFFER, texCoordId);
    checkGLError (error);
    glBufferData (GL_ARRAY_BUFFER, SF_VECTOR_SIZE*sizeof (real)*verts.size (), &(verts [0]._v[0]), GL_STATIC_DRAW);
    checkGLError (error);
    glVertexAttribPointer (texCoordLocation, SF_VECTOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray (texCoordLocation);

    // render
    glBindFramebuffer (GL_FRAMEBUFFER, fboId);
    checkGLError (error);

    glPushAttrib (GL_VIEWPORT_BIT);
    glViewport (0, 0, dim, dim);
    checkGLError (error);

    glDrawBuffer (GL_COLOR_ATTACHMENT0);

    glClearColor (0., 0., 0., 0.);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray (arrayId);
    checkGLError (error);
    glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, indexId);
    checkGLError (error);
    glDrawElements (GL_TRIANGLES, faces.size (), GL_UNSIGNED_INT, 0);
    checkGLError (error);

    glFlush ();

    // bind pixel buffer and read from framebuffer
    glReadBuffer (GL_COLOR_ATTACHMENT0);
    glReadPixels (0, 0, dim, dim, GL_RGBA, GL_FLOAT, &(rgbaData [0]));
    checkGLError (error);

    glPopAttrib ();

    glDeleteBuffers (1, &texCoordId);
    glDeleteBuffers (1, &vertexId);
    glDeleteVertexArrays (1, &arrayId);
    glDeleteBuffers (1, &indexId);
    glDeleteTextures (1, &textureId);
    glDeleteFramebuffers (1, &fboId);

    glClampColor (GL_CLAMP_VERTEX_COLOR, GL_TRUE);
    glClampColor (GL_CLAMP_READ_COLOR, GL_TRUE);
    glClampColor (GL_CLAMP_FRAGMENT_COLOR, GL_TRUE);
    glEnable (GL_CULL_FACE);

    glUseProgram (0);
  }

  // ray-trace function to reach texture boundary (returns end-position of rays)
	static inline real getAlpha (int offset1, int offset2, const Texture3D &texture, const int* intPos, const real* delta)
	{
	  int index = offset1*intPos [2] + offset2*intPos [1] + 4*intPos [0];
		real comp1 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + 3]) + delta [2] * static_cast <real> (texture._rgba [index + offset1 + 3]);
		real comp2 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + offset2 + 3]) +
      delta [2] * static_cast <real> (texture._rgba [index + offset1 + offset2 + 3]);

		real comp3 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + 7]) + delta [2] * static_cast <real> (texture._rgba [index + offset1 + 7]);
		real comp4 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + offset2 + 7]) +
      delta [2] * static_cast <real> (texture._rgba [index + offset1 + offset2 + 7]);

		return SCALE_CONSTANT * ((1. - delta [0])*( (1. - delta [1])*comp1 + delta [1]*comp2) + delta [0]*( (1. - delta [1])*comp3 + delta [1]*comp4));
	}
  void raytraceThroughVolumef (int dim, const vector <GLfloat> &coData, const vector <GLfloat> &noData, const Texture3D &texture, vector <GLfloat> &rgbaData)
  {
    real realDims [3] = {static_cast <real> (texture._dimension [0]), static_cast <real> (texture._dimension [1]), static_cast <real> (texture._dimension [2])};

    int intPos [3];
    bool outside, reverse;
    real alpha, prevAlpha, realPos [3], ray [3], delta [3];

    int offset1 = static_cast <int> (4*texture._dimension [1]*texture._dimension [0]);
    int offset2 = static_cast <int> (4*texture._dimension [0]);

    for (int i = 0; i < 4*dim*dim; i += 4){
      if (coData [i + 3] > .5){

        // see if current position is outside texture borders
        outside = false;
        for (int j = 0; j < 3; ++j){
          realPos [j] = (realDims [j] - 10.)*coData [i + j] + 5.;
        }
        if (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
            realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
            realPos [2] < 0. || realPos [2] > realDims [2] - 1.) {
          outside = true;
        }
        for (int j = 0; j < 3; ++j){
          ray [j] = -RAY_SCALE * noData [i + j];
        }

        // if point is outside, then advance along the reverse normal till it gets inside texture border
        if (outside){
          while (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
              realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
              realPos [2] < 0. || realPos [2] > realDims [2] - 1.) {
            for (int j = 0; j < 3; ++j){
              realPos [j] += ray [j];
            }
          }
        }

        // get rounded integer coordinates
        for (int j = 0; j < 3; ++j){
          intPos [j] = static_cast <int> (floor (realPos [j]));
        }
        for (int j = 0; j < 3; ++j){
          delta [j] = realPos [j] - static_cast <real> (intPos [j]);
        }
        alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

        // ray-tracing starts here
        if (! (outside && alpha > ALPHA_THRESHOLD)){

          // calculate if ray should go forward or back
          if (alpha < ALPHA_THRESHOLD){
            reverse = true;
          } else {
            reverse = false;
            for (int j = 0; j < 3; ++j){
              ray [j] = -ray [j];
            }
          }

          // advance ray till it hits texture iso-surface boundary
          if (reverse){

            while (realPos [0] > 0. && realPos [0] < realDims [0] - 1. &&
                   realPos [1] > 0. && realPos [1] < realDims [1] - 1. &&
                   realPos [2] > 0. && realPos [2] < realDims [2] - 1. && alpha < ALPHA_THRESHOLD) {
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
            }
          }
          else { // forward ray step

            while (realPos [0] > 0. && realPos [0] < realDims [0] - 1. &&
                   realPos [1] > 0. && realPos [1] < realDims [1] - 1. &&
                   realPos [2] > 0. && realPos [2] < realDims [2] - 1. && alpha > ALPHA_THRESHOLD) {
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
            }
          }

          // if out of texture border, take a step back and get the value
          if (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
              realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
              realPos [2] < 0. || realPos [2] > realDims [2] - 1.){
            for (int j = 0; j < 3; ++j){
              realPos [j] -= ray [j];
            }
            for (int j = 0; j < 3; ++j){
              intPos [j] = static_cast <int> (floor (realPos [j]));
            }
            for (int j = 0; j < 3; ++j){
              delta [j] = realPos [j] - static_cast <real> (intPos [j]);
            }
            alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
          }
          else { // inside texture border, so we can converge on isosurface boundary

            for (int j = 0; j < 3; ++j){
              ray [j] *= -.5;
            }
            for (int j = 0; j < 3; ++j){
              realPos [j] += ray [j];
            }
            for (int j = 0; j < 3; ++j){
              intPos [j] = static_cast <int> (floor (realPos [j]));
            }
            for (int j = 0; j < 3; ++j){
              delta [j] = realPos [j] - static_cast <real> (intPos [j]);
            }
            alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

            while (ABS (alpha - ALPHA_THRESHOLD) > ALPHA_DISTANCE){
              if (sqrt (ray [0]*ray [0] + ray [1]*ray [1] + ray [2]*ray [2]) < EPSILON){
                break;
              }
              prevAlpha = alpha;

              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

              if ((alpha - ALPHA_THRESHOLD)*(prevAlpha - ALPHA_THRESHOLD) < 0.){
                for (int j = 0; j < 3; ++j){
                  ray [j] *= -.5;
                }
              }
            } // while (ABS (alpha - ALPHA_THRESHOLD) > ALPHA_DISTANCE)
          }

        } // end - if (!outside && alpha > ALPHA_THRESHOLD)

        rgbaData [i + 3] = alpha;
        for (int j = 0; j < 3; ++j){
          rgbaData [i + j] = realPos [j];
        }

      } // end - if (coData [i + 3] > .5)
    } // end - for (int i = 0; i < 4*dim*dim; i += 4)

    // scale vertices to between 0 and 1, because these are 3D texture coordinates
    for (int j = 0; j < 3; ++j){
      realDims [j] = 1./ realDims [j];
    }
    for (int i = 0; i < 4*dim*dim; i += 4){
      if (rgbaData [i + 3] > EPSILON){
        for (int j = 0; j < 3; ++j){
          rgbaData [i + j] += .5;
          rgbaData [i + j] *= realDims [j];
        }
      }
    }

  }

  // ray-trace function to reach texture boundary (returns trilinearly interpolated texel color at end of ray)
  static inline void getColor (int offset1, int offset2, const Texture3D &texture, const int *intPos, const real *delta, GLubyte *rgb)
  {
    int index;
    real sum, comp1, comp2, comp3, comp4;
    for (int i = 0; i < 3; ++i){
      index = offset1*intPos [2] + offset2*intPos [1] + 4*intPos [0];
      comp1 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + i]) + delta [2] * static_cast <real> (texture._rgba [index + offset1 + i]);
      comp2 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + offset2 + i]) +
        delta [2] * static_cast <real> (texture._rgba [index + offset1 + offset2 + i]);

      comp3 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + 4 + i]) + delta [2] * static_cast <real> (texture._rgba [index + offset1 + 4 + i]);
      comp4 =  (1. - delta [2]) * static_cast <real> (texture._rgba [index + offset2 + 4 + i]) +
        delta [2] * static_cast <real> (texture._rgba [index + offset1 + offset2 + 4 + i]);

      sum = SCALE_CONSTANT * ((1. - delta [0])*( (1. - delta [1])*comp1 + delta [1]*comp2) + delta [0]*( (1. - delta [1])*comp3 + delta [1]*comp4));
      sum = sum > 1. ? 1. : sum;
      rgb [i] = static_cast <GLubyte> (floor (255. * sum));
    }
  }
  void raytraceThroughVolumeb (int dim, const vector <GLfloat> &coData, const vector <GLfloat> &noData, const Texture3D &texture, vector <GLubyte> &rgbaData)
  {
    real realDims [3] = {static_cast <real> (texture._dimension [0]), static_cast <real> (texture._dimension [1]), static_cast <real> (texture._dimension [2])};

    int intPos [3];
    bool outside, reverse;
    real alpha, prevAlpha, realPos [3], ray [3], delta [3];

    int offset1 = static_cast <int> (4*texture._dimension [1]*texture._dimension [0]);
    int offset2 = static_cast <int> (4*texture._dimension [0]);

    for (int i = 0; i < 4*dim*dim; i += 4){
      if (coData [i + 3] > .5){

        // see if current position is outside texture borders
        outside = false;
        for (int j = 0; j < 3; ++j){
          realPos [j] = (realDims [j] - 10.)*coData [i + j] + 5.;
        }
        if (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
            realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
            realPos [2] < 0. || realPos [2] > realDims [2] - 1.) {
          outside = true;
        }
        for (int j = 0; j < 3; ++j){
          ray [j] = -RAY_SCALE * noData [i + j];
        }

        // if point is outside, then advance along the reverse normal till it gets inside texture border
        if (outside){
          while (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
              realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
              realPos [2] < 0. || realPos [2] > realDims [2] - 1.) {
            for (int j = 0; j < 3; ++j){
              realPos [j] += ray [j];
            }
          }
        }

        // get rounded integer coordinates
        for (int j = 0; j < 3; ++j){
          intPos [j] = static_cast <int> (floor (realPos [j]));
        }
        for (int j = 0; j < 3; ++j){
          delta [j] = realPos [j] - static_cast <real> (intPos [j]);
        }
        alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

        // ray-tracing starts here
        if (! (outside && alpha > ALPHA_THRESHOLD)){

          // calculate if ray should go forward or back
          if (alpha < ALPHA_THRESHOLD){
            reverse = true;
          } else {
            reverse = false;
            for (int j = 0; j < 3; ++j){
              ray [j] = -ray [j];
            }
          }

          // advance ray till it hits texture iso-surface boundary
          if (reverse){

            while (realPos [0] > 0. && realPos [0] < realDims [0] - 1. &&
                   realPos [1] > 0. && realPos [1] < realDims [1] - 1. &&
                   realPos [2] > 0. && realPos [2] < realDims [2] - 1. && alpha < ALPHA_THRESHOLD) {
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
            }
          }
          else { // forward ray step

            while (realPos [0] > 0. && realPos [0] < realDims [0] - 1. &&
                   realPos [1] > 0. && realPos [1] < realDims [1] - 1. &&
                   realPos [2] > 0. && realPos [2] < realDims [2] - 1. && alpha > ALPHA_THRESHOLD) {
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
            }
          }

          // if out of texture border, take a step back and get the value
          if (realPos [0] < 0. || realPos [0] > realDims [0] - 1. ||
              realPos [1] < 0. || realPos [1] > realDims [1] - 1. ||
              realPos [2] < 0. || realPos [2] > realDims [2] - 1.){
            for (int j = 0; j < 3; ++j){
              realPos [j] -= ray [j];
            }
            for (int j = 0; j < 3; ++j){
              intPos [j] = static_cast <int> (floor (realPos [j]));
            }
            for (int j = 0; j < 3; ++j){
              delta [j] = realPos [j] - static_cast <real> (intPos [j]);
            }
            alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));
          }
          else { // inside texture border, so we can converge on isosurface boundary

            for (int j = 0; j < 3; ++j){
              ray [j] *= -.5;
            }
            for (int j = 0; j < 3; ++j){
              realPos [j] += ray [j];
            }
            for (int j = 0; j < 3; ++j){
              intPos [j] = static_cast <int> (floor (realPos [j]));
            }
            for (int j = 0; j < 3; ++j){
              delta [j] = realPos [j] - static_cast <real> (intPos [j]);
            }
            alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

            while (ABS (alpha - ALPHA_THRESHOLD) > ALPHA_DISTANCE){
              if (sqrt (ray [0]*ray [0] + ray [1]*ray [1] + ray [2]*ray [2]) < EPSILON){
                break;
              }
              prevAlpha = alpha;

              for (int j = 0; j < 3; ++j){
                realPos [j] += ray [j];
              }
              for (int j = 0; j < 3; ++j){
                intPos [j] = static_cast <int> (floor (realPos [j]));
              }
              for (int j = 0; j < 3; ++j){
                delta [j] = realPos [j] - static_cast <real> (intPos [j]);
              }
              alpha = getAlpha (offset1, offset2, texture, &(intPos [0]), &(delta [0]));

              if ((alpha - ALPHA_THRESHOLD)*(prevAlpha - ALPHA_THRESHOLD) < 0.){
                for (int j = 0; j < 3; ++j){
                  ray [j] *= -.5;
                }
              }
            } // while (ABS (alpha - ALPHA_THRESHOLD) > ALPHA_DISTANCE)
          }

        } // end - if (!outside && alpha > ALPHA_THRESHOLD)

        rgbaData [i + 3] = static_cast <GLubyte> (floor (255.*alpha));
        getColor (offset1, offset2, texture, &(intPos [0]), &(delta [0]), &(rgbaData [i]));

      } // end - if (coData [i + 3] > .5)
    } // end - for (int i = 0; i < 4*dim*dim; i += 4)
  }
}
