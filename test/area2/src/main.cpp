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
#include <cstdio>
#include <vector>
#include <string>

#include "Preprocess.h"

#include "Collide/triTriCollide.h"
#include "Collide/lineTriCollide.h"

#ifdef SF_VECTOR3_ENABLED
#include "vec3.h"
#else
#include "vec4.h"
#endif

using namespace std;
using namespace SF;
    // static function to read OFF style file
    static bool readOFFMeshFile (const string &file, vector <vec> &verts, vector <unsigned int> &indices)
    {
      assert (!file.empty ());

      FILE *fp = fopen (file.c_str (), "r");
      assert (fp);

      char header [32];
      int status = fscanf (fp, "%s\n", header);
      assert (status != 0);

      int nverts = 0, nfaces = 0, nedges = 0;
      status = fscanf (fp, "%d %d %d\n", &nverts, &nfaces, &nedges);
      assert (status != 0);
      assert (nverts > 0 && (nfaces > 0 || nedges > 0));

      // read vertex data
      real tmpr [SF_VECTOR_SIZE] = {0., 0., 0.
#ifdef SF_VECTOR4_ENABLED
        , 1.
#endif
      };
      verts.reserve (static_cast <size_t> (nverts));
      for (int i = 0; i < nverts; ++i){
#ifdef SF_DOUBLE_PRECISION
        status = fscanf (fp, "%lf %lf %lf\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#else
        status = fscanf (fp, "%f %f %f\n", &(tmpr [0]), &(tmpr [1]), &(tmpr [2]));
#endif
        assert (status != 0);
        verts.push_back (vec (tmpr));
      }

      int dummy;
      int tmpd [3];
      if (nfaces){
        for (int i = 0; i < nfaces; ++i){
          status = fscanf (fp, "%d %d %d %d\n", &dummy, &(tmpd [0]), &(tmpd [1]), &(tmpd [2]));
          assert (status != 0);
          assert (tmpd [0] >= 0 && tmpd [0] < nverts);
          assert (tmpd [1] >= 0 && tmpd [1] < nverts);
          assert (tmpd [2] >= 0 && tmpd [2] < nverts);
          for (int j = 0; j < 3; ++j){
            indices.push_back ( static_cast <unsigned int> (tmpd [j]));
          }
        }
      } else if (nedges){
        for (int i = 0; i < nedges; ++i){
          status = fscanf (fp, "%d %d %d\n", &dummy, &(tmpd [0]), &(tmpd [1]));
          assert (status != 0);
          assert (tmpd [0] >= 0 && tmpd [0] < nverts);
          assert (tmpd [1] >= 0 && tmpd [1] < nverts);
          for (int j = 0; j < 2; ++j){
            indices.push_back ( static_cast <unsigned int> (tmpd [j]));
          }
        }
      }
      fclose (fp);
      return true;
    }

int main (int argc, char **argv)
{
  FILE *fp = fopen ("/home/kish1/Data/Cube/cube.tet", "r");

  unsigned int nverts = 0, ncells = 0;
  int status = fscanf (fp, "%u %u\n", &nverts, &ncells);

  vector <vec> verts (nverts + 4);
  real tmpr1, tmpr2, tmpr3;
  for (unsigned int i = 0; i < nverts; ++i){
    status = fscanf (fp, "%f %f %f\n", &tmpr1, &tmpr2, &tmpr3);
    verts [i] = vec (tmpr1, tmpr2, tmpr3);
  }

  vector <unsigned int> indices (4*ncells + 4);
  unsigned int tmpu1, tmpu2, tmpu3, tmpu4;
  for (unsigned int i = 0; i < ncells; ++i){
    status = fscanf (fp, "%u %u %u %u\n", &tmpu1, &tmpu2, &tmpu3, &tmpu4);
    indices [4*i] = tmpu1;
    indices [4*i + 1] = tmpu2;
    indices [4*i + 2] = tmpu3;
    indices [4*i + 3] = tmpu4;
  }

  fclose (fp);

  verts [nverts] = vec (174.063, 90.8677, 155.);
  verts [nverts + 1] = vec (168.833, 85.4993, 153.169);
  verts [nverts + 2] = vec (177.065, 85.5662, 147.);
  verts [nverts + 3] = vec (166.849, 89.2717, 148.477);

  indices [4*ncells] = nverts;
  indices [4*ncells + 1] = nverts + 1;
  indices [4*ncells + 2] = nverts + 2;
  indices [4*ncells + 3] = nverts + 3;

  vector <vec> bverts1, bverts2;
  vector <unsigned int> binds;
  string bfile ("/home/kish1/Data/Scalpel/scalpel_blade.off");
  readOFFMeshFile (bfile, bverts1, binds);
  bverts2 = bverts1;

  for (unsigned int i = 0; i < bverts2.size (); ++i){
    bverts1 [i]._v [1] -= 100.;
    bverts1 [i]._v [2] += 150.;
    bverts2 [i]._v [1] -= 100.;
    bverts2 [i]._v [2] += 150.;

    bverts1 [i]._v [0] += 169.075;
    bverts2 [i]._v [0] += 169.055;
  }

  vec e1, e2;
  vector <vec> normals [2];
  normals [0].resize (binds.size ()/ 2);
  normals [1].resize (normals [0].size ());
  for (unsigned int i = 0; i < normals [0].size (); ++i){
    e1 = bverts2 [binds [2*i + 1]] - bverts2 [binds [2*i]];
    e2 = bverts1 [binds [2*i + 1]] - bverts2 [binds [2*i]];
    e1.fast_ncross (normals [0][i], e2);

    e1 = bverts1 [binds [2*i]] - bverts1 [binds [2*i + 1]];
    e2 = bverts2 [binds [2*i]] - bverts1 [binds [2*i + 1]];
    e1.fast_ncross (normals [1][i], e2);
  }


  real eu1, eu2;
  for (unsigned int i = ncells; i < ncells + 1; ++i){
    PRINT ("Cube [%u]: %u %u %u %u\n", i, indices [4*i], indices [4*i + 1], indices [4*i + 2], indices [4*i + 3]);
    for (unsigned int j = 0; j < 6; ++j){
      eu1 = 0.; eu2 = 0.;
      switch (j){
        case 0:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i]], verts [indices [4*i + 1]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i]], verts [indices [4*i + 1]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i], indices [4*i + 1], eu1, eu2);
        break;

        case 1:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i]], verts [indices [4*i + 2]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i]], verts [indices [4*i + 2]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i], indices [4*i + 2], eu1, eu2);
        break;

        case 2:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i]], verts [indices [4*i + 3]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i]], verts [indices [4*i + 3]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i], indices [4*i + 3], eu1, eu2);
        break;

        case 3:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i + 1]], verts [indices [4*i + 2]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i + 1]], verts [indices [4*i + 2]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i + 1], indices [4*i + 2], eu1, eu2);
        break;

        case 4:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i + 1]], verts [indices [4*i + 3]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i + 1]], verts [indices [4*i + 3]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i + 1], indices [4*i + 3], eu1, eu2);
        break;

        case 5:
        for (unsigned int j = 0; j < normals [0].size (); ++j){
          if (lineTriCollide (eu1, verts [indices [4*i + 2]], verts [indices [4*i + 3]], bverts2 [binds [2*j]], bverts2 [binds [2*j + 1]], bverts1 [binds [2*j + 1]], normals [0][j]) ||
              lineTriCollide (eu2, verts [indices [4*i + 2]], verts [indices [4*i + 3]], bverts1 [binds [2*j + 1]], bverts1 [binds [2*j]], bverts2 [binds [2*j]], normals [1][j])){
            break;
          }
        }
        PRINT ("Edge [%u %u]: %g %g\n", indices [4*i + 2], indices [4*i + 3], eu1, eu2);
        break;
      }
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "%g %g %g\n%g %g %g\n\n", bverts1 [0]._v [0], bverts1 [0]._v [1], bverts1 [0]._v [2], bverts2 [0]._v [0], bverts2 [0]._v [1], bverts2 [0]._v [2]);
  fprintf (stdout, "%g %g %g\n%g %g %g\n", bverts1 [24]._v [0], bverts1 [24]._v [1], bverts1 [24]._v [2], bverts2 [24]._v [0], bverts2 [24]._v [1], bverts2 [24]._v [2]);
}
