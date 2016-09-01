/**
 * @file Cut.h
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
 * The cut info-structure class for the CU_XFEM library.
 */

#pragma once

#include <cstdlib>

#include <stack>
#include <vector>

using namespace std;

namespace SF {
  namespace XFE {

    class Cut {

      public:

      unsigned int _numExVertices;
      unsigned int *_exVertices;
      unsigned int *_exUVCoords;

      unsigned int _numExFaces;
      unsigned int *_exFaces;

      unsigned int _numInVertices;
      unsigned int *_inVertices;
      unsigned int *_inUVCoords;

      unsigned int _numInFaces;
      unsigned int *_inFaces;

      inline Cut ()
      : _numExVertices (0), _exVertices (NULL), _exUVCoords (NULL), _numExFaces (0), _exFaces (NULL),
      _numInVertices (0), _inVertices (NULL), _inUVCoords (NULL), _numInFaces (0), _inFaces (NULL)
      { }

      inline ~Cut ()
      {
        if (_numExVertices){
          free (_exVertices);
          if (_exUVCoords){
            free (_exUVCoords);
          }
        }
        if (_numExFaces){
          free (_exFaces);
        }
        if (_numInVertices){
          free (_inVertices);
          if (_inUVCoords){
            free (_inUVCoords);
          }
        }
        if (_numInFaces){
          free (_inFaces);
        }
      }

      // copy constructor
      inline Cut (const Cut &c)
      : _numExVertices (c._numExVertices), _exVertices (NULL), _exUVCoords (NULL), _numExFaces (c._numExFaces), _exFaces (NULL),
      _numInVertices (c._numInVertices), _inVertices (NULL), _inUVCoords (NULL), _numInFaces (c._numInFaces), _inFaces (NULL)
      {
        if (_numExVertices){
          _exVertices = (unsigned int *)malloc (_numExVertices * sizeof (unsigned int));
          memcpy (_exVertices, c._exVertices, _numExVertices * sizeof (unsigned int));
          if (c._exUVCoords){
            _exUVCoords = (unsigned int *) malloc (_numExVertices * sizeof (unsigned int));
            memcpy (_exUVCoords, c._exUVCoords, _numExVertices * sizeof (unsigned int));
          }
        }
        if (_numExFaces){
          _exFaces = (unsigned int *)malloc (_numExFaces * sizeof (unsigned int));
          memcpy (_exFaces, c._exFaces, _numExFaces * sizeof (unsigned int));
        }
        if (_numInVertices){
          _inVertices = (unsigned int *)malloc (_numInVertices * sizeof (unsigned int));
          memcpy (_inVertices, c._inVertices, _numInVertices * sizeof (unsigned int));
          if (c._inUVCoords){
            _inUVCoords = (unsigned int *) malloc (_numInVertices * sizeof (unsigned int));
            memcpy (_inUVCoords, c._inUVCoords, _numInVertices * sizeof (unsigned int));
          }
        }
        if (_numInFaces){
          _inFaces = (unsigned int *)malloc (_numInFaces * sizeof (unsigned int));
          memcpy (_inFaces, c._inFaces, _numInFaces * sizeof (unsigned int));
        }
      }

      // assignment operator
      inline Cut& operator = (const Cut &c)
      {
        _numExVertices = c._numExVertices;
        if (_numExVertices){
          _exVertices = (unsigned int *)malloc (_numExVertices * sizeof (unsigned int));
          memcpy (_exVertices, c._exVertices, _numExVertices * sizeof (unsigned int));
          if (c._exUVCoords){
            _exUVCoords = (unsigned int *) malloc (_numExVertices * sizeof (unsigned int));
            memcpy (_exUVCoords, c._exUVCoords, _numExVertices * sizeof (unsigned int));
          }
        }
        _numExFaces = c._numExFaces;
        if (_numExFaces){
          _exFaces = (unsigned int *)malloc (_numExFaces * sizeof (unsigned int));
          memcpy (_exFaces, c._exFaces, _numExFaces * sizeof (unsigned int));
        }
        _numInVertices = c._numInVertices;
        if (_numInVertices){
          _inVertices = (unsigned int *)malloc (_numInVertices * sizeof (unsigned int));
          memcpy (_inVertices, c._inVertices, _numInVertices * sizeof (unsigned int));
          if (c._inUVCoords){
            _inUVCoords = (unsigned int *) malloc (_numInVertices * sizeof (unsigned int));
            memcpy (_inUVCoords, c._inUVCoords, _numInVertices * sizeof (unsigned int));
          }
        }
        _numInFaces = c._numInFaces;
        if (_numInFaces){
          _inFaces = (unsigned int *)malloc (_numInFaces * sizeof (unsigned int));
          memcpy (_inFaces, c._inFaces, _numInFaces * sizeof (unsigned int));
        }

        return *this;
      }

      // method to deallocate variables
      inline void
      deallocate (unsigned int nelems, unsigned int &myelem, unsigned int *arr, stack <unsigned int> &earr)
      {
        for (unsigned int i = nelems; i < myelem; ++i){
          earr.push (arr [i]);
        }
        myelem = nelems;
        arr = (unsigned int *) realloc (arr, nelems * sizeof (unsigned int));
      }

      // method to allocate internal variables
      inline void
      allocateInternalVariables (unsigned int nVerts, unsigned int nFaces, vector <vec> &verts, vector <float> &surfaceStatus,
                                 vector <vec2> &tex2D, vector <vec3> &tex3D, vector <unsigned int> &faces,
                                 stack <unsigned int> &eVerts, stack <unsigned int> &eFaces)
      {
        if (!_numInVertices){
          _inVertices = (unsigned int *) malloc (nVerts * sizeof (unsigned int));
        } else {
          _inVertices = (unsigned int *) realloc (_inVertices, nVerts * sizeof (unsigned int));
        }
        while (!eVerts.empty () && _numInVertices < nVerts){
          _inVertices [_numInVertices] = eVerts.top ();
          eVerts.pop ();
          ++_numInVertices;
        }

        unsigned int size = verts.size ();
        for (unsigned int i = _numInVertices; i < nVerts; ++i){
          _inVertices [i] = size;
          ++size;
          verts.push_back (vec ());
          surfaceStatus.push_back (0.);
          tex2D.push_back (vec2 ());
          tex3D.push_back (vec3 ());
        }
        _numInVertices = nVerts;

        if (!_numInFaces){
          _inFaces = (unsigned int *) malloc (nFaces * sizeof (unsigned int));
        } else {
          _inFaces = (unsigned int *) realloc (_inFaces, nFaces * sizeof (unsigned int));
        }
        while (!eFaces.empty () && _numInFaces < nFaces){
          _inFaces [_numInFaces] = eFaces.top ();
          eFaces.pop ();
          ++_numInFaces;
        }

        size = faces.size ();
        for (unsigned int i = _numInFaces; i < nFaces; ++i){
          _inFaces [i] = size;
          size += 3;
          for (unsigned int j = 0; j < 3; ++j){
            faces.push_back (0);
          }
        }
        _numInFaces = nFaces;
      }

      // method to allocate external variables
      inline void
      allocateExternalVariables (unsigned int nVerts, unsigned int nFaces, vector <vec> &verts, vector <vec2> &tex2D,
                                 vector <unsigned int> &faces, stack <unsigned int> &eVerts, stack <unsigned int> &eFaces)
      {
        if (_numExVertices < nVerts){
          if (!_numExVertices){
            _exVertices = (unsigned int *) malloc (nVerts * sizeof (unsigned int));
          } else {
            _exVertices = (unsigned int *) realloc (_exVertices, nVerts * sizeof (unsigned int));
          }
          while (!eVerts.empty () && _numExVertices < nVerts){
            _exVertices [_numExVertices] = eVerts.top ();
            eVerts.pop ();
            ++_numExVertices;
          }

          unsigned int size = verts.size ();
          for (unsigned int i = _numExVertices; i < nVerts; ++i){
            _exVertices [i] = size;
            ++size;
            verts.push_back (vec ());
            tex2D.push_back (vec2 ());
          }
        }
        else if (_numExVertices > nVerts){
          for (unsigned int i = nVerts; i < _numExVertices; ++i){
            eVerts.push (_exVertices [i]);
          }
          _exVertices = (unsigned int *) realloc (_exVertices, nVerts * sizeof (unsigned int));
        }
        _numExVertices = nVerts;

        if (_numExFaces < nFaces){
          if (!_numExFaces){
            _exFaces = (unsigned int *) malloc (nFaces * sizeof (unsigned int));
          } else {
            _exFaces = (unsigned int *) realloc (_exFaces, nFaces * sizeof (unsigned int));
          }
          while (!eFaces.empty () && _numExFaces < nFaces){
            _exFaces [_numExFaces] = eFaces.top ();
            eFaces.pop ();
            ++_numExFaces;
          }

          unsigned int size = faces.size ();
          for (unsigned int i = _numExFaces; i < nFaces; ++i){
            _exFaces [i] = size;
            size += 3;
            for (unsigned int j = 0; j < 3; ++j){
              faces.push_back (0);
            }
          }
        }
        else if (_numExFaces > nFaces){
          for (unsigned int i = nFaces; i < _numExFaces; ++i){
            eFaces.push (_exFaces [i]);
          }
          _exFaces = (unsigned int *) realloc (_exFaces, nFaces * sizeof (unsigned int));
        }
        _numExFaces = nFaces;
      }

      // method to allocate internal uv-coords
      inline void
      allocateInternalUVCoords (vector <vec3> &uvCoords)
      {
        _inUVCoords = (unsigned int *) malloc (_numInVertices * sizeof (unsigned int));

        unsigned int size = uvCoords.size ();
        for (unsigned int i = 0; i < _numInVertices; ++i){
          _inUVCoords [i] = size;
          ++size;
          uvCoords.push_back (vec3 ());
        }
      }

      // method to allocate internal uv-coords
      inline void
      allocateExternalUVCoords (vector <vec3> &uvCoords)
      {
        _exUVCoords = (unsigned int *) malloc (_numExVertices * sizeof (unsigned int));

        unsigned int size = uvCoords.size ();
        for (unsigned int i = 0; i < _numExVertices; ++i){
          _exUVCoords [i] = size;
          ++size;
          uvCoords.push_back (vec3 ());
        }
      }
    };

  }
}
