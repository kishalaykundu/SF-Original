/**
 * @file lineTriCollide.cpp
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
 * Functions for line segment-triangle collision detections. It also
 * returns the local co-ordinate of the point in the line-segment.
 */

#include "lineTriCollide.h"

namespace SF {

  bool pointInTriangle (vec &p, vec &t1, vec &t2, vec &t3, vec &normal, bool planeTestFlag)
  {
    // plane test
    vec w = p - t1;
    real uu = w.dot (normal);
    if (planeTestFlag && ABS (uu) > EPSILON){
      return false;
    }

    vec u = t2 - t1;
    vec v = t3 - t1;

    uu = u.dot (u);
    real uv = u.dot (v);
    real vv = v.dot (v);

    real wu = w.dot (u);
    real wv = w.dot (v);

    real D = 1./ (uv * uv - uu * vv);

    real s = (uv * wv - vv * wu)* D;
    if (s < 0. || s > 1.){
      return false;
    }
    real t = (uv * wu - uu * wv)* D;
    if (t < 0. || s + t > 1.){
      return false;
    }
    return true;
  }

  bool lineLineCollide (vec &l11, vec &l12, vec &l21, vec &l22)
  {
    vec d1 = l12 - l11;
    vec d2 = l22 - l21;

    // check for parallelism
    real ang = d1.dot (d2);
    if (ABS(ang) > (1. - EPSILON)*d1.length ()*d2.length ()){

      // check if an end point lies inside the other line
      vec dir = l11 - l21;
      real d = dir.dot (d2);
      if (ABS (d) < (1. - EPSILON)*dir.length ()*d2.length ()){
        return false;
      }
      if (ABS (d) >= (1. - EPSILON)*dir.length ()*d2.length () && dir.length () <= d2.length ()){
        return true;
      }

      dir = l12 - l21;
      d = dir.dot (d2);
      if (ABS (d) >= (1. - EPSILON)*dir.length ()*d2.length () && dir.length () <= d2.length ()){
        return true;
      }

      dir = l21 - l11;
      d = dir.dot (d1);
      if (ABS (d) >= (1. - EPSILON)*dir.length ()*d1.length () && dir.length () <= d1.length ()){
        return true;
      }

      dir = l22 - l11;
      d = dir.dot (d1);
      if (ABS (d) >= (1. - EPSILON)*dir.length ()*d1.length () && dir.length () <= d1.length ()){
        return true;
      }

      return false;
    }

    vec c = d1.cross (d2);
    real sl = 1./c.square_length ();
    mat3x3 m (l21._v [0] - l11._v [0], l21._v [1] - l11._v [1], l21._v [2] - l11._v [2], d2._v [0], d2._v [1], d2._v [2], c._v [0], c._v [1], c._v [2]);

    real s = sl * m.determinant ();
    if (s < 0. || s > d1.length ()){
      return false;
    }
    m (1, 0) = d1._v [0];
    m (1, 1) = d1._v [1];
    m (1, 2) = d1._v [2];
    s = sl * m.determinant ();
    if (s < 0. || s > d2.length ()){
      return false;
    }
    return true;
  }

  bool lineTriCollide (real &eu, vec &l1, vec &l2, vec &t1, vec &t2, vec &t3, vec &normal)
  {
    vec dir = l2 - l1;
    vec w0 = l1 - t1;
    real a = normal.dot (w0);
    real b = -normal.dot (dir);

    // line lies in triangle plane
    if (ABS (b) < EPSILON && ABS (a) < EPSILON){
      if (pointInTriangle(l1, t1, t2, t3, normal, false) || pointInTriangle(l1, t1, t2, t3, normal, false)){
        eu = 2.;
        return true;
      }
      if (lineLineCollide (l1, l2, t1, t2) || lineLineCollide (l1, l2, t2, t3) || lineLineCollide (l1, l2, t3, t1)){
        eu = 3.;
        return true;
      }
      return false;
    }

    real r = a/ b;
    if (r < 0. || r > 1.){
      return false;
    }

    vec p = l1 + dir * r;

    vec u = t2 - t1;
    vec v = t3 - t1;

    real uu = u.dot (u);
    real uv = u.dot (v);
    real vv = v.dot (v);

    vec w = p - t1;
    real wu = w.dot (u);
    real wv = w.dot (v);

    real D = 1./ (uv * uv - uu * vv);

    real s = (uv * wv - vv * wu)* D;
    if (s < -EPSILON || s > 1. + EPSILON){
      return false;
    }
    real t = (uv * wu - uu * wv)* D;
    if (t < -EPSILON || s + t > 1. + EPSILON){
      return false;
    }

    eu = r;

    return true;
  }
}
