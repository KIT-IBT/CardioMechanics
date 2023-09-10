/*! \file kaTriangle.h
   \brief Class for handling of triangles consisting of three points

   \author cw,fs,km, IBT - Universität Karlsruhe (TH)
 */

#ifndef KATRIANGLE_H
#define KATRIANGLE_H

#include <kaPoint.h>


template<class T>
class kaTriangle {
 public:
  kaPoint<T> p[3];

  kaTriangle() {}

  kaTriangle(const kaTriangle<T> &t)  {*this = t;}

  kaTriangle(const kaPoint<T> &p1, const kaPoint<T> &p2, const kaPoint<T> &p3);
  virtual ~kaTriangle() {}

  virtual inline const kaTriangle<T> & operator=(const kaTriangle &_triangle) {
    p[0] = _triangle.p[0]; p[1] = _triangle.p[1]; p[2] = _triangle.p[2]; return *this;
  }

  virtual inline bool operator==(const kaTriangle &t) const {
    return p[0] == t.p[0] && p[1] == t.p[1] && p[2] == t.p[2];
  }

  template<class U> void XMatrixN(const T &rm) {
    p[0].XMatrixN(rm);
    p[1].XMatrixN(rm);
    p[2].XMatrixN(rm);
  }

  const kaPoint<T> S() {
    return (p[0]+p[1]+p[2])/3.0;
  }

  // returns the surface of this triangle...using inner product method
  inline T GetSurface() const;

  // returns the Normal on the the plane the triangle's points define.
  inline kaPoint<T> Normal() const;

  // return the distance between a point and the plane the triangle's points define.
  inline T GetPointDistance(const kaPoint<T> &pt) const;
};  // class kaTriangle

template<class T>
kaTriangle<T>::kaTriangle(const kaPoint<T> &p1, const kaPoint<T> &p2, const kaPoint<T> &p3) {
  p[0] = p1;
  p[1] = p2;
  p[2] = p3;
}

template<class T>
kaPoint<T> kaTriangle<T>::Normal() const {
  return NormalizedCrossProduct(p[0]-p[1], p[0]-p[2]);
}

// http://en.wikipedia.org/wiki/Triangle

template<class T>
inline T kaTriangle<T>::GetSurface() const {
  kaPoint<T> L01(p[0]-p[1]);
  kaPoint<T> L02(p[0]-p[2]);

  //    T L01L02 = L01*L02;
  //    return T(0.5 * sqrt(fabs(L01.abs2()*L02.abs2()-L01L02*L01L02)));
  return 0.5*((CrossProduct(L01, L02)).abs());
}

// http://mathworld.wolfram.com/Plane.html

template<class T>
inline T kaTriangle<T>::GetPointDistance(const kaPoint<T> &pt) const {
  kaPoint<T> triangleNormal = Normal();

  // d = -a*x0 -b*y0 -x*z0;
  T d = -triangleNormal*p[0];

  // D = (a*xp+b*yp+c*zp+d)/sqrt(a*a+b*b+c*c);
  // rewritten...
  // T D = (triangleNormal*p + d)/triangleNormal.abs();
  // since triangleNormal is a normalized vector then triangleNormal.abs()==1;
  T D = (triangleNormal*pt + d);
  return fabs(D);
}

typedef kaTriangle<float>  kaTriangleFloat;
typedef kaTriangle<double> kaTriangleDouble;

#endif  // ifndef KATRIANGLE_H
