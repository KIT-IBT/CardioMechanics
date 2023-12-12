/*
 * File: kaPoint.h
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
 *
 */

#ifndef KAPOINT_H
#define KAPOINT_H

#include <kaMachineOS.h>
#include <kaExceptions.h>
#include <cmath>


//! Class for handling of points in 3D
/*!
   kaPoint includes standard mathematical and printing functions.

 */

template<class T>
class kaPoint {
 public:
  T x;  //!< x-coordinate of point
  T y;  //!< y-coordinate of point
  T z;  //!< z-coordinate of point

  //! Default constructor
  kaPoint() : x(T(0)), y(T(0)), z(T(0)) {}

  //! Constructor with coordinates
  kaPoint(const T xi, const T yi, const T zi) : x(xi), y(yi), z(zi) { /*x=xi; y=yi; z=zi;*/}

  //! Copy constructor
  kaPoint(const kaPoint &pt) {*this = pt;}

  kaPoint(T val) : x(val), y(val), z(val) {}

  //! assign new coordinates at 1 step
  inline void Set(const T &newx, const T &newy, const T &newz) {x = newx; y = newy; z = newz;}  // df+ 26.09.2008

  inline const kaPoint<T> & operator=(const T a) {x = y = z = a; return *this;}

  inline const kaPoint<T> & operator=(const kaPoint &pt) {
    x = pt.x;
    y = pt.y;
    z = pt.z;
    return *this;
  }

  inline void operator+=(const kaPoint &pt) {
    x += pt.x;
    y += pt.y;
    z += pt.z;
  }

  inline void operator-=(const kaPoint &pt) {
    x -= pt.x;
    y -= pt.y;
    z -= pt.z;
  }

  void operator*=(const T t) {
    x *= t;
    y *= t;
    z *= t;
  }

  void operator/=(const T t) {
    x /= t;
    y /= t;
    z /= t;
  }

  //! Overloading the comparison operator == to compare a template value with the point.
  /*!
      \param val a template variable
      \return true if x, y and z are equal to the passed parameter val
      \sa operator!=(T val), operator==(const kaPoint& pt)
   */
  virtual inline bool operator==(T val) const {
    return x == val && y == val && z == val;
  }

  //! Overloading the comparison operator != to compare a template value with the point.
  /*!
      \param val a template variable
      \return true if x, y and z are equal to the passed parameter val
      \sa operator==(T val), operator!=(const kaPoint& pt)
   */
  virtual inline bool operator!=(T val) const {
    return x != val || y != val || z != val;
  }

  virtual inline bool operator==(const kaPoint &pt) const {
    return x == pt.x && y == pt.y && z == pt.z;
  }

  virtual inline bool operator!=(const kaPoint &pt) const {
    return x != pt.x || y != pt.y || z != pt.z;
  }

  T abs() const {return sqrt(x*x+y*y+z*z);}

  T abs2() const {return x*x+y*y+z*z;}

  const kaPoint norm() const {return (*this)/this->abs();}

  const kaPoint operator+(const kaPoint &b) const {
    return kaPoint<T>(x+b.x, y+b.y, z+b.z);
  }

  const kaPoint operator-(const kaPoint &b) const {
    return kaPoint<T>(x-b.x, y-b.y, z-b.z);
  }

  const kaPoint operator*(const T t) const {
    return kaPoint<T>(x*t, y*t, z*t);
  }

  const kaPoint operator/(const T t) const {
    return kaPoint<T>(x/t, y/t, z/t);
  }

  T operator*(const kaPoint &b) const {
    return x*b.x + y*b.y + z*b.z;
  }

  const kaPoint operator%(const kaPoint &b) const {
    return kaPoint<T>(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
  }

  // this function transform the point's coordinate from the current coordinate system
  // to a new coordinate System according to the transformation matrix U
  // U can be a (3x3) matrix that means the transformation contains no translation (vectors transformation)
  // U can also be a (4x4) matrix and that means the transformation contains translation (points transformation)
  // for more information about that check the page below
  // http://www.euclideanspace.com/maths/geometry/affine/matrix4x4/index.htm
  template<class U> void XMatrixN(const U &rm) {
    T xa = x, ya = y, za = z, w;

    switch (rm.Dimension()) {
      case 4:
        x = rm.Get(0)*xa+rm.Get(1)*ya+rm.Get(2)*za+rm.Get(3);
        y = rm.Get(4)*xa+rm.Get(5)*ya+rm.Get(6)*za+rm.Get(7);
        z = rm.Get(8)*xa+rm.Get(9)*ya+rm.Get(10)*za+rm.Get(11);
        w = rm.Get(12)*xa+rm.Get(13)*ya+rm.Get(14)*za+rm.Get(15);

        if (w != 0.0) {
          x /= w;
          y /= w;
          z /= w;
        }
        break;

      case 3:
        x = rm.Get(0)*xa+rm.Get(1)*ya+rm.Get(2)*za;
        y = rm.Get(3)*xa+rm.Get(4)*ya+rm.Get(5)*za;
        z = rm.Get(6)*xa+rm.Get(7)*ya+rm.Get(8)*za;
        break;

      default:
        throw kaBaseException("kaPoint::XMatrixN<T> not defined for dimension !=3,4");
    }
  }  // XMatrixN

  void Print(const char *format = (char *)"< %lf %lf %lf >\n") const {
    printf(format, (double)x, (double)y, (double)z);
  }

  bool        Serialize(ostream &os);
  bool        Deserialize(istream &is);
  inline void ExportToStream(ostream &os, char del = ' ', bool bLF = true) const;  // GNU Compatible

  inline T Sum() const {return x+y+z;}
};  // class kaPoint

template<class T>
bool kaPoint<T>::Serialize(ostream &os) {
  os << scientific << x <<"\n";
  os << scientific << y <<"\n";
  os << scientific << z <<"\n";
  return true;
}

template<class T>
bool kaPoint<T>::Deserialize(istream &is) {
  is >> x;
  is >> y;
  is >> z;
  return true;
}

template<class T>
const kaPoint<T> operator*(const T t, const kaPoint<T> &pt) {
  return kaPoint<T>(t*pt.x, t*pt.y, t*pt.z);
}

template<class T>
const kaPoint<T> operator-(const kaPoint<T> &pt) {
  return kaPoint<T>(-pt.x, -pt.y, -pt.z);
}

template<class T>
inline ostream & operator<<(ostream &os, const kaPoint<T> &pt) {
  return os << "<" << pt.x << " " << pt.y << " " << pt.z << ">";
}

template<class T>
inline void kaPoint<T>::ExportToStream(ostream &os, char del, bool bLF) const// GNU Compatible
{
  os << scientific << x << del << scientific  << y <<del << scientific  << z;
  if (bLF) {
    os << "\n";
  }
}

template<class T>
inline istream & operator>>(istream &is, kaPoint<T> &pt) {
  return is >> pt.x >> pt.y >> pt.z;
}

// http://jtaylor1142001.net/calcjat/Solutions/VCrossProduct/VCPOrth.htm

template<class T>
inline kaPoint<T> CrossProduct(const kaPoint<T> &pt1, const kaPoint<T> &pt2) {
  return kaPoint<T>(pt1.y*pt2.z-pt2.y*pt1.z, pt1.z*pt2.x - pt1.x*pt2.z, pt1.x*pt2.y - pt2.x*pt1.y);
}

template<class T>
inline kaPoint<T> NormalizedCrossProduct(const kaPoint<T> &pt1, const kaPoint<T> &pt2) {
  return (kaPoint<T>(pt1.y*pt2.z-pt2.y*pt1.z, pt1.z*pt2.x - pt1.x*pt2.z, pt1.x*pt2.y - pt2.x*pt1.y)).norm();
}

template<class T>

inline T GetSinus(const kaPoint<T> &v1, const kaPoint<T> &v2) {
  return (CrossProduct(v1.norm(), v2.norm())).abs();
}

template<class T>

inline T GetRectangleSurface(const kaPoint<T> &P0, const kaPoint<T> &P1, const kaPoint<T> &P2) {
  return ((P0-P1).abs())*((P0-P2).abs());
}

template<class T>

T CalculateDistanceBetweenTwoLines(const kaPoint<T> &Line1Point, const kaPoint<T> &Line1Orientation,
                                   const kaPoint<T> &Line2Point, const kaPoint<T> &Line2Orientation) {
  T d = 0;

  kaPoint<T> U = CrossProduct(Line1Orientation, Line2Orientation);
  if (U.abs() == 0) {  // the lines are parallel
    d = CrossProduct((Line2Point - Line1Point), Line1Orientation).abs()/Line1Orientation.abs();
  } else {
    d = fabs( (Line2Point - Line1Point) * U.norm());
  }
  return d;
}

template<class T>
kaPoint<T> GetNormalizedVectorOrthogonalTo(const kaPoint<T> &V, const kaPoint<T> &HelperVector = kaPoint<T>(1.0, 1.0,
                                                                                                            1.0)) {  // Error
                                                                                                                     // of
                                                                                                                     // order
                                                                                                                     // of
                                                                                                                     // 1e-16
  // The task is to find P: <P,V> = 0
  // I can start by choosing a vector contained on the Plane defined by V....
  if (V.z != T(0.0)) {
    return (kaPoint<T>(HelperVector.x, HelperVector.y, -(V.x*HelperVector.x + V.y*HelperVector.y)/V.z)).norm();
  } else {
    if (V.y != T(0.0)) {
      return (kaPoint<T>(HelperVector.x, -(V.x*HelperVector.x + V.z*HelperVector.z)/V.y, HelperVector.z)).norm();
    } else {
      if (V.x != T(0.0)) {
        return (kaPoint<T>(-(HelperVector.y*V.y+HelperVector.z*V.z)/V.x, HelperVector.y, HelperVector.z)).norm();
      }
    }
  }
  throw(kaBaseException("GetNormalizedVectorOrthogonalTo(): could not return a Value"));
  return kaPoint<T>(T(0.0), T(0.0), T(0.0));
}

template<class T>
inline kaPoint<T> norm(const kaPoint<T> &p) {
  return p.norm();
}

template<class T>

inline T abs(const kaPoint<T> &p) {
  return p.abs();
}

template<class T>

inline T abs2(const kaPoint<T> &p) {
  return p.abs2();
}

typedef kaPoint<float>  kaPointFloat;
typedef kaPoint<double> kaPointDouble;

#endif  // ifndef KAPOINT_H
