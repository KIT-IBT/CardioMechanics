/*
 * File: FormFunctions.h
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


#ifndef FORMFUNCTIONS_H
#define FORMFUNCTIONS_H

#include <Tensor2.h>

// calculates the inverse of a square 3x3 matrix! the determinant must be passed as an argument...

template<class T>

inline void InverseMatrix3x3(const T *p, T *q, T detp) {
  const T a = p[0], b = p[1], c = p[2];
  const T d = p[3], e = p[4], f = p[5];
  const T g = p[6], h = p[7], i = p[8];

  const T detinvp = 1.0/detp;

  q[0] = (e*i-f*h)*detinvp;
  q[1] = (c*h-b*i)*detinvp;
  q[2] = (b*f-c*e)*detinvp;
  q[3] = (f*g-d*i)*detinvp;
  q[4] = (a*i-c*g)*detinvp;
  q[5] = (c*d-a*f)*detinvp;
  q[6] = (d*h-e*g)*detinvp;
  q[7] = (b*g-a*h)*detinvp;
  q[8] = (a*e-b*d)*detinvp;
}

template<class T>

inline void VecSetVal(T *x, const T v, int cnt) {
  const T *xe = x+cnt;

  while (x < xe)
    *x++ = v;
}

template<class T>

inline void Vec3SetVal(T *x, T v) {
  x[0] = v;
  x[1] = v;
  x[2] = v;
}

template<class T>

inline void Vec6SetVal(T *x, T v) {
  x[0] = v;
  x[1] = v;
  x[2] = v;
  x[3] = v;
  x[4] = v;
  x[5] = v;
}

template<class T>

inline void Vec8SetVal(T *x, T v) {
  x[0] = v;
  x[1] = v;
  x[2] = v;
  x[3] = v;
  x[4] = v;
  x[5] = v;
  x[6] = v;
  x[7] = v;
}

template<class T>

inline void VecSetVec(T *x, const T *y, int cnt) {
  const T *xe = x+cnt;

  while (x < xe)
    *x++ = *y++;
}

template<class T>

inline T VecMultVec(const T *x, const T *y, int cnt) {
  const T *xe = x+cnt;
  T rc        = 0.;

  while (x < xe)
    rc += *x++ **y++;
  return rc;
}

template<class T>

inline void VecInnerMultVec(T *z, const T *x, const T *y, int cnt) {
  const T *xe = x+cnt;

  while (x < xe)
    *z++ = *x++ **y++;
}

// Dot product of two 6x1 elements vectors

template<class T>

inline T Vec6MultVec6(const T *x, const T *y) {
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]+x[3]*y[3]+x[4]*y[4]+x[5]*y[5];
}

// Dot product of two 8x1 elements vectors

template<class T>

inline T Vec8MultVec8(const T *x, const T *y) {
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]+x[3]*y[3]+x[4]*y[4]+x[5]*y[5]+x[6]*y[6]+x[7]*y[7];
}

// Hexahedral shape funtions

// also:
// The function implimented as in page 240 of the holy book:
// Nonlinear Finite Elements for Continua and Structures
// Ted Belytschko, Wing Kam Liu, Brian Moran
// ISBN-13: 978-0471-98773-1
// http://www.amazon.com/Nonlinear-Finite-Elements-Continua-Structures/dp/0471987735
// Page 626
template<class T, const int Samples>
class HexaederGaussianQuadrature {
 public:
  T W[Samples*Samples*Samples];
  T gp[Samples*Samples*Samples][3];

  // Samples3 = Samples^3;
  int Samples3;

  // Constructor set the interpolation points and the weights assossiated with those points according to the count of
  // nodes...1, 8 or 27
  HexaederGaussianQuadrature() {
    Samples3 = Samples*Samples*Samples;

    switch (Samples) {
      case 1:
        W[0]     = 1.;
        gp[0][0] = .5;
        gp[0][1] = .5;
        gp[0][2] = .5;
        break;

      case 2:

        // the weight is 1/8 = .125
        W[0] = .125; gp[0][0] = .5-.5/sqrt(3.0); gp[0][1] = .5-.5/sqrt(3.0); gp[0][2] = .5-.5/sqrt(3.0);
        W[1] = .125; gp[1][0] = .5+.5/sqrt(3.0); gp[1][1] = .5-.5/sqrt(3.0); gp[1][2] = .5-.5/sqrt(3.0);
        W[2] = .125; gp[2][0] = .5+.5/sqrt(3.0); gp[2][1] = .5+.5/sqrt(3.0); gp[2][2] = .5-.5/sqrt(3.0);
        W[3] = .125; gp[3][0] = .5-.5/sqrt(3.0); gp[3][1] = .5+.5/sqrt(3.0); gp[3][2] = .5-.5/sqrt(3.0);
        W[4] = .125; gp[4][0] = .5-.5/sqrt(3.0); gp[4][1] = .5-.5/sqrt(3.0); gp[4][2] = .5+.5/sqrt(3.0);
        W[5] = .125; gp[5][0] = .5+.5/sqrt(3.0); gp[5][1] = .5-.5/sqrt(3.0); gp[5][2] = .5+.5/sqrt(3.0);
        W[6] = .125; gp[6][0] = .5+.5/sqrt(3.0); gp[6][1] = .5+.5/sqrt(3.0); gp[6][2] = .5+.5/sqrt(3.0);
        W[7] = .125; gp[7][0] = .5-.5/sqrt(3.0); gp[7][1] = .5+.5/sqrt(3.0); gp[7][2] = .5+.5/sqrt(3.0);
        break;

      case 3:
        for (int i = 0; i < 3; i++) {
          T posx = .5+(1-i)*.5*sqrt(3.0/5.0);
          T Wx   = (i != 1 ? 5.0/18.0 : 4.0/9.0);
          for (int j = 0; j < 3; j++) {
            T posy = .5+(1-j)*.5*sqrt(3.0/5.0);
            T Wy   = (j != 1 ? 5.0/18.0 : 4.0/9.0);
            for (int k = 0; k < 3; k++) {
              T posz    = .5+(1-k)*.5*sqrt(3.0/5.0);
              T Wz      = (k != 1 ? 5.0/18.0 : 4.0/9.0);
              int index = i+j*3+k*9;
              W[index]     = Wx*Wy*Wz;
              gp[index][0] = posx;
              gp[index][1] = posy;
              gp[index][2] = posz;
            }
          }
        }
        break;
    }  // switch
  }
};  // class HexaederGaussianQuadrature

// Abstract Base class for HexaederFormFunction classes
template<class T>
class HexaederFormFunctionBase {
 public:
  virtual ~HexaederFormFunctionBase() {}

  // Declaring a virtual function with =0 makes the class abstract. that means no instances can be created from the
  // class HexaederFormFunctionBase...
  // however it is possible to declare a pointer to HexaederFormFunctionBase, the pointer will be pointing to objects of
  // non abstract derived classes..
  // Used with virtual functions, this provides the support for object-oriented polymorphism in c++
  virtual void H(const T *p, T *q) const                                                                           = 0;
  virtual void Hx(const T *p, T *q) const                                                                          = 0;
  virtual void Hy(const T *p, T *q) const                                                                          = 0;
  virtual void Hz(const T *p, T *q) const                                                                          = 0;
  virtual T    JacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz, T *q) const = 0;
  virtual T    DetJacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz) const    = 0;
};


// Stands for Eight-Node hexahedral element See:

// check:
// http://www.colorado.edu/engineering/cas/courses.d/AFEM.d/AFEM.Ch18.d/AFEM.Ch18.index.html

// also:
// The function implimented as in page 240 of the holy book:
// Nonlinear Finite Elements for Continua and Structures
// Ted Belytschko, Wing Kam Liu, Brian Moran
// ISBN-13: 978-0471-98773-1
// http://www.amazon.com/Nonlinear-Finite-Elements-Continua-Structures/dp/0471987735
// Page 626

template<class T>
class HexaederFormFunctionTrilinear : public HexaederFormFunctionBase<T> {
 public:
  enum { DOFs = 8 };
  virtual ~HexaederFormFunctionTrilinear() {}

  inline virtual void H(const T *p, T *q) const {
    const T x1 = *p++, x2 = *p++, x3 = *p;
    const T x1m1 = 1.0-x1, x2m1 = 1.0-x2, x3m1 = 1.0-x3;

    *q++ = x1m1 *x2m1    *x3m1;
    *q++ = x1   *x2m1    *x3m1;
    *q++ = x1   *x2      *x3m1;
    *q++ = x1m1 *x2      *x3m1;
    *q++ = x1m1 *x2m1    *x3;
    *q++ = x1   *x2m1    *x3;
    *q++ = x1   *x2      *x3;
    *q   =  x1m1 *x2      *x3;
  }

  inline virtual void Hx(const T *p, T *q) const {
    // ignore first T p points to
    const T x2 = *(++p)++, x3 = *p;
    const T x2m1 = 1.0-x2, x3m1 = 1.0-x3;

    *q++ = -x2m1    *x3m1;
    *q++ = x2m1    *x3m1;
    *q++ = x2      *x3m1;
    *q++ = -x2      *x3m1;
    *q++ = -x2m1    *x3;
    *q++ = x2m1    *x3;
    *q++ = x2      *x3;
    *q   =  -x2      *x3;
  }

  inline virtual void Hy(const T *p, T *q) const {
    // ignore second T p points to
    const T x1 = *p++, x3 = *(++p);
    const T x1m1 = 1.0-x1, x3m1 = 1.0-x3;

    *q++ = -x1m1 *x3m1;
    *q++ = -x1   *x3m1;
    *q++ = x1   *x3m1;
    *q++ = x1m1 *x3m1;
    *q++ = -x1m1 *x3;
    *q++ = -x1   *x3;
    *q++ = x1   *x3;
    *q   =   x1m1 *x3;
  }

  inline virtual void Hz(const T *p, T *q) const {
    // ignore third T p points to
    const T x1 = *p++, x2 = *p;
    const T x1m1 = 1.0-x1, x2m1 = 1.0-x2;

    *q++ = -x1m1 *x2m1;
    *q++ = -x1   *x2m1;
    *q++ = -x1   *x2;
    *q++ = -x1m1 *x2;
    *q++ = x1m1 *x2m1;
    *q++ = x1   *x2m1;
    *q++ = x1   *x2;
    *q   =   x1m1 *x2;
  }

  inline T virtual JacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz, T *q) const {
    const T a = q[0] = Vec8MultVec8(x, Nx);
    const T b = q[1] = Vec8MultVec8(y, Nx);
    const T c = q[2] = Vec8MultVec8(z, Nx);

    const T d = q[3] = Vec8MultVec8(x, Ny);
    const T e = q[4] = Vec8MultVec8(y, Ny);
    const T f = q[5] = Vec8MultVec8(z, Ny);

    const T g = q[6] = Vec8MultVec8(x, Nz);
    const T h = q[7] = Vec8MultVec8(y, Nz);
    const T i = q[8] = Vec8MultVec8(z, Nz);

    return a*e*i-c*e*g+b*f*g+c*d*h-a*f*h-b*d*i;
  }

  inline T virtual DetJacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz) const {
    const T a = Vec8MultVec8(x, Nx);
    const T b = Vec8MultVec8(y, Nx);
    const T c = Vec8MultVec8(z, Nx);

    const T d = Vec8MultVec8(x, Ny);
    const T e = Vec8MultVec8(y, Ny);
    const T f = Vec8MultVec8(z, Ny);

    const T g = Vec8MultVec8(x, Nz);
    const T h = Vec8MultVec8(y, Nz);
    const T i = Vec8MultVec8(z, Nz);

    return a*e*i-c*e*g+b*f*g+c*d*h-a*f*h-b*d*i;
  }
};  // class HexaederFormFunctionTrilinear

typedef HexaederFormFunctionTrilinear<double> HexaederFormFunctionTrilinearDouble;
typedef HexaederFormFunctionTrilinear<float>  HexaederFormFunctionTrilinearFloat;

// check:
// http://www.colorado.edu/engineering/cas/courses.d/AFEM.d/AFEM.Ch18.d/AFEM.Ch18.index.html

template<class T>
class HexaederFormFunctionTriquadratic : public HexaederFormFunctionBase<T> {
 public:
  enum { DOFs = 27 };
  virtual ~HexaederFormFunctionTriquadratic() {}

  inline void H(const T *p, T *q) const {
    const T x1 = *p++, x2 = *p++, x3 = *p;
    const T x1m1 = x1-1., x2m1 = x2-1., x3m1 = x3-1.;
    const T x12m1 = x1+x1m1, x22m1 = x2+x2m1, x32m1 = x3+x3m1;

    *q++ = x12m1*x1m1*x22m1*x2m1*x32m1*x3m1;
    *q++ = x1*x12m1*x22m1*x2m1*x32m1*x3m1;
    *q++ = x1*x12m1*x2*x22m1*x32m1*x3m1;
    *q++ = x12m1*x1m1*x2*x22m1*x32m1*x3m1;
    *q++ = x12m1*x1m1*x22m1*x2m1*x3*x32m1;
    *q++ = x1*x12m1*x22m1*x2m1*x3*x32m1;
    *q++ = x1*x12m1*x2*x22m1*x3*x32m1;
    *q++ = x12m1*x1m1*x2*x22m1*x3*x32m1;
    *q++ = -4.*x1*x1m1*x22m1*x2m1*x32m1*x3m1;
    *q++ = -4.*x1*x12m1*x2*x2m1*x32m1*x3m1;
    *q++ = -4.*x1*x1m1*x2*x22m1*x32m1*x3m1;
    *q++ = -4.*x12m1*x1m1*x2*x2m1*x32m1*x3m1;
    *q++ = -4.*x12m1*x1m1*x22m1*x2m1*x3*x3m1;
    *q++ = -4.*x1*x12m1*x22m1*x2m1*x3*x3m1;
    *q++ = -4.*x1*x12m1*x2*x22m1*x3*x3m1;
    *q++ = -4.*x12m1*x1m1*x2*x22m1*x3*x3m1;
    *q++ = -4.*x1*x1m1*x22m1*x2m1*x3*x32m1;
    *q++ = -4.*x1*x12m1*x2*x2m1*x3*x32m1;
    *q++ = -4.*x1*x1m1*x2*x22m1*x3*x32m1;
    *q++ = -4.*x12m1*x1m1*x2*x2m1*x3*x32m1;
    *q++ = 16.*x1*x1m1*x2*x2m1*x32m1*x3m1;
    *q++ = 16.*x1*x1m1*x22m1*x2m1*x3*x3m1;
    *q++ = 16.*x1*x12m1*x2*x2m1*x3*x3m1;
    *q++ = 16.*x1*x1m1*x2*x22m1*x3*x3m1;
    *q++ = 16.*x12m1*x1m1*x2*x2m1*x3*x3m1;
    *q++ = 16.*x1*x1m1*x2*x2m1*x3*x32m1;
    *q   = -64.*x1*x1m1*x2*x2m1*x3*x3m1;
  }  // H

  inline void Hx(const T *p, T *q) const {
    const T x1 = *p++, x2 = *p++, x3 = *p;
    const T x1m1 = x1-1., x2m1 = x2-1., x3m1 = x3-1.;
    const T x12m1 = x1+x1m1, x22m1 = x2+x2m1, x32m1 = x3+x3m1;

    *q++ = (-3.+4.*x1)*x22m1*x2m1*x32m1*x3m1;
    *q++ = (-1.+4.*x1)*x22m1*x2m1*x32m1*x3m1;
    *q++ = (-1.+4.*x1)*x2*x22m1*x32m1*x3m1;
    *q++ = (-3.+4.*x1)*x2*x22m1*x32m1*x3m1;
    *q++ = (-3.+4.*x1)*x22m1*x2m1*x3*x32m1;
    *q++ = (-1.+4.*x1)*x22m1*x2m1*x3*x32m1;
    *q++ = (-1.+4.*x1)*x2*x22m1*x3*x32m1;
    *q++ = (-3.+4.*x1)*x2*x22m1*x3*x32m1;
    *q++ = -4.*x12m1*x22m1*x2m1*x32m1*x3m1;
    *q++ = -4.*(-1.+4.*x1)*x2*x2m1*x32m1*x3m1;
    *q++ = -4.*x12m1*x2*x22m1*x32m1*x3m1;
    *q++ = -4.*(-3.+4.*x1)*x2*x2m1*x32m1*x3m1;
    *q++ = -4.*(-3.+4.*x1)*x22m1*x2m1*x3*x3m1;
    *q++ = -4.*(-1.+4.*x1)*x22m1*x2m1*x3*x3m1;
    *q++ = -4.*(-1.+4.*x1)*x2*x22m1*x3*x3m1;
    *q++ = -4.*(-3.+4.*x1)*x2*x22m1*x3*x3m1;
    *q++ = -4.*x12m1*x22m1*x2m1*x3*x32m1;
    *q++ = -4.*(-1.+4.*x1)*x2*x2m1*x3*x32m1;
    *q++ = -4.*x12m1*x2*x22m1*x3*x32m1;
    *q++ = -4.*(-3.+4.*x1)*x2*x2m1*x3*x32m1;
    *q++ = 16.*x12m1*x2*x2m1*x32m1*x3m1;
    *q++ = 16.*x12m1*x22m1*x2m1*x3*x3m1;
    *q++ = 16.*(-1.+4.*x1)*x2*x2m1*x3*x3m1;
    *q++ = 16.*x12m1*x2*x22m1*x3*x3m1;
    *q++ = 16.*(-3.+4.*x1)*x2*x2m1*x3*x3m1;
    *q++ = 16.*x12m1*x2*x2m1*x3*x32m1;
    *q   = -64.*x12m1*x2*x2m1*x3*x3m1;
  }  // Hx

  inline void Hy(const T *p, T *q) const {
    const T x1 = *p++, x2 = *p++, x3 = *p;
    const T x1m1 = x1-1., x2m1 = x2-1., x3m1 = x3-1.;
    const T x12m1 = x1+x1m1, x22m1 = x2+x2m1, x32m1 = x3+x3m1;

    *q++ = x12m1*x1m1*(-3.+4.*x2)*x32m1*x3m1;
    *q++ = x1*x12m1*(-3.+4.*x2)*x32m1*x3m1;
    *q++ = x1*x12m1*(-1.+4.*x2)*x32m1*x3m1;
    *q++ = x12m1*x1m1*(-1.+4.*x2)*x32m1*x3m1;
    *q++ = x12m1*x1m1*(-3.+4.*x2)*x3*x32m1;
    *q++ = x1*x12m1*(-3.+4.*x2)*x3*x32m1;
    *q++ = x1*x12m1*(-1.+4.*x2)*x3*x32m1;
    *q++ = x12m1*x1m1*(-1.+4.*x2)*x3*x32m1;
    *q++ = -4.*x1*x1m1*(-3.+4.*x2)*x32m1*x3m1;
    *q++ = -4.*x1*x12m1*x22m1*x32m1*x3m1;
    *q++ = -4.*x1*x1m1*(-1.+4.*x2)*x32m1*x3m1;
    *q++ = -4.*x12m1*x1m1*x22m1*x32m1*x3m1;
    *q++ = -4.*x12m1*x1m1*(-3.+4.*x2)*x3*x3m1;
    *q++ = -4.*x1*x12m1*(-3.+4.*x2)*x3*x3m1;
    *q++ = -4.*x1*x12m1*(-1.+4.*x2)*x3*x3m1;
    *q++ = -4.*x12m1*x1m1*(-1.+4.*x2)*x3*x3m1;
    *q++ = -4.*x1*x1m1*(-3.+4.*x2)*x3*x32m1;
    *q++ = -4.*x1*x12m1*x22m1*x3*x32m1;
    *q++ = -4.*x1*x1m1*(-1.+4.*x2)*x3*x32m1;
    *q++ = -4.*x12m1*x1m1*x22m1*x3*x32m1;
    *q++ = 16.*x1*x1m1*x22m1*x32m1*x3m1;
    *q++ = 16.*x1*x1m1*(-3.+4.*x2)*x3*x3m1;
    *q++ = 16.*x1*x12m1*x22m1*x3*x3m1;
    *q++ = 16.*x1*x1m1*(-1.+4.*x2)*x3*x3m1;
    *q++ = 16.*x12m1*x1m1*x22m1*x3*x3m1;
    *q++ = 16.*x1*x1m1*x22m1*x3*x32m1;
    *q   = -64.*x1*x1m1*x22m1*x3*x3m1;
  }  // Hy

  inline void Hz(const T *p, T *q) const {
    const T x1 = *p++, x2 = *p++, x3 = *p;
    const T x1m1 = x1-1., x2m1 = x2-1., x3m1 = x3-1.;
    const T x12m1 = x1+x1m1, x22m1 = x2+x2m1, x32m1 = x3+x3m1;

    *q++ = x12m1*x1m1*x22m1*x2m1*(-3.+4.*x3);
    *q++ = x1*x12m1*x22m1*x2m1*(-3.+4.*x3);
    *q++ = x1*x12m1*x2*x22m1*(-3.+4.*x3);
    *q++ = x12m1*x1m1*x2*x22m1*(-3.+4.*x3);
    *q++ = x12m1*x1m1*x22m1*x2m1*(-1.+4.*x3);
    *q++ = x1*x12m1*x22m1*x2m1*(-1.+4.*x3);
    *q++ = x1*x12m1*x2*x22m1*(-1.+4.*x3);
    *q++ = x12m1*x1m1*x2*x22m1*(-1.+4.*x3);
    *q++ = -4.*x1*x1m1*x22m1*x2m1*(-3.+4.*x3);
    *q++ = -4.*x1*x12m1*x2*x2m1*(-3.+4.*x3);
    *q++ = -4.*x1*x1m1*x2*x22m1*(-3.+4.*x3);
    *q++ = -4.*x12m1*x1m1*x2*x2m1*(-3.+4.*x3);
    *q++ = -4.*x12m1*x1m1*x22m1*x2m1*x32m1;
    *q++ = -4.*x1*x12m1*x22m1*x2m1*x32m1;
    *q++ = -4.*x1*x12m1*x2*x22m1*x32m1;
    *q++ = -4.*x12m1*x1m1*x2*x22m1*x32m1;
    *q++ = -4.*x1*x1m1*x22m1*x2m1*(-1.+4.*x3);
    *q++ = -4.*x1*x12m1*x2*x2m1*(-1.+4.*x3);
    *q++ = -4.*x1*x1m1*x2*x22m1*(-1.+4.*x3);
    *q++ = -4.*x12m1*x1m1*x2*x2m1*(-1.+4.*x3);
    *q++ = 16.*x1*x1m1*x2*x2m1*(-3.+4.*x3);
    *q++ = 16.*x1*x1m1*x22m1*x2m1*x32m1;
    *q++ = 16.*x1*x12m1*x2*x2m1*x32m1;
    *q++ = 16.*x1*x1m1*x2*x22m1*x32m1;
    *q++ = 16.*x12m1*x1m1*x2*x2m1*x32m1;
    *q++ = 16.*x1*x1m1*x2*x2m1*(-1.+4.*x3);
    *q   = -64.*x1*x1m1*x2*x2m1*x32m1;
  }  // Hz

  inline T JacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz, T *q) const {
    const T a = q[0] = VecMultVec(x, Nx, DOFs);
    const T b = q[1] = VecMultVec(y, Nx, DOFs);
    const T c = q[2] = VecMultVec(z, Nx, DOFs);

    const T d = q[3] = VecMultVec(x, Ny, DOFs);
    const T e = q[4] = VecMultVec(y, Ny, DOFs);
    const T f = q[5] = VecMultVec(z, Ny, DOFs);

    const T g = q[6] = VecMultVec(x, Nz, DOFs);
    const T h = q[7] = VecMultVec(y, Nz, DOFs);
    const T i = q[8] = VecMultVec(z, Nz, DOFs);

    return a*e*i-c*e*g+b*f*g+c*d*h-a*f*h-b*d*i;
  }

  inline T DetJacobiMatrix(const T *x, const T *y, const T *z, const T *Nx, const T *Ny, const T *Nz) const {
    const T a = VecMultVec(x, Nx, DOFs);
    const T b = VecMultVec(y, Nx, DOFs);
    const T c = VecMultVec(z, Nx, DOFs);

    const T d = VecMultVec(x, Ny, DOFs);
    const T e = VecMultVec(y, Ny, DOFs);
    const T f = VecMultVec(z, Ny, DOFs);

    const T g = VecMultVec(x, Nz, DOFs);
    const T h = VecMultVec(y, Nz, DOFs);
    const T i = VecMultVec(z, Nz, DOFs);

    return a*e*i-c*e*g+b*f*g+c*d*h-a*f*h-b*d*i;
  }
};  // class HexaederFormFunctionTriquadratic

typedef HexaederFormFunctionTriquadratic<double> HexaederFormFunctionTriquadraticDouble;
typedef HexaederFormFunctionTriquadratic<float>  HexaederFormFunctionTriquadraticFloat;

// HexaederFormFunction is a template class that takes four template parameters
// Class T, defines the type of values represented in the class.
// const int Samples, is a number that is infact a template parameter of the parent class
// HexaederGaussianQuadrature<T,Samples>... check the HexahaederGaussianQuadrature class for more info
// class F, defines the parent class from which HexaederFormFunctions inherits beside the HexaederGaussianQuadrature
// template class. This makes it possible to adapt member functions and variables according to the object declaration.
// const int DOFs, with a default value of template class F::DOFs, an enumeration that is declared in
// HexaederFormFunctionTrilinear and HexaederFormFunctionTriquadratic

template<class T, const int Samples, class F, const int DOFs = F::DOFs>
class HexaederFormFunctions : public HexaederGaussianQuadrature<T, Samples>, public F {
 public:
  T N[Samples*Samples*Samples][DOFs];
  T Nx[Samples*Samples*Samples][DOFs];
  T Ny[Samples*Samples*Samples][DOFs];
  T Nz[Samples*Samples*Samples][DOFs];

  // constructor function
  HexaederFormFunctions() {
    for (int i = 0; i < this->Samples3; i++) {
      this->H(this->gp[i], N[i]);
      this->Hx(this->gp[i], Nx[i]);
      this->Hy(this->gp[i], Ny[i]);
      this->Hz(this->gp[i], Nz[i]);
      /*
         fprintf(stderr, "N[%d]: %lf %lf %lf %lf %lf %lf %lf %lf\n", i, N[i][0], N[i][1], N[i][2], N[i][3], N[i][4],
            N[i][5], N[i][6], N[i][7]);
         fprintf(stderr, "Nx[%d]: %lf %lf %lf %lf %lf %lf %lf %lf\n", i, Nx[i][0], Nx[i][1], Nx[i][2], Nx[i][3],
            Nx[i][4], Nx[i][5], Nx[i][6], Nx[i][7]);
         fprintf(stderr, "Ny[%d]: %lf %lf %lf %lf %lf %lf %lf %lf\n", i, Ny[i][0], Ny[i][1], Ny[i][2], Ny[i][3],
            Ny[i][4], Ny[i][5], Ny[i][6], Ny[i][7]);
         fprintf(stderr, "Nz[%d]: %lf %lf %lf %lf %lf %lf %lf %lf\n", i, Nz[i][0], Nz[i][1], Nz[i][2], Nz[i][3],
            Nz[i][4], Nz[i][5], Nz[i][6], Nz[i][7]);
       */
    }
  }

  virtual ~HexaederFormFunctions() {}

  inline void GaussianQuadratureVolume(const T *x, const T *y, const T *z, T *mass) const {
    VecSetVal(mass, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i];
      T detJ = this->DetJacobiMatrix(x, y, z, Nxi, Nyi, Nzi);

      detJ *= this->W[i];

      T *pmass = mass;
      const T *pemass = mass+DOFs, *Ni = N[i];
      for (; pmass < pemass; Ni++)
        *pmass++ += detJ **Ni;
    }
  }

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, const int p, T *stiff) const {
    T J[9], Ji[9];

    VecSetVal(stiff, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i];
      T detJ = this->JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);

      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];
      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];
      const T uxp = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);

      T *pstiff        = stiff;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += uxp*ux+uyp*uy+uzp*uz;
      }
    }
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, const T sigma, const int p,
                                        T *stiff) const {
    T J[9], Ji[9];

    VecSetVal(stiff, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i];
      T detJ = this->JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);

      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];
      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];
      const T uxp = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);

      T *pstiff        = stiff;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += uxp*ux+uyp*uy+uzp*uz;
      }
    }
    T *pstiff        = stiff;
    const T *pestiff = stiff+DOFs;
    for (; pstiff < pestiff;)
      *pstiff++ *= sigma;
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, T *s, const int p, T *stiff) const {
    //! hack for xCellent/border element
    for (int i = 0; i < DOFs; i++)
      if (s[i] < 1e-7) {
        GaussianQuadraturePoisson(x, y, z, s[i], p, stiff);
        return;
      }

    T J[9], Ji[9];

    VecSetVal(stiff, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i], *Ni = N[i];
      T detJ = this->JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);
      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];

      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];

      const T uxp = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);

      const T *ps = s;
      const T *pNi = Ni, *peNi = pNi+DOFs;
      T sigma = 0.;

      // Interpolate conductivity at sample point
      for (; pNi < peNi; ps++, pNi++)
        sigma += *pNi * (*ps);

      const T uxyzpabc = uxp*sigma;
      const T uxyzpbde = uyp*sigma;
      const T uxyzpcef = uzp*sigma;

      T *pstiff        = stiff;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += ux*uxyzpabc+uy*uxyzpbde+uz*uxyzpcef;
      }
    }
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, const int p, T *stiff, T *mass) const {
    T J[9], Ji[9];

    VecSetVal(stiff, (T)0.0, DOFs);
    VecSetVal(mass, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i];
      T detJ = JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);
      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];
      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];

      const T uxp = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);

      const T *Ni      = N[i];
      const T  NipdetJ = detJ*Ni[p];

      T *pstiff        = stiff;
      T *pmass         = mass;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += uxp*ux+uyp*uy+uzp*uz;
        *mass++   += NipdetJ **Ni++;
      }
    }
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, const SymmetricTensor2<T> &s, const int p,
                                        T *stiff) const {
    T J[9], Ji[9];

    T a = s.a;
    T b = s.b;
    T c = s.c;
    T d = s.d;
    T e = s.e;
    T f = s.f;

    VecSetVal(stiff, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i];
      T detJ = this->JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);
      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];

      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];

      const T uxp      = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp      = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp      = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);
      const T uxyzpabc = uxp*a+uyp*b+uzp*c;
      const T uxyzpbde = uxp*b+uyp*d+uzp*e;
      const T uxyzpcef = uxp*c+uyp*e+uzp*f;

      T *pstiff        = stiff;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += ux*uxyzpabc+uy*uxyzpbde+uz*uxyzpcef;
      }
    }
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, const SymmetricTensor2<T> &s, const int p,
                                        T *stiff, T *mass) const {
    T J[9], Ji[9];

    T a = s.a;
    T b = s.b;
    T c = s.c;
    T d = s.d;
    T e = s.e;
    T f = s.f;

    VecSetVal(stiff, (T)0.0, DOFs);
    VecSetVal(mass, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i], *Ni = N[i];
      T detJ = JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);
      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];

      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];

      const T uxp      = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp      = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp      = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);
      const T uxyzpabc = uxp*a+uyp*b+uzp*c;
      const T uxyzpbde = uxp*b+uyp*d+uzp*e;
      const T uxyzpcef = uxp*c+uyp*e+uzp*f;

      const T NipdetJ = detJ*Ni[p];

      T *pstiff        = stiff;
      T *pmass         = mass;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += ux*uxyzpabc+uy*uxyzpbde+uz*uxyzpcef;
        *mass++   += NipdetJ **Ni++;
      }
    }
  }  // GaussianQuadraturePoisson

  inline void GaussianQuadraturePoisson(const T *x, const T *y, const T *z, SymmetricTensor2<T> **s, const int p,
                                        T *stiff) const {
    //! hack for xCellent/border element
    for (int i = 0; i < DOFs; i++)
      if (s[i]->a < 1e-7) {
        GaussianQuadraturePoisson(x, y, z, *s[i], p, stiff);
        return;
      }

    T J[9], Ji[9];

    VecSetVal(stiff, (T)0.0, DOFs);

    for (int i = 0; i < this->Samples3; i++) {
      const T *Nxi = Nx[i], *Nyi = Ny[i], *Nzi = Nz[i], *Ni = N[i];
      T detJ = this->JacobiMatrix(x, y, z, Nxi, Nyi, Nzi, J);
      InverseMatrix3x3(J, Ji, detJ);
      detJ *= this->W[i];

      const T &Ji0 = Ji[0];
      const T &Ji1 = Ji[1];
      const T &Ji2 = Ji[2];
      const T &Ji3 = Ji[3];
      const T &Ji4 = Ji[4];
      const T &Ji5 = Ji[5];
      const T &Ji6 = Ji[6];
      const T &Ji7 = Ji[7];
      const T &Ji8 = Ji[8];

      const T Nxip = Nxi[p], Nyip = Nyi[p], Nzip = Nzi[p];

      const T uxp = detJ*(Nxip*Ji0+Nyip*Ji1+Nzip*Ji2);
      const T uyp = detJ*(Nxip*Ji3+Nyip*Ji4+Nzip*Ji5);
      const T uzp = detJ*(Nxip*Ji6+Nyip*Ji7+Nzip*Ji8);

      SymmetricTensor2<T> **ps = s;
      const T *pNi = Ni, *peNi = pNi+DOFs;
      T a, b, c, d, e, f;

      // Interpolate conductivity at sample point
      a = 0., b = 0., c = 0., d = 0., e = 0., f = 0.;
      for (; pNi < peNi; ps++, pNi++) {
        a += *pNi * (*ps)->a;
        b += *pNi * (*ps)->b;
        c += *pNi * (*ps)->c;
        d += *pNi * (*ps)->d;
        e += *pNi * (*ps)->e;
        f += *pNi * (*ps)->f;
      }

      // }
      const T uxyzpabc = uxp*a+uyp*b+uzp*c;
      const T uxyzpbde = uxp*b+uyp*d+uzp*e;
      const T uxyzpcef = uxp*c+uyp*e+uzp*f;

      T *pstiff        = stiff;
      const T *pestiff = stiff+DOFs;
      for (; pstiff < pestiff; Nxi++, Nyi++, Nzi++) {
        const T ux = *Nxi*Ji0+*Nyi*Ji1+*Nzi*Ji2;
        const T uy = *Nxi*Ji3+*Nyi*Ji4+*Nzi*Ji5;
        const T uz = *Nxi*Ji6+*Nyi*Ji7+*Nzi*Ji8;
        *pstiff++ += ux*uxyzpabc+uy*uxyzpbde+uz*uxyzpcef;
      }
    }
  }  // GaussianQuadraturePoisson
};  // class HexaederFormFunctions

// type definitions for a friendly api interface :-)
typedef HexaederFormFunctions<float, 2, HexaederFormFunctionTrilinearFloat>      HexaederFFTrilinearFloat;
typedef HexaederFormFunctions<double, 2, HexaederFormFunctionTrilinearDouble>    HexaederFFTrilinearDouble;
typedef HexaederFormFunctions<float, 2, HexaederFormFunctionTriquadraticFloat>   HexaederFFTriquadraticFloat;
typedef HexaederFormFunctions<double, 2, HexaederFormFunctionTriquadraticDouble> HexaederFFTriquadraticDouble;

#endif  // ifndef FORMFUNCTIONS_H
