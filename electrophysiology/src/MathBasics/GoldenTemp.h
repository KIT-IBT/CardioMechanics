/*
   routine for golden section search

   (template version of golden.c from NR)

 */

#ifndef _GoldenTemp_h
#define _GoldenTemp_h

#include <cmath>
#define SHFT2(a, b, c) (a)    = (b); (b) = (c);
#define SHFT3(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);

// Description:
// recommended precision
template<class T>
class GoldenTempPrecision {
 public:
  GoldenTempPrecision() {}

  ~GoldenTempPrecision() {}

  // Description:
  // returns recommended precision
  inline T operator()(void);
};

template<>
inline float GoldenTempPrecision<float>::operator()(void) {
  return 1e-3;
}

template<>
inline double GoldenTempPrecision<double>::operator()(void) {
  return 1e-7;
}

// Description:
// Golden search in on dimension. See NR
// U & u is a functional which returns the value
// it must have operator(T)

template<class T, class U>

T GoldenTemp(T ax, T bx, T cx, U &u, T tol, T &xmin) {
  const T R = 0.6180339887498949;
  const T C = 1.0-R;

  T f1, f2, x0, x1, x2, x3;

  x0 = ax;
  x3 = cx;
  if (fabs(cx-bx) > fabs(bx-ax)) {
    x1 = bx;
    x2 = bx+C*(cx-bx);
  } else {
    x2 = bx;
    x1 = bx-C*(bx-ax);
  }
  f1 = u(x1);
  f2 = u(x2);
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) {
      SHFT3(x0, x1, x2, R*x1+C*x3)
      SHFT2(f1, f2, u(x2))
    } else {
      SHFT3(x3, x2, x1, R*x2+C*x0)
      SHFT2(f2, f1, u(x1))
    }
  }
  if (f1 < f2) {
    xmin = x1;
    return f1;
  } else {
    xmin = x2;
    return f2;
  }
}  // GoldenTemp

#undef SHFT2
#undef SHFT3

#endif  // ifndef _GoldenTemp_h

/* (C) Copr. 1986-92 Numerical Recipes Software ... */
