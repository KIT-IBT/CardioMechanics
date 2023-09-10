/*
   SVDTemp.h

   contents:
   svdcmp - SVD
   svdbksb - backsubstitution for SVD


   last modified:
   os Mar 15, 2001

 */
#ifndef _SVDTemp_h
#define _SVDTemp_h


#include "kaMatrix.h"

#define IMIN(a, b) ((a) < (b) ? (a) : (b))

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

template<class T> inline T pythag(T a, T b) {
  T absa, absb;

  absa = fabs(a);
  absb = fabs(b);

  T sqrarg;

  if (absa > absb)
    return absa*sqrt(1.0+((sqrarg = absb/absa) == 0.0 ? 0.0 : sqrarg *sqrarg) );
  else
    return absb == 0.0 ? 0.0 : absb *sqrt(1.0+((sqrarg = absa/absb) == 0.0 ? 0.0 : sqrarg*sqrarg));
}

/*Description:
     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).

     THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
          T
     A=USV  OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
     GRSVD ASSUMES THAT A COPY OF THE MATRIX A IS IN THE ARRAY U. IT
     ALSO ASSUMES M .GE. N.  IF M .LT. N, THEN COMPUTE THE SINGULAR
                             T       T    T             T
     VALUE DECOMPOSITION OF A .  IF A =UWV  , THEN A=VWU  .

     GRSVD CAN ALSO BE USED TO COMPUTE THE MINIMAL LENGTH LEAST SQUARES
     SOLUTION TO THE OVERDETERMINED LINEAR SYSTEM A*X=B.

 */
template<class T>  void svdcmp(kaMatrix<T> &a, T *w, kaMatrix<T> &v) {
  const int m = a.NumRows();
  const int n = a.NumCols();

  int flag, i, its, j, jj, k, l, nm;

  T anorm, c, f, g, h, s, scale, x, y, z;

  const T one  = 1.0;
  const T zero = 0.0;

  T *rv1 = new T[n];

  g = scale = anorm = zero;
  for (i = 0; i < n; i++) {
    l      = i+1;
    rv1[i] = scale*g;
    g      = s = scale = zero;
    if (i <  m) {
      for (k = i; k < m; k++) scale += fabs(a.Get(k, i));
      if (scale) {
        for (k = i; k < m; k++) {
          a.a(k, i) /= scale;
          s         += a.Get(k, i)*a.Get(k, i);
        }
        f         = a.Get(i, i);
        g         = -SIGN(sqrt(s), f);
        h         = f*g-s;
        a.a(i, i) = f-g;
        for (j = l; j < n; j++) {
          for (s = zero, k = i; k < m; k++) s += a.Get(k, i)*a.Get(k, j);
          f = s/h;
          for (k = i; k < m; k++) a.a(k, j) += f*a.Get(k, i);
        }
        for (k = i; k < m; k++) a.a(k, i) *= scale;
      }
    }
    w[i] = scale *g;
    g    = s = scale = zero;
    if ((i < m) && (i != n-1)) {
      for (k = l; k < n; k++) scale += fabs(a.Get(i, k));
      if (scale) {
        for (k = l; k < n; k++) {
          a.a(i, k) /= scale;
          s         += a.Get(i, k)*a.Get(i, k);
        }
        f         = a.Get(i, l);
        g         = -SIGN(sqrt(s), f);
        h         = f*g-s;
        a.a(i, l) = f-g;
        for (k = l; k < n; k++) rv1[k] = a.Get(i, k)/h;
        for (j = l; j < m; j++) {
          for (s = zero, k = l; k < n; k++) s += a.Get(j, k)*a.Get(i, k);
          for (k = l; k < n; k++) a.a(j, k) += s*rv1[k];
        }
        for (k = l; k < n; k++) a.a(i, k) *= scale;
      }
    }

    T maxarg1 = anorm;
    T maxarg2 = fabs(w[i])+fabs(rv1[i]);

    anorm = maxarg1 > maxarg2 ? maxarg1 : maxarg2;
  }
  for (i = n-1; i >= 0; i--) {
    if (i < n-1) {
      if (g) {
        for (j = l; j < n; j++) v.a(j, i) = (a.Get(i, j)/a.Get(i, l))/g;
        for (j = l; j < n; j++) {
          for (s = zero, k = l; k < n; k++) s += a.Get(i, k)*v.Get(k, j);
          for (k = l; k < n; k++) v.a(k, j) += s*v.Get(k, i);
        }
      }
      for (j = l; j < n; j++) v.a(i, j) = v.a(j, i) = zero;
    }
    v.a(i, i) = one;
    g         = rv1[i];
    l         = i;
  }
  for (i = IMIN(m, n)-1; i >= 0; i--) {
    l = i+1;
    g = w[i];
    for (j = l; j < n; j++) a.a(i, j) = zero;
    if (g) {
      g = one/g;
      for (j = l; j < n; j++) {
        for (s = zero, k = l; k < m; k++) s += a.Get(k, i)*a.Get(k, j);
        f = (s/a.Get(i, i))*g;
        for (k = i; k < m; k++) a.a(k, j) += f*a.Get(k, i);
      }
      for (j = i; j < m; j++) a.a(j, i) *= g;
    } else {for (j = i; j < m; j++) a.a(j, i) = zero; }
    ++a.a(i, i);
  }
  for (k = n-1; k >= 0; k--) {
    for (its = 0; its < 30; its++) {
      flag = 1;
      for (l = k; l >= 0; l--) {
        nm = l-1;
        if ((T)(fabs(rv1[l])+anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((T)(fabs(w[nm])+anorm) == anorm)
          break;
      }
      if (flag) {
        c = zero;
        s = one;
        for (i = l; i <= k; i++) {
          f      = s*rv1[i];
          rv1[i] = c*rv1[i];
          if ((T)(fabs(f)+anorm) == anorm)
            break;
          g    = w[i];
          h    = pythag(f, g);
          w[i] = h;
          h    = one/h;
          c    = g*h;
          s    = -f*h;
          for (j = 0; j < m; j++) {
            y          = a.Get(j, nm);
            z          = a.Get(j, i);
            a.a(j, nm) = y*c+z*s;
            a.a(j, i)  = z*c-y*s;
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < zero) {
          w[k] = -z;
          for (j = 0; j < n; j++) v.a(j, k) = -v.Get(j, k);
        }
        break;
      }
      if (its == 49) {printf("no convergence in 30 svdcmp<T> iterations"); delete[] rv1; exit(1);}
      x  = w[l];
      nm = k-1;
      y  = w[nm];
      g  = rv1[nm];
      h  = rv1[k];
      f  = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g  = pythag(f, one);
      f  = ((x-z)*(x+z)+h*((y/(f+SIGN(g, f)))-h))/x;
      c  = s = one;
      for (j = l; j <= nm; j++) {
        i      = j+1;
        g      = rv1[i];
        y      = w[i];
        h      = s*g;
        g      = c*g;
        z      = pythag(f, h);
        rv1[j] = z;
        c      = f/z;
        s      = h/z;
        f      = x*c+g*s;
        g      = g*c-x*s;
        h      = y*s;
        y     *= c;
        for (jj = 0; jj < n; jj++) {
          x          = v.Get(jj, j);
          z          = v.Get(jj, i);
          v.a(jj, j) = x*c+z*s;
          v.a(jj, i) = z*c-x*s;
        }
        z    = pythag(f, h);
        w[j] = z;
        if (z) {
          z = one/z;
          c = f*z;
          s = h*z;
        }
        f = c*g+s*y;
        x = c*y-s*g;
        for (jj = 0; jj < m; jj++) {
          y          = a.Get(jj, j);
          z          = a.Get(jj, i);
          a.a(jj, j) = y*c+z*s;
          a.a(jj, i) = z*c-y*s;
        }
      }
      rv1[l] = zero;
      rv1[k] = f;
      w[k]   = x;
    }
  }

  // sort singular values and exchange columns of u and v
  // selection sort minimizes swapping of u and v

  int id;
  T   temp;

  for (i = 0; i < n-1; i++) {
    // find index of maximum singular value
    id = i;

    for (j = i+1; j < n; j++) {
      if (w[j] > w[id])
        id = j;
    }

    if (id != i) {
      // swap singular values and vectors
      temp  = w[i];
      w[i]  = w[id];
      w[id] = temp;
      a.SwapColumns(i, id);
      v.SwapColumns(i, id);
    }
  }


  delete[] rv1;
}  // svdcmp

// Description:
// Backsubstitution routine.
// solves A.X = B for vector X
// A is specified by matrix U(m x n), array w[0..n], and matrix V(n x n)
// as returned by svdcmp
//
// b[0..m] is the input right-hand side
// x[0..n] is the output solution vector

template<class T> void svdbksb(kaMatrix<T> &u, T *w, kaMatrix<T> &v, T *b, T *x) {
  const unsigned int m = u.NumRows();
  const unsigned int n = u.NumCols();

  int jj, j, i;

  const T zero = 0.0;

  T s;

  T *tmp = new T[n];

  for (j = 0; j < n; j++) {  // calcuate U'B
    s = zero;
    if (w[j]) {  // nonzero result only if w[j] is nonzero
      for (i = 0; i < m; i++) s += u.Get(i, j)*b[i];
      s /= w[j];
    }
    tmp[j] = s;
  }

  for (j = 0; j < n; j++) {  // matrix multiply by V to get answer
    s = zero;
    for (jj = 0; jj < n; jj++) s += v.Get(j, jj)*tmp[jj];
    x[j] = s;
  }

  delete[] tmp;
}  // svdbksb

#endif  // ifndef _SVDTemp_h
