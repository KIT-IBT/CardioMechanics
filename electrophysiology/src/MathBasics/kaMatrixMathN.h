/*! \file kaMatrixMathN.h
   \brief Class with mathematical functions of NxN matrices

   \author cw,fs,os, IBT - Universität Karlsruhe (TH)
 */
#ifndef KAMATRIXMATHN_H
#define KAMATRIXMATHN_H

#include <kaMatrix.h>

#define NOT_USED(x) (void)(x)

//! Class for mathematical handling of quadratic matrices
/*!
   kaSharedMathN includes variant methods for matrix mathematics.

   \author cw,fs,os, IBT - Universität Karlsruhe (TH)
 */

template<class T> class kaMatrixMathN : public kaMatrix<T> {
 public:
  //! Set matrix to identity
  inline void Identity();

  //! Transpose matrix
  inline void Transpose();

  //! Get determinant of matrix
  T Determinant() const;

  //! Get dimension of matrix
  virtual unsigned int Dimension() const = 0;

  //! Transform translation in homogeneous matrix
  void Translate(T, T, T);

  //! Transform scaling in homogeneous matrix
  void Scale(T, T, T);

  //! Transform frustum in homogeneous matrix
  void Frustum(T, T, T, T, T, T);

  //! Transform ortho in homogeneous matrix
  void Ortho(T, T, T, T, T, T);
};  // class kaMatrixMathN

template<class T>
inline void kaMatrixMathN<T>::Identity() {
  unsigned int n = Dimension();
  unsigned int i;

  kaMatrix<T>::Null();
  for (i = 0; i < n; i++)
    kaMatrix<T>::a(i, i) = 1.0;
}

template<class T>
inline void kaMatrixMathN<T>::Transpose() {
  unsigned int n = Dimension();
  unsigned int i, j;
  T t;

  for (i = 0; i < n; i++)
    for (j = i+1; j < n; j++) {
      t                    = kaMatrix<T>::Get(i, j);
      kaMatrix<T>::a(i, j) = kaMatrix<T>::Get(j, i);
      kaMatrix<T>::a(j, i) = t;
    }
}

template<class T>
T kaMatrixMathN<T>::Determinant() const {
  unsigned int n = Dimension();

  NOT_USED(n);
  assert(n == 3);

  return (T)(kaMatrix<T>::Get(0)*kaMatrix<T>::Get(4)*kaMatrix<T>::Get(8) + kaMatrix<T>::Get(1)*kaMatrix<T>::Get(5)*
             kaMatrix<T>::Get(6) + kaMatrix<T>::Get(2)*kaMatrix<T>::Get(3)*kaMatrix<T>::Get(7)
             - kaMatrix<T>::Get(2)*kaMatrix<T>::Get(4)*kaMatrix<T>::Get(6) - kaMatrix<T>::Get(5)*kaMatrix<T>::Get(7)*
             kaMatrix<T>::Get(0) - kaMatrix<T>::Get(8)*kaMatrix<T>::Get(1)*kaMatrix<T>::Get(3) );
}

template<class T>
void kaMatrixMathN<T>::Translate(T x, T y, T z) {
  unsigned int n = Dimension();

  NOT_USED(n);
  assert(n == 4);

  kaMatrixMathN<T>::Identity();
  kaMatrix<T>::a(0, 3) = x;
  kaMatrix<T>::a(1, 3) = y;
  kaMatrix<T>::a(2, 3) = z;
}

template<class T>
void kaMatrixMathN<T>::Scale(T x, T y, T z) {
  unsigned int n = Dimension();

  NOT_USED(n);
  assert(n == 4);

  kaMatrixMathN<T>::Identity();
  kaMatrix<T>::a(0, 0) = x;
  kaMatrix<T>::a(1, 1) = y;
  kaMatrix<T>::a(2, 2) = z;
}

template<class T>
void kaMatrixMathN<T>::Frustum(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  unsigned int n = Dimension();

  NOT_USED(n);
  assert(n == 4);

  kaMatrix<T>::Null();
  kaMatrix<T>::a(0, 0) = 2.0*near_arg/(right-left);
  kaMatrix<T>::a(0, 2) = (right+left)/(right-left);
  kaMatrix<T>::a(1, 1) = 2.0*near_arg/(top-bottom);
  kaMatrix<T>::a(1, 2) = (top+bottom)/(top-bottom);
  kaMatrix<T>::a(2, 2) = -(far_arg+near_arg)/(far_arg-near_arg);
  kaMatrix<T>::a(2, 3) = -2.0*far_arg*near_arg/(far_arg-near_arg);
  kaMatrix<T>::a(3, 2) = -1.0;
}

template<class T> void kaMatrixMathN<T>::Ortho(T left, T right, T bottom, T top, T near_arg, T far_arg) {
  unsigned int n = Dimension();

  NOT_USED(n);
  assert(n == 4);

  kaMatrix<T>::Null();
  kaMatrix<T>::a(0, 0) = 2.0/(right-left);
  kaMatrix<T>::a(0, 3) = -(right+left)/(right-left);
  kaMatrix<T>::a(1, 1) = 2.0/(top-bottom);
  kaMatrix<T>::a(1, 3) = -(top+bottom)/(top-bottom);
  kaMatrix<T>::a(2, 2) = -2.0/(far_arg-near_arg);
  kaMatrix<T>::a(2, 3) = -(far_arg+near_arg)/(far_arg-near_arg);
  kaMatrix<T>::a(3, 3) = 1.0;
}

#endif  // ifndef KAMATRIXMATHN_H
