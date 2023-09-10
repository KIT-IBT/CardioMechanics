/*! \file kaMatrixN.h
   \brief Class for handling of quadratic matrices

   \author fs,os, IBT - Universität Karlsruhe (TH)
 */

#ifndef KAMATRIXN_H
#define KAMATRIXN_H

#include "kaMatrixMathN.h"

template<class T, const unsigned int n> class kaMatrixN : public kaMatrixMathN<T> {
 protected:
  T pElements[n*n];

 public:
  //! Dimension of matrix
  virtual inline unsigned int Dimension() const {return n;}

  //! Number of rows (=Dimension)
  virtual inline unsigned int NumRows() const {return n;}

  //! Number of columns (=Dimension)
  virtual inline unsigned int NumCols() const {return n;}

  //! Pointer to elements
  virtual inline T *Pointer() {return pElements;}

  //! Const pointer to elements
  virtual inline const T *Pointer() const {return pElements;}

  // For efficiency - perhaps not necessary in future compiler versions
  virtual inline T & a(const unsigned int z, const unsigned int s) {return (T &)*(pElements+z*n+s);}

  virtual inline T & a(const unsigned int p) {return (T &)*(pElements+p);}

  virtual inline T *operator[](const unsigned int i) {return pElements+i*n;}

  virtual inline const T *operator[](const unsigned int i) const {return pElements+i*n;}

  virtual inline const T Get(const unsigned int z, const unsigned int s) const {return *(pElements+z*n+s);}

  virtual inline const T Get(const unsigned int p) const {return *(pElements+p);}

  inline const kaMatrixN operator*(const T) const;
  inline const kaMatrixN operator/(const T);

  inline const kaMatrixN operator*(const kaMatrixN &) const;
  inline const kaMatrixN operator+(const kaMatrixN &) const;
  void                   XMatrixN(const kaMatrixN<T, n> &);
  void                   MatrixNX(const kaMatrixN<T, n> &);
  kaMatrixN            & operator=(const kaMatrixN &);

  kaMatrixN & operator*=(const T);
  kaMatrixN & operator/=(const T);

  kaMatrixN & operator+=(const kaMatrixN &);
  kaMatrixN & operator-=(const kaMatrixN &);
  kaMatrixN & operator*=(const kaMatrixN &);
  int         Inverse();

  void Rotate(T, T, T, T);
  void XTranslate(T, T, T);
  void XScale(T, T, T);
  void XRotate(T, T, T, T);
  void XFrustum(T, T, T, T, T, T);
  void XOrtho(T, T, T, T, T, T);
  void TranslateX(T, T, T);
  void ScaleX(T, T, T);
  void FrustumX(T, T, T, T, T, T);
  void OrthoX(T, T, T, T, T, T);
  void RotateX(T, T, T, T);
  virtual inline const T Trace();

  friend class kaMatrixMathN<T>;
  T Determinant();
};  // class kaMatrixN

//! kaMatrixN = kaMatrixN.
/*!
 * Set kaMatrixN to kaMatrixN rm.
 * overload = operator
 * \param rm kaMatrixN to set
 */

template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T,
                                                                           n>::operator=(const kaMatrixN<T, n> &rm) {
  T *ps        = pElements;
  T *pe        = ps + n*n;
  const T *pss = rm.pElements;

  while (ps < pe)
    *ps++ = *pss++;
  return *this;
}

//! kaMatrixN *= scalar.
/*!
 * Multiply kaMatrixN with scalar f.
 * overload *= operator
 */

template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T, n>::operator*=(const T f) {
  T *ps = pElements, *pse = ps+n*n;

  while (ps < pse)
    *ps++ *= f;

  return *this;
}

//! kaMatrixN /= scalar.
/*!
 * Divide kaMatrixN with scalar f.
 * overload /= operator
 */

template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T, n>::operator/=(const T f) {
  return *this *= 1/f;
}

//! kaMatrixN * scalar.
/*!
 * Multiply tkaMatrixN with scalar f.
 *
 */


template<class T, const unsigned int n> inline const kaMatrixN<T, n> kaMatrixN<T, n>::operator*(const T f) const {
  kaMatrixN<T, n> temp;

  T *ps = temp.pElements, *pse = (temp.pElements+n*n);
  const T *pss = pElements;

  while (ps < pse)
    *ps++ = *pss++ *f;

  return temp;
}

//! kaMatrixN multiplication kaMatrixN.
/*!
 * Multiply the given kaMatrixN objects. overload * operator.
 * Attention: a*b<>b*a in matrix multiplication
 * \sa XMatrixN, MatrixNX, *=
 */

template<class T, const unsigned int n> inline const kaMatrixN<T, n> kaMatrixN<T, n>::operator*(const kaMatrixN<T,
                                                                                                                n> &rm)
const {
  kaMatrixN<T, n> temp;
  T *pTemp = temp.pElements;

  unsigned int i, j;

  const T *pl = pElements;
  for (i = 0; i < n; i++, pl += n) {
    const T *pr = rm.pElements;
    for (j = 0; j < n; j++, pr++) {
      T val = 0;
      const T *plk = pl, *plke = pl+n, *prk = pr;
      for (; plk < plke; plk++, prk += n)
        val += *plk **prk;
      *pTemp++ = val;
    }
  }
  return temp;
}

//! kaMatrixN + kaMatrixN.
/*! Add two kaMatrixN objects
 * overload the + operator
 */

template<class T, const unsigned int n> inline const kaMatrixN<T, n> kaMatrixN<T, n>::operator+(const kaMatrixN<T,
                                                                                                                n> &rm)
const {
  kaMatrixN<T, n> temp;

  T *ps = temp.pElements, *pse = (temp.pElements+n*n);
  const T *psa = pElements, *psb = rm.pElements;

  while (ps < pse)
    *ps++ = *psa++ + *psb++;

  return temp;
}

//! kaMatrixN / scalar.
/*!
 * Divide the kaMatrixN by scalar f
 * overload the / operator
 */

template<class T, const unsigned int n> inline const kaMatrixN<T, n> kaMatrixN<T, n>::operator/(const T f) {
  return *this*(1/f);
}

//! kaMatrixN *= kaMatrixN.
/*!
 * Multiply given kaMatrixN objects from the RIGHT side (this*rm).
 * overload the *= operator
 * \sa XMatrixN
 */


template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T,
                                                                           n>::operator*=(const kaMatrixN<T, n> &rm) {
  unsigned int i, j, k;

  kaMatrixN<T, n> mHelp(*this);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      T t = 0.0;
      for (k = 0; k < n; k++)
        t += mHelp.Get(i, k)*rm.Get(k, j);
      a(i, j) = t;
    }

  return *this;
}

//! kaMatrixN += kaMatrixN.
/*!
 * Add given kaMatrixN
 * overload the += operator
 */

template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T,
                                                                           n>::operator+=(const kaMatrixN<T, n> &rm) {
  T *ps = pElements, *pse = ps+n*n;
  const T *psm = rm.pElements;

  while (ps < pse)
    *ps++ += *psm++;

  return *this;
}

//! kaMatrixN -= kaMatrixN.
/*!
 * Subtract given kaMatrixN
 * overload the -= operator
 */

template<class T, const unsigned int n> inline kaMatrixN<T, n> & kaMatrixN<T,
                                                                           n>::operator-=(const kaMatrixN<T, n> &rm) {
  T *ps = pElements, *pse = ps+n*n;
  const T *psm = rm.pElements;

  while (ps < pse)
    *ps++ -= *psm++;

  return *this;
}

//! kaMatrixN multiplication from the RIGHT side.
/*!
 * Multiply the given kaMatrixN rm from the RIGHT side
 * (this*rm)
 * \param rm kaMatrixN to multiply
 * \sa MatrixNX *=
 */

template<class T, const unsigned int n> inline void kaMatrixN<T, n>::XMatrixN(const kaMatrixN<T, n> &rm) {
  kaMatrixN<T, n> mHelp(*this);

  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++) {
      T t = 0.0;
      for (unsigned int k = 0; k < n; k++)
        t += mHelp.Get(i, k)*rm.Get(k, j);
      a(i, j) = t;
    }
}

//! kaMatrixN multiplication from the LEFT side.
/*!
 * Multiply the given kaMatrixN rm from the LEFT side
 * (this*rm)
 * \param rm kaMatrixN to multiply
 * \sa MatrixNX *=
 */

template<class T, const unsigned int n> inline void kaMatrixN<T, n>::MatrixNX(const kaMatrixN<T, n> &rm) {
  kaMatrixN<T, n> mHelp(*this);

  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++) {
      T t = 0.0;
      for (unsigned int k = 0; k < n; k++)
        t += rm.Get(i, k)*mHelp.Get(k, j);
      a(i, j) = t;
    }
}

//! Invert kaMatrixN.
/*!
 * Invert an kaMatrixN. Return true if inversion was possible
 * else false
 * \sa Transpose
 * \return true if inversion was possible, else false
 */

template<class T, const unsigned int n>
int kaMatrixN<T, n>::Inverse() {
  unsigned int i, j, k, m;
  T t;

  kaMatrixN<T, n> InverseMatrix;
  InverseMatrix.Identity();


  for (i = 0; i < n; i++) {
    m = i;
    for (j = i+1; j < n; j++)
      if (fabs(Get(j, i)) > fabs(Get(m, i)))
        m = j;

    for (k = i; k < n; k++) {
      t       = Get(i, k);
      a(i, k) = Get(m, k);
      a(m, k) = t;
    }

    for (k = 0; k < n; k++) {
      t                     = InverseMatrix.Get(i, k);
      InverseMatrix.a(i, k) = InverseMatrix.Get(m, k);
      InverseMatrix.a(m, k) = t;
    }

    for (j = i+1; j < n; j++) {
      if (Get(i, i) == 0.0)
        return 0;

      t = Get(j, i)/Get(i, i);

      for (k = i; k < n; k++)
        a(j, k) -= Get(i, k)*t;

      for (k = 0; k < n; k++)
        InverseMatrix.a(j, k) -= InverseMatrix.Get(i, k)*t;
    }
  }

  i = n;
  do {
    i--;
    if (Get(i, i) == 0.0)
      return 0;

    for (k = 0; k < n; k++)
      InverseMatrix.a(i, k) /= Get(i, i);
    a(i, i) = 1.0;

    for (j = 0; j < i; j++) {
      t       = Get(j, i);
      a(j, i) = 0.0;

      for (k = 0; k < n; k++)
        InverseMatrix.a(j, k) -= InverseMatrix.Get(i, k)*t;
    }
  } while (i);

  *this = InverseMatrix;
  return 1;
}  // >::Inverse

//! Invert kaMatrixN<double,3>.
/*!
 * Invert a kaMatrixN<double,3>. Return true if inversion was possible
 * else false
 * \sa Transpose
 * \return true if inversion was possible, else false
 */

// *
template<class T, const unsigned int n> T kaMatrixN<T, n>::Determinant() {
  return kaMatrixMathN<T>::Determinant();
}

template<>
inline double kaMatrixN<double, 4>::Determinant() {
  // Do the calculations here!
  // according to the first row + - + -
  // the first element =
  double det = 0;

  det += kaMatrix<double>::Get(0)*
    (kaMatrix<double>::Get(5)*kaMatrix<double>::Get(10)*kaMatrix<double>::Get(15) + kaMatrix<double>::Get(6)*
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(13) + kaMatrix<double>::Get(7)*kaMatrix<double>::Get(9)*
     kaMatrix<double>::Get(14)
     - kaMatrix<double>::Get(7)*kaMatrix<double>::Get(10)*kaMatrix<double>::Get(13) -
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(14)*kaMatrix<double>::Get(5) - kaMatrix<double>::Get(15)*
     kaMatrix<double>::Get(6)*kaMatrix<double>::Get(9) );

  det -= kaMatrix<double>::Get(1)*
    (kaMatrix<double>::Get(4)*kaMatrix<double>::Get(10)*kaMatrix<double>::Get(15) + kaMatrix<double>::Get(6)*
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(12) + kaMatrix<double>::Get(7)*kaMatrix<double>::Get(8)*
     kaMatrix<double>::Get(14)
     - kaMatrix<double>::Get(7)*kaMatrix<double>::Get(10)*kaMatrix<double>::Get(12) -
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(14)*kaMatrix<double>::Get(4) - kaMatrix<double>::Get(15)*
     kaMatrix<double>::Get(6)*kaMatrix<double>::Get(8) );

  det += kaMatrix<double>::Get(2)*
    (kaMatrix<double>::Get(4)*kaMatrix<double>::Get(9)*kaMatrix<double>::Get(15) + kaMatrix<double>::Get(5)*
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(12) + kaMatrix<double>::Get(7)*kaMatrix<double>::Get(8)*
     kaMatrix<double>::Get(13)
     - kaMatrix<double>::Get(7)*kaMatrix<double>::Get(9)*kaMatrix<double>::Get(12) -
     kaMatrix<double>::Get(11)*kaMatrix<double>::Get(13)*kaMatrix<double>::Get(4) - kaMatrix<double>::Get(15)*
     kaMatrix<double>::Get(5)*kaMatrix<double>::Get(8) );

  det -= kaMatrix<double>::Get(3)*
    (kaMatrix<double>::Get(4)*kaMatrix<double>::Get(9)*kaMatrix<double>::Get(14) + kaMatrix<double>::Get(5)*
     kaMatrix<double>::Get(10)*kaMatrix<double>::Get(12) + kaMatrix<double>::Get(6)*kaMatrix<double>::Get(8)*
     kaMatrix<double>::Get(13)
     - kaMatrix<double>::Get(6)*kaMatrix<double>::Get(9)*kaMatrix<double>::Get(12) -
     kaMatrix<double>::Get(10)*kaMatrix<double>::Get(13)*kaMatrix<double>::Get(4) - kaMatrix<double>::Get(14)*
     kaMatrix<double>::Get(5)*kaMatrix<double>::Get(8) );
  return det;
} // >::Determinant

// */

template<>
inline int kaMatrixN<double, 3>::Inverse() {
  double d = Determinant();

  if (!d)
    return 0;

  d = 1/d;

  kaMatrixN<double, 3> M;
  M.a(0) = -Get(5)*Get(7) + Get(4)*Get(8);
  M.a(1) = Get(2)*Get(7) - Get(1)*Get(8);
  M.a(2) = -Get(2)*Get(4) + Get(1)*Get(5);
  M.a(3) = Get(5)*Get(6) - Get(3)*Get(8);
  M.a(4) = -Get(2)*Get(6) + Get(0)*Get(8);
  M.a(5) = Get(2)*Get(3) - Get(0)*Get(5);
  M.a(6) = -Get(4)*Get(6) + Get(3)*Get(7);
  M.a(7) = Get(1)*Get(6) - Get(0)*Get(7);
  M.a(8) = -Get(1)*Get(3) + Get(0)*Get(4);

  *this = M*d;
  return 1;
}

//! Translate kaMatrixN (right side).
/*!
 * The vector (x,y,z) is used to translate the kaMatrixN (from the RIGHT side).
 * \param x translate in x direction
 * \param y translate in y direction
 * \param z translate in z direction
 */

template<class T, const unsigned int n> void kaMatrixN<T, n>::XTranslate(T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp.Translate(x, y, z);
  *this *= temp;
}

//! Scale kaMatrixN (right side).
/*!
 * The vector (x,y,z) is used to scale the kaMatrixN (from the RIGHT side).
 * \param x scale in x direction
 * \param y scale in y direction
 * \param z scale in z direction
 */

template<class T, const unsigned int n> void kaMatrixN<T, n>::XScale(T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp.Scale(x, y, z);
  *this *= temp;
}

//! Rotate kaMatrixN (right side).

/*!
 * The kaMatrixN is rotated via the rotation axis (x,y,z) and the angel a.
 * \param a angel to rotate
 * \param x rotation axis component
 * \param y rotation axis component
 * \param z rotation axis component
 */
template<class T, const unsigned int n> void kaMatrixN<T, n>::XRotate(T a, T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp.Rotate(a, x, y, z);
  *this *= temp;
}

template<class T, const unsigned int n> void kaMatrixN<T, n>::XFrustum(T left, T right, T bottom, T top, T near_arg,
                                                                       T far_arg) {
  kaMatrixN<T, 4> temp;
  temp.Frustum(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T, const unsigned int n> void kaMatrixN<T, n>::XOrtho(T left, T right, T bottom, T top, T near_arg,
                                                                     T far_arg) {
  kaMatrixN<T, 4> temp;
  temp.Ortho(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

//! Translate kaMatrixN (left side).
/*!
 * The vector (x,y,z) is used to translate the kaMatrixN (from the LEFT side).
 * \param x translate in x direction
 * \param y translate in y direction
 * \param z translate in z direction
 */

template<class T, const unsigned int n> void kaMatrixN<T, n>::TranslateX(T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp = *this;
  Translate(x, y, z);
  *this *= temp;
}

//! Scale kaMatrixN (left side).

/*!
 * The vector (x,y,z) is used to scale the kaMatrixN (from the LEFT side).
 * \param x scale in x direction
 * \param y scale in y direction
 * \param z scale in z direction
 */
template<class T, const unsigned int n> void kaMatrixN<T, n>::ScaleX(T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp = *this;
  Scale(x, y, z);
  *this *= temp;
}

//! Rotate matrix (left side).
/*!
 * The kaMatrixN is rotated via the rotation axis (x,y,z) and the angel a.
 * \param a angel to rotate
 * \param x rotation axis component
 * \param y rotation axis component
 * \param z rotation axis component
 */

template<class T, const unsigned int n> void kaMatrixN<T, n>::RotateX(T a, T x, T y, T z) {
  kaMatrixN<T, 4> temp;
  temp = *this;
  Rotate(a, x, y, z);
  *this *= temp;
}

template<class T, const unsigned int n> void kaMatrixN<T, n>::FrustumX(T left, T right, T bottom, T top, T near_arg,
                                                                       T far_arg) {
  kaMatrixN<T, 4> temp;
  temp = *this;
  Frustum(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

template<class T, const unsigned int n> void kaMatrixN<T, n>::OrthoX(T left, T right, T bottom, T top, T near_arg,
                                                                     T far_arg) {
  kaMatrixN<T, 4> temp;
  temp = *this;
  Ortho(left, right, bottom, top, near_arg, far_arg);
  *this *= temp;
}

//! Rotate kaMatrixN.
/*!
 * The kaMatrixN is rotated via the rotation axis (x,y,z) and the angel phi.
 * \param phi angle to rotate
 * \param x rotation axis component
 * \param y rotation axis component
 * \param z rotation axis component
 */

template<class T, const unsigned int n>
void kaMatrixN<T, n>::Rotate(T phi, T x, T y, T z) {
  assert(n == 4);

  T iBetrag = 1.0/sqrt(x*x + y*y + z*z);
  x *= iBetrag;
  y *= iBetrag;
  z *= iBetrag;

  kaMatrixN<T, 3> uut;
  uut.a(0, 0) = x*x; uut.a(0, 1) = x*y; uut.a(0, 2) = x*z;
  uut.a(1, 0) = y*x; uut.a(1, 1) = y*y; uut.a(1, 2) = y*z;
  uut.a(2, 0) = z*x; uut.a(2, 1) = z*y; uut.a(2, 2) = z*z;

  kaMatrixN<T, 3> s;
  s.a(0, 0) = 0;  s.a(0, 1) = -z; s.a(0, 2) = y;
  s.a(1, 0) = z;  s.a(1, 1) = 0;  s.a(1, 2) = -x;
  s.a(2, 0) = -y; s.a(2, 1) = x;  s.a(2, 2) = 0;
  s        *= sin(phi);

  kaMatrixN<T, 3> I;
  I.Identity();
  I -= uut;
  I *= cos(phi);

  uut += s;
  uut += I;

  kaMatrixMathN<T>::Null();
  a(0, 0) = uut.Get(0, 0); a(0, 1) = uut.Get(0, 1); a(0, 2) = uut.Get(0, 2);
  a(1, 0) = uut.Get(1, 0); a(1, 1) = uut.Get(1, 1); a(1, 2) = uut.Get(1, 2);
  a(2, 0) = uut.Get(2, 0); a(2, 1) = uut.Get(2, 1); a(2, 2) = uut.Get(2, 2);
}  // >::Rotate

//! Trace
/*!
 * The trace of the matrix is calculated and returned
 */

template<class T, const unsigned int n>
const T kaMatrixN<T, n>::Trace() {
  T *ps = pElements, *pse = ps+n*n;
  T  result = T(0.0);

  while (ps < pse) {
    result += *ps;
    ps     += (n+1);
  }
  return result;
}

#endif  // ifndef KAMATRIXN_H
