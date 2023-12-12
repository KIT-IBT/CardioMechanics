/*
 * File: kaMatrixMxN.h
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


#ifndef M_MATRIXMXN_H
#define M_MATRIXMXN_H


#include "kaMatrix.h"

template<class T> class MatrixFastSymBandMxN;

template<class T, const int m, const int n>
class kaMatrixMxN : public kaMatrix<T> {
 protected:
  T pElements[m*n];

 public:
  inline kaMatrixMxN() {}

  inline kaMatrixMxN(const kaMatrix<T> &);

  inline ~kaMatrixMxN() {}

  inline T *Pointer() {return pElements;}  // the first element

  inline const T *Pointer() const {return pElements;}

  inline unsigned int NumRows() const {return m;}

  inline unsigned int NumCols() const {return n;}

  inline kaMatrixMxN<T, m, n> & operator=(const kaMatrixMxN<T, m, n> &);

  inline kaMatrixMxN<T, m, n> & operator*=(const T);
  inline kaMatrixMxN<T, m, n> & operator/=(const T);

  inline kaMatrixMxN<T, m, n> & operator+=(const kaMatrixMxN<T, m, n> &);
  inline kaMatrixMxN<T, m, n> & operator-=(const kaMatrixMxN<T, m, n> &);

  friend class MatrixFastSymBandMxN<T>;
};  // class kaMatrixMxN


template<class T, const int m, const int n> inline kaMatrixMxN<T, m, n> & kaMatrixMxN<T, m, n>::operator=(
  const kaMatrixMxN<T, m,
                    n> &rm)
{
  T *ps        = Pointer();
  T *pe        = Pointer() + NumRows()*NumCols();
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ = *pss++;
  return *this;
}

template<class T, const int m, const int n> inline kaMatrixMxN<T, m, n> & kaMatrixMxN<T, m, n>::operator*=(const T f) {
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());

  while (ps < pse)
    *ps++ *= f;

  return *this;
}

template<class T, const int m, const int n> inline kaMatrixMxN<T, m, n> & kaMatrixMxN<T, m, n>::operator/=(const T f) {
  *this *= 1/f;
}

template<class T, const int m, const int n> inline kaMatrixMxN<T, m, n> & kaMatrixMxN<T, m, n>::operator+=(
  const kaMatrixMxN<T, m,
                    n> &rm)
{
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());
  T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ += *psm++;
}

template<class T, const int m, const int n> inline kaMatrixMxN<T, m, n> & kaMatrixMxN<T, m, n>::operator-=(
  const kaMatrixMxN<T, m,
                    n> &rm)
{
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());
  T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ -= *psm++;
}

// .SECTION Operators
// Matrix multiplication
template<class T, const int m, const int p, const int n>
const inline kaMatrixMxN<T, m, n> operator*(const kaMatrixMxN<T, m, p> &a, const kaMatrixMxN<T, p, n> &b) {
  kaMatrixMxN<T, m, n> temp;
  int i, j, k;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      T val = 0;
      for (k = 0; k < p; k++)
        val += a.Get(i, k)*b.Get(k, j);
      temp.a(i, j) = val;
    }

  return temp;
}

/*-------------------------------------------------------------------------------*/

// kaMatrixMxN implementation


template<class T, const int m, const int n>
inline kaMatrixMxN<T, m, n>::kaMatrixMxN(const kaMatrix<T> &rm) {
  T *ps = pElements, *pse = (ps+n*n);
  const T *pss = rm.Pointer();

  while (ps < pse)
    *ps++ = *pss++;
}

#endif  // ifndef M_MATRIXMXN_H
