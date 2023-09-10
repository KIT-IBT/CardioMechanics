/*
   MatrixFastSymBandMxN.h

   derived form kaMatrix

   os Feb 28, 2001
   fs Jun 29, 2001

 */

#ifndef M_MATRIXFASTSYMBANDMXN_H
#define M_MATRIXFASTSYMBANDMXN_H


#include "kaMatrix.h"

template<class T, const int m, const int n>
class MatrixFastSymBandMxN : public kaMatrixMxN<T> {
 protected:
  int pIndex[n];

 public:
  inline MatrixFastSymBandMxN(const MatrixFastSymBandMxN<T> &);

  inline int *Index() {return pIndex;}

  inline const int *ConstIndex() const {return pIndex;}

  MatrixFastSymBandMxN<T, m, n> & operator=(const MatrixFastSymBandMxN<T, m, n> &);
  MatrixFastSymBandMxN<T, m, n> & operator=(const T);

  MatrixFastSymBandMxN<T, m, n> & operator*=(const T);
  MatrixFastSymBandMxN<T, m, n> & operator/=(const T);

  MatrixFastSymBandMxN<T, m, n>        & operator+=(const MatrixFastSymBandMxN<T, m, n> &);
  MatrixFastSymBandMxN<T, m, n>        & operator-=(const MatrixFastSymBandMxN<T, m, n> &);
  inline MatrixFastSymBandMxN<T, m, n> & operator-=(const MatrixFastSymBandMxN<T, m, n> &);
  virtual void                           Save(FILE *);
  virtual void                           Restore(FILE *);
};

/*-------------------------------------------------------------------------------*/

// MatrixFastSymBandMxN implementation


template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator=(
  const MatrixFastSymBandMxN<T, m,
                             n> &rm)
{
  T *ps        = Pointer();
  T *pe        = Pointer() + NumRows()*NumCols();
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ = *pss++;
  return *this;
}

template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator=(
  const T a) {
  T *ps = Pointer();
  T *pe = Pointer() + NumRows()*NumCols();

  while (ps < pe)
    *ps++ = a;
  return *this;
}

template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator*=(
  const T f) {
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());

  while (ps < pse)
    *ps++ *= f;

  return *this;
}

template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator/=(
  const T f) {
  *this *= 1/f;
}

template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator+=(
  const MatrixFastSymBandMxN<T, m,
                             n> &rm)
{
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());
  T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ += *psm++;
}

template<class T, const int m, const int n> inline MatrixFastSymBandMxN<T, m, n> & MatrixFastSymBandMxN<T, m,
                                                                                                        n>::operator-=(
  const MatrixFastSymBandMxN<T, m,
                             n> &rm)
{
  T *ps = Pointer(), *pse = (Pointer()+NumRows()*NumCols());
  T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ -= *psm++;
}

template<class T, const int m, const int n>
inline MatrixFastSymBandMxN<T, m, n>::MatrixFastSymBandMxN(const MatrixFastSymBandMxN<T> &rm) {
  int *pi = pIndex, *pie = pi+n;
  const int *pis = rm.ConstIndex();

  while (pi < pie)
    *pi++ = *pis++;

  T *ps = pElements, *pse = ps+n*n;
  const T *pss = rm.Pointer();
  while (ps < pse)
    *ps++ = *pss++;
}

template<class T, const int m, const int n> inline void MatrixFastSymBandMxN<T, m, n>::Save(FILE *hFile) {
  fwrite(Index(), sizeof(int), NumCols(), hFile);  // elements
  fwrite(Pointer(), sizeof(T), Size(), hFile);  // elements
}

template<class T, const int m, const int n> inline void MatrixFastSymBandMxN<T, m, n>::Restore(FILE *hFile) {
  fread(Index(), sizeof(int), NumCols(), hFile);  // elements
  fread(Pointer(), sizeof(T), Size(), hFile);    // elements
}

#endif  // ifndef M_MATRIXFASTSYMBANDMXN_H
