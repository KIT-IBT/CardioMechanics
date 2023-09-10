/*
   MatrixSymBandMxN.h

   derived form kaMatrix

   os Feb 28, 2001
   fs Jun 29, 2001

 */

#ifndef M_MATRIXSYMBANDMXN_H
#define M_MATRIXSYMBANDMXN_H


#include "kaSharedMatrixMxN.h"


template<class T>
class MatrixSymBandMxN : public kaSharedMatrixMxN<T> {
 protected:
  unsigned int *pIndex;

 public:
  inline MatrixSymBandMxN();
  inline MatrixSymBandMxN(const MatrixSymBandMxN<T> &);

  inline virtual ~MatrixSymBandMxN();

  inline unsigned int *Index() {return pIndex;}

  inline const unsigned int *ConstIndex() const {return pIndex;}

  inline unsigned int NumRows() const {return *pm;}

  inline unsigned int NumCols() const {return *pn;}

  MatrixSymBandMxN & operator=(T);
  MatrixSymBandMxN & operator=(const MatrixSymBandMxN &);

  MatrixSymBandMxN & operator*=(const T);
  MatrixSymBandMxN & operator/=(const T);

  MatrixSymBandMxN & operator+=(const MatrixSymBandMxN &);
  MatrixSymBandMxN & operator-=(const MatrixSymBandMxN &);
  void               VectorMult(const T *, T *) const;
  void               AddColumn(unsigned int i, T factor, T *b) const;
  void               SetColumn(unsigned int i, T val);
  void               SetRow(unsigned int i, T val);
  virtual void       Save(FILE *) const;
  virtual void       Restore(FILE *);
  void               New(unsigned int, unsigned int);
  void               Delete();
};  // class MatrixSymBandMxN


/*-------------------------------------------------------------------------------*/

// MatrixSymBandMxN implementation


template<class T>
inline MatrixSymBandMxN<T>::MatrixSymBandMxN() {
  pm        = pn = NULL;
  pIndex    = NULL;
  pElements = NULL;
}

template<class T>
inline MatrixSymBandMxN<T>::MatrixSymBandMxN(const MatrixSymBandMxN<T> &rm) {
  pm        = pn = NULL;
  pIndex    = NULL;
  pElements = NULL;
  New(rm.NumRows(), rm.NumCols());

  int *pi = pIndex, *pie = pi+n;
  const int *pis = rm.ConstIndex();
  while (pi < pie)
    *pi++ = *pis++;

  T *ps = pElements, *pse = ps+m*n;
  const T *pss = rm.Pointer();
  while (ps < pse)
    *ps++ = *pss++;
}

template<class T>
inline MatrixSymBandMxN<T>::~MatrixSymBandMxN() {
  delete pm;
  delete pn;
  delete pIndex;
  delete[] pElements;

  // printf("MatrixSymBandMxN<T>::~MatrixSymBandMxN() %ld %ld %ld %ld\n", pm, pn, pIndex, pElements);
}

template<class T>
inline void MatrixSymBandMxN<T>::Delete() {
  delete pm;
  pm = NULL;
  delete pn;
  pn = NULL;
  delete pIndex;
  pIndex = NULL;
  delete[] pElements;
  pElements = NULL;
}

template<class T>
inline void MatrixSymBandMxN<T>::New(unsigned int m, unsigned int n) {
  Delete();
  pm        = new unsigned int;
  pn        = new unsigned int;
  *pm       = m;
  *pn       = n;
  pIndex    = new unsigned int[n];
  pElements = new T[m*n];

  // printf("MatrixSymBandMxN<T>::New %ld %ld %ld %ld\n", pm, pn, pIndex, pElements);
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator=(T a) {
  T *ps = Pointer(), *pe = Pointer()+*pm **pn;

  while (ps < pe)
    *ps++ = a;
  return *this;
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator=(const MatrixSymBandMxN<T> &rm) {
  T *ps = Pointer(), *pe = Pointer()+*pm **pn;
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ = *pss++;
  return *this;
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator*=(const T f) {
  T *ps = Pointer(), *pse = Pointer()+*pm **pn;

  while (ps < pse)
    *ps++ *= f;

  return *this;
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator/=(const T f) {
  *this *= 1/f;
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator+=(const MatrixSymBandMxN<T> &rm) {
  T *ps = Pointer(), *pse = Pointer()+*pm **pn;
  const T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ += *psm++;
}

template<class T> inline MatrixSymBandMxN<T> & MatrixSymBandMxN<T>::operator-=(const MatrixSymBandMxN<T> &rm) {
  T *ps = Pointer(), *pse = Pointer()+*pm **pn;
  const T *psm = rm.Pointer();

  while (ps < pse)
    *ps++ -= *psm++;
}

template<class T>
inline void MatrixSymBandMxN<T>::VectorMult(const T *x, T *b) const {
  unsigned int i, j;
  const unsigned int m = *pm, n = *pn;
  const T *psmi0 = pElements;

  for (i = 0; i < m; i++, x++, b++) {
    *b = *psmi0++ *x[0];
    for (j = 1; j < n; j++, psmi0++) {
      unsigned int pIndexj = pIndex[j];

      if (i+pIndexj < m)
        *b += *psmi0*x[pIndexj];
      if (i >= pIndexj)
        *b += *(psmi0-n*pIndexj)*x[-pIndexj];
    }
  }
}

template<class T>
inline void MatrixSymBandMxN<T>::AddColumn(unsigned int i, T factor, T *b) const {
  unsigned int j = 0;

  b[i] += factor*Get(i, j);
  for (j = 1; j < NumCols(); j++) {
    if (i+pIndex[j] < NumRows())
      b[i+pIndex[j]] += factor*Get(i, j);
    if (i >= pIndex[j])
      b[i-pIndex[j]] += factor*Get(i-pIndex[j], j);
  }
}

template<class T>
inline void MatrixSymBandMxN<T>::SetColumn(unsigned int i, T val) {
  unsigned int j;

  for (j = 0; j < NumCols(); j++)
    if (i >= pIndex[j])
      a(i-pIndex[j], j) = val;
}

template<class T>
inline void MatrixSymBandMxN<T>::SetRow(unsigned int i, T val) {
  unsigned int j;

  for (j = 0; j < NumCols(); j++)
    a(i, j) = val;
}

template<class T> inline void MatrixSymBandMxN<T>::Save(FILE *hFile) const {
  fwrite(pm, sizeof(unsigned int), 1, hFile);
  fwrite(pn, sizeof(unsigned int), 1, hFile);
  fwrite(ConstIndex(), sizeof(T), NumCols(), hFile);
  fwrite(Pointer(), sizeof(T), Size(), hFile);
}

template<class T> inline void MatrixSymBandMxN<T>::Restore(FILE *hFile) {
  unsigned int m, n;

  fread(&m, sizeof(unsigned int), 1, hFile);
  fread(&n, sizeof(unsigned int), 1, hFile);
  New(m, n);
  fread(Index(), sizeof(T), NumCols(), hFile);
  fread(Pointer(), sizeof(T), Size(), hFile);
}

#endif  // ifndef M_MATRIXSYMBANDMXN_H
