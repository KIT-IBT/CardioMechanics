/*! \file kaSharedMatrixMxN.h
   \brief Class for handling of MxN matrices in shared memory

   \author cw,fs,os, IBT - Universität Karlsruhe (TH)
 */

#ifndef KASHAREDMATRIXMXN_H
#define KASHAREDMATRIXMXN_H

#include "kaMatrix.h"

template<class T> class MatrixSymBandMxN;


// .NAME kaSharedMatrixMxN - one more matrix class
// .SECTION Description
// Matrix with dinamicaly allocated memory.

#define __MATRIXMxNFILEID  "kaSharedMatrixNxNFile"
template<class T>
class kaSharedMatrixMxN : public kaMatrix<T> {
 protected:
  unsigned int *pm;  //!< Pointer to number of rows
  unsigned int *pn;  //!< Pointer to number of columns
  T *pElements;  //!< Pointer to elements

 public:
  //! Default constructor
  inline kaSharedMatrixMxN();

  //! Constructor with rows und colums
  inline kaSharedMatrixMxN(unsigned int, unsigned int);

  //! Constructor with kaMatrix
  inline kaSharedMatrixMxN(const kaMatrix<T> &);

  //! Constructor with kaSharedMatrixMxN
  inline kaSharedMatrixMxN(const kaSharedMatrixMxN<T> &);

  //! Destructor
  inline ~kaSharedMatrixMxN();

  //! Pointer to elements
  inline T *Pointer() {return pElements;}

  //! Const pointer to elements
  inline const T *Pointer() const {return pElements;}

  // (by ouss)
  //! Const pointer to first element of row c
  inline const T *PointerToRow(const unsigned int c) const {return pElements+*pm*c;}

  //! Number of rows
  inline unsigned int NumRows() const {return *pm;}

  //! Number of columns
  inline unsigned int NumCols() const {return *pn;}

  virtual inline T & a(unsigned int z, unsigned int s) {return (T &)*(pElements+z **pn+s);}

  virtual inline T & a(unsigned int p) {return (T &)*(pElements+p);}

  virtual inline T *operator[](unsigned int i) {return pElements+i **pn;}

  virtual inline const T *operator[](unsigned int i) const {return pElements+i **pn;}

  virtual inline const T Get(unsigned int z, unsigned int s) const {return *(pElements+z **pn+s);}

  virtual inline const T Get(unsigned int p) const {return *(pElements+p);}

  inline kaSharedMatrixMxN & operator=(T);
  inline kaSharedMatrixMxN & operator=(const kaSharedMatrixMxN<T> &rm);

  inline void operator+=(const kaSharedMatrixMxN<T> &rm);
  inline void operator-=(const kaSharedMatrixMxN<T> &rm);

  // Description:
  // Returns pointer to transpose of this matrix
  // kaSharedMatrixMxN<T> * Transpose();
  inline kaSharedMatrixMxN<T> Transpose();

  // Description:
  // Read from binary file. Matrix must be empty (created with default constructor).
  int ReadBinary(const char *);
  int SaveBinary(const char *) const;  // write to binary file
  // Description:
  // Returns the Sum of All elements of the Matrix.
  inline T ElementsSum();
  void     New(const unsigned int m, const unsigned int n);
  void     Delete();
};  // class kaSharedMatrixMxN


// .SECTION Operators
// Matrix multiplication
template<class T>
const inline kaSharedMatrixMxN<T> operator*(const kaMatrix<T> &a, const kaMatrix<T> &b) {
  kaSharedMatrixMxN<T> temp(a.NumRows(), b.NumCols());
  int i, j, k;

  const int m = temp.NumRows();
  const int n = temp.NumCols();
  const int p = a.NumCols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      T val = 0;
      for (k = 0; k < p; k++)
        val += a.Get(i, k)*b.Get(k, j);
      temp.a(i, j) = val;
    }

  return temp;
}

// Plus matrix. Matrices must be of the same size.

template<class T>
void kaSharedMatrixMxN<T>::operator+=(const kaSharedMatrixMxN<T> &rm) {
  T *ps        = Pointer();
  T *pe        = Pointer() + this->Size();
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ += *pss++;
}

// Minus matrix. Matrices must be of the same size.
template<class T>
void kaSharedMatrixMxN<T>::operator-=(const kaSharedMatrixMxN<T> &rm) {
  T *ps        = Pointer();
  T *pe        = Pointer() + this->Size();
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ -= *pss++;
}

/*-------------------------------------------------------------------------------*/

template<class T>
inline kaSharedMatrixMxN<T> kaSharedMatrixMxN<T>::Transpose() {
  kaSharedMatrixMxN<T> temp(*pn, *pm);
  unsigned long int i, j;
  for (i = 0; i < *pn; i++) {
    for (j = i+1; j < *pn; j++) {
      temp.a(i, j) = kaMatrix<T>::Get(j, i);
      temp.a(j, i) = kaMatrix<T>::Get(i, j);
    }
  }
  return temp;
}

/*-------------------------------------------------------------------------------*/

// kaSharedMatrixMxN implementation


template<class T>
inline kaSharedMatrixMxN<T>::kaSharedMatrixMxN() {
  pm        = NULL;
  pn        = NULL;
  pElements = NULL;
}

template<class T>
inline kaSharedMatrixMxN<T>::kaSharedMatrixMxN(unsigned int m, unsigned int n) {
  pm        = NULL;
  pn        = NULL;
  pElements = NULL;
  New(m, n);
}

template<class T>
inline kaSharedMatrixMxN<T>::kaSharedMatrixMxN(const kaMatrix<T> &rm) {
  pm        = NULL;
  pn        = NULL;
  pElements = NULL;
  New(rm.NumRows(), rm.NumCols());
  *this = rm;
}

template<class T>
inline kaSharedMatrixMxN<T>::kaSharedMatrixMxN(const kaSharedMatrixMxN<T> &rm) {
  pm        = NULL;
  pn        = NULL;
  pElements = NULL;
  New(rm.NumRows(), rm.NumCols());
  *this = rm;
}

template<class T>
inline kaSharedMatrixMxN<T>::~kaSharedMatrixMxN() {
  Delete();
}

template<class T>
inline void kaSharedMatrixMxN<T>::Delete() {
  if (pm)
    delete pm;
  pm = NULL;
  if (pn)
    delete pn;
  pn = NULL;
  if (pElements)
    delete[] pElements;
  pElements = NULL;
}

template<class T>
inline void kaSharedMatrixMxN<T>::New(const unsigned int m, const unsigned int n) {
  Delete();
  pm        = new unsigned int;
  pn        = new unsigned int;
  *pm       = m;
  *pn       = n;
  pElements = new T[m*n];
}

template<class T>
inline kaSharedMatrixMxN<T> & kaSharedMatrixMxN<T>::operator=(const T a) {
  T *ps = Pointer(), *pe = Pointer() + NumRows()*NumCols();

  while (ps < pe)
    *ps++ = a;

  return *this;
}

template<class T>
inline kaSharedMatrixMxN<T> & kaSharedMatrixMxN<T>::operator=(const kaSharedMatrixMxN<T> &rm) {
  T *ps = Pointer(), *pe = Pointer() + NumRows()*NumCols();
  const T *pss = rm.Pointer();

  while (ps < pe)
    *ps++ = *pss++;

  return *this;
}

template<class T> int kaSharedMatrixMxN<T>::SaveBinary(const char *name) const {
  FILE *to = fopen(name, "w");

  if (!to) {
    fprintf(stderr, "kaSharedMatrixMxN<T>::Save: cannot open %s \n", name);
    return 0;
  }
  char d[] = __MATRIXMxNFILEID;
  fwrite(d, sizeof(__MATRIXMxNFILEID), 1, to);
  fwrite(pm, sizeof(unsigned int), 1, to);
  fwrite(pn, sizeof(unsigned int), 1, to);
  fwrite(pElements, sizeof(T), NumRows()*NumCols(), to);
  fclose(to);
  return 1;
}

template<class T> int kaSharedMatrixMxN<T>::ReadBinary(const char *name) {
  if (pElements)

    if (!((pm == NULL) && (pn == NULL) && (pElements == NULL))) {
      fprintf(stderr, "kaSharedMatrixMxN<T>::Read: cannot read into already initialized matrix\n");
      return 0;
    }

  FILE *from = fopen(name, "r");
  if (!from) {
    fprintf(stderr, "kaSharedMatrixMxN<T>::Read: cannot open %s \n", name);
    return 0;
  }
  char d[] = __MATRIXMxNFILEID;  // bug found 13.07.2008
  fread(d, sizeof(__MATRIXMxNFILEID), 1, from);
  if (strcmp(d, __MATRIXMxNFILEID) != 0) {
    fprintf(stderr, "kaSharedMatrixMxN<T>::Read: %s  is not a binary matrix file \n", name);
    return 0;
  }

  unsigned int m, n;
  fread(&m, sizeof(unsigned int), 1, from);
  fread(&n, sizeof(unsigned int), 1, from);
  New(m, n);
  fread(pElements, sizeof(T), NumRows()*NumCols(), from);
  return 1;
}  // >::ReadBinary

template<class T, class V>
inline const kaSharedMatrixMxN<V> operator*(const kaMatrix<T> &a, const kaMatrix<V *> &b) {
  kaSharedMatrixMxN<V> temp(a.NumRows(), b.NumCols());
  unsigned long int i, j, k;
  const int m = temp.NumRows();
  const int n = temp.NumCols();

  const int p = a.NumCols();

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      V val = V(0);
      for (k = 0; k < p; k++) {
        val += (operator*(a.Get(i, k), *b.Get(k, j)));
      }
      temp.a(i, j) = val;
    }
  }
  return temp;
}

/*
   //Another Matrix multiplication
   template<class T,class V>
   inline const kaSharedMatrixMxN<V> operator*( const kaMatrix<T>& a, const kaMatrix<V>& b)
   {
        const int m = a.NumRows();
        const int p = a.NumCols();
        const int n = b.NumCols();

        kaSharedMatrixMxN<V> temp(m, n);
        unsigned long int i, j, k;

        for(i=0; i<m; i++)
        {
                for(j=0; j<n; j++)
                {
                        V val=V(0);
                        for(k=0; k<p; k++)
                        {
                                val += (operator*(a.Get(i,k),b.Get(k,j)));
                        }
                        temp.a(i,j)=val;
                }
        }
        return temp;
   }
 */

template<class T, class V>
inline const kaSharedMatrixMxN<V> MutiplyMatrixWithOneColumn(const kaMatrix<T> &a, const kaMatrix<V> &b,
                                                             const unsigned int j) {
  const int m = a.NumRows();
  const int p = a.NumCols();

  kaSharedMatrixMxN<V> temp(m, 1);
  unsigned long int i, k;

  for (i = 0; i < m; i++) {
    V val = V(0);
    for (k = 0; k < p; k++) {
      val += (operator*(a.Get(i, k), b.Get(k, j)));
    }
    temp.a(i) = val;
  }
  return temp;
}

template<class T>
inline const kaSharedMatrixMxN<T> ElementsNorm(kaSharedMatrixMxN<T> a) {
  const int m = a.NumRows();
  const int n = a.NumCols();

  kaSharedMatrixMxN<T> temp(m, n);
  unsigned long int i, j;

  const T *ps = a.Pointer();
  T *pss      = temp.Pointer();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      *pss = norm(*ps);
      pss++; ps++;
    }
  return temp;
}

template<class T, class V>
inline const kaSharedMatrixMxN<T> ElementsAbs(kaSharedMatrixMxN<V> a) {
  const int m = a.NumRows();
  const int n = a.NumCols();

  kaSharedMatrixMxN<T> temp(m, n);
  unsigned long int i, j;

  const V *ps = a.Pointer();
  T *pss      = temp.Pointer();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      *pss = abs(*ps);
      pss++; ps++;
    }
  return temp;
}

template<class T>
inline T kaSharedMatrixMxN<T>::ElementsSum() {
  T sum = 0.0;
  T *p  = Pointer();

  for (int i = 0; i < kaMatrix<T>::Size(); i++, p++)
    sum += *p;
  return sum;
}

#endif  // ifndef KASHAREDMATRIXMXN_H
