/*
 * File: UPSymMatrix.h
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

#ifndef __UPSymMatrix_h
#define __UPSymMatrix_h

#include <kaMatrixBase.h>

// .SECTION Name UPSymMatrix
// .SECTION Description

//  The matrix is packed column by column in n(n+1)/2 elements of a
// one-dimensional array. To calculate the location of each element aij of matrix A in an array AP using the upper
// triangular
// packed technique, use the following formula:
//
// AP(i+((j+1)j)/2)) = aij    where j >= i
//
// AP[0] = a00 (start the first column)
// AP[1] = a01 (start the second column)
// AP[2] = a11
// AP[3] = a02 (start the third column)
// AP[4] = a12
// AP[5] = a22
// AP[6] = a03 (start the fourth column)
// .       .
// .       .
// and so on
//
// folowing is the example of how to transform your symmetric matrix
// to upper-packed storage mode:
//
//   int k = 0;
//   for(int j=0; j<N; j++){
//    for(int i=0; i<=j; i++){
//     AP[k]=A[i][j]
//     k++;
//    }
//   }
//
//
// WARNING: to access elements ALWAYS use j>=i
//
//
template<class T, int N> class UPSymMatrix : public kaMatrixBase<T> {
  T pElements[N*(N+1)/2];

 public:
  UPSymMatrix() {}

  ~UPSymMatrix() {}

  inline T *Pointer() {return pElements;}

  inline const T *Pointer() const {return pElements;}

  inline unsigned int NumCols() const {return N;}

  inline unsigned int NumRows() const {return N;}

  // Description:
  // R/W element access
  inline T & a(int p) {return (T &)*(pElements+p);}

  inline T & a(int i, int j) {return (T &)*(pElements+i+((j+1)*j)/2);}

  // Description:
  // Use this if there is no need to change the element.
  inline const T Get(int i, int j) const {return *(pElements+i+((j+1)*j)/2);}

  inline const T Get(int p) const {return *(pElements+p);}

  inline UPSymMatrix<T, N> & operator=(const UPSymMatrix<T, N> &);

  void Print(const char *format = "%7.3lf\t") const {
    for (int i = 0; i < N*(N+1)/2; i++) printf(format, Get(i)); printf("\n");
  }

  void Save(FILE *) const;
  void Restore(FILE *);

 private:
  T *operator[](int i) {return NULL;}

  const T *operator[](int i) const {return NULL;}
};  // class UPSymMatrix

// === Some Operators ===

// Description:
// Output operator. Writes matrix in text form
template<class T, int N> ostream & operator<<(ostream &os, const UPSymMatrix<T, N> &rm) {
  unsigned int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (j >= i)
        os<<rm.Get(i, j)<<' ';
      else
        os<<rm.Get(j, i)<<' ';
    }
    os<<'\n';
  }
  return os;
}

// === implementation

template<class T, int N>
UPSymMatrix<T, N> & UPSymMatrix<T, N>::operator=(const UPSymMatrix<T, N> &rm) {
  T *p        = this->Pointer();
  const T *rp = rm.Pointer();

  int i = 0;

  while ((i++) != N*(N+1)/2) {*(p++) = *(rp++);}
}

template<class T, int N> inline void  UPSymMatrix<T, N>::Save(FILE *hFile) const {
  fwrite(Pointer(), sizeof(T), N*(N+1)/2, hFile);  // elements
}

template<class T, int N> inline void  UPSymMatrix<T, N>::Restore(FILE *hFile) {
  fread(Pointer(), sizeof(T), N*(N+1)/2, hFile);  // elements
}

#endif  // ifndef __UPSymMatrix_h
