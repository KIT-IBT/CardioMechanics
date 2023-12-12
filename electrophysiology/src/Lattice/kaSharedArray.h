/*
 * File: kaSharedArray.h
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


#ifndef KASHAREDARRAY_H
#define KASHAREDARRAY_H

#include <kaMachineOS.h>
#include <kaBasicIO.h>

using namespace nskaGlobal;

//! Class for handling of array with template elements in shared memory
/*!
   The array's size and data is addressed via pointers.
   The class serves a basis for kaLattice.

 */

template<class T> class kaSharedArray {
  uint32_t *n;  //!< Pointer to size of array
  T *pArray;  //!< Pointer to data of array

  void Alloc(uint32_t);
  void Free();

 public:
  //! Default constructor
  inline kaSharedArray() {
#if KALATTICEDEBUG
    cerr << "kaSharedArray<X>::kaSharedArray" << endl;
#endif  // if KALATTICEDEBUG
    n      = NULL;
    pArray = NULL;
  }

  //! Constructor with number of elements
  inline kaSharedArray(uint32_t num) {
#if KALATTICEDEBUG
    cerr << "kaSharedArray<X>::kaSharedArray(uint32_t num)" << endl;
#endif  // if KALATTICEDEBUG
    Alloc(num);
  }

  //! Destructor
  inline ~kaSharedArray() {
#if KALATTICEDEBUG
    cerr << "kaSharedArray<X>::~kaSharedArray" << endl;
#endif  // if KALATTICEDEBUG
    Free();
  }

  //! Add matrix memory size to pointer
  void GetSize(uint32_t &, uint32_t) const;

  //! Attach matrix to memory block
  void Attach(char * &, uint32_t);

  //! Detach matrix from memory block
  void Detach();

  //! Get array element
  inline T & operator[](uint32_t i) {return (T &)(*(pArray+i));}

  //! Get const array element
  inline const T operator[](uint32_t i) const {return *(pArray+i);}

  //! Fill array with given value
  void FillWith(T);

  //! Write dimension and array to file
  void Save(FILE *);

  //! Read dimension and array from file
  void Restore(FILE *);
};  // class kaSharedArray

template<class T> void kaSharedArray<T>::Alloc(uint32_t num) {
  n      = new uint32_t;
  pArray = new T[num];
  *n     = num;
}

template<class T> void kaSharedArray<T>::Free() {
#if KALATTICEDEBUG
  cerr << "kaSharedArray<X>::Free()" << endl;
#endif  // if KALATTICEDEBUG
  delete[] pArray;
  delete n;
}

template<class T> void kaSharedArray<T>::GetSize(uint32_t &s, uint32_t a) const {
  s += sizeof(uint32_t);
  s += a*sizeof(T);
}

template<class T> void kaSharedArray<T>::Attach(char * &rp, uint32_t a) {
#if KALATTICEDEBUG
  cerr << "kaSharedArray<X>::Attach" << endl;
#endif  // if KALATTICEDEBUG
  n      = (uint32_t *)rp;
  *n     = a;
  rp    += sizeof(uint32_t);
  pArray = (T *)rp;
  rp    += *n*sizeof(T);
}

template<class T> void kaSharedArray<T>::Detach() {
  n      = (uint32_t *)NULL;
  pArray = (T *)NULL;
}

template<class T> void kaSharedArray<T>::FillWith(T t) {
  T *p = pArray, *pe = p+*n;

  while (p < pe)
    *p++ = t;
}

template<class T> void kaSharedArray<T>::Save(FILE *hFile) {
  kaWrite(n, 1, hFile);
  kaWrite(pArray, *n, hFile);
}

template<class T> void kaSharedArray<T>::Restore(FILE *hFile) {
  uint32_t num;

  kaRead(&num, 1, hFile);

  if (num != *n) {
    Free();
    Alloc(num);
  }
  kaRead(pArray, *n, hFile);
}

#endif  // ifndef KASHAREDARRAY_H
