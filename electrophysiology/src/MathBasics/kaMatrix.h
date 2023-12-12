/*
 * File: kaMatrix.h
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


#ifndef KAMATRIX_H
#define KAMATRIX_H

#include <kaMachineOS.h>

//! Base class for matrix classes.
/*!
   kaMatrix is used to define uniform interfaces to classes of matrices of different types.

 */

template<class T> class kaMatrix  {
 public:
  //! Default constructor
  inline kaMatrix() {}

  //! Destructor
  virtual ~kaMatrix() {}

  //! Virtual pointer to elements
  virtual T *Pointer() = 0;

  //! Virtual const pointer to elements
  virtual const T *Pointer() const = 0;

  //! Number of columns
  virtual unsigned int NumCols() const = 0;

  //! Number of rows
  virtual unsigned int NumRows() const = 0;

  //! Get matrix element for r/w access by coordinate
  virtual inline T & a(unsigned int z, unsigned int s) {return (T &)*(Pointer()+z*NumCols()+s);}

  //! Get matrix element for r/w access by number
  virtual inline T & a(unsigned int p) {return (T &)*(Pointer()+p);}

  //! Get matrix line for r/w access by number
  virtual inline T *operator[](unsigned int i) {return Pointer()+i*NumCols();}

  //! Get matrix line for read access by number
  virtual inline const T *operator[](unsigned int i) const {return Pointer()+i*NumCols();}

  //! Get matrix element for read access by coordinates
  virtual inline const T Get(unsigned int z, unsigned int s) const {return *(Pointer()+z*NumCols()+s);}

  //! Get matrix element for read access by number
  virtual inline const T Get(unsigned int p) const {return *(Pointer()+p);}

  //! Write matrix elements in flat array
  inline kaMatrix & PutIntoArray(T *);

  //! Read matrix elements from flat array
  inline kaMatrix & SetFromArray(const T *);

  //! Swap columns of matrix by exchanging element's values
  virtual void SwapColumns(const unsigned int k, const unsigned int l);

  // by ouss
  //! Copy column k to column l by copying element's values
  virtual void CopyColumn(const unsigned int k, const unsigned int l);
  virtual void CopyColumn(kaMatrix<T> *ext, const unsigned int k, const unsigned int l);

  //! Shift the columns n_col to the right, discard the right most n_col and fill the left most n_col with fill_val
  virtual void RightShift(const unsigned int n_col, const T fill_val);
  virtual void RightShift(const unsigned int n_col);
  virtual void RightShiftRow(const unsigned int i_row, const unsigned int n_col, T fill_val);
  virtual void RightShiftRow(const unsigned int i_row, const unsigned int n_col);
  virtual void LeftShift(const unsigned int n_col, const T fill_val);
  virtual void LeftShift(const unsigned int n_col);  // don't fill
  virtual void LeftShiftRow(const unsigned int i_row, unsigned long int n_col, T fill_val);
  virtual void LeftShiftRow(const unsigned int i_row, unsigned long int n_col);  // no fill
  virtual void RepeatFirstColumn();
  virtual void RepeatFirstValueInRow(const unsigned int i_row);

  //! Set matrix to 0
  inline void Null();

  // by ouss
  //! Set matrix to val
  inline void Fill(const T val);
  inline void FillRow(const unsigned int i_row, const T val);
  inline void FillColumn(const unsigned int j_column, const T val);

  //! Calculate Frobenius norm of matrix
  T FrobeniusNorm();

  //! Print matrix elements to stdout via printf with given format
  virtual void Print(const char *format = "%7.3lg\t") const;

  //! Write matrix elements in file
  virtual void Save(FILE *) const;

  //! Read matrix elements from file
  virtual void Restore(FILE *);

  //! Calculate number of elements
  inline unsigned int Size() const {return NumCols()*NumRows();}
};  // class kaMatrix


/*-----------------------------------------------------------------------------*/
/* kaMatrix implementation */

template<class T> inline kaMatrix<T> & kaMatrix<T>::PutIntoArray(T *d) {
  T *ps = Pointer(), *pe = Pointer()+Size();

  while (ps < pe)
    *d++ = *ps++;

  return *this;
}

template<class T> inline kaMatrix<T> & kaMatrix<T>::SetFromArray(const T *d) {
  T *ps = Pointer(), *pe = Pointer()+Size();

  while (ps < pe)
    *ps++ = *d++;

  return *this;
}

template<class T> void kaMatrix<T>::SwapColumns(const unsigned int k, const unsigned int l) {
  T *sx = Pointer()+k, *sxe = sx+NumRows()*NumCols();
  T *sy = Pointer()+l;
  T  temp;

  while (sx < sxe) {
    temp = *sx;
    *sx  = *sy;
    *sy  = temp;
    sx  += NumCols();
    sy  += NumCols();
  }
}

template<class T> inline void kaMatrix<T>::RightShiftRow(const unsigned int i_row, const unsigned int n_col,
                                                         T fill_val) {
  if (n_col == 0)
    return;

  if (i_row >= NumRows())
    return;

  unsigned int n = NumCols();
  if (n_col >= n) {
    T *sy = Pointer()+i_row*n, *sye = Pointer()+i_row*n + n;
    while (sy < sye)
      *sy++ = fill_val;
    return;
  }
  T *sy = Pointer()+i_row*n + n -1;
  T *sx = sy - n_col, *ssxe = sy-n;
  while (sx > ssxe)
    *sy-- = *sx--;
  while (sy > ssxe)
    *sy-- = fill_val;
}

template<class T> inline void kaMatrix<T>::RightShiftRow(const unsigned int i_row, const unsigned int n_col) {
  if (n_col == 0)
    return;

  if (i_row >= NumRows())
    return;

  unsigned int n = NumCols();
  if (n_col >= n)
    return;

  T *sy = Pointer()+i_row*n + n -1;
  T *sx = sy - n_col, *ssxe = sy-n;
  while (sx > ssxe)
    *sy-- = *sx--;
}

template<class T> inline void kaMatrix<T>::RightShift(const unsigned int n_col, T fill_val) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n) {
    Fill(fill_val);
    return;
  }
  T *sy = Pointer()+n-1, *sye = Pointer()+Size();
  while (sy < sye) {
    T *ssy = sy;
    T *sx = ssy - n_col, *ssxe = sy-n;
    while (sx > ssxe)
      *ssy-- = *sx--;
    while (ssy > ssxe)
      *ssy-- = fill_val;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::RightShift(const unsigned int n_col) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n)
    return;

  T *sy = Pointer()+n-1, *sye = Pointer()+Size();
  while (sy < sye) {
    T *ssy = sy;
    T *sx = sy - n_col, *ssxe = sy-n;
    while (sx > ssxe)
      *ssy-- = *sx--;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::LeftShift(const unsigned int n_col, const T fill_val) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n) {
    Fill(fill_val);
    return;
  }
  T *sy = Pointer(), *sye = Pointer()+Size();
  while (sy < sye) {
    T *ssy = sy;
    T *sx = ssy + n_col, *ssxe = sy + n;
    while (sx < ssxe)
      *ssy++ = *sx++;
    while (ssy < ssxe)
      *ssy++ = fill_val;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::LeftShift(const unsigned int n_col) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n)
    return;

  T *sy = Pointer(), *sye = Pointer()+Size();
  while (sy < sye) {
    T *ssy = sy;
    T *sx = ssy + n_col, *ssxe = sy + n;
    while (sx < ssxe)
      *ssy++ = *sx++;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::LeftShiftRow(const unsigned int i_row, unsigned long int n_col,
                                                        T fill_val) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n) {
    // fill
    T *sy = Pointer()+i_row*n, *sye = sy+n;
    while (sy < sye)
      *sy++ = fill_val;
    return;
  }
  T *sy = Pointer()+i_row*n;
  T *sx = sy + n_col, *sxe = sy + n;
  while (sx < sxe)
    *sy++ = *sx++;
  while (sy < sxe)
    *sy++ = fill_val;
}

template<class T> inline void kaMatrix<T>::LeftShiftRow(const unsigned int i_row, unsigned long int n_col) {
  if (n_col == 0)
    return;

  unsigned int n = NumCols();
  if (n_col >= n)
    return;

  T *sy = Pointer()+i_row*n;
  T *sx = sy + n_col, *sxe = sy + n;
  while (sx < sxe)
    *sy++ = *sx++;
}

template<class T> inline void kaMatrix<T>::RepeatFirstColumn() {
  unsigned int n = NumCols();

  if (n == 1)
    return;

  T *sy = Pointer() + 1, *sx = sy - 1, *ssye = sx + n, *sye = Pointer()+Size();
  while (sx < sye) {
    while (sy < ssye)
      *sy++ = *sx;
    sx   += n;
    ssye += n;
    sy++;
  }
}

template<class T> inline void kaMatrix<T>::RepeatFirstValueInRow(const unsigned int i_row) {
  if (i_row > NumRows())
    return;

  unsigned int n = NumCols();
  T *sy = Pointer() +i_row*n + 1, *sx =  sy - 1, *sye = sy + n -1;
  while (sy < sye)
    *sy++ = *sx;
}

template<class T> inline void kaMatrix<T>::CopyColumn(const unsigned int k, const unsigned int l) {
  T *sx = Pointer()+k, *sxe = sx + Size();
  T *sy               = Pointer()+l;
  unsigned long int n = NumCols();

  while (sx < sxe) {
    *sy = *sx;
    sx += n;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::CopyColumn(kaMatrix<T> *ext, const unsigned int k, const unsigned int l) {
  T *sx = ext->Pointer() + k;
  T *sy = Pointer()+l;

  // T *sxe = sx + ext->Size();
  // T *sye = sy + Size();
  T *se;

  ext->NumRows() > NumRows() ? se = sy+Size() : se = sx+ext->Size();
  unsigned long int n = NumCols(), n_ext = ext->NumCols();
  while (sx < se) {
    *sy = *sx;
    sx += n_ext;
    sy += n;
  }
}

template<class T> inline void kaMatrix<T>::Null() {
  T *ps = Pointer(), *pse = (Pointer()+Size());

  while (ps < pse)
    *ps++ = T(0.0);
}

template<class T> inline void kaMatrix<T>::Fill(const T val) {
  T *ps = Pointer(), *pse = (Pointer()+Size());

  while (ps < pse)
    *ps++ = val;
}

template<class T> inline void kaMatrix<T>::FillRow(const unsigned int i_row, const T val) {
  if (i_row > NumRows())
    return;

  unsigned int n = NumCols();
  T *ps = (Pointer() + i_row*n), *pse = (ps+n);
  while (ps < pse)
    *ps++ = val;
}

template<class T> inline void kaMatrix<T>::FillColumn(const unsigned int j_column, const T val) {
  unsigned int n = NumCols();

  if (j_column > n)
    return;

  T *ps = Pointer() + j_column, *pse = (ps+Size());
  while (ps < pse) {
    *ps = val;
    ps += n;
  }
}

template<class T> T kaMatrix<T>::FrobeniusNorm() {
  T *ps = Pointer(), *pse = (Pointer()+Size());
  T  norm = 0.;

  while (ps < pse) {
    norm += *ps **ps;
    ps++;
  }
  return norm;
}

template<class T> void kaMatrix<T>::Print(const char *format) const {
  unsigned int i, j;

  for (i = 0; i < NumRows(); i++) {
    for (j = 0; j < NumCols(); j++) {
      cerr << Get(i, j) << " ";
    }
    printf("\n");
  }
  printf("\n");
}

template<class T> inline void kaMatrix<T>::Save(FILE *hFile) const {
  fwrite(Pointer(), sizeof(T), Size(), hFile);  // elements
}

template<class T> inline void kaMatrix<T>::Restore(FILE *hFile) {
  fread(Pointer(), sizeof(T), Size(), hFile);  // elements
}

//! Output operator. Writes matrix in text form
template<class T> ostream & operator<<(ostream &os, const kaMatrix<T> &rm) {
  unsigned int i, j;

  for (i = 0; i < rm.NumRows(); i++) {
    for (j = 0; j < rm.NumCols(); j++)
      os<<rm.Get(i, j)<<' ';

    os<<'\n';
  }

  return os;
}

#endif  // ifndef KAMATRIX_H
