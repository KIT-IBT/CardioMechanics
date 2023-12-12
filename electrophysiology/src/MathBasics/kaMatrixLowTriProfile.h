/*
 * File: kaMatrixLowTriProfile.h
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

#ifndef _kaMatrixLowTriProfile_h
#define _kaMatrixLowTriProfile_h


#include <kaMatrix.h>
#include <valarray>

#ifdef  __APPLE__
# include <Accelerate/Accelerate.h>

inline float _cblas_dot(unsigned int N, float *X, float *Y) {
  return cblas_sdot(N, X, 1, Y, 1);
}

inline double _cblas_dot(unsigned int N, double *X, double *Y) {
  return cblas_ddot(N, X, 1, Y, 1);
}

inline void _cblas_axpy(unsigned int N, float A, float *X, float *Y) {
  return cblas_saxpy(N, A, X, 1, Y, 1);
}

inline void _cblas_axpy(unsigned int N, double A, double *X, double *Y) {
  return cblas_daxpy(N, A, X, 1, Y, 1);
}

template<class T> class kaMatrixSparseCol;

#endif  // __APPLE__

/*!
   \class kaMatrixLowTriProfile

   Lower triangular matrix in skyline (profile) format.
   The rows are stored beginning from the leftmost nonzero element.

 */


template<class T> class kaMatrixLowTriProfile : public kaMatrixBase<T> {
  T *pElements;
  unsigned int *pRow;
  unsigned int n;
  unsigned int nelem;  // number of elements

 public:
  //! default constructor;
  kaMatrixLowTriProfile();
  /*! \param n is the dimension \\
     \param prow is the array of integers, prow[i] is the beginning of the row i in the matrix.

   */
  kaMatrixLowTriProfile(const unsigned int n, const unsigned int *prow);
  kaMatrixLowTriProfile(const kaMatrixSparseCol<T> &m);

  ~kaMatrixLowTriProfile() {Delete();}

  //! Access the storage array.
  inline T *Pointer() {return pElements;}

  //! Access the storage array.
  inline const T *Pointer() const {return pElements;}

  //! Access the internal array of row pointers

  /*!
     The i_th element of the array is the index of the first element of the row i in the storage array.
     The n_th element of teh array equals the number of the elements stored in the matrix.
   */
  inline unsigned int *RowPointer() {return pRow;}

  //! Access the internal array of row pointers
  inline const  unsigned int *RowPointer() const {return pRow;}

  inline unsigned int NumCols() const {return n;}

  inline unsigned int NumRows() const {return n;}

  inline unsigned int Dimension() const {return n;}

  //! R/W element access to the storage array
  inline T & a(const int p) {return (T &)*(pElements+p);}

  //! R/W element access by the element indices.
  /*! NOTE: i>=j !
     access to a non-existent element may cause a segementation fault
   */
  inline T & a(const int i, const int j) {return (T &)*(pElements+pRow[i+1]-i+j-1);}


  //! Use this if there is no need to change the elements.

  /*! NOTE: i>=j !
     access to a non-existent element may cause a segementation fault
   */
  inline const T Get(const int i, const int j) const {return *(pElements+pRow[i+1]-i+j-1);}

  //! Use this if there is no need to change the elements.
  inline const T Get(const int p) const {return *(pElements+p);}

  //! return index of the element in the storage array
  // returns -1 if the element does not exist
  inline int GetIndex(const int i, const int j) const {
    if (i < j)
      return -1;

    int p;
    if (pRow[i] > (p = pRow[i+1]-i+j-1))
      return -1;

    return p;
  }

  inline void Null();

  //! Access the row of the matrix (only valid indices will work!)
  /*!This operator can be used to access the elements of the matrix by the
     usual notation , e. g. a[i][j]
   */
  inline T *operator[](int i) {return pElements+pRow[i+1]-1-i;}

  //! Access the row of the matrix (only valid indices will work!)
  /*!This operator can be used to access the elements of the matrix by the
     usual notation , e. g. a[i][j]
   */
  inline const T *operator[](int i) const {return pElements+pRow[i+1]-1-i;}

  //! assess the size of the matrix in the skyline format given the profile info
  // \param prow and \param n are  the same as in the constructor
  static unsigned int MemorySize(const unsigned int n, const unsigned int *prow);

  //! size in bytes
  unsigned int MemorySize() const {return nelem*sizeof(T)+(n+1+2)*sizeof(unsigned int);}

  //! number of elements stored in the skyline format
  inline unsigned int Size() const {return nelem;}

  //! assignement operator
  template<class U>
  inline kaMatrixLowTriProfile<T> & operator=(const kaMatrixLowTriProfile<U> &rm) {
    if ((n != rm.NumCols()) || (nelem != rm.Size())) {
      Delete();
      n         = rm.NumCols();
      nelem     = rm.Size();
      pElements = new T[nelem];
      pRow      = new unsigned int[n+1];
    }

    const unsigned int *pcrm = rm.RowPointer();
    unsigned int *pc         = pRow;
    while (pcrm != rm.RowPointer()+n+1)
      *(pc++) = *(pcrm++);

    const U *perm = rm.Pointer();
    T *pe         = pElements;
    while (perm != rm.Pointer()+rm.Size())
      *(pe++) = *(perm++);

    return *this;
  }

  //! Conversion from kaMatrix, the lower triangle of the input matrix is taken.
  template<class U>
  inline kaMatrixLowTriProfile<T> & operator=(const kaMatrix<U> &rm) {
    Delete();
    n    = rm.NumCols() < rm.NumRows() ? rm.NumCols() : rm.NumRows();
    pRow = new unsigned int[n+1];


    int i, j;
    pRow[0] = 0;
    int jprev = 0;
    for (i = 1; i < n; i++) {
      j = 0;
      while (0.0 == rm.Get(i, j) && j < i) j++;
      pRow[i] = pRow[i-1]+i-jprev;
      jprev   = j;
    }
    nelem = pRow[n] = pRow[n-1]+n-jprev;

    pElements = new T[nelem];

    int k = 0;
    for (i = 0; i < n; i++) {
      j = 0;
      while (0.0 == rm.Get(i, j) && j < i) j++;
      while (j <= i) pElements[k++] = rm.Get(i, j++);
    }

    return *this;
  }  // =

  //! Assignement of kaMatrixLowTriProfile to kaMatrixSparseCol. Only overlapping part of the profile will be assigned.

  inline  kaMatrixLowTriProfile<T> & operator=(const kaMatrixSparseCol<T> &rm);

  void Print(const char *format = "%7.3lf\t") const {
    unsigned int i, j;

    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++)
        printf(format, pRow[i+1]-i+j-1 < pRow[i] ? 0 : Get(i, j) ); // access in the "right" order
      printf("\n");
    }
  }

  /*!
     compute Cholesky decomposition of the matrix
     returns 1 if OK
     returns 0 if failed.
   */
  int cholesky();
  /*!
     Assuming that this is the lower Cholesky factor of an equation system matrix,
     solve the equation system by forward and back substitution.
     Vector \param y is overwritten by the solution.
   */
  void chol_solve(T *y);

  //! Scale matrix to get a unit diagonal.
  /*!Scale matrix get a unit diagonal.\n
     Result: scale(i) = 1/A(i,i)
   */
  void Scale(valarray<T> &scale);

  //! Scale matrix to get a unit diagonal.
  /*! Scale equation sysytem to get ones on the diagonal of the
     system matrix. \n
     Only upper or lower part of the system matrix is needed.\n
     \param B is the right-hand vector.
     \param scale scaling factors: scale(i) = 1/sqrt(A(i,i)).
     For equation Ax=b the scaling works as follows:\n
     diag(scale)*A*diag(scale)*inv(diag(scale))*x=diag(scale)*b \n
     => we get z=inv(diag(scale))*x as a solution.\n
     Therefore the solution must be 'unscaled':  x = diag(scale)*x
   */
  void SymmScale(valarray<T> &scale);

  //! Set rows and colums of matrix to zero if conditions are specified on correspinding nodes.
  /*!\param c is the vector of conditions, size of  c equals the dimension of the matrix. If
     c[i]!=0, the row and column of the matrix are set to 0, the diagonal element i is set to 1.0.
   */
  void ModifyForConditions(const unsigned int *cind, const unsigned int size);

  //! use  SymmMatrixVectorProduct if only the lower (upper) triangle of A is stored
  void SymmMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const;

  //! not implemented
  void Save(FILE *) const {cerr <<" kaMatrixLowTriProfile::Save is not implemented \n";}

  //! not implemented
  void Restore(FILE *) {cerr <<" kaMatrixLowTriProfile::Restore is not implemented \n";}

 private:
  void Delete() {
    delete[] pElements;
    delete[] pRow;
  }

  template<class U> friend ostream & operator<<(ostream &os, const  kaMatrixLowTriProfile<U> &rm);
};  // class kaMatrixLowTriProfile

//! Output operator. Writes matrix in text form
template<class U> ostream & operator<<(ostream &os, const  kaMatrixLowTriProfile<U> &rm) {
  unsigned int i, j;

  for (i = 0; i < rm.n; i++) {
    for (j = 0; j < rm.n; j++)
      if (i >= j)
        os<< (rm.pRow[i+1]-i+j-1 < rm.pRow[i] ? 0 : rm.Get(i, j) )<<' ';
      else
        os<<0.0<<' ';

    os<<'\n';
  }
  return os;
}

// === implementation
template<class T>
kaMatrixLowTriProfile<T>::kaMatrixLowTriProfile() {
  n         = 0;
  nelem     = 0;
  pRow      = NULL;
  pElements = NULL;
}

template<class T>
kaMatrixLowTriProfile<T>::kaMatrixLowTriProfile(const unsigned int nn, const unsigned int *prow) {
  n    = nn;
  pRow = new unsigned int[n+1];

  pRow[0] = 0;

  for (int i = 1; i < n+1; i++)
    pRow[i] = pRow[i-1]+i-prow[i-1];

  nelem = pRow[n];

  pElements = new T[nelem];
}

template<class T>
kaMatrixLowTriProfile<T>::kaMatrixLowTriProfile(const kaMatrixSparseCol<T> &rm) {
  // 1) reserve space
  n       = rm.Dimension();
  pRow    = new unsigned int[n+1];
  pRow[0] = 0;

  // get profile of the matrix
  unsigned int *profile = new unsigned int[n];
  for (int i = 0; i < n; i++) profile[i] = i;

  for (int j = 0; j < rm.Dimension(); j++)
    for (int k = rm.CCLPTR[j]; k < rm.CCLPTR[j+1]; k++) {
      int row = rm.CRNUM[k];
      if (profile[row] > j)
        profile[row] = j;
    }

  for (int i = 1; i < n+1; i++)
    pRow[i] = pRow[i-1]+i-profile[i-1];
  nelem     = pRow[n];
  pElements = new T[nelem];

  delete[] profile;

  // 2) perform assignment
  *this = rm;
}

template<class T>
void kaMatrixLowTriProfile<T>::Null() {
  T *pe  = pElements;
  T *ppe = pElements + nelem;

  while (pe != ppe) *(pe++) = 0.0;
}

template<class T>
unsigned int kaMatrixLowTriProfile<T>::MemorySize(const unsigned int n, const unsigned int *prow) {
  unsigned int prev = 0, next = 0;

  for (int i = 1; i < n+1; i++) {
    next = prev + i - prow[i-1];
    prev = next;
  }
  return next*sizeof(T)+(n+1+2)*sizeof(unsigned int);
}

template<class T>
kaMatrixLowTriProfile<T> &  kaMatrixLowTriProfile<T>::operator=(const kaMatrixSparseCol<T> &rm) {
  int j;

  for (j = 0; j < rm.NumCols(); j++) {
    a(j, j) = rm.CDIAG[j];
    int clptrbeg = rm.CCLPTR[j];
    int clptrend = rm.CCLPTR[j+1];
    int i, rnum;
    int p;
    for (i = clptrbeg; i < clptrend; i++) {
      rnum = rm.CRNUM[i];
      if (-1 != (p = GetIndex(rnum, j)) )
        a(p) =  rm.COFDIAG[i];
    }
  }
  return *this;
}

template<class T>
int kaMatrixLowTriProfile<T>::cholesky() {
  int i, j;
  T   adiag;

  int ibeg, iend, jbeg, jend;

  unsigned int li, lj;

  T *psum, *pdiag;
  double sum;
  T *pi, *pj;

  for (j = 0; j < n; j++) {
    if (!(j%1000) ) {cout<<"j = "<<j<<"  \r"; cout.flush();}
    jbeg = pRow[j];
    jend = pRow[j+1]-1;  // the needed part of the row j (all except the diag element)
    // this is the diagonal element jj

    lj = jend - jbeg;  // working length of the row j

    // proceed i==j
    pdiag = pElements+jend;
    sum   = *pdiag;
    pj    = pElements+jbeg;
    while (pj != pdiag) {
      sum -= *pj **pj;
      pj++;
    }

    if (sum <= 0.0)
      return 0;

    *pdiag = adiag = sqrt(sum);

#pragma omp parallel for private(ibeg,iend,psum,sum,li,pi,pj)
    for (i = j+1; i < n; i++) {
      ibeg = pRow[i];
      iend = pRow[i+1]-1-i+j;  // from the beginning to j of the row i
      // iend is also the index of the element i,j

      if (iend >= ibeg) {
        psum = pElements+iend;
        sum  = *psum;
        li   = iend - ibeg; // working length of the row i

        unsigned int length = lj < li ? lj : li;

        pi = psum  - length;
        pj = pdiag - length;

#ifndef __APPLE__
        while (pi != psum)
          sum -= *(pi++) **(pj++);
#else  // ifndef __APPLE__
        sum -= _cblas_dot(length, pi, pj);
#endif  // ifndef __APPLE__
        *psum = sum/adiag;
      }
    }
  }
  return 1;
}  // >::cholesky

template<class T>
void kaMatrixLowTriProfile<T>::chol_solve(T *y) {
  // forward substitution (row-oriented)
  *y /= *pElements;

  // unsigned int i,j;
  //  unsigned int kbeg, kend;//, py;
  T *py, *pkbeg, *pkend, *py_out, *py_out_end;
  double sum;

  unsigned int *ppRow;

  ppRow = pRow+1;

  py_out     = y + 1;
  py_out_end = y+n;

  while (py_out != py_out_end) {  //  for( i = 1; i<n; i++)
    sum   = *py_out;
    pkbeg = pElements + *(ppRow++);  // beginning of the row i
    pkend = pElements+ *ppRow -1;  // diag element in the row i

    py = py_out+(pkbeg-pkend);

#ifndef __APPLE__
    while (pkbeg != pkend)
      sum -= *(pkbeg++) **(py++);
#else  // ifndef __APPLE__
    sum -= _cblas_dot(pkend-pkbeg, pkbeg, py);
#endif  // __APPLE__

    sum /= *pkend;

    *py_out = sum;
    py_out++;
  }

  // back substitution (column-oriented)

  ppRow  = pRow + n-1;
  py_out = y+n-1;

  // for( j = n-1; j>0; j--)
  while (py_out != y) {
    pkbeg    = pElements + *ppRow;
    pkend    = pElements + *(ppRow+1)-1;
    py       = py_out+(pkbeg-pkend);
    *py_out /= *pkend;

#ifndef __APPLE__
    while (pkbeg != pkend)
      *(py++) -= (*py_out) **(pkbeg++);
#else  // ifndef __APPLE__
    _cblas_axpy(pkend-pkbeg, -(*py_out), pkbeg, py);
#endif  // __APPLE__

    py_out--;
    ppRow--;
  }
  *y /= *pElements;
}  // >::chol_solve

template<class T>
void kaMatrixLowTriProfile<T>::Scale(valarray<T> &scale) {
  scale[0]     = 1.0/pElements[0];
  pElements[0] = 1.0;

  for (int i = 1; i < n; i++) {
    int pdiag  = pRow[i+1]-1;
    T   factor = 1.0/pElements[pdiag];
    pElements[pdiag] = 1.0;
    scale[i]         = factor;
    T *pj     = pElements+pRow[i];
    T *pj_end = pElements + pdiag;
    while (pj != pj_end)
      *(pj++) *= factor;
  }
}

template<class T>
void kaMatrixLowTriProfile<T>::SymmScale(valarray<T> &scale) {
  scale[0]     = 1.0/sqrt(pElements[0]);
  pElements[0] = 1.0;

  for (int i = 1; i < n; i++) {
    int pdiag  = pRow[i+1]-1;
    T   factor = 1.0/sqrt(pElements[pdiag]);
    pElements[pdiag] = 1.0;
    scale[i]         = factor;
    for (int pj = pRow[i]; pj < pdiag; pj++)
      pElements[pj] *= factor*scale[i - pdiag + pj];
  }
}

template<class T>
void kaMatrixLowTriProfile<T>::ModifyForConditions(const unsigned int *cind, const unsigned int size) {
  valarray<unsigned char> flags((unsigned char)0, Dimension());
  for (int i = 0; i < size; i++)
    flags[cind[i]] = 1;

  if (flags[0])
    pElements[0] = 1.0;

  for (int i = 1; i < n; i++) {
    int pdiag = pRow[i+1]-1;
    if (flags[i])
      pElements[pdiag] = 1.0;
    for (int pj = pRow[i]; pj < pdiag; pj++)
      if (flags[i] || flags[i-pdiag+pj])
        pElements[pj] = 0.0;
  }
}

template<class T>
inline void kaMatrixLowTriProfile<T>::SymmMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const {
  r = x;
  for (int i = 0; i < n; i++)
    r *= pElements[pRow[i+1]-1];

  for (int i = 1; i < n; i++) {
    int pdiag = pRow[i+1]-1;
    T   xi    = x[i];
    T   ri    = 0; // df: should it be here?
    for (int pj = pRow[i]; pj <= pdiag; pj++) {
      //        T ri = 0;   // df: was standing here
      int j = i-pdiag+pj;
      ri += pElements[pj]*(x[j]+xi);
    }
    r[i] += ri;
  }
}

#endif  // ifndef _kaMatrixLowTriProfile_h
