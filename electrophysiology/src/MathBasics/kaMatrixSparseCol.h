/*
 * File: kaMatrixSparseCol.h
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


/*! \class kaMatrixSparseCol
   Sparse matrix. Compressed column storage with distinct diagonal.
 */

#ifndef _kaMatrixSparseCol_h
#define _kaMatrixSparseCol_h


#include <kaExceptions.h>
#include <kaMatrixBase.h>
#include <algorithm>
#include <valarray>
#include <map>

/* Enables saving in PETSc format
 * This should be defined in Imakefile! (because it changes includes and librarys)
 * Exmaple: ../FESolver/Imakefile (NOT only line 5!) */
#ifdef WITH_PETSC
# include <petscmat.h>
# include <PETScLSE.h>
#endif  // ifdef WITH_PETSC

using namespace std;

template<class T> class kaMatrixLowTriProfile;

template<class T>
class kaMatrixSparseCol : public kaMatrixBase<T> {
  unsigned int MAXN;  // max number of rows and columns
  unsigned int NZ;   // number of off-diagonal non-zeroes

 public:
  valarray<T> CDIAG;  // size MAXN
  valarray<T> COFDIAG;  // size NZ
  valarray<int> CCLPTR;  // size MAXN+1
  valarray<int> CRNUM;  // size NZ


  //! Number of rows and columns and number of non-zeros are supplied.
  kaMatrixSparseCol(const unsigned int maxn, const unsigned int nz) : MAXN(maxn), NZ(nz), CDIAG(maxn), COFDIAG(nz),
    CCLPTR(maxn+1), CRNUM(nz) {}

  //! Copy constructor
  kaMatrixSparseCol(const kaMatrixSparseCol<T> &a) : MAXN(a.MAXN), NZ(a.NZ), CDIAG(a.CDIAG), COFDIAG(a.COFDIAG), CCLPTR(
      a.CCLPTR), CRNUM(a.CRNUM) {}

  //! Constructor creating a submatrix
  kaMatrixSparseCol(kaMatrixSparseCol<T> &a, const std::valarray<size_t> &indices);

  ~kaMatrixSparseCol() {}

  kaMatrixSparseCol<T> & operator=(const kaMatrixSparseCol<T> &);

  //! matrices must have at least the same size!
  // only reserved elements in the matrix will be copied
  kaMatrixSparseCol<T> & operator=(const kaMatrixLowTriProfile<T> &);
  //! replace the whole structure with that of <i>matrix</i>
  void assignMatrix(const kaMatrixSparseCol<T> &matrix);

  inline unsigned int Dimension() const {return MAXN;}

  inline unsigned int NumRows() const {return MAXN;}

  inline unsigned int NumCols() const {return MAXN;}

  //! Number of offdiagonal non-zeros
  inline unsigned int NumNZ() const {return NZ;}
  //! Read-only access (NOT an optimal way to acces elements. Involves search procedure).
  inline const T Get(int, int) const;

  //! Read/write access (NOT an optimal way to acces elements. Involves search procedure.)
  inline T & a(int, int);
  //! Set the element of the matrix to zero and return its previous value.
  inline const T zero(const int i, const int j);

  inline void Null() {
    CDIAG   = 0;
    COFDIAG = 0;
  }

  //! not implemented
  const T Get(int) const {
    throw kaBaseException("kaMatrixSparseCol::Get(int) not implemented");
    return 0;
  }

  //! not implemented
  T & a(int) {
    throw kaBaseException("kaMatrixSparseCol::a(int) not implemented");
    return (T &)CDIAG[0];
  }

  //! not implemented
  T *operator[](int i) {
    throw kaBaseException("kaMatrixSparseCol::operator[] not implemented");
    return 0;
  }

  //! not implemented
  const T *operator[](int i) const {
    throw kaBaseException("kaMatrixSparseCol::operator[] not implemented");
    return 0;
  }

  inline T *Pointer() {
    throw kaBaseException("kaMatrixSparseCol::Pointer() not implemented");
    return NULL;
  }

  inline const T *Pointer() const {
    throw kaBaseException("kaMatrixSparseCol::Pointer() not implemented");
    return NULL;
  }

  //! Size in bytes.
  inline unsigned int MemorySize() const {
    return sizeof(T)*(CDIAG.size()+COFDIAG.size())+sizeof(int)*(CCLPTR.size()+CRNUM.size());
  }

  //! Not an optimal way to save the matrix.
  // For debug purpose only.
  void Print(const char *format = "%7.3lf\t") const {
    unsigned int i, j;

    for (i = 0; i < MAXN; i++) {
      for (j = 0; j < MAXN; j++)
        printf(format, Get(i, j));
      printf("\n");
    }
  }

#ifdef WITH_PETSC
  void SavePETSc(const char *, int, char **, IOType = ot_bin, int = 1);
#endif  // ifdef WITH_PETSC
  void Save(FILE *) const;
  void Restore(FILE *);
  /*!Save in ASCII format, which can be read by MATLAB.
     Indices start with 1.
     Set symmetric to 1 to save the full matrix except of only lower (upper) triangular part.
   */
  int         save_ascii(const char *name, const int symmetric = 0);
  inline void MatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const;
  inline void TranspMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const;
  //! use  SymmMatrixVectorProduct if only the lower (upper) triangle of A is stored
  inline void SymmMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const;
  inline T    FrobeniusNorm(const int symmetric = 0) const;
  //! Standard incomplete Cholesky decomposition.
  int istdic();

  //! Jones/Plassman incomplete Cholesky: column oriented.
  /*! Jones/Plassman incomplete Cholesky: column oriented
     ACM Transactions on Mathematical Software, Vol 21, No 1, March 1995, pp. 5-17
     FORTRAN version is available via netlib
   */
  int jpicc();

  //! Cholesky-type factorization solve.
  /*! Cholesky-type factorization solve. \n
     The matrix is a lower triangular matrix. \n
     The procedure solves A*A'*x=y   to find x. \n
   */
  void chol_solve(const valarray<T> &y, valarray<T> &x);

  //! Scale equation system to get a unit diagonal.
  /*!Scale matrix get a unit diagonal.\n
     Result: scale(i) = 1/A(i,i)
   */
  void Scale(valarray<T> &scale);

  //! Scale equation system to get a unit diagonal.
  /*! Scale equation system to get ones on the diagonal of the
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

  /*!\param cind is the list of nodes where the conditions are set.
     The corresponding row and column of the matrix are set to 0, the diagonal element is set to 1.0.
   */
  inline void ModifyForConditions(const unsigned int *cind, const unsigned int size) {
    valarray<unsigned char> flags((unsigned char)0, Dimension());
    for (int i = 0; i < size; i++)
      flags[cind[i]] = 1;

    for (int i = 0; i < MAXN; i++)
      if (flags[i])
        CDIAG[i] = 1.0;

    for (int j = 0; j < MAXN; j++)
      for (int k = CCLPTR[j]; k < CCLPTR[j+1]; k++)
        if (flags[CRNUM[k]] || flags[j])
          COFDIAG[k] = 0.0;
  }

 private:
  // get pointer to a element
  // inline T * pa(int,int);
  // bisection search for i in CCLPTR
  inline int search(const int i, const int j) const;
};  // class kaMatrixSparseCol

// class to implement isort for jpicc
template<class T>
class isort {
  const T *keys;

 public:
  inline int operator()(const int &a1, const int &a2) const {
    return fabs(*(keys+a1)) > fabs(*(keys+a2));
  }

  isort(const T *k) : keys(k) {}

  inline static void sort(const int &n, const int &k, valarray<T> &akeys, valarray<int> &indvec) {
    isort<T> Cmp(&akeys[0]);
    int *start = &indvec[0];
    int *end   = &indvec[n-1];
    int *mid   = &indvec[k-1];
    partial_sort(start, mid+1, end+1, Cmp);
  }
};


//! Constructor creating a submatrix
template<class T>
kaMatrixSparseCol<T>::kaMatrixSparseCol(kaMatrixSparseCol<T> &a, const std::valarray<size_t> &indices) {
  MAXN = indices.size();
  if (!MAXN)
    throw kaBaseException("kaMatrixSparseCol::kaMatrixSparseCol: empty list of indices");

  // fill up diagonal elements
  CDIAG.resize(MAXN);
  for (size_t i = 0; i < MAXN; i++) CDIAG[i] = a.CDIAG[indices[i]];

  // create backward mapper
  std::valarray<int> idxBack(-1, a.Dimension());
  for (size_t i = 0; i < MAXN; i++) idxBack[indices[i]] = (int)i;

  // compute the number of non-diagonal elements: overall and in each row
  NZ = 0;
  std::valarray<size_t> elementsInCol((size_t)0, MAXN);
  std::valarray<std::map<size_t, T>> matrixStructure(MAXN);  // structure of the future matrix

  for (size_t i = 0; i < a.Dimension(); i++) {
    for (int j = a.CCLPTR[i]; j < a.CCLPTR[i+1]; j++) {
      int idx1 = idxBack[i];
      int idx2 = idxBack[a.CRNUM[j]];
      if ((idx1 < 0) || (idx2 < 0) )
        continue; // skip this element of a
      NZ++;
      if (idx1 < idx2) {  // idx1 = column index, idx2 = row index
        elementsInCol[idx1]++;
        matrixStructure[idx1][idx2] = a.COFDIAG[j];
      } else {           // idx1 = row index, idx2 = column index
        elementsInCol[idx2]++;
        matrixStructure[idx2][idx1] = a.COFDIAG[j];
      }
    }
  }

  // create the data
  COFDIAG.resize(NZ);
  CCLPTR.resize(MAXN+1);
  CRNUM.resize(NZ);
  CCLPTR[0] = 0;

  for (size_t col = 0; col < MAXN; col++) {
    CCLPTR[col+1] = CCLPTR[col] + elementsInCol[col];
    size_t counter = 0;
    for (typename std::map<size_t, T>::iterator mit = matrixStructure[col].begin(); mit != matrixStructure[col].end();
         mit++) {
      size_t offset = CCLPTR[col]+(counter++);
      CRNUM[offset]   = mit->first;
      COFDIAG[offset] = mit->second;
    }
  }
}

template<class T>
int kaMatrixSparseCol<T>::search(const int i, const int j) const {
  int nb = CCLPTR[j];
  int ne = CCLPTR[j+1];

  if (ne == nb)
    return -1;

  while (ne-nb > 1) {
    int pos = (nb+ne)/2;
    if (CRNUM[pos] > i)
      ne = pos; else
      nb = pos;
  }
  if (CRNUM[nb] != i)
    return -1;

  return nb;
}

template<class T>
void kaMatrixSparseCol<T>::assignMatrix(const kaMatrixSparseCol<T> &matrix) {
  MAXN = matrix.MAXN;
  NZ   = matrix.NZ;

  CDIAG.resize(MAXN);
  COFDIAG.resize(NZ);
  CCLPTR.resize(MAXN+1);
  CRNUM.resize(NZ);

  CDIAG   = matrix.CDIAG;
  COFDIAG = matrix.COFDIAG;
  CCLPTR  = matrix.CCLPTR;
  CRNUM   = matrix.CRNUM;
}

template<class T>
kaMatrixSparseCol<T> & kaMatrixSparseCol<T>::operator=(const kaMatrixSparseCol<T> &a) {
  if (MAXN != a.MAXN) {
    CDIAG.resize(a.MAXN);
    CCLPTR.resize(a.MAXN+1);
  }
  if (NZ != a.NZ) {
    COFDIAG.resize(NZ);
    CRNUM.resize(NZ);
  }

  CDIAG   = a.CDIAG;
  COFDIAG = a.COFDIAG;
  CCLPTR  = a.CCLPTR;
  CRNUM   = a.CRNUM;
  return *this;
}

template<class T>
kaMatrixSparseCol<T> & kaMatrixSparseCol<T>::operator=(const kaMatrixLowTriProfile<T> &rm) {
  int j;

  for (j = 0; j < MAXN; j++) {
    CDIAG[j] = rm.Get(j, j);
    int clptrbeg = CCLPTR[j];
    int clptrend = CCLPTR[j+1];
    int i, rnum;
    int p;
    for (i = clptrbeg; i < clptrend; i++) {
      rnum = CRNUM[i];
      if (-1 != (p = rm.GetIndex(rnum, j)) )
        COFDIAG[i] = rm.Get(p);
    }
  }
  return *this;
}

template<class T>
const T kaMatrixSparseCol<T>::Get(int i, int j) const {
  int p;

  if (i == j)
    return CDIAG[i];

  if ( (p = search(i, j)) < 0)
    return 0;

  return COFDIAG[p];
}

template<class T>
T &  kaMatrixSparseCol<T>::a(int i, int j) {
  int p;

  if (i == j)
    return CDIAG[i];

  if ( (p = search(i, j) ) < 0)
    throw kaBaseException("kaMatrixSparseCol::a(int,int): the element does not exist");
  return COFDIAG[p];
}

template<class T>
const T kaMatrixSparseCol<T>::zero(const int i, const int j) {
  T *p;
  int k;
  T   temp;

  if (i == j)
    p = &CDIAG[i];
  else if ( (k = search(i, j) ) < 0)
    return 0;
  else
    p = &COFDIAG[k];

  temp = *p;
  *p   = 0;
  return temp;
}

template<class T>
void kaMatrixSparseCol<T>::MatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const {
  r = 0.0;
  int j;
  T   xj;

#pragma omp parallel for private(j, xj)
  for (j = 0; j < MAXN; j++) {
    xj    = x[j];
    r[j] += CDIAG[j]*xj;
    int clptrbeg = CCLPTR[j];
    int clptrend = CCLPTR[j+1];
    int i, rnum;
    T   a;
    for (i = clptrbeg; i < clptrend; i++) {
      rnum = CRNUM[i];
      a    = COFDIAG[i]*xj;
#pragma omp atomic
      r[rnum] += a;
    }
  }
}

template<class T>
void kaMatrixSparseCol<T>::TranspMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const {
  r = 0.0;
  int j;
  T   rj;
#pragma omp parallel for private(j, rj)
  for (j = 0; j < MAXN; j++) {
    rj = CDIAG[j]*x[j];
    for (int i = CCLPTR[j]; i < CCLPTR[j+1]; i++)
      rj += COFDIAG[i]*x[CRNUM[i]];
    r[j] += rj;
  }
}

template<class T>
void kaMatrixSparseCol<T>::SymmMatrixVectorProduct(const valarray<T> &x, valarray<T> &r) const {
  r = 0.0;
  int j1;
  T   xj1;
#pragma omp parallel for private(xj1)
  for (j1 = 0; j1 < MAXN; j1++) {
    xj1 = x[j1];
    int clptrbeg = CCLPTR[j1];
    int clptrend = CCLPTR[j1+1];
    int i, rnum;
    T   a;
    for (i = clptrbeg; i < clptrend; i++) {
      rnum = CRNUM[i];
      a    = COFDIAG[i]*xj1;
#pragma omp atomic
      r[rnum] += a;
    }
  }

  int j2;
  T   rj2;
#pragma omp parallel for private(rj2)
  for (j2 = 0; j2 < MAXN; j2++) {
    rj2 = CDIAG[j2]*x[j2];
    for (int i = CCLPTR[j2]; i < CCLPTR[j2+1]; i++)
      rj2 += COFDIAG[i]*x[CRNUM[i]]; // r[j]
    r[j2] += rj2;
  }
}  // >::SymmMatrixVectorProduct

template<class T>
T kaMatrixSparseCol<T>::FrobeniusNorm(const int symmetric) const {
  T norm = 0, offdiagnorm = 0;
  int i;

#pragma omp parallel for private(i) reduction(+:norm)
  for (i = 0; i < MAXN; i++) norm += CDIAG[i]*CDIAG[i];

  int k;
#pragma omp parallel for private(k) reduction(+:offdiagnorm)
  for (k = 0; k < NZ; k++) offdiagnorm += COFDIAG[k]*COFDIAG[k];
  if (symmetric)
    return sqrt(norm+2.0*offdiagnorm);

  return sqrt(norm+offdiagnorm);
}

template<class T>
int kaMatrixSparseCol<T>::save_ascii(const char *name, const int symmetric) {
  ofstream os(name);

  if (!os)
    return 0;

  os<<setprecision(15);

  for (int j = 0; j < MAXN; j++) {
    if (CDIAG[j])
      os<<j+1<<' '<<j+1<<' '<<CDIAG[j]<<'\n';

    for (int i = CCLPTR[j]; i < CCLPTR[j+1]; i++) {
      os<<CRNUM[i]+1<<' '<<j+1<<' '<<COFDIAG[i]<<'\n';
      if (symmetric)
        os<<j+1<<' '<<CRNUM[i]+1<<' '<<COFDIAG[i]<<'\n';
    }
  }

  return 1;
}

template<class T>
void kaMatrixSparseCol<T>::chol_solve(const valarray<T> &y, valarray<T> &x) {
  valarray<T> temp_y(y);

  /* first solve  L.z=y where z=Transp(L).x  here z overwrites x*/
  for (int i = 0; i < MAXN; i++) {
    x[i] = temp_y[i]/CDIAG[i];
    for (int k = CCLPTR[i]; k < CCLPTR[i+1]; k++)
      temp_y[CRNUM[k]] -= x[i]*COFDIAG[k];
  }

  /* now solve Transp(L).x=z */
  for (int i = MAXN-1; i >= 0; i--) {
    T sum = 0;
    for (int k = CCLPTR[i]; k < CCLPTR[i+1]; k++)
      sum += x[CRNUM[k]]*COFDIAG[k];
    x[i] = (x[i]-sum)/CDIAG[i];
  }
}

template<class T>
int kaMatrixSparseCol<T>::istdic() {
  int isk, iek, isj, iej;
  int i, j, k;
  int row;
  T   lval, t;
  int iptr;

  const int N = MAXN;  // dimensionality of input matrix

  valarray<T> ta(MAXN);  // a temporary work vector of length N to keep the current column
  valarray<int> ifirst(-1, MAXN);  /* ifirst(j) points to the next value in  column j to use.
                                      ifirst also has a dual use. At step K, only  the first K-1 elements are used for
                                         the
                                      above purpose. For the last N-K elements, ifirst(j) indicates if a nonzero value
                                      exists in the position j of column K */
  valarray<int> List(-1, MAXN);  /* List(j) points to the linked list of columns which will update column j */

  for (k = 0; k < N; k++) {  // over all columns
    // load column k into ta
    isk = CCLPTR[k];
    iek = CCLPTR[k+1];
    for (j = isk; j < iek; j++) {
      row         = CRNUM[j];
      ta[row]     = COFDIAG[j];
      ifirst[row] = 0;
    }

    // make sure the diagonal of k is okay and then take the sqrt
    if (CDIAG[k] < 0)
      return -k;
    else
      CDIAG[k] = sqrt(CDIAG[k]);

    // update column k using previous columns
    j = List[k];
    while (j != -1) {
      isj  = ifirst[j];
      iej  = CCLPTR[j+1]-1;
      lval = COFDIAG[isj];
      isj++;
      if (isj < iej) {
        ifirst[j]        = isj;
        iptr             = j;
        j                = List[j];
        List[iptr]       = List[CRNUM[isj]];
        List[CRNUM[isj]] = iptr;
      } else {
        j = List[j];
      }
      for (i = isj; i <= iej; i++) {
        row = CRNUM[i];
        if (ifirst[row] != -1)
          ta[row] -= lval*COFDIAG[i];
      }
    }

    // ifirst and List keep track of where in column k we are
    if (isk < (iek-1)) {
      iptr       = CRNUM[isk];
      List[k]    = List[iptr];
      List[iptr] = k;
      ifirst[k]  = isk;
    }

    // update remaining diagonals using column k
    lval = CDIAG[k];
    for (j = isk; j < iek; j++) {
      row         = CRNUM[j];
      t           = ta[row];
      ifirst[row] = -1;
      t          /= lval;
      CDIAG[row]  = CDIAG[row]-t*t;
      COFDIAG[j]  = t;
    }
  }
  return 0;
}  // >::istdic

template<class T>
int kaMatrixSparseCol<T>::jpicc() {
  const int N = MAXN;  // dimensionality of input matrix

  valarray<T> ta(N);  // a vector to keep the current column
  int talen = 0;  // actual size of ta

  valarray<int> itcol(N);  // keeps the row values of the current column

  valarray<int> ifirst(-1, N);  /* ifirst(j) points to the next value in  column j to use.
                                   ifirst also has a dual use. At step K, only  the first K-1 elements are used for the
                                   above purpose. For the last N-K elements, ifirst(j) indicates if a nonzero value
                                   exists in the position j of column K */
  valarray<int> List(-1, N);  /* List(j) points to the linked list of columns which will update column j */

  /* here begins ... */
  for (int k = 0; k < N; k++) {  // over all columns
    // load column k into ta
    talen = 0;
    int isk = CCLPTR[k];
    int iek = CCLPTR[k+1];
    for (int j = isk; j < iek; j++) {
      int row = CRNUM[j];
      ta[row]      = COFDIAG[j];
      itcol[talen] = row;
      ifirst[row]  = 0;
      talen++;
    }

    // take the sqrt of the diagonal of column k
    if (CDIAG[k] < 0)
      return -k; // factorization failed
    else
      CDIAG[k] = sqrt(CDIAG[k]);

    // update column k using the previous columns
    int j = List[k];
    while (j != -1) {
      int isj  = ifirst[j];
      int iej  = CCLPTR[j+1]-1;
      T   lval = COFDIAG[isj];
      isj++;
      if (isj < iej) {
        ifirst[j] = isj;
        int iptr = j;
        j                = List[j];
        List[iptr]       = List[CRNUM[isj]];
        List[CRNUM[isj]] = iptr;
      } else {j = List[j];}

      for (int i = isj; i <= iej; i++) {
        int row = CRNUM[i];
        if (ifirst[row] != -1) {ta[row] -= lval*COFDIAG[i];} else {
          ifirst[row]  = 0;
          itcol[talen] = row;
          ta[row]      = (-1.0)*lval*COFDIAG[i];
          talen++;
        }
      }
    }


    // update remaining diagonals using column k
    for (int j = 0; j < talen; j++) {
      int row = itcol[j];
      ta[row]    /= CDIAG[k];
      CDIAG[row] -= ta[row]*ta[row];
    }


    // find the lagest elements in column k
    int klen  = iek-isk;
    int count = (klen < talen) ? klen : talen;  // k=min(klen,talen)
    isort<T>::sort(talen, count, ta, itcol);

    // sort itcol
    int *first = &itcol[0];
    int *last  = first+count;
    sort(first, last);

    // put the lagest elements back into sparse data structure
    count = 0;
    for (int j = isk; j < iek; j++) {
      COFDIAG[j] = ta[itcol[count]];
      CRNUM[j]   = itcol[count];
      count++;
    }

    // ifirst and List keep track of where in column k we are
    if (isk < (iek-1)) {
      int iptr = CRNUM[isk];
      List[k]    = List[iptr];
      List[iptr] = k;
      ifirst[k]  = isk;
    }
    for (int j = 0; j < talen; j++) {
      ifirst[itcol[j]] = -1;
    }
  }
  return 0;
}  // >::jpicc

template<class T>
void kaMatrixSparseCol<T>::Scale(valarray<T> &scale) {
  for (int i = 0; i < MAXN; i++) {
    if (CDIAG[i] != 0) {
      T factor = 1.0/CDIAG[i];
      CDIAG[i] = factor;

      //      B[i]*=factor;
    } else {throw kaBaseException("Scale: zero diagonal element encountered\n");}
  }
  for (int j = 0; j < MAXN; j++)
    for (int k = CCLPTR[j]; k < CCLPTR[j+1]; k++) {
      int row = CRNUM[k];
      COFDIAG[k] *= CDIAG[row];
    }
  scale = CDIAG;
  CDIAG = 1.0;
}

template<class T>
void kaMatrixSparseCol<T>::SymmScale(valarray<T> &scale) {
  for (int i = 0; i < MAXN; i++) {
    if (CDIAG[i] != 0.0) {
      T factor = 1.0/sqrt(CDIAG[i]);
      CDIAG[i] = factor;
    } else {throw kaBaseException("SymmScale: zero diagonal element encountered\n");}
  }
  for (int j = 0; j < MAXN; j++)
    for (int k = CCLPTR[j]; k < CCLPTR[j+1]; k++) {
      int row = CRNUM[k];
      COFDIAG[k] *= CDIAG[row]*CDIAG[j];
    }
  scale = CDIAG;
  CDIAG = 1.0;
}

#ifdef WITH_PETSC

template<class T>
void kaMatrixSparseCol<T>::SavePETSc(const char *file, int argc, char **argv, IOType ot, int scale) {
  /* define variables*/
  PetscErrorCode ierr;
  Mat pmat;
  PetscViewer fd;
  MPI_Comm comm            = MPI_COMM_SELF;
  PetscInt *nonzerosPerRow = new PetscInt[MAXN];

  kaSharedMatrixN<double> dummy(4);

  /* initialize PETSc */
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  /* No multi-CPU support in derived classes! */
  int numCPU;
  MPI_Comm_size(comm, &numCPU);
  if (numCPU > 1) {
    ierr = PetscFinalize(); CHKERRQ(ierr);
    throw kaBaseException("kaMatrixSparseCol::SavePETSc: Can only save with single CPU. Number of CPUs is %i.", numCPU);
  }

  /* Create matrix */
  ierr = MatCreate(comm, &pmat); CHKERRQ(ierr);
  ierr = MatSetSizes(pmat, PETSC_DECIDE, PETSC_DECIDE, MAXN, MAXN); CHKERRQ(ierr);

  /* ierr = MatSetType(pmat, MATSEQAIJ); CHKERRQ(ierr); MATAIJ is default */
  ierr = MatSetFromOptions(pmat); CHKERRQ(ierr);  /* Allows eg. using -mat_type parameter */

  /* Better way to initialize this without creating it on the stack instead of the heap? */
  for (unsigned int j = 0; j < MAXN; j++) {
    nonzerosPerRow[j] = 0;
  }

  /* Count non-zeros per row */
  for (unsigned int j = 0; j < MAXN; j++) {
    if (CDIAG[j]) {
      nonzerosPerRow[j]++;
    }

    for (unsigned int i = CCLPTR[j]; i < CCLPTR[j+1]; i++) {
      if (COFDIAG[i]) {
        nonzerosPerRow[CRNUM[i]]++;
        nonzerosPerRow[j]++;
      }
    }
  }

  /* Preallocate matrix */
  ierr = MatSeqAIJSetPreallocation(pmat, 0, nonzerosPerRow); CHKERRQ(ierr);

  /* alternatively instead of ll. 815-817, 839:
   * ierr = MatCreateSeqAIJ(comm, MAXN, MAXN, 0, nonzerosPerRow, &pmat)
   * But: PETSc documentation says:
   * "It is recommended that one use the MatCreate(), MatSetType() and/or MatSetFromOptions(), MatXXXXSetPreallocation()
   * paradgm instead of this routine directly." */

  /* Set all entries to 0 */
  ierr = MatZeroEntries(pmat); CHKERRQ(ierr);

  /* Set the values (cached) */
  for (unsigned int j = 0; j < MAXN; j++) {
    /* Set diagonal element j (if not 0) */
    if (CDIAG[j]) {ierr = MatSetValue(pmat, j, j, CDIAG[j], INSERT_VALUES); CHKERRQ(ierr);}

    /* Add non-zero values of row j */
    for (unsigned int i = CCLPTR[j]; i < CCLPTR[j+1]; i++) {
      /*Set off-diagonal numbers (symmetric) if not 0 */
      if (COFDIAG[i]) {
        ierr = MatSetValue(pmat, CRNUM[i], j, COFDIAG[i], INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(pmat, j, CRNUM[i], COFDIAG[i], INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  /* finalize the matrix */
  ierr = MatAssemblyBegin(pmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(pmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /* Apply scale */
  MatScale(pmat, scale);

  /* Open file, use specified output format */
  switch (ot) {
    case ot_bin:
      ierr = PetscViewerBinaryOpen(comm, file, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
      break;
    case ot_ascii:
      ierr = PetscViewerASCIIOpen(comm, file, &fd); CHKERRQ(ierr);
      break;
    case ot_matlab:
      ierr = PetscViewerASCIIOpen(comm, file, &fd); CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(fd, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
      break;
    default:
      throw kaBaseException("kaMatrixSparseCol::SavePETSc: Can not save to file %s: Output format %i not implemented.",
                            file, ot);
      break;
  }

  /* Prepare dummy transition matrix */
  for (char i = 0; i < 4; i++)
    for (char j = 0; j < 4; j++)
      dummy.a(i, j) = (i == j ? 1.0 : 0.0);

  /* Write to file */
  ierr = MatView(pmat, fd); CHKERRQ(ierr);

  SaveLatInfo(file, dummy, MAXN, 1, 1);

  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

  /* Free resources */
  ierr = MatDestroy(&pmat); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}  // >::SavePETSc

#endif  // ifdef WITH_PETSC

template<class T>
void kaMatrixSparseCol<T>::Save(FILE *to) const {
  const char d[] = "kaMatrixSparseCol";

  fwrite(d, sizeof("kaMatrixSparseCol"), 1, to);
  fwrite(&MAXN, sizeof(int), 1, to);
  fwrite(&NZ, sizeof(int), 1, to);

  int i;
  for (i = 0; i < MAXN; i++) {
    T a = CDIAG[i];
    fwrite(&a, sizeof(T), 1, to);
  }
  for (i = 0; i < NZ; i++) {
    T a = COFDIAG[i];
    fwrite(&a, sizeof(T), 1, to);
  }
  for (i = 0; i < MAXN+1; i++) {
    int b = CCLPTR[i];
    fwrite(&b, sizeof(int), 1, to);
  }
  for (i = 0; i < NZ; i++) {
    int b = CRNUM[i];
    fwrite(&b, sizeof(int), 1, to);
  }
}

template<class T>
void kaMatrixSparseCol<T>::Restore(FILE *from) {
  char d[20];

  fread(d, sizeof("kaMatrixSparseCol"), 1, from);
  if (strcmp(d, "kaMatrixSparseCol"))
    throw kaBaseException("kaMatrixSparseCol<T>::Restore: not a matrix file");

  fread(&MAXN, sizeof(int), 1, from);
  fread(&NZ, sizeof(int), 1, from);

  CDIAG.resize(MAXN);
  COFDIAG.resize(NZ);
  CCLPTR.resize(MAXN+1);
  CRNUM.resize(NZ);

  int i;
  for (i = 0; i < MAXN; i++)
    fread(&CDIAG[i], sizeof(T), 1, from);
  for (i = 0; i < NZ; i++)
    fread(&COFDIAG[i], sizeof(T), 1, from);
  for (i = 0; i < MAXN+1; i++)
    fread(&CCLPTR[i], sizeof(int), 1, from);
  for (i = 0; i < NZ; i++)
    fread(&CRNUM[i], sizeof(int), 1, from);
}

#endif  // ifndef _kaMatrixSparseCol_h
