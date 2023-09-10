/**@file MatrixTools.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#ifndef _MATRIXTOOLS_H
#define _MATRIXTOOLS_H

#include <valarray>
#include <set>


#include <kaMachineOS.h>
#include <kaExceptions.h>

using namespace std;

class mt_excp : public kaBaseException {
 public:
  mt_excp(const char *msg) : kaBaseException("exception in MatrixTools: %s", msg) {}
};


struct CMatrixStructure {
  int nz;  // number of non - zeroes
  int MAXN;
  vector<set<int>> mstruct;
  void add(const int &, const int &);  // reserve a place in the matrix
  int  size() const; // size in bytes
  CMatrixStructure(const int &);  // MAXN is supplied
  void       Print(ostream &);
  set<int> & operator[](const int &i) {return mstruct[i];}
};

class CMatrix {
 public:
  /* STORAGE FOR THE MATRIX IN COLUMN MAJOR ORDER
     (DIAGONALS, OFF-DIAGONALS, COLUMN PTR'S, ROW NUMBERS) */
  int MAXN;  // max number of rows and columns
  int NZ;                 // number of off-diagonal non-zeroes
  valarray<double> CDIAG;  // MAXN
  valarray<double> COFDIAG;  // NZ
  valarray<int> CCLPTR;  // MAXN+1
  valarray<int> CRNUM;  // NZ
  /*you should supply at least MAXN and NZ
     to create an empty CMatrix which can be filled
     by the user in any convinient way*/
  CMatrix(const int &, const int &);
  /*you can also supply a CMatrixStructure
     This results into the creating of an empty
     CDIAG and COFDIAG, but the CCLPTR and CRNUM are
     filled basing on the CMatrixStructure */
  CMatrix(CMatrixStructure &);
  CMatrix(const CMatrix &a);  // copy of a
  ~CMatrix() {}

  CMatrix     & operator=(const CMatrix &);
  inline void   add(const int &, const int &, const double &);
  inline void   set(const int &, const int &, const double &);
  inline double zero(const int &, const int &);  // set to zero and return prev value

  inline int Dimension() const {return MAXN;}

  inline int NumRows() const {return MAXN;}

  inline int NumCols() const {return MAXN;}

  inline double   Get(const int &, const int &) const;
  inline double & a(const int &, const int &);

  inline int size() {return sizeof(double)*(CDIAG.size()+COFDIAG.size())+sizeof(int)*(CCLPTR.size()+CRNUM.size());}

  int save(const char *name);
  int read(const char *name);

  // Description:
  // Save in ASCII format, which can be read by MATLAB
  // Indices start with 1.
  // Set symmetric to 1 to save the full matrix except of only lower (upper) triangular part.
  int save_ascii(const char *name, const int symmetric = 0);

  //! use  SymmMatrixVectorProduct if only the lower (upper) triangle of A is stored
  inline void MatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const;
  inline void TranspMatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const;

  //! use  SymmMatrixVectorProduct if only the lower (upper) triangle of A is stored
  inline void   SymmMatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const;
  inline double FrobeniusNorm(const int symmetric = 0) const;
};  // class CMatrix


/*** Some operators ***/
extern ostream & operator<<(ostream &, CMatrix &);

/*** Some Useful Functions ***/
extern  void SetColToZero(CMatrix &, const int &);
extern  void SetRowToZero(CMatrix &, const int &);
/*! Scale equation sysytem to get ones on the diagonal of the
   system matrix. B is the right-hand vector.
 */
extern void Scale(CMatrix &A, valarray<double> &B);
/*! Scale equation sysytem to get ones on the diagonal of the
   system matrix. B is the right-hand vector.
   A is a symmetric matrix (only upper or lower part is stored)
   in scale the scaling factors will be written to compute the solution
   Equation: Ax=b
   Scaled inv(S)Ainv(S)Sx=inv(S)b => we get Sx => must do inv(S)Sx
 */
extern void SymmScale(CMatrix &A, valarray<double> &B, valarray<double> &scale);
/*
   Jones/Plassman incomplete Cholesky: column oriented
   ACM Transactions on Mathematical Software, Vol 21, No 1, March 1995, pp. 5-17
   FORTRAN version is available via netlib

   Important: ONLY A STRICTLY LOWER TRIANGLE OF A IS NEEDED!
 */
extern int JPICC(CMatrix &A);
/* standard incomplete Cholesky: column oriented
   FORTRAN version is available via netlib together with the JPICC */
extern int ISTDIC(CMatrix &A);
/* Cholesky-type factorization solve. A-lower triangular matrix
   the procedure solves A.Transpose[A].x=y   to find x
   Parameters : A, y, x */
extern void CholFactSolve(const CMatrix &A, const valarray<double> &y, valarray<double> &x);
/*  ibsort
   get the k lagest nonzeroes in akeys indirectly addressed by indvec
   upon exit the first k elements in indvec will contain the indices of
   the k largest elements in akeys
   n is the length of the index vector;

   performance test(compare to isort): 1000000 full sortings of array of 10 elements, SGI workstation: 0:57.42sec
   (includes re-initialisations of arrays)
   after compiling with -O3 "ibsort" works ~10% faster as "isort" (full sort);
   partial (50 of 100) sort with ibsort works 30% slower as isort .....
   further tests show approx. equal performance of both procedures

 */
extern void ibsort(const int &n, const int &k, const valarray<double> &akeys, valarray<int> &indvec);
/*  ihsort
   sorts a double precision vector indirectly addressed by an integer vector
   (uses heap sort)
   the indvec is rearranged such that indvec[0] addresses the largest
   element in akeys, indvec[2] addresses the next largest ...
 */
extern void ihsort(const int &len, valarray<int> &indvec, const valarray<double> &akeys);
/*  isort - analog of ibsort,
    implementation based on the c++ standard library

    performance test(compare to ibsort): 1000000 full sortings of
    array of 10 elements, SGI workstation: 0:47.53sec(includes re-initialisations of arrays)
    after compiling with -O3 "ibsort" works ~10% as "isort"
    partial (50 of 100) sort for ibsort works 30% slower as isort ....
    further tests show approx. equal performance of both procedures
 */
extern void isort(const int &n, const int &k, valarray<double> &akeys, valarray<int> &indvec);


// inline function implementation

void CMatrix::add(const int &i, const int &j, const double &value) {
  if (j > MAXN-1)
    throw mt_excp("CMatrix::add: index j greater than MAXN-1");
  if (i != j) {
    int nb = CCLPTR[j];
    int ne = CCLPTR[j+1];
    if (ne == nb)
      throw mt_excp("CMatrix::set cannot set value");
    while (ne-nb > 1) {  // bisection to find i
      int pos = (nb+ne)/2;
      if (CRNUM[pos] > i)
        ne = pos; else
        nb = pos;
    }
    if (CRNUM[nb] != i)
      throw mt_excp("CMatrix::add cannot add value");
    COFDIAG[nb] += value;
  } else {
    CDIAG[j] += value;
  }
}

double CMatrix::zero(const int &i, const int &j) {
  double retval;

  if (j > MAXN-1)
    throw mt_excp("CMatrix::zero: index j greater than MAXN-1");
  if (i != j) {
    int nb = CCLPTR[j];
    int ne = CCLPTR[j+1];
    if (ne == nb)
      return 0;

    while (ne-nb > 1) {  // bisection to find i
      int pos = (nb+ne)/2;
      if (CRNUM[pos] > i)
        ne = pos; else
        nb = pos;
    }
    if (CRNUM[nb] != i)
      return 0;

    retval      = COFDIAG[nb];
    COFDIAG[nb] = 0;
  } else {
    retval   = CDIAG[j];
    CDIAG[j] = 0;
  }
  return retval;
}  // CMatrix::zero

void CMatrix::set(const int &i, const int &j, const double &value) {
  if (j > MAXN-1)
    throw mt_excp("CMatrix::set: index j greater than MAXN-1");
  if (i != j) {
    int nb = CCLPTR[j];
    int ne = CCLPTR[j+1];
    if (ne == nb)
      throw mt_excp("CMatrix::set cannot set value");
    while (ne-nb > 1) {  // bisection to find i
      int pos = (nb+ne)/2;
      if (pos == NZ)
        throw mt_excp("Uoooi!!!");
      if (CRNUM[pos] > i)
        ne = pos; else
        nb = pos;
    }
    if (CRNUM[nb] != i) {
      cout<<"CMatrix:: ("<<i<<','<<j<<")\n";
      throw mt_excp("CMatrix::set cannot set value");
    }
    COFDIAG[nb] = value;
  } else {
    CDIAG[j] = value;
  }
}

double CMatrix::Get(const int &i, const int &j) const {
  if (j > MAXN-1)
    throw mt_excp("CMatrix::at: index j greater than MAXN-1");
  if (i != j) {
    int nb = CCLPTR[j];
    int ne = CCLPTR[j+1];
    if (ne == nb)
      return 0.0;

    while (ne-nb > 1) {  // bisection to find i
      int pos = (nb+ne)/2;
      if (CRNUM[pos] > i)
        ne = pos; else
        nb = pos;
    }
    if (CRNUM[nb] != i)
      return 0.0;

    return COFDIAG[nb];
  } else {
    return CDIAG[j];
  }
}

double & CMatrix::a(const int &i, const int &j) {
  if (j > MAXN-1)
    throw mt_excp("CMatrix::at: index j greater than MAXN-1");
  if (i != j) {
    int nb = CCLPTR[j];
    int ne = CCLPTR[j+1];
    if (ne == nb)
      throw mt_excp("the element does not exist");
    while (ne-nb > 1) {  // bisection to find i
      int pos = (nb+ne)/2;
      if (CRNUM[pos] > i)
        ne = pos; else
        nb = pos;
    }
    if (CRNUM[nb] != i)
      throw mt_excp("the element does not exist");
    return COFDIAG[nb];
  } else {
    return CDIAG[j];
  }
}

void CMatrix::MatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const {
  if (x.size() != MAXN)
    throw kaBaseException("MatrixVectorProduct: size of input vector != number of columns\n");
  r = 0.0;
  int j;
  double xj;

#pragma omp parallel for private(j, xj)
  for (j = 0; j < MAXN; j++) {
    xj    = x[j];
    r[j] += CDIAG[j]*xj;
    int clptrbeg = CCLPTR[j];
    int clptrend = CCLPTR[j+1];
    int i, rnum;
    double a;
    for (i = clptrbeg; i < clptrend; i++) {  // no check on the row length is performed!
      rnum = CRNUM[i];
      a    = COFDIAG[i]*xj;
#pragma omp atomic
      r[rnum] += a;
    }
  }
}

void CMatrix::TranspMatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const {
  if (r.size() != MAXN)
    throw kaBaseException("TranspMatrixVectorProduct: size of output vector != number of columns\n!");
  r = 0.0;

  int j;
  double rj;
#pragma omp parallel for private(j, rj)
  for (j = 0; j < MAXN; j++) {
    rj = CDIAG[j]*x[j];
    for (int i = CCLPTR[j]; i < CCLPTR[j+1]; i++)
      rj += COFDIAG[i]*x[CRNUM[i]];
    r[j] += rj;
  }
}

void CMatrix::SymmMatrixVectorProduct(const valarray<double> &x, valarray<double> &r) const {
  r = 0.0;

  int j1;
  double xj1;
#pragma omp parallel for private(xj1)
  for (j1 = 0; j1 < MAXN; j1++) {
    xj1 = x[j1];
    int clptrbeg = CCLPTR[j1];
    int clptrend = CCLPTR[j1+1];
    int i, rnum;
    double a;
    for (i = clptrbeg; i < clptrend; i++) {  // no check on the row length is performed!
      rnum = CRNUM[i];
      a    = COFDIAG[i]*xj1;
#pragma omp atomic
      r[rnum] += a;
    }
  }

  int j2;
  double rj2;
#pragma omp parallel for private(rj2)
  for (j2 = 0; j2 < MAXN; j2++) {
    rj2 = CDIAG[j2]*x[j2];
    for (int i = CCLPTR[j2]; i < CCLPTR[j2+1]; i++)
      rj2 += COFDIAG[i]*x[CRNUM[i]]; // r[j]
    r[j2] += rj2;
  }
}  // CMatrix::SymmMatrixVectorProduct

inline double CMatrix::FrobeniusNorm(const int symmetric) const {
  double norm = 0, offdiagnorm = 0;
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

#endif  // ifndef _MATRIXTOOLS_H
