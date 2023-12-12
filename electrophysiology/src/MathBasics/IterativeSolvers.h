/*
 * File: IterativeSolvers.h
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


#ifndef _ITERATIVESOLVERS_H
#define _ITERATIVESOLVERS_H

#include <kaMatrixSparseCol.h>
#include <IBTDebug.h>
#include <limits>
#include <numeric>

#include <sys/types.h>
#include <sys/times.h>

// linux_jn
// #ifndef CLK_TCK
// #define CLK_TCK CLOCKS_PER_SEC
// #endif

template<class T>
class IterativeSolvers : public IBTDebug {
 public:
  IterativeSolvers() : tolerance(numeric_limits<T>::epsilon()*10.0), max_iter(1000), alfa(0.0005) {}

  ~IterativeSolvers() {}

  inline T Tolerance() const {return tolerance;}

  inline int MaxIter() const {return max_iter;}

  inline T Alfa() const {return alfa;}

  inline void SetTolerance(const T t) {
    tolerance = t;
  }

  void SetMaxIter(const int maxi) {
    max_iter = maxi;
  }

  void SetAlfa(const T a) {
    alfa = a;
  }

  /*! Conjugate Gradient (CG) solver;
     A-matrix; b-right vector; x-solution vector;
     x can be also an initial approximation to the solution
     pm- precond. method
     pm=0 - diagonal preconditioning
     pm=1 - improved incomplete cholesky
     pm=2 - standard no-fill incomplete cholesky
     returns 0 -tolerance satisfied; 1-maxiter reached
   */
  int CG(const kaMatrixSparseCol<T> &A, const valarray<T> &b, valarray<T> &x, const int pm = 0);

  // Biconjugate gradient solver
  // int BiCG(CMatrix &A, const std::valarray<double> &b, std::valarray<double> &x, char pm=0);

 private:
  T tolerance;
  int max_iter;
  T alfa;  // needed if IC doesn't succeed
};  // class IterativeSolvers

template<class In, class T>

inline void update_x(In first, In last, const T val) {  // x[i]=x[i]+a*x[i];
  while (first != last) *first += *first++ *val;
}

template<class In1, class In2, class T>

inline void Saypx(In1 first1, In1 last1, In2 first2, const T val) {  // y[i]=a*y[i]+x[i]
  while (first1 != last1) {*first1 = *first1 * val + *first2; first1++; first2++;}
}

template<class In1, class In2, class T>

inline void Saxpy(In1 first1, In1 last1, In2 first2, const T val) {  // y[i]=a*x[i]+y[i]
  while (first1 != last1) {(*first1) += *first2 * val; first1++; first2++;}
}

template<class T>
int IterativeSolvers<T>::CG(const kaMatrixSparseCol<T> &A, const valarray<T> &b, valarray<T> &x, const int pm) {
  struct tms tmsbuf;
  double usrtime;

  int retcode;  // return code
  int MAXN = A.Dimension();

  kaMatrixSparseCol<T> *CH = NULL;

  // incomplete cholesky
  if ((1 == pm) || (2 == pm)) {
    CH = new kaMatrixSparseCol<T>(A);
    double factor = 1.0;
    double a      = 0.0;

    int retcode = -1;


    times(&tmsbuf);
    usrtime = -(double)tmsbuf.tms_utime/(double)CLK_TCK;

    if (1 == pm)
      retcode = CH->jpicc(); // improved incomplete chol
    else if (2 == pm)
      retcode = CH->istdic(); // standard incomplete chol

    while (retcode) {
      info_stream(2)<<"CG:Cholesky failed for alfa="<<a<<", return code:"<<retcode<<'\n';
      factor *= 2.0;
      a       = alfa*factor;

      // restore CH and try again
      *CH        = A;
      CH->CDIAG += a;

      if (1 == pm)
        retcode = CH->jpicc(); // improved incomplete chol
      else if (2 == pm)
        retcode = CH->istdic(); // standard incomplete chol
    }


    info_stream(2)<<"CG: Cholesky OK, alfa="<<a<<'\n';

    times(&tmsbuf);
    usrtime += (double)tmsbuf.tms_utime/(double)CLK_TCK;
    cout<<"Time for incomplete Cholesky [sec]:"<<setprecision(2)<<fixed<<usrtime<<'\n';
  }


  valarray<T> res(0.0, MAXN);
  valarray<T> p(0.0, MAXN);
  valarray<T> q(0.0, MAXN);
  valarray<T> z(0.0, MAXN);

  T ro0  = 0;
  T ro1  = 0;
  T beta = 0;

  A.SymmMatrixVectorProduct(x, p);  // p is used to save memory

  res = b - p;

  p = b;  // copy of b to get access through pointers; p will be initialized later anyway
  T norm_b = sqrt(inner_product(&p[0], &p[MAXN], &p[0], 0.0) );

  T rlimit = tolerance*norm_b;

  if (norm_b == 0.0) {
    info_stream(2) << " CG: Norm of B is zero! Solution is set to zero.";
    x = 0.0;
    return 0;
  }


  int iter = 1;

  for (;;) {
    if (pm == -1)
      z = res;
    else if (pm == 0)
      z = res/A.CDIAG;
    else
      CH->chol_solve(res, z);

    ro0 = ro1;
    ro1 = 0;

    // for(int i=0; i<MAXN; i++)  ro1=ro1+res[i]*z[i];
    ro1 = inner_product(&res[0], &res[MAXN], &z[0], 0.0);

    if (iter == 1) {
      copy(&z[0], &z[MAXN], &p[0]);  // for(int i=0; i<MAXN; i++) p[i]=z[i];
    } else {
      beta = ro1/ro0;
      Saypx(&p[0], &p[MAXN], &z[0], beta);

      // for(int i=0; i<MAXN; i++) p[i]=z[i]+beta*p[i];
    }


    A.SymmMatrixVectorProduct(p, q);
    double denom = 0;

    // for(int i=0; i<MAXN; i++) denom+=p[i]*q[i];
    denom = inner_product(&p[0], &p[MAXN], &q[0], 0.0);
    double aa = ro1/denom;

    Saxpy(&x[0], &x[MAXN], &p[0], aa);
    Saxpy(&res[0], &res[MAXN], &q[0], -aa);

    double norm_r = 0;
    norm_r = inner_product(&res[0], &res[MAXN], &res[0], 0.0);
    norm_r = sqrt(norm_r);


    info_stream(2)<<setprecision(8)<<scientific<<" CG: Iteration = "<< iter <<" ||res|| = "<<norm_r;
    info_stream(2)<<setprecision(8)<<scientific<<" rlimit = "<<rlimit<<'\n';

    if (norm_r < rlimit) {
      retcode = 0;
      break;
    }

    if (iter == max_iter) {
      retcode = 1;
      break;
    }
    iter++;
  }

  if (CH)
    delete CH;
  return retcode;
}  // >::CG

#endif  // ifndef _ITERATIVESOLVERS_H
