/*
 * File: ExpFit.h
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


#ifndef EXPFIT_H
#define EXPFIT_H


#include <kaExceptions.h>
#include <nrutil.h>

template<class T>

void covsrt(T **covar, int ma, int ia[], int mfit) {
  int i, j, k;
  T   swap;

  for (i = mfit+1; i <= ma; i++)
    for (j = 1; j <= i; j++)
      covar[i][j] = covar[j][i] = 0.0;

  k = mfit;
  for (j = ma; j >= 1; j--)
    if (ia[j]) {
      for (i = 1; i <= ma; i++) {
        swap        = covar[i][k];
        covar[i][k] = covar[i][j];
        covar[i][j] = swap;
      }
      for (i = 1; i <= ma; i++) {
        swap        = covar[k][i];
        covar[k][i] = covar[i][j];
        covar[i][j] = swap;
      }
      k--;
    }
}

template<class T>

void gaussj(T **a, int n, T **b, int m) {
  int *indxc, *indxr, *ipiv;
  int  i, icol = 1, irow = 1, j, k, l, ll;
  T big, dum, pivinv, swap;

  indxc = ivector(1, n);
  indxr = ivector(1, n);
  ipiv  = ivector(1, n);
  for (j = 1; j <= n; j++)
    ipiv[j] = 0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if (ipiv[j] != 1)
        for (k = 1; k <= n; k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
              big  = fabs(a[j][k]);
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1) {
            throw kaBaseException("gaussj: Singular Matrix-1");
          }
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l = 1; l <= n; l++) {
        swap       = a[irow][l];
        a[irow][l] = a[icol][l];
        a[icol][l] = swap;
      }
      for (l = 1; l <= m; l++) {
        swap       = b[irow][l];
        b[irow][l] = b[icol][l];
        b[icol][l] = swap;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0)
      throw kaBaseException("gaussj: Singular Matrix-2");

    pivinv        = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 1; l <= n; l++)
      a[icol][l] *= pivinv;
    for (l = 1; l <= m; l++)
      b[icol][l] *= pivinv;
    for (ll = 1; ll <= n; ll++)
      if (ll != icol) {
        dum         = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l = 1; l <= n; l++)
          a[ll][l] -= a[icol][l]*dum;
        for (l = 1; l <= m; l++)
          b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l = n; l >= 1; l--) {
    if (indxr[l] != indxc[l])
      for (k = 1; k <= n; k++) {
        swap           = a[k][indxr[l]];
        a[k][indxr[l]] = a[k][indxc[l]];
        a[k][indxc[l]] = swap;
      }
  }
  free_ivector(ipiv, 1, n);
  free_ivector(indxr, 1, n);
  free_ivector(indxc, 1, n);
}  // gaussj

template<class T>

void mrqcof(const T x[], const T y[], const T sig[], const int ndata, T a[], int ia[], const int ma, T **alpha,
            T beta[], int funcNum, T *chisq) {
  int i, j, k, l, m, mfit = 0;
  T   ymod, wt, sig2i, dy, *dyda;

  dyda = vecto(1, ma);
  for (j = 1; j <= ma; j++)
    if (ia[j])
      mfit++;
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= j; k++)
      alpha[j][k] = 0.0;
    beta[j] = 0.0;
  }
  *chisq = 0.0;
  double tmp;
  for (i = 1; i <= ndata; i++) {
    switch (funcNum) {
      case 0:

        // Boltzmann equation
        tmp     = exp((x[i]-a[2])/a[3]);
        ymod    = a[1]/(1.+tmp);
        dyda[1] = 1./(1.+tmp);
        dyda[2] = a[1]/((1.+tmp)*(1.+tmp))/a[3]*tmp;
        dyda[3] = a[1]/((1.+tmp)*(1.+tmp))*(x[i]-a[2])/a[3]/a[3]*tmp;
        break;

      case 1:

        // Exponential equation
        tmp     = exp(-a[2]*(x[i]-a[3]));
        ymod    = a[1]*(1.-tmp)+a[4];
        dyda[1] = 1.-tmp;
        dyda[2] = -a[1]*(x[i]-a[3])*tmp;
        dyda[3] = -a[1]*a[2]*tmp;
        dyda[4] = 1;
        break;

      case 2:

        // Hill equation
        ymod    = a[3]*pow(x[i], a[1]) / (pow(x[i], a[1])+pow(a[2], a[1]));
        dyda[1] = a[3]*pow(a[2], a[1])*pow(x[i], a[1])*(log(x[i])-log(a[2]))/pow((pow(x[i], a[1])+pow(a[2], a[1])), 2);
        dyda[2] = -a[3]*pow(x[i], a[1])*a[1]*pow(a[2], a[1]-1)/pow((pow(x[i], a[1])+pow(a[2], a[1])), 2);
        dyda[3] = pow(x[i], a[1]) / (pow(x[i], a[1])+pow(a[2], a[1]));
        break;

      case 3:

        // Hill equation with offset
        ymod    = a[3]*pow(x[i], a[1]) / (pow(x[i], a[1])+pow(a[2], a[1]))+a[4];
        dyda[1] = a[3]*pow(a[2], a[1])*pow(x[i], a[1])*(log(x[i])-log(a[2]))/pow((pow(x[i], a[1])+pow(a[2], a[1])), 2);
        dyda[2] = -a[3]*pow(x[i], a[1])*a[1]*pow(a[2], a[1]-1)/pow((pow(x[i], a[1])+pow(a[2], a[1])), 2);
        dyda[3] = pow(x[i], a[1]) / (pow(x[i], a[1])+pow(a[2], a[1]));
        dyda[4] = 1;
        break;

      default:
        throw kaBaseException("mrqcof: no usable function number.");
    }  // switch

    sig2i = 1.0/(sig[i]*sig[i]);
    dy    = y[i]-ymod;
    for (j = 0, l = 1; l <= ma; l++) {
      if (ia[l]) {
        wt = dyda[l]*sig2i;
        for (j++, k = 0, m = 1; m <= l; m++)
          if (ia[m])
            alpha[j][++k] += wt*dyda[m];
        beta[j] += dy*wt;
      }
    }
    *chisq += dy*dy*sig2i;
  }
  for (j = 2; j <= mfit; j++)
    for (k = 1; k < j; k++)
      alpha[k][j] = alpha[j][k];
  free_vector(dyda, 1, ma);
}  // mrqcof

template<class T>

void mrqmin(T x[], T y[], T sig[], int ndata, T a[], int ia[], int ma, T **covar, T **alpha, T *chisq, int funcNum,
            T *alamda) {
  int j, k, l;
  static int mfit;
  static T   ochisq, *atry, *beta, *da, **oneda;

  if (*alamda < 0.0) {
    atry = vecto(1, ma);
    beta = vecto(1, ma);
    da   = vecto(1, ma);
    for (mfit = 0, j = 1; j <= ma; j++)
      if (ia[j])
        mfit++;
    oneda   = matrix(1, mfit, 1, 1);
    *alamda = 0.001;
    mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, funcNum, chisq);
    ochisq = (*chisq);
    for (j = 1; j <= ma; j++)
      atry[j] = a[j];
  }
  for (j = 1; j <= mfit; j++) {
    for (k = 1; k <= mfit; k++)
      covar[j][k] = alpha[j][k];
    covar[j][j] = alpha[j][j]*(1.0+(*alamda));
    oneda[j][1] = beta[j];
  }
  gaussj(covar, mfit, oneda, 1);
  for (j = 1; j <= mfit; j++)
    da[j] = oneda[j][1];
  if (*alamda == 0.0) {
    covsrt(covar, ma, ia, mfit);
    free_matrix(oneda, 1, mfit, 1, 1);
    free_vector(da, 1, ma);
    free_vector(beta, 1, ma);
    free_vector(atry, 1, ma);
    return;
  }
  for (j = 0, l = 1; l <= ma; l++)
    if (ia[l])
      atry[l] = a[l]+da[++j];
  mrqcof(x, y, sig, ndata, atry, ia, ma, covar, da, funcNum, chisq);
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq   = (*chisq);
    for (j = 1; j <= mfit; j++) {
      for (k = 1; k <= mfit; k++)
        alpha[j][k] = covar[j][k];
      beta[j] = da[j];
    }
    for (l = 1; l <= ma; l++)
      a[l] = atry[l];
  } else {
    *alamda *= 10.0;
    *chisq   = ochisq;
  }
}  // mrqmin

#endif  // ifndef EXPFIT_H
