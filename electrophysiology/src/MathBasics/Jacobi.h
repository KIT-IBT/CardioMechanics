/**@file Jacobi.h
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

#ifndef JACOBI_H
#define JACOBI_H

#include <kaMachineOS.h>

#define ROTATE(a, i, j, k, l) {g = a[i][j]; h = a[k][l]; a[i][j] = g-s*(h+g*tau); a[k][l] = h+s*(g-h*tau);}

template<class T, const int n>

bool jacobi(T a[n][n], T d[n], T v[n][n]) {
  int j, iq, ip, i, nrot = 0;
  double tresh, theta, tau, t, sm, s, h, g, c;

  double b[n], z[n];

  for (ip = 0; ip < n; ip++) {
    for (iq = 0; iq < n; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
    b[ip]     = d[ip] = a[ip][ip];
    z[ip]     = 0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0;
    for (ip = 0; ip < n-1; ip++)
      for (iq = ip+1; iq < n; iq++)
        sm += fabs(a[ip][iq]);
    if (sm == 0)
      return true;

    if (i < 4)
      tresh = .2*sm/(n*n);
    else
      tresh = 0;
    for (ip = 0; ip < n-1; ip++) {
      for (iq = ip+1; iq < n; iq++) {
        g = 100.0 * fabs(a[ip][iq]);
        if ((i > 4) && ( (T)(fabs(d[ip])+g) == (T)fabs(d[ip])) && ( (T)(fabs(d[iq])+g) == (T)fabs(d[iq])) ) {
          a[ip][iq] = 0;
        } else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
          if ((T)(fabs(h)+g) == (T)fabs(h)) {
            t = a[ip][iq] / h;
          } else {
            theta = .5*h/a[ip][iq];
            t     = 1.0 / (fabs(theta)+sqrt(1+theta*theta));
            if (theta < 0)
              t = -t;
          }
          c         = 1 / sqrt(1+t*t);
          s         = t * c;
          tau       = s / (1+c);
          h         = t * a[ip][iq];
          z[ip]    -= h; z[iq] += h; d[ip] -= h; d[iq] += h;
          a[ip][iq] = 0;
          for (j = 0; j <= ip-1; j++)
            ROTATE(a, j, ip, j, iq);
          for (j = ip+1; j <= iq-1; j++)
            ROTATE(a, ip, j, j, iq);
          for (j = iq+1; j < n; j++)
            ROTATE(a, ip, j, iq, j);
          for (j = 0; j < n; j++)
            ROTATE(v, j, ip, j, iq);
          nrot++;
        }
      }
    }
    for (ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      d[ip]  = b[ip];
      z[ip]  = 0;
    }
  }
  return false;
}  // jacobi

#endif  // ifndef JACOBI_H
