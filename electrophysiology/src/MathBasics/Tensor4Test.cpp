/*
 * File: Tensor4Test.cpp
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


#include <Tensor4.h>

int main(int argc, char *argv[]) {
  kaMatrixN<double, 3> R;
  int cphi   = 0;
  int ctheta = 0;
  R.a(0) = SinTable.val[ctheta]*CosTable.val[cphi]; R.a(1) = -SinTable.val[cphi];
  R.a(2) = CosTable.val[ctheta]*CosTable.val[cphi];
  R.a(3) = SinTable.val[ctheta]*SinTable.val[cphi]; R.a(4) = CosTable.val[cphi];
  R.a(5) = CosTable.val[ctheta]*SinTable.val[cphi];
  R.a(6) = CosTable.val[ctheta];                    R.a(7) = 0.;                  R.a(8) = -SinTable.val[ctheta];


  cerr << "R:\n" << R <<  "\n";

  Tensor4<double> T, U;

  T.Null();
  int i, j, k, l;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          T[i][j][k][l] = i+j*3+k*9+l*27;

  U = T;

  // cerr << "T:\n" << T << "\n";
  // cerr << "U:\n" << U << "\n";

  //    for(i=0; i<100000; i++)
  //  T.CoordinateTransform(R);
  //  T.CoordinateTransform2(R);
  cerr << "T:\n" << T << "\n";

  //  U.CoordinateTransform(R);
  //     U.CoordinateTransform2(R);
  U.CoordinateTransform(R);
  int rc = R.Inverse();
  assert(rc);
  cerr << "R:\n" << R <<  "\n";
  U.CoordinateTransform(R);

  cerr << "U:\n" << U << "\n";

  cerr << "U-T:\n" << U-T << "\n";
}  // main
