/*
 * File: FormFunctionsTest.cpp
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



#include <kaPoint.h>
#include <kaSharedMatrixN.h>
#include <Tensor2.h>
#include "FormFunctions.h"

int main(int argc, char *argv[]) {
  typedef HexaederFormFunctionTrilinear<double> F;
  HexaederFormFunctions<double, 2, F> A;

  double x[F::DOFs], y[F::DOFs], z[F::DOFs];
  double o = 0, h = (argc == 2 ? atof(argv[1]) : 1);
  double h2 = h/2;

  x[0] = o;   y[0] = o;   z[0] = o;
  x[1] = o+h; y[1] = o;   z[1] = o;
  x[2] = o+h; y[2] = o+h; z[2] = o;
  x[3] = o;   y[3] = o+h; z[3] = o;
  x[4] = o;   y[4] = o;   z[4] = o+h;
  x[5] = o+h; y[5] = o;   z[5] = o+h;
  x[6] = o+h; y[6] = o+h; z[6] = o+h;
  x[7] = o;   y[7] = o+h; z[7] = o+h;

  if (F::DOFs > 8) {
    x[8]  = o+h2;  y[8] = o;     z[8] = o;
    x[9]  = o+h;   y[9] = o+h2;  z[9] = o;
    x[10] = o+h2; y[10] = o+h;  z[10] = o;
    x[11] = o;    y[11] = o+h2; z[11] = o;
    x[12] = o;    y[12] = o;    z[12] = o+h2;
    x[13] = o+h;  y[13] = o;    z[13] = o+h2;
    x[14] = o+h;  y[14] = o+h;  z[14] = o+h2;
    x[15] = o;    y[15] = o+h;  z[15] = o+h2;
    x[16] = o+h2; y[16] = o;    z[16] = o+h;
    x[17] = o+h;  y[17] = o+h2; z[17] = o+h;
    x[18] = o+h2; y[18] = o+h;  z[18] = o+h;
    x[19] = o;    y[19] = o+h2; z[19] = o+h;
  }

  if (F::DOFs > 20) {
    x[20] = o+h2; y[20] = o+h2; z[20] = o;
    x[21] = o+h2; y[21] = o;    z[21] = o+h2;
    x[22] = o+h;  y[22] = o+h2; z[22] = o+h2;
    x[23] = o+h2; y[23] = o+h;  z[23] = o+h2;
    x[24] = o;    y[24] = o+h2; z[24] = o+h2;
    x[25] = o+h2; y[25] = o+h2; z[25] = o+h;
    x[26] = o+h2; y[26] = o+h2; z[26] = o+h2;
  }

  if (argc > 2) {
    kaSharedMatrixN<double> M(4);
    FILE *fp = fopen(argv[2], "r");
    if (fp) {
      M.Restore(fp);
      fclose(fp);
    }
    printf("Transform matrix:\n");
    M.Print();
    for (int i = 0; i < F::DOFs; i++) {
      kaPointDouble O(x[i], y[i], z[i]);
      O.XMatrixN(M);
      x[i] = O.x;
      y[i] = O.y;
      z[i] = O.z;
    }
  }
  for (int i = 0; i < F::DOFs; i++)
    printf("%d %lf %lf %lf\n", i, x[i], y[i], z[i]);

  printf("Stiffness matrix:\n");
  double stiff[F::DOFs *F::DOFs];
  for (int i = 0; i < F::DOFs; i++) {
    A.GaussianQuadraturePoisson(x, y, z, SymmetricTensor2<double>(1, 0, 0, 1, 0, 1), i, stiff+i*F::DOFs);
    for (int j = 0; j < F::DOFs; j++)
      printf("%+8.5le ", stiff[i*F::DOFs+j]);
    printf("\n");
  }
  double pot[F::DOFs];
  printf("Pot:");
  for (int i = 0; i < F::DOFs; i++) {
    //    pot[i]=(x[i]-o)/h*(y[i]-o)/h*(z[i]-o)/h;
    pot[i] = (i > 3 ? 1 : -1);
    printf("\t%lf", pot[i]);
  }
  printf("\n");

  double tmp[F::DOFs];
  for (int i = 0; i < F::DOFs; i++) {
    tmp[i] = VecMultVec(stiff+i*F::DOFs, pot, F::DOFs);
    printf("\t%lf", tmp[i]);
  }
  printf("\n");

  double p = VecMultVec(tmp, pot, F::DOFs);

  printf("p: %lf\n", 0.5*p);
  exit(0);

  SymmetricTensor2<double> *s[8];
  s[0] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[1] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[2] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[3] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[4] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[5] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[6] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);
  s[7] = new SymmetricTensor2<double>(1, 0, 0, 1, 0, 1);

  printf("N:\n");
  for (int i = 0; i < F::DOFs; i++) {
    for (int j = 0; j < F::DOFs; j++)
      printf("%+8.5le ", A.N[i][j]);
    printf("\n");
  }

  printf("Stiffness matrix:\n");
  for (int i = 0; i < F::DOFs; i++) {
    double stiff[F::DOFs];
    A.GaussianQuadraturePoisson(x, y, z, s, i, stiff);
    for (int j = 0; j < F::DOFs; j++)
      printf("%+8.5le ", stiff[j]);
    printf("\n");
  }

  for (int i = 0; i < 4; i++)
    s[i] = new SymmetricTensor2<double>(1e-9, 0, 0, 1e-9, 0, 1e-9);

  printf("Stiffness matrix:\n");
  for (int i = 0; i < F::DOFs; i++) {
    double stiff[F::DOFs];
    A.GaussianQuadraturePoisson(x, y, z, s, i, stiff);
    for (int j = 0; j < F::DOFs; j++)
      printf("%+8.5le ", stiff[j]);
    printf("\n");
  }


  //   printf("Mass element matrix:\n");
  //   for(i=0; i<8; i++) {
  //     double stiff[8], mass[8];
  //     A.GaussianQuadraturePoisson(x, y, z, i, stiff, mass);
  //     for(int j=0; j<8; j++)
  //       printf("%+lf ", mass[j]);
  //     printf("\n");
  //   }
}  // main
