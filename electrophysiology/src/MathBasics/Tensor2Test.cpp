/*! \file Tensor2Test.cpp
   \brief Function for testing class tensor2

   \author fs, CVRTI - University of Utah
 */

#include <Tensor2.h>

int main(int argc, char *argv[]) {
  Tensor2<double> T, R, UT;
  SymmetricTensor2<double> U;
  T.Null(); R.Null(); U.Null();

  R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
  R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
  R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;

  U.a = 2; U.b = 1.4; U.c = 1.3; U.d = 2; U.e = 1.1; U.f = 1;
  UT  = U;
  T   = R*UT;
  cerr << "T:\n" << T << "\n";

  cerr << "R:\n" << R << "\n";
  cerr << "U:\n" << U << "\n";
  R.Null(); U.Null();
  for (int i = 0; i < 1000000; i++)
    T.PolarDecomposition(R, U);
  cerr << "R:\n" << R << "\n";
  cerr << "U:\n" << U << "\n";
  UT = U;
  T  = R*UT;
  cerr << "T:\n" << T << "\n";
}
