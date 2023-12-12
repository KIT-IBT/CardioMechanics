/*
 * File: MatrixNOpMatrixN.cpp
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


#include <kaSharedMatrixN.h>
#include <kaExceptions.h>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "%s <Matrix1> <Matrix2> <MatrixResult>\n", argv[0]);
    fprintf(stderr, "\t\t[-inv1] [-inv2] [-mult] [-inv]\n");
    exit(-1);
  }
  try {
    // Lese Referenz-Transformationsmatrix
    kaSharedMatrixN<double> TransformMatrix1(4);
    FILE *fp = fopen(argv[1], "r");
    TransformMatrix1.Restore(fp);
    fclose(fp);

    // Lese Modell-Transformationsmatrix
    kaSharedMatrixN<double> TransformMatrix2(4);
    fp = fopen(argv[2], "r");
    TransformMatrix2.Restore(fp);
    fclose(fp);

    kaSharedMatrixN<double> TransformMatrix3(4);
    bool verbose = false;
    for (int i = 4; i < argc; i++)
      if (!strcmp(argv[i], "-verbose")) {
        verbose = true;
      } else if (!strcmp(argv[i], "-inv1")) {
        TransformMatrix1.Inverse();
        if (verbose)
          fprintf(stderr, "%s\n", argv[i]);
      } else if (!strcmp(argv[i], "-inv2")) {
        TransformMatrix2.Inverse();
        if (verbose)
          fprintf(stderr, "%s\n", argv[i]);
      } else if (!strcmp(argv[i], "-mult")) {
        TransformMatrix3 = TransformMatrix1*TransformMatrix2;
        if (verbose)
          fprintf(stderr, "%s\n", argv[i]);
      } else if (!strcmp(argv[i], "-inv")) {
        TransformMatrix3.Inverse();
        if (verbose)
          fprintf(stderr, "%s\n", argv[i]);
      } else {
        throw kaBaseException("Unknown option %s", argv[i]);
      }

    if (verbose) {
      fprintf(stderr, "Matrix result\n");
      TransformMatrix3.Print();
    }

    // Ergebnismatrix in File schreiben
    fp = fopen(argv[3], "w");
    TransformMatrix3.Save(fp);
    fclose(fp);
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
