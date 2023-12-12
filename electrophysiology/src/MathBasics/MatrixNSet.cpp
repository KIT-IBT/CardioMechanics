/*
 * File: MatrixNSet.cpp
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
    fprintf(stderr, "%s <kaMatrixN> [-verbose] [-unit]\n", argv[0]);
    fprintf(stderr, "\t\t[-trans val val val] [-transX val] [-transY val] [-transZ val]\n");
    fprintf(stderr, "\t\t[-rot val val val] [-rotX val] [-rotY val] [-rotZ val] [-rotA angle x y z (axis)] \n");
    fprintf(stderr, "\t\t[-scale val val val] [-scaleX val] [-scaleY val] [-scaleZ val]\n");
    fprintf(stderr, "\t\t[-set row col val] [-add row col val]\n");
    fprintf(stderr, "\t\t[-inv]\n");
    exit(-1);
  }
  try {
    kaSharedMatrixN<double> M(4);
    FILE *fp = fopen(argv[1], "r");
    if (fp) {
      M.Restore(fp);
      fclose(fp);
    }

    bool verbose = false;
    for (int i = 2; i < argc; i++)
      if (!strcasecmp(argv[i], "-verbose")) {
        verbose = true;
      } else if (!strcasecmp(argv[i], "-unit")) {
        M.Identity();
      } else if (!strcasecmp(argv[i], "-trans") && (i+3 < argc) ) {
        double x = atof(argv[++i]);
        double y = atof(argv[++i]);
        double z = atof(argv[++i]);
        M.TranslateX(x, y, z);
        if (verbose)
          fprintf(stderr, "Trans %lf %lf %lf\n", x, y, z);
      } else if (!strcasecmp(argv[i], "-transX") && (i+1 < argc)) {
        double x = atof(argv[++i]);
        M.TranslateX(x, 0, 0);
        if (verbose)
          fprintf(stderr, "TransX %lf\n", x);
      } else if (!strcasecmp(argv[i], "-transY") && (i+1 < argc)) {
        double y = atof(argv[++i]);
        M.TranslateX(0, y, 0);
        if (verbose)
          fprintf(stderr, "TransY %lf\n", y);
      } else if (!strcasecmp(argv[i], "-transZ") && (i+1 < argc)) {
        double z = atof(argv[++i]);
        M.TranslateX(0, 0, z);
        if (verbose)
          fprintf(stderr, "TransZ %lf\n", z);
      } else if (!strcasecmp(argv[i], "-rot") && (i+3 < argc)) {
        double x = atof(argv[++i]);
        double y = atof(argv[++i]);
        double z = atof(argv[++i]);
        M.RotateX(x, 1, 0, 0);
        M.RotateX(y, 0, 1, 0);
        M.RotateX(z, 0, 0, 1);
        if (verbose)
          fprintf(stderr, "Rot %lf %lf %lf\n", x, y, z);
      } else if (!strcasecmp(argv[i], "-rotA") && (i+4 < argc)) {
        double a = atof(argv[++i]);
        double x = atof(argv[++i]);
        double y = atof(argv[++i]);
        double z = atof(argv[++i]);
        M.RotateX(a, x, y, z);
        if (verbose)
          fprintf(stderr, "RotA %lf %lf %lf %lf\n", a, x, y, z);
      } else if (!strcasecmp(argv[i], "-rotX") && (i+1 < argc)) {
        double x = atof(argv[++i]);
        M.RotateX(x, 1, 0, 0);
        if (verbose)
          fprintf(stderr, "RotX %lf\n", x);
      } else if (!strcasecmp(argv[i], "-rotY") && (i+1 < argc)) {
        double y = atof(argv[++i]);
        M.RotateX(y, 0, 1, 0);
        if (verbose)
          fprintf(stderr, "RotY %lf\n", y);
      } else if (!strcasecmp(argv[i], "-rotZ") && (i+1 < argc)) {
        double z = atof(argv[++i]);
        M.RotateX(z, 0, 0, 1);
        if (verbose)
          fprintf(stderr, "RotZ %lf\n", z);
      } else if (!strcasecmp(argv[i], "-scale") && (i+3 < argc)) {
        double x = atof(argv[++i]);
        double y = atof(argv[++i]);
        double z = atof(argv[++i]);
        M.ScaleX(x, y, z);
        if (verbose)
          fprintf(stderr, "Scale %lf %lf %lf\n", x, y, z);
      } else if (!strcasecmp(argv[i], "-scaleX") && (i+1 < argc)) {
        double x = atof(argv[++i]);
        M.ScaleX(x, 1, 1);
        if (verbose)
          fprintf(stderr, "ScaleX %lf\n", x);
      } else if (!strcasecmp(argv[i], "-scaleY") && (i+1 < argc)) {
        double y = atof(argv[++i]);
        M.ScaleX(1, y, 1);
        if (verbose)
          fprintf(stderr, "ScaleY %lf\n", y);
      } else if (!strcasecmp(argv[i], "-scaleZ") && (i+1 < argc)) {
        double z = atof(argv[++i]);
        M.ScaleX(1, 1, z);
        if (verbose)
          fprintf(stderr, "ScaleZ %lf\n", z);
      } else if (!strcasecmp(argv[i], "-inv") && (i < argc)) {
        M.Inverse();
        if (verbose)
          fprintf(stderr, "Inv\n");
      } else if ((!strcasecmp(argv[i], "-set") || !strcasecmp(argv[i], "-member")) && (i+3 < argc)) {
        int r      = atoi(argv[++i]);
        int c      = atoi(argv[++i]);
        double val = atof(argv[++i]);
        if ((r < 0) || (c < 0) || (r > 3) || (c > 3))
          throw kaBaseException("%s: element(%d,%d) not in matrix", argv[0], r, c);
        M.a(r, c) = val;
        if (verbose)
          fprintf(stderr, "member %d %d %lf\n", r, c, val);
      } else if (!strcasecmp(argv[i], "-add") && (i+3 < argc)) {
        int r      = atoi(argv[++i]);
        int c      = atoi(argv[++i]);
        double val = atof(argv[++i]);
        if ((r < 0) || (c < 0) || (r > 3) || (c > 3))
          throw kaBaseException("%s: element(%d,%d) not in matrix", argv[0], r, c);
        M.a(r, c) += val;
        if (verbose)
          fprintf(stderr, "member %d %d %lf\n", r, c, val);
      } else {
        throw kaBaseException("Unknown option %s", argv[i]);
      }

    if (verbose) {
      fprintf(stderr, "Matrix result\n");
      M.Print();
    }

    // Write matrix to file
    fp = fopen(argv[1], "w");
    if (fp) {
      M.Save(fp);
      fclose(fp);
    } else {
      throw kaBaseException("Matrix %s writing failed", argv[1]);
    }
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
