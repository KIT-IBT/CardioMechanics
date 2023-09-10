/*------------------------------------------------------------------------------
   --   Name
   --           kaSharedMatrixNOpkaSharedMatrixN.cpp
   --
   --   Remark
   --           Bestimmung einer Transformationsmatrix zwischen zwei homogenen
   --           Tranformationsmatrizen.
   --           Da die affinen Abbildungen in beiden Faellen durch eine
   --           homogene Transformation beschrieben werden, sollte es reichen
   --           die Abbildungsmatrix des Templates zu invertieren und von rechts
   --           an die Abbildungsmatrix der Referenz zu multiplizieren, um die
   --           direkte Transformation zwischen beiden Koordinatensystemen zu erhalten.
   --           In Weltkoordinaten gilt:
   --                   P1(Welt|Template)=A(Template)*P1
   --                   P2(Welt|Referenz)=A(Referenz)*P2
   --                   P1(Welt|Referenz)=A(Referenz)*A(Template)^-1*A(Template)*P1
   --           Als Transformationsmatrix T ergibt sich damit:
   --                   T = A(Referenz)*A(Template)^-1
   --           Als Argumente bei Programmaufruf werden in Files
   --           gespeicherte Matrizen eingelesen.
   --
   --
   --
   --   History
   --           Erstellung      01.04.98        fs
   --           Erweiterung     22.07.98        gvw
   --           last modified   02.02.99        gvw
   --
   --   created at IBT - Universitaet Karlsruhe
   ------------------------------------------------------------------------------*/

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
