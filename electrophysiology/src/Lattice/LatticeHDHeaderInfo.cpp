/*! \file LatticeHeaderInfo.cpp
   \brief Functions for program LatticeHeaderInfo

   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
 */

#include <kaLatticeHeader.h>

//! Main function

/*! The function handles the printing of a help message, the parsing of command line arguments, the printing of lattice
   headers and provides error handler via exception catching
 */
int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << argv[0] << " <lattice>" << endl;
    cerr << "\t[-shmsize]" << endl;
    cerr << "\t[-width]" << endl;
    cerr << "\t[-height]" << endl;
    cerr << "\t[-depth]" << endl;
    cerr << "\t[-matrix]" << endl;
    cerr << "\t[-type]" << endl;
    cerr << "\t[-elemsize]" << endl;
    cerr << "\t[-numvox]" << endl;
    cerr << "\t[-compression]" << endl;
    return -1;
  }

  try {
    kaLatticeCreate cls;
    strcpy(cls.Name, argv[1]);
    kaLatticeHeader sh(cls);
    sh.headerRead();

    bool all = true;
    for (int a = 2; a < argc; a++) {
      if (!strcmp(argv[a], "-shmsize")) {
        all = false; cout << sh.lhSMSize() << endl;
      } else if (!strcmp(argv[a], "-width")) {
        all = false; cout << sh.lhWidth() << endl;
      } else if (!strcmp(argv[a], "-height")) {
        all = false; cout << sh.lhHeight() << endl;
      } else if (!strcmp(argv[a], "-depth")) {
        all = false; cout << sh.lhDepth() << endl;
      } else if (!strcmp(argv[a],
                         "-matrix")) {
        cout.setf(ios::fixed, ios::floatfield);
        all = false;
        kaSharedMatrixN<double> &M = sh.lhMatrix();
        cout << M;
      } else if (!strcmp(argv[a], "-type")) {
        all = false;
        DataTypeInfo DTI;
        cout << DTI.TypeName(sh.lhDtype()) << endl;
      } else if (!strcmp(argv[a], "-elemsize")) {
        all = false; cout << sh.lhElemSize() << endl;
      } else if (!strcmp(argv[a], "-numvox")) {
        all = false; cout << sh.lhDsize() << endl;
      } else if (!strcmp(argv[a],
                         "-compression")) {
        all = false;
        cout << compTypeName[sh.lhDtype()] << endl;
      } else {
        throw kaBaseException("Unknown option %s", argv[a]);
      }
    }

    if (all)
      sh.FileHeaderOut(stdout);
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
}  // main
