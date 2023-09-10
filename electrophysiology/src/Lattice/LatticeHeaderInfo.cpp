/*! \file LatticeHeaderInfo.cpp
   \brief Functions for program LatticeHeaderInfo
   \date created   ?/??/?? ??\n
        man   02/26/03 mm
   \author cw,fs,idb, IBT - Universität Karlsruhe (TH)
   \sa Synopsis \ref LatticeHeaderInfo
 */

// doxygen manpage LatticeHeaderInfo
/*! \page LatticeHeaderInfo LatticeHeaderInfo
    displays the lattice header

   \section SYNOPSIS_LatticeHeaderInfo SYNOPSIS
   LatticeHeaderInfo \<kaLattice\> [-shmsize] [-width] [-height] [-depth] [-matrix] [-type] [-elemsize] [-numvox]
      [-compression]

   \section OPTIONS_LatticeHeaderInfo OPTIONS
   \param  "<kaLattice>"  name of lattice
   \param  "-shmsiz" display only size of shared memory block
   \param  "-width" display only width of data block
   \param  "-height" display only of data block
   \param  "-depth" display only of data block
   \param  "-matrix" display only ^ation matrix
   \param  "-type" display only  type of data stored in lattice
   \param  "-elemsize" display only size of one voxel in bytes
   \param  "-numvox" display only  number of voxels (data size)
   \param  "-compression" dispaly only compression type

   \section DESCRIPTION_LatticeHeaderInfo DESCRIPTION
   A lattice data structure consists of a header and a data block.
   LatticeHeaderInfo displays the header of the lattice data structure.

   The data stored in the header is described in the class kaLatticeHeader.

   The load flag . (??)

   \section SOURCE_LatticeHeaderInfo SOURCE
   LatticeHeaderInfo.cpp

   \section SEEALSO_LatticeHeaderInfo SEE ALSO
   \ref LatticeInfo \ref LatticeCreate \ref LatticeSet \ref LatticeSize
 *
 */
#include <kaLattice.h>

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
    kaLatticeSharedHeader sh(cls);

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
        M.Print("%9.6lg\t");
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
        cout << compTypeName[sh.lhCompression()] << endl;
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
