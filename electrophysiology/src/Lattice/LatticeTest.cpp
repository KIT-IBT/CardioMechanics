/*
 * File: LatticeTest.cpp
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



#include <kaLattice.h>


//! Template function for testing a lattice file and illustrating the access
/*!
   The arguments includes the lattice name.
   The lattices can be of different type.
 */

template<class ValTyp>

void kaLatticeTest(int argc, char *argv[]) {
  kaLatticeCreate cls;

  strcpy(cls.Name, argv[1]);
  cls.xLattice = 16;
  cls.yLattice = 16;
  cls.zLattice = 16;

  kaLattice<ValTyp> l(cls);

  printf("dim %d %d %d\n", l.xLattice, l.yLattice, l.zLattice);

  printf("\nMatrix:\n");
  l.m.Print();

  int istop = (l.xLattice < 16) ? l.xLattice : 16;
  int jstop = (l.yLattice < 16) ? l.yLattice : 16;
  int i, j;

  printf("First %d voxels:\n", istop*jstop);
  for (i = 0; i < istop; i++) {
    printf("\t");
    for (j = 0; j < jstop; j++)
      printf("%4.4lf ", (double)l.lat[i*30+j]);
    printf("\n");
  }

  IndexType steps  = 4;
  IndexType xsteps = max(1, l.xLattice / steps);
  IndexType ysteps = max(1, l.yLattice / steps);
  IndexType zsteps = max(1, l.zLattice / steps);
  double vsteps    = 256.0 / (double)(steps*steps*steps);
  double val       = 1;

  IndexType x, y, z;
  for (z = 0; z < l.zLattice; z += zsteps)
    for (y = 0; y < l.yLattice; y += ysteps)
      for (x = 0; x < l.xLattice; x += xsteps) {
        cerr << "Voxel(" << x << ", " << y << ", " << z << ") = " << l.val(x, y, z) << endl;
        l.val(x, y, z) = (ValTyp)val;
        cerr << "Voxel(" << x << ", " << y << ", " << z << ") = " << l.val(x, y, z) << endl;
        val += vsteps;
      }
}  // kaLatticeTest

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "%s Lattice\n", argv[0]);
    exit(-1);
  }

  try {
    kaLatticeCreate cls;
    strcpy(cls.Name, argv[1]);
    kaLatticeSharedHeader sh(cls);
    DataType dt = sh.lhDtype();

    switch (dt) {
      case dtUchar:  {
        kaLatticeTest<uint8_t>(argc, argv);
      }
      break;
      case dtSchar:  {
        kaLatticeTest<int8_t>(argc, argv);
      }
      break;
      case dtShort:  {
        kaLatticeTest<int16_t>(argc, argv);
      }
      break;
      case dtUshort: {
        kaLatticeTest<uint16_t>(argc, argv);
      }
      break;
      case dtInt:    {
        kaLatticeTest<int32_t>(argc, argv);
      }
      break;
      case dtUint:   {
        kaLatticeTest<uint32_t>(argc, argv);
      }
      break;
      case dtLong:   {
        kaLatticeTest<int64_t>(argc, argv);
      }
      break;
      case dtUlong:  {
        kaLatticeTest<uint64_t>(argc, argv);
      }
      break;
      case dtFloat:  {
        kaLatticeTest<float>(argc, argv);
      }
      break;
      case dtDouble: {
        kaLatticeTest<double>(argc, argv);
      }
      break;
      case dtUnknown:
      default:
        throw kaBaseException("Unknown or unsupported type of lattice %s", argv[1]);
    }  // switch
  } catch (kaBaseException &e) {
    cerr << argv[0] << " Error: " << e << endl;
    exit(-1);
  }
  return 0;
}  // main
