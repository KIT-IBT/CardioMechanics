/*
 * File: vtkLatticeSG.h
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


#ifndef VTKLATTICESG_H
#define VTKLATTICESG_H

#include "vtkUnsignedCharArray.h"
#include "vtkIntArray.h"
#include "vtkCharArray.h"
#include "vtkShortArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMatrix4x4.h"
#include "vtkCellData.h"

#include "vtkUniformGrid.h"
#include "vtkUniformGrid.h"
#include <kaLattice.h>
#include <kaLatticeHeader.h>

class vtkLatticeSG : public vtkUniformGrid {
 public:
  static vtkLatticeSG *New();
  vtkTypeMacro(vtkLatticeSG, vtkUniformGrid);
  int  SetLattice(const char *name);
  void ReleaseLattice();

  //! Set compression type (before save, if needed)
  void SetCompressionType(CompressionType cmp) {
    if (this->header)
      this->header->lhCompression(cmp);
  }

  //! Copy the lattice matrix into supplied vtkMatrix4x4
  void GetMatrix(vtkMatrix4x4 *a) const;

  //! Set the lattice matrix for the supplied vtkMatrix4x4
  void SetMatrix(const vtkMatrix4x4 *a);

  //! Write lattice to file
  void Save(char *name = NULL);

 protected:
  vtkLatticeSG();

  ~vtkLatticeSG() {this->ReleaseLattice();}

  void setScalarType(int);

  vtkLatticeSG(const vtkLatticeSG &) {}

  void operator=(const vtkLatticeSG &) {}

  vtkDataArray *latticeScalars;

  kaLatticeSharedHeader *header;

  kaLattice<unsigned char> *uclat;
  vtkUnsignedCharArray *uchararray;

  kaLattice<short> *shlat;
  vtkShortArray *sharray;

  kaLattice<int> *ilat;
  vtkIntArray *intarray;

  kaLattice<float> *flat;
  vtkFloatArray *floatarray;

  kaLattice<double> *dlat;
  vtkDoubleArray *doublearray;
};  // class vtkLatticeSG

#endif  // ifndef VTKLATTICESG_H
