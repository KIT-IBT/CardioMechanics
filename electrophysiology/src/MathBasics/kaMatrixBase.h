/*
 * File: kaMatrixBase.h
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

#ifndef KAMATRIXBASE_H
#define KAMATRIXBASE_H

#include <kaMachineOS.h>

/*!
   \class  kaMatrixBase
   \brief Abstract matrix class
 */
template<class T> class kaMatrixBase {
 public:
  virtual T       *Pointer()       = 0;
  virtual const T *Pointer() const = 0;

  virtual T          & a(int, int)     = 0;
  virtual T          & a(int)          = 0;
  virtual unsigned int NumCols() const = 0;
  virtual unsigned int NumRows() const = 0;

  virtual T       *operator[](int i)        = 0;
  virtual const T *operator[](int i) const  = 0;
  virtual const T  Get(int z, int s) const  = 0;
  virtual const T  Get(int p) const         = 0;
  virtual void     Print(const char *format = "%7.3lf\t") const = 0;
  virtual void     Save(FILE *) const       = 0;
  virtual void     Restore(FILE *)          = 0;
};


#endif  // ifndef KAMATRIXBASE_H
