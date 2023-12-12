/*
 * File: ForceDummy.h
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


#ifndef FORCE_DUMMY_H
#define FORCE_DUMMY_H

#include <ForceModelBasis.h>

template<class T>
class ForceDummy : virtual public vbForceModel<T> {
 public:
  ForceDummy() {}

  ~ForceDummy() {}

  virtual void Init() {}

  virtual T Calc(double tinc, T stretch, T velocity, T &Ca, int euler = 1) {return 0.0;}

  virtual void Print(ostream &tempstr) {}

  virtual int GetSize(void) {return 0;}

  virtual T *GetBase(void) {return NULL;}

  virtual int GetNumParameters() {return 0;}

  virtual void WriteForceStatus(FILE *writeto) {}

  virtual void ReadForceStatus(FILE *readfrom) {}

  virtual void GetParameterNames(vector<string> &getpara) {}
};
#endif  // ifndef FORCE_DUMMY_H
