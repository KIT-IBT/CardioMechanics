/*
 * File: ElphyDummy.h
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


#ifndef ELPHY_DUMMY_H
#define ELPHY_DUMMY_H

#include <ElphyModelBasis.h>

template<class T>
class ElphyDummy : public vbElphyModel<T> {
 public:
  ElphyDummy() {}

  ~ElphyDummy() {}

  virtual void Init() {}

  virtual inline bool AddHeteroValue(string desc, double val) {return false;}

  virtual inline T Volume() {return 2.064403e-13;}

  virtual inline T GetVm() {return -.085;}

  virtual inline T GetCai() {return .0;}

  virtual inline T GetCao() {return .0;}

  virtual inline T GetNai() {return .0;}

  virtual inline T GetNao() {return .0;}

  virtual inline T GetKi() {return .0;}

  virtual inline T GetKo() {return .0;}

  virtual inline int GetSize(void) {return 0;}

  virtual inline T *GetBase(void) {return 0;}

  virtual inline T GetSpeedupMax(void) {return .0;}

  virtual T GetAmplitude(void) {return .0;}

  virtual inline T GetStimTime() {return 0.003;}

  virtual inline unsigned char getSpeed(T adVm) {
    return (unsigned char)5;
  }

  virtual T Calc(double tinc, T V, T i_extern = 0., T stretch = 1., int euler = 1) {return .0;}

  virtual void Print(ostream &tempstr, double t, T V) {}

  virtual void LongPrint(ostream &tempstr, double t, T V) {}

  virtual void GetParameterNames(vector<string> &getpara) {}

  virtual void GetLongParameterNames(vector<string> &getpara) {}
};  // class ElphyDummy
#endif  // ifndef ELPHY_DUMMY_H
