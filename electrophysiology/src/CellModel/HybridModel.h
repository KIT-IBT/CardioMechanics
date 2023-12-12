/*
 * File: HybridModel.h
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


#ifndef HYBRIDMODEL_H
#define HYBRIDMODEL_H

#include <HybridModelParameters.h>

#undef v
#define v(a) hmp->P[NS_HybridModelParameters::a].value

class HybridModel : public vbForceModel<ML_CalcType> {
 public:
  HybridModelParameters *hmp;

  ML_CalcType AMATP, MATP, MADPP, AwMADPP, AsMADPP, AsMADP, AMADP, TCa, TMon, MADP, M;

  HybridModel(HybridModelParameters *);
  ~HybridModel() {}

  virtual inline int GetSize(void) {
    return sizeof(HybridModel)-sizeof(vbForceModel<ML_CalcType>)-sizeof(HybridModelParameters *);
  }

  virtual inline ML_CalcType *GetBase(void) {return &AMATP;}

  void Init();

  virtual inline bool InIsTCa(void) {return true;}

  virtual ML_CalcType Calc(double, ML_CalcType, ML_CalcType, ML_CalcType &, int);

  virtual ML_CalcType CalcTrop(double, ML_CalcType, ML_CalcType, ML_CalcType, int);

  virtual inline ML_CalcType ForceEulerNEuler(int, ML_CalcType, ML_CalcType, ML_CalcType, ML_CalcType, ML_CalcType);
  virtual void Print(ostream &);
  virtual void GetParameterNames(vector<string> &);
}; // class HybridModel
#endif  // ifndef HYBRIDMODEL_H
