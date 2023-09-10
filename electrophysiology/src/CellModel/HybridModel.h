/**@file HybridModel.h
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
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
