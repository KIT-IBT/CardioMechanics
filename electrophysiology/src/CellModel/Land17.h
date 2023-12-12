/*
 * File: Land17.h
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

// based on the available Matlab code from cemrg.co.uk

#ifndef LAND17
#define LAND17

#include <Land17Parameters.h>

#undef v
#define v(a) L17p->P[NS_Land17Parameters::a].value

class Land17 : public vbForceModel<ML_CalcType> {
 public:
  Land17Parameters *L17p;
  ML_CalcType XS;
  ML_CalcType XW;
  ML_CalcType TRPN;
  ML_CalcType TmBlocked;
  ML_CalcType ZETAS;
  ML_CalcType ZETAW;
  ML_CalcType Ta;
  ML_CalcType T;
  ML_CalcType Cd;

  Land17(Land17Parameters *);
  ~Land17() {}

  virtual inline ML_CalcType *GetBase(void) {return &XS; }

  virtual inline int GetSize(void) {
    return sizeof(Land17)-sizeof(vbForceModel<ML_CalcType>)-sizeof(Land17Parameters *);
  }

  virtual inline bool InIsTCa(void) {return false; }

  void Init();
  virtual ML_CalcType Calc(double, ML_CalcType, ML_CalcType, ML_CalcType &, int);
  virtual void Print(ostream &);
  virtual void GetParameterNames(vector<string> &);
};  // class Land17
#endif  // ifndef LAND17
