/*
 * File: Kurata.h
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


#ifndef KURATA
#define KURATA

#include <KurataParameters.h>

#undef HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_KurataParameters::a)
#else  // ifdef HETERO
# define v(a) pCmP->P[NS_KurataParameters::a].value
#endif  // ifdef HETERO

class Kurata : public vbElphyModel<ML_CalcType> {
 public:
  KurataParameters *pCmP;
  ML_CalcType gdL, gfL, gpa;
  ML_CalcType gn, gq, gh;
  ML_CalcType Cai, Carel, Caup;
  ML_CalcType gdR, gfR, Rtc;
  ML_CalcType Nai;
  ML_CalcType Ki;
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  Kurata(KurataParameters *pp);
  ~Kurata();
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType Volume() {return 2.89529179*10e-13;}

  virtual inline  ML_CalcType GetAmplitude() {return v(VT_Amp); /*80;*/}

  virtual inline  ML_CalcType GetStimTime() {return 0.001;}

  virtual inline  ML_CalcType GetVm();

  virtual inline  ML_CalcType GetCai() {return Cai;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Cao);}

  virtual inline  ML_CalcType GetNai() {return Nai;}

  virtual inline  ML_CalcType GetNao() {return v(VT_Nao);}

  virtual inline  ML_CalcType GetKi() {return Ki;}

  virtual inline  ML_CalcType GetKo() {return v(VT_Ko);}

  virtual inline int GetSize(void);

  virtual inline ML_CalcType *GetBase(void) {return &gdL;}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                    int euler);
  virtual void Print(ostream &tempstr, double t,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double t,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class Kurata
#endif  // ifndef KURATA
