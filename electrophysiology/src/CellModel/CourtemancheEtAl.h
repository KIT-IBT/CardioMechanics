/* -------------------------------------------------------

   CourtemancheEtAl.h

   Ver. 1.1.0

   Created:       dw (27.02.2007)
   Last modified: Tobias Gerach (30.06.2022)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#ifndef COURTEMANCHE
#define COURTEMANCHE

#include <CourtemancheParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_CourtemancheParameters::a)
#else  // ifdef HETERO
# define v(a) pCmP->P[NS_CourtemancheParameters::a].value
#endif  // ifdef HETERO

class Courtemanche : public vbElphyModel<ML_CalcType> {
 public:
  CourtemancheParameters *pCmP;
  ML_CalcType Ca_i;
  ML_CalcType Ca_up, Ca_rel;
  ML_CalcType h, j, oi, ui, Xr, Xs, u, v, w;
  ML_CalcType Na_i, K_i;
  ML_CalcType m, oa, ua, d, f, f_Ca;
#ifdef ISAC
  ML_CalcType I_SAC;
#endif // ifdef ISAC
#ifdef TRPN
  ML_CalcType Ca_TRPN;
#endif // ifdef TRPN
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  Courtemanche(CourtemancheParameters *pp);
  ~Courtemanche();
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType Volume() {return v(VT_Vcell_outside)*1e-18; }

  virtual inline  ML_CalcType GetAmplitude() {return 20.0; } /// 2nA = 1e3* 2/100 pA/pF = 20 pA/pF

  virtual inline  ML_CalcType GetStimTime() {return 0.002; }

  virtual inline  ML_CalcType GetVm() {return v(VT_Init_Vm); }

  virtual inline  ML_CalcType GetCai() {return Ca_i; }

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o); }

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o); }

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o); }

  virtual inline  ML_CalcType GetNai() {return Na_i; }

  virtual inline  ML_CalcType GetKi() {return K_i; }

  virtual inline int GetSize(void);

  virtual inline ML_CalcType *GetBase(void) {return &Ca_i; }

  virtual inline void SetCai(ML_CalcType val) {Ca_i = val; }

  virtual inline void SetNai(ML_CalcType val) {Na_i = val; }

  virtual inline void SetKi(ML_CalcType val) {K_i = val; }

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                    int euler);
  virtual void Print(ostream &tempstr, double t,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double t,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class Courtemanche
#endif  // ifndef COURTEMANCHE
