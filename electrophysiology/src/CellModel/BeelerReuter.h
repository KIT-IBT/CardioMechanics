/*
 * File: BeelerReuter.h
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


#ifndef BEELER_REUTER
#define BEELER_REUTER

#include <BeelerReuterParameters.h>

#undef HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_BeelerReuterParameters::a)
#else  // ifdef HETERO
# define v(a) pBRP->P[NS_BeelerReuterParameters::a].value
#endif  // ifdef HETERO

class BeelerReuter : public vbElphyModel<ML_CalcType> {
 public:
  BeelerReuterParameters *pBRP;
  ML_CalcType Ca_i;
  ML_CalcType m, h, j, d, f, x1;

#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO

  BeelerReuter(BeelerReuterParameters *);
  ~BeelerReuter() {}

  virtual void Init();

  virtual inline int GetSize(void) {
    return sizeof(BeelerReuter)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(BeelerReuterParameters *)
#ifdef HETERO
           -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
    ;
  }

  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline ML_CalcType *GetBase(void) {return &Ca_i;}

  virtual inline ML_CalcType Volume() {return v(VT_Vol);}

  virtual inline ML_CalcType GetAmplitude() {return v(VT_Amp); /*30.0;*/}

  virtual inline ML_CalcType GetStimTime() {return 0.003;}

  virtual inline ML_CalcType GetVm() {return v(VT_Init_Vm);}

  virtual inline ML_CalcType GetCai() {return Ca_i;}

  virtual inline ML_CalcType GetCao() {return 0.0;}

  virtual inline ML_CalcType GetNai() {return 0.0;}

  virtual inline ML_CalcType GetNao() {return 0.0;}

  virtual inline ML_CalcType GetKi() {return 0.0;}

  virtual inline ML_CalcType GetKo() {return 0.0;}

  virtual inline void SetCai(ML_CalcType val) {Ca_i = val;}

  virtual void Print(ostream &, double, ML_CalcType);
  virtual void LongPrint(ostream &, double, ML_CalcType);
  virtual void GetParameterNames(vector<string> &);
  virtual void GetLongParameterNames(vector<string> &);

  virtual inline unsigned char getSpeed(ML_CalcType adVm) {
    return (unsigned char)(adVm < .15e-6 ? 3 : (adVm < .3e-6 ? 2 : 1));
  }

  virtual ML_CalcType Calc(double, ML_CalcType, ML_CalcType, ML_CalcType, int);

  virtual int GetNumStatus() {return 7;}

  virtual void GetStatus(double *) const;
  virtual void SetStatus(const double *);
};  // class BeelerReuter

#endif  // ifndef BEELER_REUTER
