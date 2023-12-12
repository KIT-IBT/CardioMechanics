/*
 * File: SachseEtAl.h
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



#include <SachseEtAlParameters.h>

#undef HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_SachseEtAlParameters::a)
#else  // ifdef HETERO
# define v(a) pParam->P[NS_SachseEtAlParameters::a].value
#endif  // ifdef HETERO

class SachseEtAl : public vbElphyModel<ML_CalcType> {
 public:
  double Ki;

  double C0Shaker;
  double C1Shaker;
  double C2Shaker;
  double C3Shaker;
  double C4Shaker;
  double OShaker;

  SachseEtAlParameters *pParam;

  SachseEtAl(SachseEtAlParameters *pParamArg) {
    pParam = pParamArg;
    Init();
  }

  ~SachseEtAl() {}

  virtual void Init();

  virtual inline ML_CalcType Volume() {return v(VT_Vfibro);}

  virtual inline ML_CalcType GetAmplitude() {return v(VT_Amp);}

  virtual inline ML_CalcType GetStimTime() {return 0.;}

  virtual inline ML_CalcType GetVm() {return v(VT_Vm);}

  virtual inline ML_CalcType GetCai() {return v(VT_Cai);}

  virtual inline ML_CalcType GetCao() {return v(VT_Cao);}

  virtual inline ML_CalcType GetNai() {return 0;}

  virtual inline ML_CalcType GetNao() {return 0;}

  virtual inline ML_CalcType GetKi() {return Ki;}

  virtual inline ML_CalcType GetKo() {return v(VT_Ko);}

  virtual inline int GetSize(void) {
    return sizeof(SachseEtAl)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(SachseEtAlParameters *);
  }

  virtual inline ML_CalcType *GetBase(void) {return (ML_CalcType *)&Ki;}

  virtual inline void SetCai(ML_CalcType val) {}

  virtual void        Print(ostream &tempstr, double t, ML_CalcType V);
  virtual void        LongPrint(ostream &tempstr, double t, ML_CalcType V);
  virtual void        GetParameterNames(vector<string> &getpara);
  virtual void        GetLongParameterNames(vector<string> &getpara);
  virtual ML_CalcType Calc(double tinc, ML_CalcType V, ML_CalcType I_Stim = .0, ML_CalcType stretch = 1.,
                           int euler                                      = 1);

  virtual int GetNumStatus() {return 7;}

  virtual void GetStatus(double *p) const;
  virtual void SetStatus(const double *p);
};  // class SachseEtAl
