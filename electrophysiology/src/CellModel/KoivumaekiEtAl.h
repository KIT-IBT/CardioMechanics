/*
 * File: KoivumaekiEtAl.h
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


#ifndef KOIVUMAEKIETAL
#define KOIVUMAEKIETAL

#include <KoivumaekiEtAlParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_KoivumaekiEtAlParameters::a)
#else  // ifdef HETERO
# define v(a) ptKeaP->P[NS_KoivumaekiEtAlParameters::a].value
#endif  // ifdef HETERO

class KoivumaekiEtAl : public vbElphyModel<ML_CalcType> {
 public:
  KoivumaekiEtAlParameters *ptKeaP;
  ML_CalcType m;
  ML_CalcType h1;
  ML_CalcType h2;
  ML_CalcType d;
  ML_CalcType f1;
  ML_CalcType fca;
  ML_CalcType r;
  ML_CalcType s;
  ML_CalcType sus_r;
  ML_CalcType sus_s;
  ML_CalcType n;
  ML_CalcType pa;
  ML_CalcType y;
  ML_CalcType RyR_oss;
  ML_CalcType RyR_css;
  ML_CalcType RyR_ass;
  ML_CalcType RyR_o1;
  ML_CalcType RyR_c1;
  ML_CalcType RyR_a1;
  ML_CalcType RyR_o2;
  ML_CalcType RyR_c2;
  ML_CalcType RyR_a2;
  ML_CalcType RyR_o3;
  ML_CalcType RyR_c3;
  ML_CalcType RyR_a3;
  ML_CalcType SERCA_Ca1;
  ML_CalcType SERCA_Ca2;
  ML_CalcType SERCA_Ca3;
  ML_CalcType SERCA_Cass;
  ML_CalcType Na_ss;
  ML_CalcType Na_i;
  ML_CalcType K_i;
  ML_CalcType HTRPNCa1;
  ML_CalcType HTRPNCa2;
  ML_CalcType HTRPNCa3;
  ML_CalcType HTRPNCa4;
  ML_CalcType LTRPNCa1;
  ML_CalcType LTRPNCa2;
  ML_CalcType LTRPNCa3;
  ML_CalcType LTRPNCa4;
  ML_CalcType Ca_ss;
  ML_CalcType Ca_i1;
  ML_CalcType Ca_i2;
  ML_CalcType Ca_i3;
  ML_CalcType Ca_i4;
  ML_CalcType Ca_SR1;
  ML_CalcType Ca_SR2;
  ML_CalcType Ca_SR3;
  ML_CalcType Ca_SR4;

  KoivumaekiEtAl(KoivumaekiEtAlParameters *pp);
  ~KoivumaekiEtAl();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return (v(VT_V_cytosol))*(v(VT_V_corrcell))*1e-12;}

  virtual inline  ML_CalcType GetVm() {return v(VT_V_init);}

  virtual inline  ML_CalcType GetCai() {return (Ca_i1*0.0035+Ca_i2*0.0025+Ca_i3*0.0015+Ca_i4*0.0005)/0.008;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o);}

  virtual inline  ML_CalcType GetNai() {return Na_i;}

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o);}

  virtual inline  ML_CalcType GetKi() {return K_i;}

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o);}

  virtual inline ML_CalcType *GetBase(void) {return &m;}

  virtual inline int GetSize(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_Amp);}

  virtual inline  ML_CalcType GetStimTime() {return v(VT_stim_duration);}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch                                  = 1., int euler = 2);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class KoivumaekiEtAl
#endif  // ifndef KOIVUMAEKIETAL
