/*
 * File: GrandiEtAlVentricle.h
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


#ifndef GRANDIETALVENTRICLE
#define GRANDIETALVENTRICLE

#include <GrandiEtAlVentricleParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_GrandiEtAlVentricleParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_GrandiEtAlVentricleParameters::a].value
#endif  // ifdef HETERO

class GrandiEtAlVentricle : public vbElphyModel<ML_CalcType> {
 public:
  GrandiEtAlVentricleParameters *ptTeaP;
  ML_CalcType Na_j;
  ML_CalcType Na_sl;
  ML_CalcType K_i;
  ML_CalcType Ca_j;
  ML_CalcType Ca_sl;
  ML_CalcType m;
  ML_CalcType h;
  ML_CalcType j;
  ML_CalcType x_kr;
  ML_CalcType x_ks;
  ML_CalcType Na_i;
  ML_CalcType x_to_s;
  ML_CalcType y_to_s;
  ML_CalcType x_to_f;
  ML_CalcType y_to_f;
  ML_CalcType d;
  ML_CalcType f;
  ML_CalcType f_Ca_Bj;
  ML_CalcType f_Ca_Bsl;
  ML_CalcType Ry_Rr;
  ML_CalcType Ry_Ro;
  ML_CalcType Ry_Ri;
  ML_CalcType Ca_sr;
  ML_CalcType Ca_i;
  ML_CalcType Na_Bj;
  ML_CalcType Na_Bsl;
  ML_CalcType Tn_CL;
  ML_CalcType Tn_CHc;
  ML_CalcType Tn_CHm;
  ML_CalcType CaM;
  ML_CalcType Myo_c;
  ML_CalcType Myo_m;
  ML_CalcType SRB;
  ML_CalcType SLL_j;
  ML_CalcType SLL_sl;
  ML_CalcType SLH_j;
  ML_CalcType SLH_sl;
  ML_CalcType Csqn_b;

  GrandiEtAlVentricle(GrandiEtAlVentricleParameters *pp);
  ~GrandiEtAlVentricle();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0;}

  virtual inline  ML_CalcType Volume() {return (v(VT_Vcell))*1000.0;}

  virtual inline  ML_CalcType GetVm() {return v(VT_V_init);}

  virtual inline  ML_CalcType GetCai() {return Ca_i;}

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o);}

  virtual inline  ML_CalcType GetNai() {return Na_i;}

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o);}

  virtual inline  ML_CalcType GetKi() {return K_i;}

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o);}

  virtual inline  ML_CalcType *GetBase(void) {return &Na_j;}

  virtual inline  ML_CalcType GetIto() {return 0.0;}

  virtual inline  ML_CalcType GetIKr() {return 0.0;}

  virtual inline  ML_CalcType GetIKs() {return 0.0;}

  virtual inline int GetSize(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0;}

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_stim_amplitude);}

  virtual inline  ML_CalcType GetStimTime() {return v(VT_stim_duration);}

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch                                  = 1., int euler = 2);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class GrandiEtAlVentricle
#endif  // ifndef GRANDIETALVENTRICLE
