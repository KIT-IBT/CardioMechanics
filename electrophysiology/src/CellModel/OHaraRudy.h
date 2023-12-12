/*
 * File: OHaraRudy.h
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


#ifndef OHARARUDY
#define OHARARUDY

#include <OHaraRudyParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_OHaraRudyParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_OHaraRudyParameters::a].value
#endif  // ifdef HETERO


class OHaraRudy : public vbElphyModel<ML_CalcType> {
 public:
  OHaraRudyParameters *ptTeaP;

  /// state variables
  ML_CalcType m, h_fast, h_slow, j, h_CaMK_fast, h_CaMK_slow, j_CaMK, m_L, h_L, h_L_CaMK;  // I_Na
  ML_CalcType a, i_fast, i_slow, a_CaMK, i_CaMK_fast, i_CaMK_slow;  // I_to
  ML_CalcType d, f_fast, f_slow, f_Ca_fast, f_Ca_slow, j_Ca, n, f_CaMK_fast, f_Ca_CaMK_fast;  // I_CaL
  ML_CalcType x_r_fast, x_r_slow;  // I_Kr
  ML_CalcType x_s1, x_s2;  // I_Ks
  ML_CalcType x_K1;  // I_K1
  ML_CalcType Na_i, Na_ss, K_i, K_ss, Ca_i, Ca_ss, Ca_nsr, Ca_jsr;  // ion concentrations and buffers
  ML_CalcType CaMK_trap;  // Calcium/Calmodulin-dependent protein kinase
  ML_CalcType J_rel_NP, J_rel_CaMK;  // SR calcium release flux via RyRs
  #ifdef TRPN
  ML_CalcType Ca_TRPN;  // TROPONIN
  #endif  // ifdef TRPN

  /// currents
  ML_CalcType I_Na_fast, I_Na_late, I_to, I_CaL, I_CaNa, I_CaK, I_Kr, I_Ks, I_K1, I_NaCa_i, I_NaCa_ss;
  ML_CalcType I_NaK, I_Nab, I_Cab, I_Kb, I_pCa;
  #ifdef ISAC
  ML_CalcType I_SAC;
  #endif // ifdef ISAC


  OHaraRudy(OHaraRudyParameters *pp);
  ~OHaraRudy();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0; }

  virtual inline  ML_CalcType Volume() {return 1.15606568e-12; }

  virtual inline  ML_CalcType GetVm() {return v(VT_Init_Vm); }

  virtual inline  ML_CalcType GetCai() {return Ca_i; }

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o); }

  virtual inline  ML_CalcType GetNai() {return Na_i; }

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o); }

  virtual inline  ML_CalcType GetKi() {return K_i; }

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o); }

  virtual inline  ML_CalcType *GetBase(void) {return &m; }

  virtual inline  ML_CalcType GetIto() {return 0.0; }

  virtual inline  ML_CalcType GetIKr() {return 0.0; }

  virtual inline  ML_CalcType GetIKs() {return 0.0; }

  virtual inline int GetSize(void);

  virtual inline  ML_CalcType GetSpeedupMax(void) {return .0; }

  virtual  ML_CalcType GetAmplitude(void) {return v(VT_amplitude); }

  virtual inline  ML_CalcType GetStimTime() {return v(VT_duration); }

  virtual inline unsigned char getSpeed(ML_CalcType adVm);
  virtual void                 Init();
  virtual  ML_CalcType         Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch = 1.,
                                    int euler = 2);
  virtual void Print(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void LongPrint(ostream &tempstr, double tArg,  ML_CalcType V);
  virtual void GetParameterNames(vector<string> &getpara);
  virtual void GetLongParameterNames(vector<string> &getpara);
};  // class OHaraRudy
#endif  // ifndef OHARARUDY
