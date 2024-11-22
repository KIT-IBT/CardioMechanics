/* -------------------------------------------------------

   Tomek.h

   Ver. 1.0.0

   Created:       Jan-Erik Duhme (12.2023)
   Last modified: Jan-Erik Duhme (04.12.2023)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#ifndef TOMEK
#define TOMEK

#include <TomekParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_TomekParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_TomekParameters::a].value
#endif  // ifdef HETERO


class Tomek : public vbElphyModel<ML_CalcType> {
 public:
  TomekParameters *ptTeaP;

  /// state variables
  ML_CalcType m, A_h, B_h, h, A_j, B_j, j, h_p, j_p, m_L, h_L, h_L_CaMK;  // I_Na
  ML_CalcType a, i_fast, i_slow, a_CaMK, i_CaMK_fast, i_CaMK_slow;  // I_to
  ML_CalcType d, f_fast, f_slow, f_Ca_fast, f_Ca_slow, j_Ca, f_CaMK_fast, f_Ca_CaMK_fast, n_ss, n_i;  // I_CaL
  ML_CalcType C_0, C_1, C_2, O, I;  // I_Kr
  ML_CalcType x_s1, x_s2;  // I_Ks
  ML_CalcType Na_i, Na_ss, K_i, K_ss, Ca_i, Ca_ss, Ca_nsr, Ca_jsr, Cl_i;  // ion concentrations and buffers
  ML_CalcType CaMK_trap;  // Calcium/Calmodulin-dependent protein kinase
  ML_CalcType J_rel_NP, J_rel_CaMK;  // SR calcium release flux via RyRs
  #ifdef TRPN
  ML_CalcType Ca_TRPN;  // TROPONIN
  #endif  // ifdef TRPN

  /// currents
  ML_CalcType I_Na, I_Na_late, I_to, I_CaL, I_CaNa, I_CaK, I_CaL_i, I_CaL_ss, I_CaNa_i, I_CaNa_ss, I_CaK_i, I_CaK_ss, I_Kr, I_Ks, I_K1, I_NaCa_i, I_NaCa_ss;
  ML_CalcType I_NaK, I_CaCl_junc, I_CaCl_sl, I_CaCl, I_Nab, I_Cab, I_Kb, I_Clb, I_pCa;
  ML_CalcType J_diff_Na, J_diff_Ca, J_diff_K, J_leak, J_rel, J_tr, J_up;
  ML_CalcType cur_I_tot, I_tot;
  #ifdef ISAC
  ML_CalcType I_SAC;
  #endif // ifdef ISAC


  Tomek(TomekParameters *pp);
  ~Tomek();
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

  virtual inline  ML_CalcType GetCli() {return Cl_i; }

  virtual inline  ML_CalcType GetClo() {return v(VT_Cl_o); }

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
};  // class Tomek
#endif  // ifndef TOMEK
