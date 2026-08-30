/*
 * File: TWorld.h
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


#ifndef TWORLD
#define TWORLD

#include <TWorldParameters.h>

#define HETERO
#undef v

#ifdef HETERO
# define v(a) PS->getValue(NS_TWorldParameters::a)
#else  // ifdef HETERO
# define v(a) ptTeaP->P[NS_TWorldParameters::a].value
#endif  // ifdef HETERO


class TWorld : public vbElphyModel<ML_CalcType> {
 public:
  TWorldParameters *ptTeaP;

  /// state variables currents
    ML_CalcType m, A_h, B_h, A_j, B_j, h, j, h_p, j_p, m_PKA, h_PKA, j_PKA, h_both, j_both, m_L, h_L, h_L_p; // Sodium current (INa, INaL)
    ML_CalcType d, f_fast, f_slow, f_Ca_fast, f_Ca_slow, j_Ca, f_p_fast, f_Ca_p_fast, d_PKA, f_PKA_fast, f_PKA_slow, f_Ca_PKA_fast, f_Ca_PKA_slow, f_both_fast, f_Ca_both_fast, n_Ca_dyad, n_Ca_sl; // L-type calcium current (I_CaL, I_CaNa, I_CaK)
    ML_CalcType a_slow, a_fast, i_slow, i_fast, a_p_slow, a_p_fast, i_p_slow, i_p_fast; // Transient outward current (Ito)
    ML_CalcType C_0, C_1, C_2, O, I; // Rapid delayed rectifier current (IKr)
    ML_CalcType xs_dyad, xs_sl; // Slow delayed rectifier current (IKs)
    ML_CalcType ryr_R, ryr_O, ryr_I, ryr_CaRI, ryr_R_p, ryr_O_p, ryr_I_p, ryr_CaRI_p; // Calcium release from SR (Jrel, Jleak)

  /// currents
    ML_CalcType I_Na_Base_NP, I_Na_Base_CaMK, I_Na_Base_PKA, I_Na_Base_Both, I_Na_Base, I_NaFast_dyad, I_NaFast_sl, I_NaL_dyad, I_NaL_sl, I_Na_dyad, I_Na_sl, I_NaFast, I_NaL; // Sodium current (INa, INaL)
    ML_CalcType I_CaL_pureCDI_dyad, I_CaL_pureCDI_sl, I_CaL_dyad_NP, I_CaL_dyad_CaMK, I_CaL_dyad_PKA, I_CaL_dyad_Both, I_CaL_sl_NP, I_CaL_sl_CaMK, I_CaL_sl_PKA, I_CaL_sl_Both, I_CaNa_dyad_NP, I_CaNa_dyad_CaMK, I_CaNa_dyad_PKA, I_CaNa_dyad_Both, I_CaNa_sl_NP, I_CaNa_sl_CaMK, I_CaNa_sl_PKA, I_CaNa_sl_Both, I_CaK_dyad_NP, I_CaK_dyad_CaMK, I_CaK_dyad_PKA, I_CaK_dyad_Both, I_CaK_sl_NP, I_CaK_sl_CaMK, I_CaK_sl_PKA, I_CaK_sl_Both, I_CaL_dyad, I_CaNa_dyad, I_CaK_dyad, I_CaL_sl, I_CaNa_sl, I_CaK_sl; // L-type calcium current (I_CaL, I_CaNa, I_CaK)
    ML_CalcType I_to_slow, I_to_fast, I_to; // Rapid delayed rectifier current (IKr)
    ML_CalcType I_Kr; // Rapid delayed rectifier current (IKr)
    ML_CalcType I_Ks_dyad, I_Ks_sl, I_Ks; // Slow delayed rectifier current (IKs)
    ML_CalcType I_K1; // Inward rectifier current (IK1)
    ML_CalcType I_NaCa_sl, I_NaCa_dyad, I_NaCa; // Sodium-calcium exchanger (INaCa)
    ML_CalcType I_NaK_dyad_noPKA, I_NaK_dyad_PKA, I_NaK_sl_noPKA, I_NaK_sl_PKA, I_NaK_dyad, I_NaK_sl, I_NaK; // Sodium-potassium pump (INaK)
    ML_CalcType I_CaCl_dyad, I_CaCl_sl, I_CaCl, I_Clb; // Chloride currents (ICaCl, IClb)
    ML_CalcType I_Nab_dyad, I_Nab_sl, I_Nab, I_Kb, I_Cab_dyad, I_Cab_sl, I_Cab; // Background currents (INab, ICab, IKb)
    ML_CalcType J_rel_ICaLdep_act, J_rel_ICaLdep_f1, J_rel_ICaLdep_f2, J_rel_ICaLdep, J_SR_Carel_NP, J_SR_Carel_CaMK, J_SR_Carel, J_SR_leak; // Calcium release from SR (Jrel, Jleak)
    ML_CalcType J_up_NP, J_up_CaMK, J_up; // Calcium reuptake to the SR (Jup)
    ML_CalcType I_pCa_dyad, I_pCa_sl, I_pCa; // Sarcolemmal calcium pump (pCa)
    
  /// ion concentrations
    ML_CalcType I_Na_tot_dyad, I_Na_tot_sl, I_Na_tot, Na_o, Na_dyad, Na_sl, Na_myo; // Sodium Concentration
    ML_CalcType I_K_tot, K_myo; // Potassium Concentration
    ML_CalcType I_Cl_tot, Cl_myo; // Cloride Concentration
    ML_CalcType I_Ca_tot_dyad, I_Ca_tot_sl, I_Ca_tot, Ca_dyad, Ca_sl, Ca_i, Ca_SR; // Calcium Concentration
    ML_CalcType I_tot;

  /// Signalling and Buffering
    ML_CalcType CaMK_trap, CaMK_f_ICaL, CaMK_f_RyR, CaMK_f_PLB, casig_serca_trap; // CaMK and Ca Signalling
    ML_CalcType Buffer_NaBj, Buffer_NaBsl, Buffer_TnClow, Buffer_TnCHc, Buffer_TnCHm, Buffer_CaM, Buffer_Myosin_ca, Buffer_Myosin_mg, Buffer_SRB, Buffer_SLLj, Buffer_SLLsl, Buffer_SLHj, Buffer_SLHsl, Buffer_Csqn; // Buffering
    
  /// Land-Niederer model of contraction
    ML_CalcType XS;
    ML_CalcType XW;
    ML_CalcType CaTRPN;
    ML_CalcType TmBlocked;
    ML_CalcType ZETAS;
    ML_CalcType ZETAW;
    ML_CalcType Ta;
    ML_CalcType Tension;
    ML_CalcType Cd;
    

  TWorld(TWorldParameters *pp);
  ~TWorld();
#ifdef HETERO
  ParameterSwitch *PS;
#endif  // ifdef HETERO
  virtual inline bool AddHeteroValue(string desc, double val);

  virtual inline  ML_CalcType SurfaceToVolumeRatio() {return 1.0; }

  virtual inline  ML_CalcType Volume() {return 1.15606568e-12; }

  virtual inline  ML_CalcType GetVm() {return v(VT_Init_Vm); }

  virtual inline  ML_CalcType GetCai() {return Ca_i; }

  virtual inline  ML_CalcType GetCao() {return v(VT_Ca_o); }

  virtual inline  ML_CalcType GetNai() {return Na_myo; }

  virtual inline  ML_CalcType GetNao() {return v(VT_Na_o); }

  virtual inline  ML_CalcType GetKi() {return K_myo; }

  virtual inline  ML_CalcType GetKo() {return v(VT_K_o); }

  virtual inline  ML_CalcType GetCli() {return Cl_myo; }

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
};  // class TWorld
#endif  // ifndef TWORLD
