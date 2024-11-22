/*
 * File: OHaraRudyIso.cpp
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


#include <OHaraRudyIso.h>

OHaraRudyIso::OHaraRudyIso(OHaraRudyIsoParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_OHaraRudyIsoParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

OHaraRudyIso::~OHaraRudyIso() {}

#ifdef HETERO

inline bool OHaraRudyIso::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool OHaraRudyIso::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int OHaraRudyIso::GetSize(void) {
  return sizeof(OHaraRudyIso)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(OHaraRudyIsoParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char OHaraRudyIso::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void OHaraRudyIso::Init() {
#if KADEBUG
  cerr << "#initializing Class: OHaraRudyIso ... " << endl;
        #endif  // if KADEBUG

  m              = v(VT_Init_m);
  h_fast         = v(VT_Init_h_fast);
  h_slow         = v(VT_Init_h_slow);
  j              = v(VT_Init_j);
  h_CaMK_slow    = v(VT_Init_h_CaMK_slow);
  j_CaMK         = v(VT_Init_j_CaMK);
  h_PKA_fast     = v(VT_Init_h_PKA_fast);
  h_PKA_slow     = v(VT_Init_h_PKA_slow);
  j_PKA          = v(VT_Init_j_PKA);
  h_both_fast    = v(VT_Init_h_both_fast);
  h_both_slow    = v(VT_Init_h_both_slow);
  j_both         = v(VT_Init_j_both);
  m_L            = v(VT_Init_m_L);
  h_L            = v(VT_Init_h_L);
  h_L_CaMK       = v(VT_Init_h_L_CaMK);
  a              = v(VT_Init_a);
  i_fast         = v(VT_Init_i_fast);
  i_slow         = v(VT_Init_i_slow);
  a_CaMK         = v(VT_Init_a_CaMK);
  i_CaMK_fast    = v(VT_Init_i_CaMK_fast);
  i_CaMK_slow    = v(VT_Init_i_CaMK_slow);
  d              = v(VT_Init_d);
  f_fast         = v(VT_Init_f_fast);
  f_slow         = v(VT_Init_f_slow);
  f_Ca_fast      = v(VT_Init_f_Ca_fast);
  f_Ca_slow      = v(VT_Init_f_Ca_slow);
  j_Ca           = v(VT_Init_j_Ca);
  n              = v(VT_Init_n);
  f_CaMK_fast    = v(VT_Init_f_CaMK_fast);
  f_Ca_CaMK_fast = v(VT_Init_f_Ca_CaMK_fast);
  d_PKA          = v(VT_Init_d_PKA);
  f_PKA_fast     = v(VT_Init_f_PKA_fast);
  f_PKA_slow     = v(VT_Init_f_PKA_slow);
  f_Ca_PKA_fast  = v(VT_Init_f_Ca_PKA_fast);
  f_Ca_PKA_slow  = v(VT_Init_f_Ca_PKA_slow);
  f_both_fast    = v(VT_Init_f_both_fast);
  f_Ca_both_fast = v(VT_Init_f_Ca_both_fast);
  x_r_fast       = v(VT_Init_x_r_fast);
  x_r_slow       = v(VT_Init_x_r_slow);
  x_s1           = v(VT_Init_x_s1);
  x_s2           = v(VT_Init_x_s2);
  x_s1_PKA       = v(VT_Init_x_s1_PKA);
  x_s2_PKA       = v(VT_Init_x_s2_PKA);
  x_K1           = v(VT_Init_x_K1);
  Na_i           = v(VT_Init_Na_i);
  Na_ss          = v(VT_Init_Na_ss);
  K_i            = v(VT_Init_K_i);
  K_ss           = v(VT_Init_K_ss);
  Ca_i           = v(VT_Init_Ca_i);
  Ca_ss          = v(VT_Init_Ca_ss);
  Ca_nsr         = v(VT_Init_Ca_nsr);
  Ca_jsr         = v(VT_Init_Ca_jsr);
  CaMK_trap      = v(VT_Init_CaMK_trap);
  J_rel_NP       = v(VT_Init_J_rel_NP);
  J_rel_CaMK     = v(VT_Init_J_rel_CaMK);
  J_rel_PKA      = v(VT_Init_J_rel_PKA);
  J_rel_both     = v(VT_Init_J_rel_both);
  #ifdef TRPN
  Ca_TRPN = v(VT_Init_Ca_TRPN)/0.07;
  #endif  // ifdef TRPN
}  // OHaraRudyIso::Init

ML_CalcType OHaraRudyIso::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
  tinc *= 1000.0;  // second to millisecond conversion
  ML_CalcType V_m = V * 1000.0;

  #if defined(TRPN) || defined(ISAC)
  ML_CalcType lambda_m = std::min(1.2, stretch);
  #endif  // if defined(TRPN) || defined(ISAC)

  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_m)+.5);  // array position

  /// CaMK constants
  double KmCaMK = 0.15;
  double aCaMK  = 0.05;
  double bCaMK  = 0.00068;
  double CaMKo  = 0.05;
  double KmCaM  = 0.0015;

  /// update CaMK
  const ML_CalcType CaMK_bound  = CaMKo * ((1.0 - CaMK_trap) / (1.0 + (KmCaM / Ca_ss)));
  const ML_CalcType CaMK_active = (CaMK_bound + CaMK_trap)
  #ifdef HF
    * v(VT_HF_Gomez_1) // 1.5
  #endif  // ifdef HF
  ;
  const ML_CalcType dCaMK_trap = aCaMK * CaMK_bound * (CaMK_bound + CaMK_trap) -
    (bCaMK * CaMK_trap);
  CaMK_trap += tinc * dCaMK_trap;

  /// reversal potentials
  const ML_CalcType E_Na = v(VT_RToverF) * log(v(VT_Na_o) / Na_i);
  const ML_CalcType E_K  = v(VT_RToverF) * log(v(VT_K_o) / K_i);
  double PKNa            = 0.01833;
  const ML_CalcType E_Ks = v(VT_RToverF) * log((v(VT_K_o) + PKNa * v(VT_Na_o)) / (K_i + PKNa * Na_i));

  // convenient shorthand calculations
  const ML_CalcType VFFoverRT = (V_m * v(VT_F) * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType VFoverRT  = (V_m * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType util_1    = v(VT_A_cap) / (v(VT_F) * v(VT_v_myo));
  const ML_CalcType util_2    = v(VT_v_ss) / v(VT_v_myo);
  const ML_CalcType util_3    = v(VT_A_cap) / (v(VT_F) * v(VT_v_ss));

  /////////////////////////////////////////////////////////
  /// calculate I_Na
  /////////////////////////////////////////////////////////
    
  // gating NP
  const ML_CalcType m_inf = 1.0 / (1.0 + exp(-(V_m + 39.57)/9.871));
  const ML_CalcType tau_m = 1.0 / (6.765 * exp((V_m + 11.64)/34.77) + 8.552 * exp(-(V_m + 77.42)/5.955));
  m = m_inf - (m_inf - m) * exp(-tinc / tau_m);
#ifdef PASSINI
  const ML_CalcType h_inf = 1.0 / (1.0 + exp((V_m + 78.5)/6.22));  // Passini et al. 2016
#else  // ifdef PASSINI
  const ML_CalcType h_inf = 1.0 / (1.0 + exp((V_m + 82.9)/6.086));  // ORd
#endif  // ifdef PASSINI
#ifdef DUTTA
  const ML_CalcType tau_h_fast = 1.0 / (3.686e-6 * exp(-(V_m + 3.8875)/7.8579) + 16.0 * exp((V_m - 0.4963)/9.1843));  // Dutta et al.2017
#else  // ifdef DUTTA
  const ML_CalcType tau_h_fast = 1.0 / (1.432e-05 * exp(-(V_m + 1.196)/6.285) + 6.149 * exp((V_m + 0.5096)/20.27));  // ORd
#endif  // ifdef DUTTA
  const ML_CalcType tau_h_slow = 1.0 / (0.009794 * exp(-(V_m + 17.95)/28.05) + 0.3343 * exp((V_m + 5.730)/56.66));
  h_fast = h_inf - (h_inf - h_fast) * exp(-tinc / tau_h_fast);
  h_slow = h_inf - (h_inf - h_slow) * exp(-tinc / tau_h_slow);
  const ML_CalcType h     = v(VT_A_h_fast) * h_fast + v(VT_A_h_slow) * h_slow;
  const ML_CalcType j_inf = h_inf;
#ifdef DUTTA
  const ML_CalcType tau_j = 4.8590 + 1.0 / (0.8628 * exp(-(V_m + 116.7258)/7.6005) + 1.1096 * exp((V_m + 6.2719)/9.0358));  // Dutta 2017
#else  // ifdef DUTTA
  const ML_CalcType tau_j = 2.038 + 1.0 / (0.02136 * exp(-(V_m + 100.6)/8.281) + 0.3052 * exp((V_m + 0.9941)/38.45));  // ORd
#endif  // ifdef DUTTA
  j = j_inf - (j_inf - j) * exp(-tinc / tau_j);

  // gating CaMK-P
#ifdef PASSINI
  const ML_CalcType h_CaMK_inf = 1.0 / (1.0 + exp((V_m + 84.7) / 6.22));  // Passini et al. 2016
#else  // ifdef PASSINI
  const ML_CalcType h_CaMK_inf = 1.0 / (1.0 + exp((V_m + 89.1) / 6.086));   // ORd
#endif  // ifdef PASSINI
  const ML_CalcType tau_h_CaMK_slow = 3.0 * tau_h_slow;
  h_CaMK_slow = h_CaMK_inf - (h_CaMK_inf - h_CaMK_slow) * exp(-tinc / tau_h_CaMK_slow);
  const ML_CalcType h_CaMK     = v(VT_A_h_fast) * h_fast + v(VT_A_h_slow) * h_CaMK_slow;
  const ML_CalcType j_CaMK_inf = j_inf;
  const ML_CalcType tau_j_CaMK = 1.46 * tau_j;
  j_CaMK = j_CaMK_inf - (j_CaMK_inf - j_CaMK) * exp(-tinc / tau_j_CaMK);

  // gating PKA-P
  const ML_CalcType h_PKA_inf = 1.0 / (1.0 + exp((V_m + 82.9 + 5.0) / 6.086));
  h_PKA_fast = h_PKA_inf - (h_PKA_inf - h_PKA_fast) * exp(-tinc / tau_h_fast);
  h_PKA_slow = h_PKA_inf - (h_PKA_inf - h_PKA_slow) * exp(-tinc / tau_h_slow);
  const ML_CalcType h_PKA = v(VT_A_h_fast) * h_PKA_fast + v(VT_A_h_slow) * h_PKA_slow;
  const ML_CalcType j_PKA_inf = h_PKA_inf;
  j_PKA = j_PKA_inf - (j_PKA_inf - j_PKA) * exp(-tinc / tau_j);
    
  // gating PKA-P & CaMK-P
  const ML_CalcType h_both_inf = 1.0 / (1.0 + exp((V_m + 89.1 + 5.0) / 6.086));
  const ML_CalcType tau_h_both_fast = tau_h_fast;
  const ML_CalcType tau_h_both_slow = tau_h_CaMK_slow;
  h_both_fast = h_both_inf - (h_both_inf - h_both_fast) * exp(-tinc / tau_h_both_fast);
  h_both_slow = h_both_inf - (h_both_inf - h_both_slow) * exp(-tinc / tau_h_both_slow);
  const ML_CalcType h_both     = v(VT_A_h_fast) * h_both_fast + v(VT_A_h_slow) * h_both_slow;
  const ML_CalcType j_both_inf = h_both_inf;
  const ML_CalcType tau_j_both = tau_j_CaMK;
  j_both = j_both_inf - (j_both_inf - j_both) * exp(-tinc / tau_j_both);
    
  // fraction of phosphorylation
  const ML_CalcType phi_INa_CaMK = 1.0 / (1.0 + KmCaMK / CaMK_active);
  const ML_CalcType phi_INa_PKA = v(VT_INa_P_frac);
  const ML_CalcType phi_INa_both = phi_INa_CaMK * phi_INa_PKA;
  const ML_CalcType phi_INa_CaMKonly = phi_INa_CaMK - phi_INa_both;
  const ML_CalcType phi_INa_PKAonly = phi_INa_PKA - phi_INa_both;
    
  // currents
  double G_Na_fast = 75.0;
  double G_Na_fast_PKA = 2.7 * G_Na_fast;
  I_Na_fast_NP = G_Na_fast * (V_m - E_Na) *m*m*m* h * j;
  I_Na_fast_CaMK = G_Na_fast * (V_m - E_Na) *m*m*m* h_CaMK * j_CaMK;
  I_Na_fast_PKA = G_Na_fast_PKA * (V_m - E_Na) *m*m*m* h_PKA * j_PKA;
  I_Na_fast_both =G_Na_fast_PKA * (V_m - E_Na) *m*m*m* h_both * j_both;
    
  I_Na_fast = v(VT_INaF_Multiplier) * ((1-phi_INa_CaMKonly-phi_INa_PKAonly-phi_INa_both)*I_Na_fast_NP + phi_INa_CaMKonly*I_Na_fast_CaMK + phi_INa_PKAonly*I_Na_fast_PKA + phi_INa_both*I_Na_fast_both);

    
  /// calculate I_NaL
  ML_CalcType thL_mod = 1.0;
  #ifdef HF
  thL_mod *= v(VT_HF_Gomez_2);
  #endif // ifdef HF
  const ML_CalcType m_L_inf = 1.0 / (1.0 + exp(-(V_m + 42.85) / 5.264));
  const ML_CalcType tau_m_L = tau_m;
  m_L = m_L_inf - (m_L_inf - m_L) * exp(-tinc / tau_m_L);
  const ML_CalcType h_L_inf = 1.0 / (1.0 + exp((V_m + 87.61) / 7.488));
  h_L = h_L_inf - (h_L_inf - h_L) * exp(-tinc / (thL_mod * v(VT_tau_h_L)));
  const ML_CalcType h_L_CaMK_inf = 1.0 / (1.0 + exp((V_m + 93.81) / 7.488));
  const ML_CalcType tau_h_L_CaMK = 3.0 * v(VT_tau_h_L) * thL_mod;
  h_L_CaMK = h_L_CaMK_inf - (h_L_CaMK_inf - h_L_CaMK) * exp(-tinc / tau_h_L_CaMK);
  const ML_CalcType phi_INaL_CaMK = phi_INa_CaMK;
  double G_Na_L                   = 0.0075;
  I_Na_late = v(VT_INaL_Multiplier)
  #ifdef HF
    * v(VT_HF_Gomez_3)
  #endif // ifdef HF
    * G_Na_L * (V_m - E_Na) * m_L *
    ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_CaMK);

  /// calculate I_to
  const ML_CalcType a_inf = 1.0 / (1.0 + exp(-(V_m - 14.34) / 14.82));
  const ML_CalcType tau_a = 1.0515 /
    (1.0/(1.2089*(1.0 + exp(-(V_m - 18.4099) / 29.3814))) + 3.5/(1.0 + exp((V_m + 100.0) / 29.3814)));
  a = a_inf - (a_inf - a) * exp(-tinc / tau_a);
  const ML_CalcType i_inf           = 1.0 / (1.0 + exp((V_m + 43.94) / 5.711));
  const ML_CalcType delta_epi       = v(VT_celltype) == 1.0 ? (1.0 - 0.95 / (1.0 + exp((V_m + 70.0) / 5.0))) : 1.0;
  const ML_CalcType tau_i_fast_base =  4.562 + 1.0 /
    (0.3933 * exp(-(V_m + 100.0) / 100.0) + 0.08004 * exp((V_m + 50.0) / 16.59));
  const ML_CalcType tau_i_slow_base = 23.62 + 1.0 /
    (0.001416 * exp(-(V_m + 96.52) / 59.05) + 1.7808e-08 * exp((V_m + 114.1) / 8.079));
  const ML_CalcType tau_i_fast = tau_i_fast_base * delta_epi;
  const ML_CalcType tau_i_slow = tau_i_slow_base * delta_epi;
  i_fast = i_inf - (i_inf - i_fast) * exp(-tinc / tau_i_fast);
  i_slow = i_inf - (i_inf - i_slow) * exp(-tinc / tau_i_slow);
  const ML_CalcType A_i_fast   = 1.0 / (1.0 + exp((V_m - 213.6) / 151.2));
  const ML_CalcType A_i_slow   = 1.0 - A_i_fast;
  const ML_CalcType i          = A_i_fast * i_fast + A_i_slow * i_slow;
  const ML_CalcType a_CaMK_inf = 1.0 / (1.0 + exp(-(V_m - 24.34) / 14.82));
  const ML_CalcType tau_a_CaMK = tau_a;
  a_CaMK = a_CaMK_inf - (a_CaMK_inf - a_CaMK) * exp(-tinc / tau_a_CaMK);
  const ML_CalcType i_CaMK_inf         = i_inf;
  const ML_CalcType delta_CaMK_develop = 1.354 + 1.0e-04 / (exp((V_m - 167.4) / 15.89) + exp(-(V_m - 12.23) / 0.2154));
  const ML_CalcType delta_CaMK_recover = 1.0 - 0.5 / (1.0 + exp((V_m + 70.0) / 20.0));
  const ML_CalcType tau_i_CaMK_fast    = tau_i_fast * delta_CaMK_develop * delta_CaMK_recover;
  const ML_CalcType tau_i_CaMK_slow    = tau_i_slow * delta_CaMK_develop * delta_CaMK_recover;
  i_CaMK_fast = i_CaMK_inf - (i_CaMK_inf - i_CaMK_fast) * exp(-tinc / tau_i_CaMK_fast);
  i_CaMK_slow = i_CaMK_inf - (i_CaMK_inf - i_CaMK_slow) * exp(-tinc / tau_i_CaMK_slow);
  const ML_CalcType i_CaMK       = A_i_fast * i_CaMK_fast + A_i_slow * i_CaMK_slow;
  const ML_CalcType phi_Ito_CaMK = phi_INa_CaMK;
  double G_to                    = 0.02;
  I_to = v(VT_Ito_Multiplier)
  #ifdef HF
    * v(VT_HF_Gomez_4)
  #endif // ifdef HF
    * G_to * (V_m - E_K) *
    ((1.0 - phi_Ito_CaMK) * a * i + phi_Ito_CaMK * a_CaMK * i_CaMK);

  /////////////////////////////////////////////////////////
  /// calculate I_CaL, I_CaNa, I_CaK
  /////////////////////////////////////////////////////////

  // gating NP
  const ML_CalcType d_inf = 1.0 / (1.0 + exp(-(V_m + 3.940) / 4.230));
  const ML_CalcType tau_d = 0.6 + 1.0 / (exp(-0.05 * (V_m + 6.0)) + exp(0.09 * (V_m + 14.0)));
  d = d_inf - (d_inf -d) * exp(-tinc / tau_d);
  const ML_CalcType f_inf      = 1.0 / (1.0 + exp((V_m + 19.58) / 3.696));
  const ML_CalcType tau_f_fast = 7.0 + 1.0 / (0.0045 * exp(-(V_m + 20.0) / 10.0) + 0.0045 * exp((V_m + 20.0) / 10.0));
  const ML_CalcType tau_f_slow = 1000.0 + 1.0 / (3.5e-05 * exp(-(V_m + 5.0) / 4.0) + 3.5e-05 * exp((V_m + 5.0) / 6.0));
  f_fast = f_inf - (f_inf - f_fast) * exp(-tinc / tau_f_fast);
  f_slow = f_inf - (f_inf - f_slow) * exp(-tinc / tau_f_slow);
  const ML_CalcType f             = v(VT_A_f_fast) * f_fast + v(VT_A_f_slow) * f_slow;
  const ML_CalcType f_Ca_inf      = f_inf;
  const ML_CalcType tau_f_Ca_fast = 7.0 + 1.0 / (0.04 * exp(-(V_m - 4.0) / 7.0) + 0.04 * exp((V_m - 4.0) / 7.0));
  const ML_CalcType tau_f_Ca_slow = 100.0 + 1.0 / (0.00012 * exp(-V_m / 3.0) + 0.00012 * exp(V_m /7.0));
  f_Ca_fast = f_Ca_inf - (f_Ca_inf - f_Ca_fast) * exp(-tinc / tau_f_Ca_fast);
  f_Ca_slow = f_Ca_inf - (f_Ca_inf - f_Ca_slow) * exp(-tinc / tau_f_Ca_slow);
  const ML_CalcType A_f_Ca_fast = 0.3 + 0.6 / (1.0 + exp((V_m - 10.0) / 10.0));
  const ML_CalcType A_f_Ca_slow = 1.0 - A_f_Ca_fast;
  const ML_CalcType f_Ca        = A_f_Ca_fast * f_Ca_fast + A_f_Ca_slow * f_Ca_slow;
  const ML_CalcType j_Ca_inf    = f_Ca_inf;
  double tau_j_Ca               = 75.0;
  j_Ca = j_Ca_inf - (j_Ca_inf - j_Ca) * exp(-tinc / tau_j_Ca);
  const ML_CalcType k_m2_n    = j_Ca * 1.0;
  double Kmn                  = 0.002;
  double k2n                  = 1000.0;
  const ML_CalcType alpha_n   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_ss)), 4.0));
  n = alpha_n * (k2n / k_m2_n) - (alpha_n * (k2n / k_m2_n) - n) * exp(-k_m2_n * tinc);
  
  // gating CaMK-P
  const ML_CalcType f_CaMK_inf      = f_inf;
  const ML_CalcType tau_f_CaMK_fast = 2.5 * tau_f_fast;
  f_CaMK_fast = f_CaMK_inf - (f_CaMK_inf - f_CaMK_fast) * exp(-tinc / tau_f_CaMK_fast);
  const ML_CalcType f_CaMK             = v(VT_A_f_fast) * f_CaMK_fast + v(VT_A_f_slow) * f_slow;
  const ML_CalcType f_Ca_CaMK_inf      = f_inf;
  const ML_CalcType tau_f_Ca_CaMK_fast = 2.5 * tau_f_Ca_fast;
  f_Ca_CaMK_fast = f_Ca_CaMK_inf - (f_Ca_CaMK_inf - f_Ca_CaMK_fast) * exp(-tinc / tau_f_Ca_CaMK_fast);
  const ML_CalcType f_Ca_CaMK = A_f_Ca_fast * f_Ca_CaMK_fast + A_f_Ca_slow * f_Ca_slow;

  // gating PKA-P
  const ML_CalcType d_PKA_inf = 1.0/(1.0+exp((-(V_m+3.940+12.0))/4.230));
  d_PKA = d_PKA_inf - (d_PKA_inf - d_PKA) * exp(-tinc / tau_d);
  const ML_CalcType f_PKA_inf = 1.0/(1.0+exp((V_m+19.58+8.0)/3.696));
  f_PKA_fast = f_PKA_inf - (f_PKA_inf - f_PKA_fast) * exp(-tinc / tau_f_fast);
  f_PKA_slow = f_PKA_inf - (f_PKA_inf - f_PKA_slow) * exp(-tinc / tau_f_slow);
  const ML_CalcType f_Ca_PKA_inf = f_PKA_inf;
  f_Ca_PKA_fast = f_Ca_PKA_inf - (f_Ca_PKA_inf - f_Ca_PKA_fast) * exp(-tinc / tau_f_Ca_fast);
  f_Ca_PKA_slow = f_Ca_PKA_inf - (f_Ca_PKA_inf - f_Ca_PKA_slow) * exp(-tinc / tau_f_Ca_slow);
  const ML_CalcType f_PKA = v(VT_A_f_fast);
  const ML_CalcType f_Ca_PKA = A_f_Ca_fast * f_Ca_PKA_fast + A_f_Ca_slow * f_Ca_PKA_slow;

  // gating PKA-P & CaMK-P
  const ML_CalcType f_both_inf = f_PKA_inf;
  f_both_fast = f_both_inf - (f_both_inf - f_both_fast) * exp(-tinc / tau_f_Ca_CaMK_fast);
  const ML_CalcType f_both = v(VT_A_f_fast) * f_both_fast + v(VT_A_f_slow) * f_PKA_slow;
  const ML_CalcType f_Ca_both_inf = f_Ca_PKA_inf;
  f_Ca_both_fast = f_Ca_both_inf - (f_Ca_both_inf - f_Ca_both_fast) * exp(-tinc / tau_f_Ca_CaMK_fast);
  const ML_CalcType f_Ca_both = A_f_Ca_fast * f_Ca_both_fast + A_f_Ca_slow * f_Ca_PKA_slow;

  //
  double z_Na                = 1.0;
  double gamma_Nai           = 0.75;
  const ML_CalcType exp_z_Na = exp(z_Na * VFoverRT);
  double gamma_Nao           = 0.75;
  const ML_CalcType Psi_CaNa = z_Na * z_Na * VFFoverRT * (gamma_Nai * Na_ss * exp_z_Na - gamma_Nao * v(VT_Na_o)) / (exp_z_Na - 1.0);
  
  double z_K                 = 1.0;
  double gamma_Ki            = 0.75;
  const ML_CalcType exp_z_K  = exp(z_K * VFoverRT);
  double gamma_Ko            = 0.75;
  const ML_CalcType Psi_CaK = z_K * z_K * VFFoverRT * (gamma_Ki * K_ss * exp_z_K - gamma_Ko * v(VT_K_o)) / (exp_z_K - 1.0);

  double z_Ca                = 2.0;
  double gamma_Cai           = 1.0;
  const ML_CalcType exp_z_Ca = exp(z_Ca * VFoverRT);
  double gamma_Cao           = 0.341;
  double Ca_ss_ceiling = 0.025e-3; //TODO: Checken ob hier die Einheit stimmt
  double Ca_ss_value = Ca_ss_ceiling;
  if (Ca_ss <= Ca_ss_ceiling) {
      Ca_ss_value = Ca_ss;
  }
  const ML_CalcType Psi_Ca   = z_Ca * z_Ca * VFFoverRT * (gamma_Cai * Ca_ss_value * exp_z_Ca - gamma_Cao * v(VT_Ca_o)) / (exp_z_Ca - 1.0);

  // fraction of phosphorylation
  const ML_CalcType phi_ICaL_CaMK = phi_INa_CaMK;
  const ML_CalcType phi_ICaL_PKA = v(VT_ICaL_P_frac);
  const ML_CalcType phi_ICaL_both = phi_ICaL_CaMK * phi_ICaL_PKA;
  const ML_CalcType phi_ICaL_CaMKonly = phi_ICaL_CaMK - phi_ICaL_both;
  const ML_CalcType phi_ICaL_PKAonly = phi_ICaL_PKA - phi_ICaL_both;

  // currents
  I_CaL_NP = v(VT_P_Ca) * Psi_Ca * d * (f * (1.0 - n) + j_Ca * f_Ca * n);
  I_CaL_CaMK = v(VT_P_Ca_CaMK) * Psi_Ca * d * (f_CaMK * (1.0 - n) + f_Ca_CaMK * n * j_Ca);
  I_CaL_PKA = v(VT_P_Ca_PKA) * Psi_Ca * d_PKA * (f_PKA * (1.0 - n) + f_Ca_PKA * n * j_Ca);
  I_CaL_both = v(VT_P_Ca_PKA) * Psi_Ca * d_PKA * (f_both * (1.0 - n) + f_Ca_both * n * j_Ca);
  I_CaL = (1-phi_ICaL_CaMKonly-phi_ICaL_PKAonly-phi_ICaL_both)*I_CaL_NP + phi_ICaL_CaMKonly*I_CaL_CaMK + phi_ICaL_PKAonly*I_CaL_PKA + phi_ICaL_both*I_CaL_both;

  I_CaNa_NP = v(VT_P_CaNa) * Psi_CaNa * d * (f * (1.0 - n) + f_Ca * n * j_Ca);
  I_CaNa_CaMK = v(VT_P_CaNa_CaMK) * Psi_CaNa * d * (f_CaMK * (1.0 - n) + f_Ca_CaMK * n * j_Ca);
  I_CaNa_PKA = v(VT_P_CaNa_PKA) * Psi_CaNa * d_PKA * (f_PKA * (1.0 - n) + f_Ca_PKA * n * j_Ca);
  I_CaNa_both = v(VT_P_CaNa_PKA) * Psi_CaNa * d_PKA * (f_both * (1.0 - n) + f_Ca_both * n * j_Ca);
  I_CaNa = (1-phi_ICaL_CaMKonly-phi_ICaL_PKAonly-phi_ICaL_both)*I_CaNa_NP + phi_ICaL_CaMKonly*I_CaNa_CaMK + phi_ICaL_PKAonly*I_CaNa_PKA + phi_ICaL_both*I_CaNa_both;

  I_CaK_NP = v(VT_P_CaK) * Psi_CaK * d * (f * (1.0 - n) + f_Ca * n * j_Ca);
  I_CaK_CaMK = v(VT_P_CaK_CaMK) * Psi_CaK * d * (f_CaMK * (1.0 - n) + f_Ca_CaMK * n * j_Ca);
  I_CaK_PKA = v(VT_P_CaK_PKA) * Psi_CaK * d_PKA * (f_PKA * (1.0 - n) + f_Ca_PKA * n * j_Ca);
  I_CaK_both = v(VT_P_CaK_PKA) * Psi_CaK * d_PKA * (f_both * (1.0 - n) + f_Ca_both * n * j_Ca);
  I_CaK = (1-phi_ICaL_CaMKonly-phi_ICaL_PKAonly-phi_ICaL_both)*I_CaK_NP + phi_ICaL_CaMKonly*I_CaK_CaMK + phi_ICaL_PKAonly*I_CaK_PKA + phi_ICaL_both*I_CaK_both;


  /// calculate I_Kr
  const ML_CalcType x_r_inf      = 1.0 / (1.0 + exp((-(V_m + 8.337)) / 6.789));
  const ML_CalcType tau_x_r_fast = 12.98 + 1.0 /
    (0.3652 * exp((V_m - 31.66) / 3.869) + 4.123e-05 * exp((-(V_m - 47.78)) / 20.38));
  const ML_CalcType tau_x_r_slow = 1.865 + 1.0 /
    (0.06629 * exp((V_m - 34.7) / 7.355) + 1.128e-05 * exp((-(V_m - 29.74)) / 25.94));
  const ML_CalcType A_x_r_fast = 1.0 / (1.0 + exp((V_m + 54.81) / 38.21));
  const ML_CalcType A_x_r_slow = 1.0 - A_x_r_fast;
  x_r_fast = x_r_inf - (x_r_inf - x_r_fast) * exp(-tinc / tau_x_r_fast);
  x_r_slow = x_r_inf - (x_r_inf - x_r_slow) * exp(-tinc / tau_x_r_slow);
  const ML_CalcType x_r  = A_x_r_fast * x_r_fast + A_x_r_slow * x_r_slow;
  const ML_CalcType R_Kr = 1.0 / ((1.0 + exp((V_m + 55.0) / 75.0)) * (1.0 + exp((V_m - 10.0) / 30.0)));
  double G_Kr            = 0.046;
  I_Kr = v(VT_IKr_Multiplier) * G_Kr * sqrt(v(VT_K_o) / 5.4) * x_r * (V_m - E_K) * R_Kr;

  /////////////////////////////////////////////////////////
  /// calculate I_Ks
  /////////////////////////////////////////////////////////

  // gating NP
  const ML_CalcType x_s1_inf = 1.0/(1.0 + exp(-(V_m + 11.6) / 8.932));
  const ML_CalcType tau_x_s1 = 817.3 + 1.0 /(2.326e-04 * exp((V_m + 48.28) / 17.80)+ 0.001292 * exp(-(V_m + 210.0) / 230.0));
  x_s1 = x_s1_inf - (x_s1_inf - x_s1) * exp(-tinc / tau_x_s1);
  const ML_CalcType x_s2_inf = x_s1_inf;
  const ML_CalcType tau_x_s2 = 1.0 / (0.01 * exp((V_m - 50.0) / 20.0) + 0.0193 * exp(-(V_m + 66.54) / 31.0));
  x_s2 = x_s2_inf - (x_s2_inf - x_s2) * exp(-tinc / tau_x_s2);

  // gating PKA-P
  const ML_CalcType x_s1_PKA_inf = x_s1_inf;
  const ML_CalcType tau_x_s1_PKA = (-1.75)*1.0/(2.326e-4*exp((10+48.28)/17.80)+0.001292*exp((-(10+210.0))/230.0)) + 817.3 + 2.75*1.0/((2.326*10e-4*exp((V_m+48.28)/17.8)) + 0.001292*exp(-(V_m+210.0)/230.0));
  x_s1_PKA = x_s1_PKA_inf - (x_s1_PKA_inf - x_s1_PKA) * exp(-tinc / tau_x_s1_PKA);
  const ML_CalcType x_s2_PKA_inf = x_s1_PKA_inf;
  const ML_CalcType tau_x_s2_PKA = tau_x_s2;
  x_s2_PKA = x_s2_PKA_inf - (x_s2_PKA_inf - x_s2_PKA) * exp(-tinc / tau_x_s2_PKA);

  // fraction of phosphorylation
  const ML_CalcType phi_IKs_PKA = v(VT_IKs_P_frac);

  // currents
  double G_Ks = 0.0034;
  double G_Ks_PKA = 8.0 * G_Ks;
  I_Ks_NP = G_Ks * (1.0 + 0.6 / (1.0 + pow((3.8e-05 / Ca_i), 1.4))) * x_s1 * x_s2 * (V_m - E_Ks);
  I_Ks_PKA = G_Ks_PKA * (1.0 + 0.6 / (1.0 + pow((3.8e-05 / Ca_i), 1.4))) * x_s1_PKA * x_s2_PKA * (V_m - E_Ks);
  I_Ks = v(VT_IKs_Multiplier) * ((1 - phi_IKs_PKA) * I_Ks_NP + phi_IKs_PKA * I_Ks_PKA);

  /// calculating I_K1
  const ML_CalcType x_K1_inf = 1.0 / (1.0 + exp(-(V_m + 2.5538 * v(VT_K_o) + 144.59) / (1.5692 * v(VT_K_o) + 3.8115)));
  const ML_CalcType tau_x_K1 = 122.2 / (exp(-(V_m + 127.2) / 20.36) + exp((V_m + 236.8) / 69.33));
  x_K1 = x_K1_inf - (x_K1_inf - x_K1) * exp(-tinc / tau_x_K1);
  const ML_CalcType R_K1 = 1.0 / (1.0 + exp((V_m + 105.8 - 2.6 * v(VT_K_o)) / 9.493));
  double G_K1            = 0.1908;
  I_K1 = v(VT_IK1_Multiplier)
  #ifdef HF
    * v(VT_HF_Gomez_5)
  #endif // ifdef HF
    * G_K1 * sqrt(v(VT_K_o)) * x_K1 * (V_m - E_K) * R_K1;

  /// calculate I_NaCa_i
  double kna3                   = 88.12;
  double wna                    = 6.0e4;
  double wca                    = 6.0e4;
  double wnaca                  = 5.0e3;
  double qna                    = 0.5224;
  double qca                    = 0.1670;
  const ML_CalcType h_Ca        = exp(qca * VFoverRT);
  const ML_CalcType h_Na        = exp(qna * VFoverRT);
  const ML_CalcType h_1_i       = 1.0 + (Na_i * (1.0 + h_Na)) / kna3;
  const ML_CalcType h_2_i       = (Na_i * h_Na) / (kna3 * h_1_i);
  const ML_CalcType h_3_i       = 1.0 / h_1_i;
  const ML_CalcType h_4_i       = 1.0 + (Na_i * (1.0 + Na_i / v(VT_k_Na2))) / v(VT_k_Na1);
  const ML_CalcType h_5_i       = (Na_i * Na_i) / (h_4_i * v(VT_k_Na1) * v(VT_k_Na2));
  const ML_CalcType h_6_i       = 1.0 / h_4_i;
  const ML_CalcType h_7         = 1.0 + (v(VT_Na_o) * (1.0 + 1.0 / h_Na)) / kna3;
  const ML_CalcType h_8         = v(VT_Na_o) / (kna3 * h_Na * h_7);
  const ML_CalcType h_9         = 1.0 / h_7;
  const ML_CalcType k_3_d       = h_9 * wca;
  const ML_CalcType k_3_dd      = h_8 * wnaca;
  const ML_CalcType k_3         = k_3_d + k_3_dd;
  const ML_CalcType k_4_d_i     = (h_3_i * wca) / h_Ca;
  const ML_CalcType k_4_dd_i    = h_2_i * wnaca;
  const ML_CalcType k_4_i       = k_4_d_i + k_4_dd_i;
  const ML_CalcType k_6_i       = h_6_i * Ca_i * v(VT_k_Ca_on);
  const ML_CalcType k_7_i       = h_5_i * h_2_i * wna;
  const ML_CalcType k_8         = h_8 * v(VT_h_11) * wna;
  const ML_CalcType x_1_i       = v(VT_k_2) * k_4_i * (k_7_i + k_6_i) + v(VT_k_5) * k_7_i * (v(VT_k_2) + k_3);
  const ML_CalcType x_2_i       = v(VT_k_1) * k_7_i * (k_4_i + v(VT_k_5)) + k_4_i * k_6_i * (v(VT_k_1) + k_8);
  const ML_CalcType x_3_i       = v(VT_k_1) * k_3 * (k_7_i + k_6_i) + k_8 * k_6_i * (v(VT_k_2) + k_3);
  const ML_CalcType x_4_i       = v(VT_k_2) * k_8 * (k_4_i + v(VT_k_5)) + k_3 * v(VT_k_5) * (v(VT_k_1) + k_8);
  const ML_CalcType E_1_i       = x_1_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_2_i       = x_2_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_3_i       = x_3_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_4_i       = x_4_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  double KmCaAct                = 150.0e-6;
  const ML_CalcType allo_i      = 1.0 / (1.0 + (KmCaAct / Ca_i) * (KmCaAct / Ca_i));
  const ML_CalcType J_NaCa_Na_i = 3.0 * (E_4_i * k_7_i - E_1_i * k_8) + E_3_i * k_4_dd_i - E_2_i * k_3_dd;
  const ML_CalcType J_NaCa_Ca_i = E_2_i * v(VT_k_2) - E_1_i * v(VT_k_1);
  double G_NaCa                 = 0.0008;
  I_NaCa_i = v(VT_INaCa_Multiplier)
  #ifdef HF
    * v(VT_HF_Gomez_6)
  #endif // ifdef HF
    * G_NaCa * 0.8 * allo_i * (z_Na * J_NaCa_Na_i + z_Ca * J_NaCa_Ca_i);

  /// calculate I_NaCa_ss
  const ML_CalcType h_1_ss       = 1.0 + (Na_ss * (1.0 + h_Na)) / kna3;
  const ML_CalcType h_2_ss       = (Na_ss * h_Na) / (kna3 * h_1_ss);
  const ML_CalcType h_3_ss       = 1.0 / h_1_ss;
  const ML_CalcType h_4_ss       = 1.0 + (Na_ss * (1.0 + Na_ss / v(VT_k_Na2))) / v(VT_k_Na1);
  const ML_CalcType h_5_ss       = (Na_ss * Na_ss) / (h_4_ss * v(VT_k_Na1) * v(VT_k_Na2));
  const ML_CalcType h_6_ss       = 1.0 / h_4_ss;
  const ML_CalcType k_4_d_ss     = (h_3_ss * wca) / h_Ca;
  const ML_CalcType k_4_dd_ss    = h_2_ss * wnaca;
  const ML_CalcType k_4_ss       = k_4_d_ss + k_4_dd_ss;
  const ML_CalcType k_6_ss       = h_6_ss * Ca_ss * v(VT_k_Ca_on);
  const ML_CalcType k_7_ss       = h_5_ss * h_2_ss * wna;
  const ML_CalcType x_1_ss       = v(VT_k_2) * k_4_ss * (k_7_ss + k_6_ss) + v(VT_k_5) * k_7_ss * (v(VT_k_2) + k_3);
  const ML_CalcType x_2_ss       = v(VT_k_1) * k_7_ss * (k_4_ss + v(VT_k_5)) + k_4_ss * k_6_ss * (v(VT_k_1) + k_8);
  const ML_CalcType x_3_ss       = v(VT_k_1) * k_3 * (k_7_ss + k_6_ss) + k_8 * k_6_ss * (v(VT_k_2) + k_3);
  const ML_CalcType x_4_ss       = v(VT_k_2) * k_8 * (k_4_ss + v(VT_k_5)) + k_3 * v(VT_k_5) * (v(VT_k_1) + k_8);
  const ML_CalcType E_1_ss       = x_1_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_2_ss       = x_2_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_3_ss       = x_3_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_4_ss       = x_4_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType allo_ss      = 1.0 / (1.0 + (KmCaAct / Ca_ss) * (KmCaAct / Ca_ss));
  const ML_CalcType J_NaCa_Na_ss = 3.0 * (E_4_ss * k_7_ss - E_1_ss * k_8) + E_3_ss * k_4_dd_ss - E_2_ss * k_3_dd;
  const ML_CalcType J_NaCa_Ca_ss = E_2_ss * v(VT_k_2) - E_1_ss * v(VT_k_1);
  I_NaCa_ss = v(VT_INaCa_Multiplier)
#ifdef HF
    * v(VT_HF_Gomez_6)
#endif // ifdef HF
    * G_NaCa * 0.2 * allo_ss *
    (z_Na * J_NaCa_Na_ss + z_Ca * J_NaCa_Ca_ss);

  /////////////////////////////////////////////////////////
  /// calculate I_NaK
  /////////////////////////////////////////////////////////
  double k1p                = 949.5;
  double k2m                = 39.4;
  double k3p                = 1899.0;
  double k3m                = 79300.0;
  double k4m                = 40.0;

  // PKA-P
  double Knai0_NP           = 9.073;
  double Knai0_PKA          = 0.7 * 9.073;
  double Knai0              = 9.073;
  const ML_CalcType phi_INaK_PKA = v(VT_INaK_P_frac);
  if (phi_INaK_PKA == 0) {
      double Knai0 = 9.073;
  } else {
      double Knai0 = (1-phi_INaK_PKA)*Knai0_NP + phi_INaK_PKA*Knai0_PKA;
  }

  double Knao0              = 27.78;
  double delta              = -0.1550;
  const ML_CalcType K_Nai   = Knai0 * exp(delta * VFoverRT * 1/3);
  const ML_CalcType K_Nao   = Knao0 * exp((1.0 - delta) * VFoverRT * 1/3);
  double Kki                = 0.5;
  double Kko                = 0.3582;
  double H                  = 1.0e-7;
  double eP                 = 4.2;
  double Khp                = 1.698e-7;
  double Knap               = 224.0;
  double Kxkur              = 292.0;
  const ML_CalcType P       = eP / (1.0 + H / Khp + Na_i / Knap + K_i / Kxkur);
  const ML_CalcType alpha_1 = (k1p * pow((Na_i / K_Nai), 3.0)) / (pow((1.0 + Na_i / K_Nai), 3.0) + pow((1.0 + K_i / Kki), 2.0) - 1.0);
  const ML_CalcType beta_2 = (k2m * pow((v(VT_Na_o) / K_Nao), 3.0)) / (pow((1.0 + v(VT_Na_o) / K_Nao), 3.0) + pow((1.0 + v(VT_K_o) / Kko), 2.0) - 1.0);
  const ML_CalcType alpha_3 = (k3p * pow((v(VT_K_o) / Kko), 2.0)) / (pow((1.0 + v(VT_Na_o) / K_Nao), 3.0) + pow((1.0 + v(VT_K_o) / Kko), 2.0) - 1.0);
  const ML_CalcType beta_3 = (k3m * P * H) / (1.0 + v(VT_MgATP));
  const ML_CalcType beta_4 = (k4m * pow((K_i / Kki), 2.0)) / (pow((1.0 + Na_i / K_Nai), 3.0) + pow((1.0 + K_i / Kki), 2.0) - 1.0);
  const ML_CalcType x_1 = v(VT_alpha_4) * alpha_1 * v(VT_alpha_2) + beta_2 * beta_4 * beta_3 + v(VT_alpha_2) * beta_4 * beta_3 + beta_3 * alpha_1 * v(VT_alpha_2);
  const ML_CalcType x_2 = beta_2 * v(VT_beta_1) * beta_4 + alpha_1 * v(VT_alpha_2) * alpha_3 + alpha_3 * v(VT_beta_1) * beta_4 + v(VT_alpha_2) * alpha_3 * beta_4;
  const ML_CalcType x_3 = v(VT_alpha_2) * alpha_3 * v(VT_alpha_4) + beta_3 * beta_2 * v(VT_beta_1) + beta_2 * v(VT_beta_1) * v(VT_alpha_4) + alpha_3 * v(VT_alpha_4) * v(VT_beta_1);
  const ML_CalcType x_4 = beta_4 * beta_3 * beta_2 + alpha_3 * v(VT_alpha_4) * alpha_1 + beta_2 * v(VT_alpha_4) * alpha_1 + beta_3 * beta_2 * alpha_1;
  const ML_CalcType E_1      = x_1 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_2      = x_2 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_3      = x_3 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_4      = x_4 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType J_NaK_Na = 3.0 * (E_1 * alpha_3 - E_2 * beta_3);
  const ML_CalcType J_NaK_K  = 2.0 * (E_4 * v(VT_beta_1) - E_3 * alpha_1);
  double P_NaK               = 30.0;
  I_NaK = v(VT_INaK_Multiplier)
#ifdef HF
    * v(VT_HF_Gomez_7)
#endif // ifdef HF
    * P_NaK * (z_Na * J_NaK_Na + z_K * J_NaK_K);

  /////////////////////////////////////////////////////////
  /// calculate I_Kb
  /////////////////////////////////////////////////////////
  const ML_CalcType x_Kb = 1.0 / (1.0 + exp(-(V_m - 14.48) / (18.34)));
  double G_Kb            = 0.003;
  double G_Kb_PKA = 2.5 * G_Kb;
  const ML_CalcType phi_IKb_PKA = v(VT_IKb_P_frac);
  I_Kb_NP = G_Kb * x_Kb * (V_m - E_K);
  I_Kb_PKA = G_Kb_PKA * x_Kb * (V_m - E_K);
  I_Kb = v(VT_IKb_Multiplier) * ((1 - phi_IKb_PKA) * I_Kb_NP + phi_IKb_PKA * I_Kb_PKA);

  /// calculate I_Nab
  double P_Nab = 3.75e-10;
  I_Nab = v(VT_INab_Multiplier) * P_Nab * z_Na * z_Na * VFFoverRT *
    ((Na_i * exp_z_Na - v(VT_Na_o)) / (exp_z_Na - 1.0));

  /// calculate I_Cab
  double P_Cab = 2.5e-8;
  I_Cab = v(VT_ICab_Multiplier) * P_Cab * z_Ca * z_Ca * VFFoverRT *
    ((gamma_Cai * Ca_i * exp_z_Ca - gamma_Cao * v(VT_Ca_o)) / (exp_z_Ca - 1.0));

  /// calculate I_pCa
  double G_pCa = 0.0005;
  I_pCa = v(VT_IpCa_Multiplier) * G_pCa * (Ca_i / (0.0005 + Ca_i));

  /// stretch activated channel
 #ifdef ISAC
  I_SAC = v(VT_G_sac) * (V_m - v(VT_E_sac)) / (1. + v(VT_K_sac) * exp(-v(VT_alpha) * (stretch - 1.)));
  I_SAC *= v(VT_ISAC_SWITCH);
 #endif  // ifdef ISAC

  /// I_tot
  const ML_CalcType cur_I_tot =
    I_Na_fast + I_Na_late + I_to + I_CaL + I_CaNa + I_CaK + I_Kr + I_Ks + I_K1 + I_NaCa_i + I_NaCa_ss + I_NaK +
    I_Nab + I_Cab + I_Kb + I_pCa
  #ifdef ISAC
    + I_SAC
  #endif  // ifdef ISAC
    + i_external;
  const ML_CalcType I_tot = -cur_I_tot;

  /// diffusion fluxes J_diff_Na, J_diff_Ca, J_diff_K
  double tau_diff_Na          = 2.0;
  double tau_diff_Ca          = 0.2;
  double tau_diff_K           = 2.0;
  const ML_CalcType J_diff_Na = (Na_ss - Na_i) / tau_diff_Na;
  const ML_CalcType J_diff_Ca = (Ca_ss - Ca_i) / tau_diff_Ca;
  const ML_CalcType J_diff_K  = (K_ss - K_i) / tau_diff_K;

  /////////////////////////////////////////////////////////
  /// SR Calcuim release flux J_rel
  /////////////////////////////////////////////////////////
  ML_CalcType J_rel_NP_inf = (v(VT_alpha_rel) * (-I_CaL)) / (1.0 + pow((1.5
  #ifdef HF
                                                                        * v(VT_HF_Gomez_8) // 0.8
  #endif  // ifdef HF
                                                                        / Ca_jsr), 8.0));
  ML_CalcType J_rel_CaMK_inf = (v(VT_alpha_rel_CaMK) * (-I_CaL)) / (1.0 + pow((1.5
  #ifdef HF
                                                                               * v(VT_HF_Gomez_8) // 0.8
  #endif  // ifdef HF
                                                                               / Ca_jsr), 8.0));
  ML_CalcType J_rel_PKA_inf = (v(VT_alpha_rel_PKA) * (-I_CaL)) / (1.0 + pow((1.5
  #ifdef HF
                                                                               * v(VT_HF_Gomez_8) // 0.8
  #endif  // ifdef HF
                                                                               / Ca_jsr), 8.0));
  ML_CalcType J_rel_both_inf = (v(VT_alpha_rel_both) * (-I_CaL)) / (1.0 + pow((1.5
  #ifdef HF
                                                                               * v(VT_HF_Gomez_8) // 0.8
  #endif  // ifdef HF
                                                                                / Ca_jsr), 8.0));
    
  if (v(VT_celltype) == 2.0) {
    J_rel_NP_inf   *= 1.7;
    J_rel_CaMK_inf *= 1.7;
    J_rel_PKA_inf *= 1.7;
    J_rel_both_inf *= 1.7;
  }
  const ML_CalcType tau_rel_NP_b = v(VT_beta_tau) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_NP   = tau_rel_NP_b < 0.001 ? 0.001 : tau_rel_NP_b;
  J_rel_NP = J_rel_NP_inf - (J_rel_NP_inf - J_rel_NP) * exp(-tinc / tau_rel_NP);
    
  const ML_CalcType tau_rel_CaMK_b = v(VT_beta_tau_CaMK) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_CaMK   = tau_rel_CaMK_b < 0.001 ? 0.001 : tau_rel_CaMK_b;
  J_rel_CaMK = J_rel_CaMK_inf - (J_rel_CaMK_inf - J_rel_CaMK) * exp(-tinc / tau_rel_CaMK);
  
  const ML_CalcType tau_rel_PKA_b = 0.75 * tau_rel_NP;
  const ML_CalcType tau_rel_PKA   = tau_rel_PKA_b < 0.001 ? 0.001 : tau_rel_PKA_b;
  J_rel_PKA = J_rel_PKA_inf - (J_rel_PKA_inf - J_rel_PKA) * exp(-tinc / tau_rel_PKA);
      
  const ML_CalcType tau_rel_both_b = 0.75 * tau_rel_CaMK;
  const ML_CalcType tau_rel_both   = tau_rel_both_b < 0.001 ? 0.001 : tau_rel_both_b;
  J_rel_both = J_rel_both_inf - (J_rel_both_inf - J_rel_both) * exp(-tinc / tau_rel_both);
  
  // fraction of phosphorylation
  const ML_CalcType phi_rel_CaMK = phi_INa_CaMK;
  const ML_CalcType phi_rel_PKA = v(VT_RyR_P_frac);
  const ML_CalcType phi_rel_both = phi_rel_CaMK * phi_rel_PKA;
  const ML_CalcType phi_rel_CaMKonly = phi_rel_CaMK - phi_rel_both;
  const ML_CalcType phi_rel_PKAonly = phi_rel_PKA - phi_rel_both;
    
  // currents
  const ML_CalcType J_rel = ((1.0-phi_rel_CaMKonly-phi_rel_PKAonly-phi_rel_both)*J_rel_NP + phi_rel_CaMKonly*J_rel_CaMK + phi_rel_PKAonly*J_rel_PKA + phi_rel_both*J_rel_both);
  
  /////////////////////////////////////////////////////////
  /// J_up
  /////////////////////////////////////////////////////////
  ML_CalcType J_up_NP   = (0.004375 * Ca_i) / (0.00092 + Ca_i);
  double DJ_up_CaMK     = 1.75;
  double DK_m_PLB       = 0.00017;
  ML_CalcType J_up_CaMK = (1.0 + DJ_up_CaMK) * (0.004375 * Ca_i) / (0.00092 - DK_m_PLB + Ca_i);
  ML_CalcType J_up_PKA = (0.004375 * Ca_i) / (0.54 * 0.00092 + Ca_i);;
  ML_CalcType J_up_both = (1.0 + DJ_up_CaMK) * (0.004375 * Ca_i) / (0.54 * (0.00092 - DK_m_PLB) + Ca_i);;
    
  if (v(VT_celltype) == 1.0) {
    J_up_NP   *= 1.3;
    J_up_CaMK *= 1.3;
    J_up_PKA   *= 1.3;
    J_up_both *= 1.3;
  }
    
  // fraction of phosphorylation
  const ML_CalcType phi_up_CaMK = phi_INa_CaMK;
  const ML_CalcType phi_up_PKA = v(VT_SERCA_P_frac);
  const ML_CalcType phi_up_both = phi_up_CaMK * phi_up_PKA;
  const ML_CalcType phi_up_CaMKonly = phi_up_CaMK - phi_up_both;
  const ML_CalcType phi_up_PKAonly = phi_up_PKA - phi_up_both;
    
  const ML_CalcType J_leak      =
  #ifdef HF
    v(VT_HF_Gomez_9) * // 1.3
  #endif  // ifdef HF
    (0.0039375 * Ca_nsr) / 15.0;

  const ML_CalcType J_up =
#ifdef HF
    v(VT_HF_Gomez_0) *
#endif // ifdef HF
    ((1.0-phi_up_CaMKonly-phi_up_PKAonly-phi_up_both)*J_up_NP + phi_up_CaMKonly*J_up_CaMK + phi_up_PKAonly*J_up_PKA + phi_up_both*J_up_both) - J_leak;


  /// J_tr
  double tau_tr          = 100.0;
  const ML_CalcType J_tr = (Ca_nsr - Ca_jsr) / tau_tr;

  /// Concentrations
  const ML_CalcType cur_Nai = I_Na_fast + I_Na_late + 3.0 * I_NaCa_i + 3.0 * I_NaK + I_Nab;
  const ML_CalcType dNa_i = -(cur_Nai) * util_1 + J_diff_Na * util_2;
  Na_i += tinc * dNa_i;
  const ML_CalcType dNa_ss = -(I_CaNa + 3.0 * I_NaCa_ss) * util_3 - J_diff_Na;
  Na_ss += tinc * dNa_ss;

  const ML_CalcType cur_Ki = I_to + I_Kr + I_Ks + I_K1 + I_Kb + i_external - 2.0 * I_NaK;
  const ML_CalcType dK_i = -(cur_Ki) * util_1 + J_diff_K * util_2;
  K_i += tinc * dK_i;
  const ML_CalcType dK_ss = -(I_CaK * util_3) - J_diff_K;
  K_ss += tinc * dK_ss;

  /// Ca buffer constants
  double cmdnmax = 0.05;
  if (v(VT_celltype) == 1) {
    cmdnmax *= 1.3;
  }
  double kmcmdn  = 0.00238;
  double trpnmax = 0.07;
  double kmtrpn_NP  = 0.0005;
  double kmtrpn_PKA = 1.6 * kmtrpn_NP;
  double BSRmax  = 0.047;
  double KmBSR   = 0.00087;
  double BSLmax  = 1.124;
  double KmBSL   = 0.0087;
  double csqnmax = 10.0;
  double kmcsqn  = 0.8;
    
  // PKA-P
  double kmtrpn  = 0.0005;
  const ML_CalcType phi_TRPN_PKA = v(VT_TnI_P_frac);
  if (phi_TRPN_PKA == 0) {
      double kmtrpn  = 0.0005;
  } else {
      double kmtrpn  = (1 - phi_TRPN_PKA) * kmtrpn_NP + phi_TRPN_PKA * kmtrpn_PKA;
  }
    
    
  #ifdef TRPN
  const ML_CalcType Ca_T50 = (v(VT_Ca_T50) + v(VT_beta1) * (lambda_m - 1.0))
    # ifdef HFM
    * v(VT_HFM_Multiplier)  // Change values for HF 0.6
    # endif  // ifdef HF
  ;
  const ML_CalcType dCaTRPN = v(VT_k_TRPN)*(pow(1000*Ca_i/Ca_T50, v(VT_n_TRPN))*(1.-Ca_TRPN)-Ca_TRPN);// needs muM with Land parameters
  Ca_TRPN += tinc * dCaTRPN;

  const ML_CalcType I_Trpn = trpnmax * dCaTRPN;
  #endif  // ifdef TRPN

  const ML_CalcType beta_Cai = 1.0 / (1.0 + ((cmdnmax * kmcmdn) / pow((kmcmdn + Ca_i), 2.0))
  #ifndef TRPN
                                      + ((trpnmax * kmtrpn) / pow((kmtrpn + Ca_i), 2.0))
  #endif  // ifndef TRPN
                                      );
  const ML_CalcType cur_Cai = I_pCa + I_Cab - 2.0 * I_NaCa_i;

  const ML_CalcType dCa_i = beta_Cai *
    (-(cur_Cai) * 0.5 * util_1 - J_up * (v(VT_v_nsr) / v(VT_v_myo)) + J_diff_Ca * util_2
  #ifdef TRPN
     - I_Trpn
  #endif  // ifdef TRPN
    );
  Ca_i += tinc * dCa_i;

  const ML_CalcType beta_Cass = 1.0 /
    (1.0 +
     ((BSRmax * KmBSR) /
      pow((KmBSR + Ca_ss), 2.0)) + ((BSLmax * KmBSL) / pow((KmBSL + Ca_ss), 2.0)));
  const ML_CalcType dCa_ss = beta_Cass *
    (-(I_CaL - 2.0 * I_NaCa_ss) * 0.5 * util_3 + J_rel * (v(VT_v_jsr) / v(VT_v_ss)) - J_diff_Ca);
  Ca_ss += tinc * dCa_ss;
  const ML_CalcType dCa_nsr = J_up - (J_tr * (v(VT_v_jsr) / v(VT_v_nsr)));
  Ca_nsr += tinc * dCa_nsr;
  const ML_CalcType beta_Cajsr = 1.0 / (1.0 + ((csqnmax * kmcsqn) / pow((kmcsqn + Ca_jsr), 2.0)));
  const ML_CalcType dCa_jsr    = beta_Cajsr * (J_tr - J_rel);
  Ca_jsr += tinc * dCa_jsr;


  return 0.001 * tinc * I_tot;
}  // OHaraRudyIso::Calc

void OHaraRudyIso::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr << tArg << ' ' << V << ' ' <<
    m << ' '  << h_fast << ' ' << h_slow << ' ' << j << ' ' << h_CaMK_slow << ' ' << j_CaMK << ' ' << h_PKA_fast << ' ' << h_PKA_slow << ' ' << j_PKA << ' ' << h_both_fast << ' ' << h_both_slow << ' ' << j_both << ' ' << m_L << ' ' <<
    h_L << ' ' << h_L_CaMK << ' ' << a << ' ' <<
    i_fast << ' ' << i_slow << ' ' << a_CaMK << ' ' << i_CaMK_fast << ' ' << i_CaMK_slow << ' ' << d << ' ' << f_fast <<
    ' ' << f_slow << ' ' << f_Ca_fast << ' ' << f_Ca_slow << ' ' <<
    j_Ca << ' ' << f_CaMK_fast << ' ' << f_Ca_CaMK_fast << ' ' << d_PKA << ' ' << f_PKA_fast << ' ' << f_PKA_slow << ' ' << f_Ca_PKA_fast << ' ' << f_Ca_PKA_slow << ' ' << f_both_fast << ' ' << f_Ca_both_fast << ' ' << n << ' ' << x_r_fast << ' ' << x_r_slow << ' ' <<
    x_s1 << ' ' << x_s2 << ' ' << x_s1_PKA << ' ' << x_s2_PKA << ' ' << x_K1 << ' ' <<
    Na_i << ' ' << Na_ss << ' ' << K_i << ' ' << K_ss << ' ' << Ca_i << ' ' << Ca_ss << ' ' << Ca_nsr << ' ' <<
    Ca_jsr << ' ' << CaMK_trap << ' ' << J_rel_NP << ' ' << J_rel_CaMK << ' ' << J_rel_PKA << ' ' << J_rel_both << ' '
    #ifdef TRPN
    << 0.07*Ca_TRPN << ' ' // export in mM
    #endif  // ifdef TRPN
  ;
}

void OHaraRudyIso::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);

  tempstr << I_Na_fast << ' ' << I_Na_late << ' ' << I_to << ' ' << I_CaL << ' ' << I_CaNa << ' ' << I_CaK << ' ' <<
    I_Kr << ' ' << I_Ks << ' ' << I_K1 << ' ' << I_NaCa_i << ' ' << I_NaCa_ss << ' ' << I_NaK << ' ' << I_Nab << ' ' <<
    I_Cab << ' ' << I_Kb << ' ' << I_pCa << ' '
  #ifdef ISAC
    << I_SAC << ' '
  #endif // ifdef ISAC
  ;
}  // OHaraRudyIso::LongPrint

void OHaraRudyIso::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {
    "m",           "h_fast",             "h_slow",                  "j",                            "h_CaMK_slow",
    "j_CaMK", "h_PKA_fast", "h_PKA_slow", "j_PKA", "h_both_fast", "h_both_slow", "j_both",
    "m_L",         "h_L",
    "h_L_CaMK",    "a",                  "i_fast",                  "i_slow",                       "a_CaMK",
    "i_CaMK_fast", "i_CaMK_slow",        "d",                       "f_fast",                       "f_slow",
    "f_Ca_fast",
    "f_Ca_slow",
    "j_Ca",        "f_CaMK_fast",        "f_Ca_CaMK_fast", "d_PKA", "f_PKA_fast", "f_PKA_slow", "f_Ca_PKA_fast", "f_Ca_PKA_slow", "f_both_fast", "f_Ca_both_fast",          "n",                            "x_r_fast",
    "x_r_slow",
    "x_s1",
    "x_s2", "x_s1_PKA", "x_s2_PKA",        "x_K1",               "Na_i",                    "Na_ss",                        "K_i",
    "K_ss",
    "Ca_i",        "Ca_ss",
    "Ca_nsr",      "Ca_jsr",             "CaMK_trap",               "J_rel_NP",
    "J_rel_CaMK", "J_rel_PKA", "J_rel_both"
   #ifdef TRPN
    ,              "Ca_TRPN"
   #endif  // ifdef TRPN
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}  // OHaraRudyIso::GetParameterNames

void OHaraRudyIso::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {
    "I_Na_fast", "I_Na_late", "I_to",      "I_CaL",      "I_CaNa",      "I_CaK",      "I_Kr",      "I_Ks",
    "I_K1",
    "I_NaCa_i",
    "I_NaCa_ss",
    "I_NaK",     "I_Nab",     "I_Cab",     "I_Kb",       "I_pCa"
#ifdef ISAC
    , "I_SAC"
#endif // ifdef ISAC
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
