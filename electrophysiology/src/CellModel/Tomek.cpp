/* -------------------------------------------------------

   Tomek.cpp

   Ver. 1.0.0

   Created:       Jan-Erik Duhme (12.2023)
   Last modified: Jan-Erik Duhme (04.12.2023)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#include <Tomek.h>

Tomek::Tomek(TomekParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_TomekParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

Tomek::~Tomek() {}

#ifdef HETERO

inline bool Tomek::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool Tomek::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int Tomek::GetSize(void) {
  return sizeof(Tomek)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(TomekParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char Tomek::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void Tomek::Init() {
#if KADEBUG
  cerr << "#initializing Class: Tomek ... " << endl;
        #endif  // if KADEBUG

  m              = v(VT_Init_m);
  A_h            = v(VT_Init_A_h);
  B_h            = v(VT_Init_B_h);
  h              = v(VT_Init_h);
  j              = v(VT_Init_j);
  h_p            = v(VT_Init_h_p);
  j_p            = v(VT_Init_j_p);
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
  n_ss           = v(VT_Init_n_ss);
  n_i            = v(VT_Init_n_i);
  f_CaMK_fast    = v(VT_Init_f_CaMK_fast);
  f_Ca_CaMK_fast = v(VT_Init_f_Ca_CaMK_fast);
  C_0            = v(VT_Init_C_0);
  C_1            = v(VT_Init_C_1);
  C_2            = v(VT_Init_C_2);
  O              = v(VT_Init_O);
  I              = v(VT_Init_I);
  x_s1           = v(VT_Init_x_s1);
  x_s2           = v(VT_Init_x_s2);
  Na_i           = v(VT_Init_Na_i);
  Na_ss          = v(VT_Init_Na_ss);
  K_i            = v(VT_Init_K_i);
  K_ss           = v(VT_Init_K_ss);
  Ca_i           = v(VT_Init_Ca_i);
  Ca_ss          = v(VT_Init_Ca_ss);
  Ca_nsr         = v(VT_Init_Ca_nsr);
  Ca_jsr         = v(VT_Init_Ca_jsr);
  Cl_i           = v(VT_Init_Cl_i);
  CaMK_trap      = v(VT_Init_CaMK_trap);
  J_rel_NP       = v(VT_Init_J_rel_NP);
  J_rel_CaMK     = v(VT_Init_J_rel_CaMK);
  #ifdef TRPN
  Ca_TRPN = v(VT_Init_Ca_TRPN)/0.07;
  #endif  // ifdef TRPN
}  // Tomek::Init

ML_CalcType Tomek::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
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
  const ML_CalcType CaMK_bound  = ((CaMKo * (1.0 - CaMK_trap)) / (1.0 + (KmCaM / Ca_ss)));
  const ML_CalcType CaMK_active = (CaMK_bound + CaMK_trap);
  const ML_CalcType dCaMK_trap = ((aCaMK * CaMK_bound * (CaMK_bound + CaMK_trap)) -
    (bCaMK * CaMK_trap));
  CaMK_trap += tinc * dCaMK_trap;

  /// reversal potentials
  const ML_CalcType E_Na = v(VT_RToverF) * log(v(VT_Na_o) / Na_i);
  const ML_CalcType E_K  = v(VT_RToverF) * log(v(VT_K_o) / K_i);
  double PKNa            = 0.01833;
  const ML_CalcType E_Ks = v(VT_RToverF) * log((v(VT_K_o) + PKNa * v(VT_Na_o)) / (K_i + PKNa * Na_i));
  double zcl = 1.0;
  const ML_CalcType E_Cl = ((v(VT_R) * v(VT_T)) / (zcl * v(VT_F))) * log((Cl_i/v(VT_Cl_o)));

  // convenient shorthand calculations
  const ML_CalcType VFFoverRT = (V_m * v(VT_F) * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType VFoverRT  = (V_m * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType util_1    = v(VT_A_cap) / (v(VT_F) * v(VT_v_myo));
  const ML_CalcType util_2    = v(VT_v_ss) / v(VT_v_myo);
  const ML_CalcType util_3    = v(VT_A_cap) / (v(VT_F) * v(VT_v_ss));

  /// calculate I_Na
  const ML_CalcType m_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86)/9.03))*(1.0 + exp(-(V_m + 56.86)/9.03)));
  const ML_CalcType tau_m = 0.1292 * exp(-(((V_m + 45.79)/15.54) * ((V_m + 45.79)/15.54))) + 0.06487 * exp(-(((V_m - 4.823)/51.12) * ((V_m - 4.823)/51.12)));
  m = m_inf - (m_inf - m) * exp(-tinc / tau_m);
    
  if (V_m < -40.0) {
    A_h = (0.057*exp((-(V_m+80.0)/6.8)));
    B_h = ((2.7*exp((0.079*V_m)))+(3.1e5*exp((0.3485*V_m))));
    A_j = ((((-2.5428e4*exp((0.2444*V_m))) - (6.948e-6*exp((-0.04391*V_m))))*(V_m+37.78))/(1.0+exp((0.311*(V_m+79.23)))));
    B_j = ((0.02424*exp((-0.01052*V_m)))/(1.0+exp((-0.1378*(V_m+40.14)))));
  } else {
    A_h = 0.0;
    B_h = (0.77/(0.13*(1.+exp((-(V_m+10.66)/11.1)))));
    A_j = 0.0;
    B_j =((0.6*exp((0.057*V_m)))/(1.0+exp((-0.1*(V_m+32.0))))) ;
  }

  const ML_CalcType tau_h = 1.0/(A_h+B_h);
  const ML_CalcType h_inf = (1.0/((1.0+exp(((V_m+71.55)/7.43)))*(1.0+exp(((V_m+71.55)/7.43)))));
  h = h_inf - (h_inf - h) * exp(-tinc / tau_h);

  const ML_CalcType tau_j = 1.0/(A_j+B_j);
  const ML_CalcType j_inf = (1.0/((1.0+exp(((V_m+71.55)/7.43)))*(1.0+exp(((V_m+71.55)/7.43)))));;
  j = j_inf - (j_inf - j) * exp(-tinc / tau_j);

  const ML_CalcType h_p_inf = (1.0/((1.0+exp(((V_m+77.55)/7.43)))*(1.0+exp(((V_m+77.55)/7.43)))));
  h_p = h_p_inf - (h_p_inf - h_p) * exp(-tinc / tau_h);

  const ML_CalcType tau_j_p = 1.46 * tau_j;
  j_p = j_inf - (j_inf - j_p) * exp(-tinc / tau_j_p);

  const ML_CalcType phi_INa_CaMK = 1.0 / (1.0 + KmCaMK / CaMK_active);

  double G_Na = 11.7802;
  I_Na = v(VT_INaF_Multiplier) * G_Na * (V_m - E_Na) *m*m*m* ((1.0 - phi_INa_CaMK) * h * j + (phi_INa_CaMK * h_p * j_p));


  /// calculate I_NaL
  ML_CalcType thL_mod = 1.0;
  const ML_CalcType m_L_inf = 1.0 / (1.0 + exp(-(V_m + 42.85) / 5.264));
  const ML_CalcType tau_m_L = tau_m;
  m_L = m_L_inf - (m_L_inf - m_L) * exp(-tinc / tau_m_L);
  const ML_CalcType h_L_inf = 1.0 / (1.0 + exp((V_m + 87.61) / 7.488));
  h_L = h_L_inf - (h_L_inf - h_L) * exp(-tinc / (thL_mod * v(VT_tau_h_L)));
  const ML_CalcType h_L_CaMK_inf = 1.0 / (1.0 + exp((V_m + 93.81) / 7.488));
  const ML_CalcType tau_h_L_CaMK = 3.0 * v(VT_tau_h_L) * thL_mod;
  h_L_CaMK = h_L_CaMK_inf - (h_L_CaMK_inf - h_L_CaMK) * exp(-tinc / tau_h_L_CaMK);
  const ML_CalcType phi_INaL_CaMK = phi_INa_CaMK;
  double G_Na_L                   = 0.0279;
  I_Na_late = v(VT_INaL_Multiplier) * G_Na_L * (V_m - E_Na) * m_L * ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_CaMK);

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
    (0.001416 * exp(-(V_m + 96.52) / 59.05) + 1.780e-08 * exp((V_m + 114.1) / 8.079));
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
  double G_to                    = 0.16;
  I_to = v(VT_Ito_Multiplier) * G_to * (V_m - E_K) * ((1.0 - phi_Ito_CaMK) * a * i + phi_Ito_CaMK * a_CaMK * i_CaMK);

  /// calculate I_CaL, I_CaNa, I_CaK
  ML_CalcType d_inf = 1.0763*exp((-1.007*exp((-0.0829*V_m))));
  if (V_m > 31.4978)
  {
    d_inf = 1.0;
  }
  const ML_CalcType tau_d = 0.6 + 1.0 / (exp(-0.05 * (V_m + 6.0)) + exp(0.09 * (V_m + 14.0)));
  d = d_inf - (d_inf -d) * exp(-tinc / tau_d);
  const ML_CalcType f_inf = 1.0 / (1.0 + exp((V_m + 19.58) / 3.696));
  const ML_CalcType tau_f_fast = 7.0 + 1.0 / (0.0045 * exp(-(V_m + 20.0) / 10.0) + 0.0045 * exp((V_m + 20.0) / 10.0));
  const ML_CalcType tau_f_slow = 1000.0 + 1.0 / (3.5e-05 * exp(-(V_m + 5.0) / 4.0) + 3.5e-05 * exp((V_m + 5.0) / 6.0));
  const ML_CalcType A_f_Ca_fast = 0.3 + 0.6 / (1.0 + exp((V_m - 10.0) / 10.0));
  const ML_CalcType A_f_Ca_slow = 1.0 - A_f_Ca_fast;
  f_fast = f_inf - (f_inf - f_fast) * exp(-tinc / tau_f_fast);
  f_slow = f_inf - (f_inf - f_slow) * exp(-tinc / tau_f_slow);
  const ML_CalcType f             = v(VT_A_f_fast) * f_fast + v(VT_A_f_slow) * f_slow;
  const ML_CalcType f_Ca_inf      = f_inf;
  const ML_CalcType tau_f_Ca_fast = 7.0 + 1.0 / (0.04 * exp(-(V_m - 4.0) / 7.0) + 0.04 * exp((V_m - 4.0) / 7.0));
  const ML_CalcType tau_f_Ca_slow = 100.0 + 1.0 / (0.00012 * exp(-V_m / 3.0) + 0.00012 * exp(V_m /7.0));
  
  f_Ca_fast = f_Ca_inf - (f_Ca_inf - f_Ca_fast) * exp(-tinc / tau_f_Ca_fast);
  f_Ca_slow = f_Ca_inf - (f_Ca_inf - f_Ca_slow) * exp(-tinc / tau_f_Ca_slow);
  const ML_CalcType f_Ca        = A_f_Ca_fast * f_Ca_fast + A_f_Ca_slow * f_Ca_slow;
  const ML_CalcType j_Ca_inf = 1. / (1 + exp((V_m + 18.08) / (2.7916)));
  double tau_j_Ca               = 75.0;
  j_Ca = j_Ca_inf - (j_Ca_inf - j_Ca) * exp(-tinc / tau_j_Ca);

  const ML_CalcType tau_f_CaMK_fast = 2.5 * tau_f_fast;
  const ML_CalcType f_CaMK_inf      = f_inf;
  f_CaMK_fast = f_CaMK_inf - (f_CaMK_inf - f_CaMK_fast) * exp(-tinc / tau_f_CaMK_fast);

  const ML_CalcType f_CaMK             = v(VT_A_f_fast) * f_CaMK_fast + v(VT_A_f_slow) * f_slow;

  const ML_CalcType tau_f_Ca_CaMK_fast = 2.5 * tau_f_Ca_fast;
  const ML_CalcType f_Ca_CaMK_inf      = f_inf;
  f_Ca_CaMK_fast = f_Ca_CaMK_inf - (f_Ca_CaMK_inf - f_Ca_CaMK_fast) * exp(-tinc / tau_f_Ca_CaMK_fast);

  const ML_CalcType f_Ca_CaMK = A_f_Ca_fast * f_Ca_CaMK_fast + A_f_Ca_slow * f_Ca_slow;
  const ML_CalcType k_m2_n    = j_Ca * 1.0;
  double Kmn                  = 0.002;
  double k2n                  = 500.0;
  const ML_CalcType alpha_n_Ca_ss   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_ss)), 4.0));
  n_ss = alpha_n_Ca_ss * (k2n / k_m2_n) - (alpha_n_Ca_ss * (k2n / k_m2_n) - n_ss) * exp(-k_m2_n * tinc);

  const ML_CalcType I_o = (0.5*(v(VT_Na_o)+v(VT_K_o)+v(VT_Cl_o)+(4.*v(VT_Ca_o)))/1000.);
  const ML_CalcType I_ss = ((0.5*(Na_ss+K_ss+Cl_i+(4.*Ca_ss)))/1000.);
  double constA = (1.82e6*pow((74.*310.),-1.5)); // Diel constant and temperature as constants

  const ML_CalcType gamma_Ca_ss = exp((-constA*4.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
  const ML_CalcType gamma_Ca_o = exp((-constA*4.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
  const ML_CalcType gamma_Na_ss = exp((-constA*1.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
  const ML_CalcType gamma_Na_o = exp((-constA*1.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
  const ML_CalcType gamma_K_ss = exp((-constA*1.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
  const ML_CalcType gamma_K_o = exp((-constA*1.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
  const ML_CalcType Psi_Ca_ss = ((4.*VFFoverRT*((gamma_Ca_ss*Ca_ss*exp((2.*VFoverRT))) - (gamma_Ca_o*v(VT_Ca_o))))/(exp((2.*VFoverRT)) - 1.));
  const ML_CalcType Psi_CaNa_ss = ((1.*VFFoverRT*((gamma_Na_ss*Na_ss*exp((1.*VFoverRT))) - (gamma_Na_o)*v(VT_Na_o)))/(exp((1.*VFoverRT)) - 1.));
  const ML_CalcType Psi_CaK_ss = ((1.*VFFoverRT*((gamma_K_ss*K_ss*exp((1.*VFoverRT))) - (gamma_K_o*v(VT_K_o))))/(exp((1.*VFoverRT)) - 1.));

  double PCa_b = 8.3757e-05;
  const ML_CalcType PCa = v(VT_IpCa_Multiplier) * PCa_b;
  const ML_CalcType PCaNa = (0.00125*PCa);
  const ML_CalcType PCaK = (3.574e-4*PCa);
  const ML_CalcType PCap = (1.1*PCa);
  const ML_CalcType PCaNap = (0.00125*PCap);
  const ML_CalcType PCaKp = (3.574e-4*PCap);

  double ICaL_fractionSS = 0.8;
  const ML_CalcType phi_ICaL_CaMK = 1.0 / (1.0 + KmCaMK / CaMK_active);

  I_CaL_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCa*Psi_Ca_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCap*Psi_Ca_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));
  I_CaNa_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCaNa*Psi_CaNa_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCaNap*Psi_CaNa_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));
  I_CaK_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCaK*Psi_CaK_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCaKp*Psi_CaK_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));


  const ML_CalcType alpha_n_Ca_i   = 1.0 / ((k2n / k_m2_n) + ((1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))));
  n_i = alpha_n_Ca_i * (k2n / k_m2_n) - (alpha_n_Ca_i * (k2n / k_m2_n) - n_i) * exp(-k_m2_n * tinc);


  const ML_CalcType I_i = ((0.5*(Na_i+K_i+Cl_i+(4.*Ca_i)))/1000.);
  const ML_CalcType gamma_Ca_i = exp((-constA*4.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
  const ML_CalcType gamma_Na_i = exp((-constA*1.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
  const ML_CalcType gamma_K_i = exp((-constA*1.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
  const ML_CalcType Psi_Ca_i = ((4.*VFFoverRT*((gamma_Ca_i*Ca_i*exp((2.*VFoverRT))) - (gamma_Ca_o*v(VT_Ca_o))))/(exp((2.*VFoverRT)) - 1.));
  const ML_CalcType Psi_CaNa_i = (1.*VFFoverRT*((gamma_Na_i*Na_i*exp((1.*VFoverRT))) - (gamma_Na_o*v(VT_Na_o)))/(exp((1.*VFoverRT)) - 1.));
  const ML_CalcType Psi_CaK_i = ((1.*VFFoverRT*((gamma_K_i*K_i*exp((1.*VFoverRT))) - (gamma_K_o*v(VT_K_o))))/(exp((1.*VFoverRT)) - 1.));

  I_CaL_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCa*Psi_Ca_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCap*Psi_Ca_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));
  I_CaNa_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCaNa*Psi_CaNa_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCaNap*Psi_CaNa_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));
  I_CaK_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCaK*Psi_CaK_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCaKp*Psi_CaK_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));


  I_CaL = I_CaL_ss + I_CaL_i;
  I_CaNa = I_CaNa_ss + I_CaNa_i;
  I_CaK = I_CaK_ss + I_CaK_i;


  /// calculate I_Kr
  const ML_CalcType alpha_Kr = 0.1161 * exp(0.299 * VFoverRT);
  const ML_CalcType beta_Kr = 0.2442 * exp(-1.604 * VFoverRT);
  double alpha_Kr_1 = 1.25 * 0.1235;
  double beta_Kr_1 = 0.1911;
  const ML_CalcType alpha_Kr_2 = 0.0578 * exp(0.971 * VFoverRT);
  const ML_CalcType beta_Kr_2 = 0.349e-3 * exp(-1.062 * VFoverRT);
  const ML_CalcType alpha_Kr_i = 0.2533 * exp(0.5953 * VFoverRT);
  const ML_CalcType beta_Kr_i = 0.06525 * exp(-0.8209 * VFoverRT);
  const ML_CalcType alpha_C2_to_I = 0.52e-4 * exp(1.525 * VFoverRT);
  const ML_CalcType beta_I_to_C2 = (beta_Kr_2 * beta_Kr_i * alpha_C2_to_I) / (alpha_Kr_2 * alpha_Kr_i);
  const ML_CalcType dC_0 = C_1 * beta_Kr - C_0 * alpha_Kr;
  C_0 += tinc * dC_0;
  const ML_CalcType dC_1 = C_0 * alpha_Kr + C_2 * beta_Kr_1 - C_1 * (beta_Kr + alpha_Kr_1);
  C_1 += tinc * dC_1;
  const ML_CalcType dC_2 = C_1 * alpha_Kr_1 + O * beta_Kr_2 + I * beta_I_to_C2 - C_2 * (beta_Kr_1 + alpha_Kr_2 + alpha_C2_to_I);
  C_2 += tinc * dC_2;
  const ML_CalcType dO = (((alpha_Kr_2 * C_2) + (beta_Kr_i * I)) - ((beta_Kr_2 + alpha_Kr_i) * O));
  O += tinc * dO;
  const ML_CalcType dI = (((alpha_C2_to_I * C_2) + (alpha_Kr_i * O)) - ((beta_I_to_C2 + beta_Kr_i) * I));
  I += tinc * dI;
  double G_Kr = 0.0321;

  I_Kr = v(VT_IKr_Multiplier) * G_Kr * sqrt(v(VT_K_o) / 5.0) * O * (V_m - E_K);



  /// calculate I_Ks
  const ML_CalcType x_s1_inf = 1.0/(1.0 + exp(-(V_m + 11.6) / 8.932));
  const ML_CalcType tau_x_s1 = 817.3 + 1.0 /
    (2.326e-04 * exp((V_m + 48.28) / 17.80)+ 0.001292 * exp(-(V_m + 210.0) / 230.0));
  x_s1 = x_s1_inf - (x_s1_inf - x_s1) * exp(-tinc / tau_x_s1);
  const ML_CalcType x_s2_inf = x_s1_inf;
  const ML_CalcType tau_x_s2 = 1.0 / (0.01 * exp((V_m - 50.0) / 20.0) + 0.0193 * exp(-(V_m + 66.54) / 31.0));
  x_s2 = x_s2_inf - (x_s2_inf - x_s2) * exp(-tinc / tau_x_s2);
  double G_Ks = 0.0011;
  I_Ks = v(VT_IKs_Multiplier) * G_Ks * (1.0 + (0.6 / (1.0 + pow((3.8e-05 / (Ca_i)), 1.4)))) * x_s1 * x_s2 * (V_m - E_Ks);

  /// calculating I_K1
  const ML_CalcType a_K1 = 4.094 / (1.0 + exp(0.1217 * (V_m - E_K - 49.934)));
  const ML_CalcType b_K1 = (15.72 * exp(0.0674 * (V_m - E_K - 3.257)) + exp(0.0618 * (V_m - E_K - 594.31))) / (1.0 + exp(-0.1629 * (V_m - E_K + 14.207)));
  const ML_CalcType K1_SS = (a_K1) / (a_K1 + b_K1);
  double G_K1            = 0.6992;
  I_K1 = v(VT_IK1_Multiplier) * G_K1 * sqrt(v(VT_K_o)/5.0) * K1_SS * (V_m - E_K);

  /// calculate I_NaCa_i
  double kna3                   = 88.12;
  double wna                    = 6.0e4;
  double wca                    = 6.0e4;
  double wnaca                  = 5.0e3;
  double qna                    = 0.5224;
  double qca                    = 0.1670;
  double z_Ca                   = 2.0;
  double z_Na                   = 1.0;
  double z_K                    = 1.0;
  const ML_CalcType exp_z_Ca    = exp(z_Ca * VFoverRT);
  const ML_CalcType exp_z_Na    = exp(z_Na * VFoverRT);
  const ML_CalcType exp_z_K     = exp(z_K * VFoverRT);
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
  const ML_CalcType allo_i      = 1.0 / (1.0 + ((KmCaAct / Ca_i) * (KmCaAct / Ca_i)));
  const ML_CalcType J_NaCa_Na_i = 3.0 * (E_4_i * k_7_i - E_1_i * k_8) + E_3_i * k_4_dd_i - E_2_i * k_3_dd;
  const ML_CalcType J_NaCa_Ca_i = E_2_i * v(VT_k_2) - E_1_i * v(VT_k_1);
  double G_NaCa                 = 0.0034;
  I_NaCa_i = v(VT_INaCa_Multiplier)
    * G_NaCa * 0.65 * allo_i * (z_Na * J_NaCa_Na_i + z_Ca * J_NaCa_Ca_i);

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
    * G_NaCa * 0.35 * allo_ss *
    (z_Na * J_NaCa_Na_ss + z_Ca * J_NaCa_Ca_ss);

  /// calculate I_NaK
  double k1p                = 949.5;
  double k2m                = 39.4;
  double k3p                = 1899.0;
  double k3m                = 79300.0;
  double k4m                = 40.0;
  double Knai0              = 9.073;
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
  const ML_CalcType alpha_1 =
    (k1p *
     pow((Na_i / K_Nai), 3.0)) / (pow((1.0 + Na_i / K_Nai), 3.0) + pow((1.0 + K_i / Kki), 2.0) - 1.0);
  const ML_CalcType beta_2 =
    (k2m *
     pow((v(VT_Na_o) / K_Nao),
         3.0)) / (pow((1.0 + v(VT_Na_o) / K_Nao), 3.0) + pow((1.0 + v(VT_K_o) / Kko), 2.0) - 1.0);
  const ML_CalcType alpha_3 =
    (k3p *
     pow((v(VT_K_o) / Kko),
         2.0)) / (pow((1.0 + v(VT_Na_o) / K_Nao), 3.0) + pow((1.0 + v(VT_K_o) / Kko), 2.0) - 1.0);
  const ML_CalcType beta_3 = (k3m * P * H) / (1.0 + v(VT_MgATP));
  const ML_CalcType beta_4 =
    (k4m *
     pow((K_i / Kki), 2.0)) / (pow((1.0 + Na_i / K_Nai), 3.0) + pow((1.0 + K_i / Kki), 2.0) - 1.0);
  const ML_CalcType x_1 = v(VT_alpha_4) * alpha_1 * v(VT_alpha_2) + beta_2 * beta_4 * beta_3 + v(VT_alpha_2) * beta_4 *
    beta_3 + beta_3 * alpha_1 * v(VT_alpha_2);
  const ML_CalcType x_2 = beta_2 * v(VT_beta_1) * beta_4 + alpha_1 * v(VT_alpha_2) * alpha_3 + alpha_3 * v(VT_beta_1) *
    beta_4 + v(VT_alpha_2) * alpha_3 * beta_4;
  const ML_CalcType x_3 = v(VT_alpha_2) * alpha_3 * v(VT_alpha_4) + beta_3 * beta_2 * v(VT_beta_1) + beta_2 * v(
    VT_beta_1) * v(VT_alpha_4) + alpha_3 * v(VT_alpha_4) * v(VT_beta_1);
  const ML_CalcType x_4 = beta_4 * beta_3 * beta_2 + alpha_3 * v(VT_alpha_4) * alpha_1 + beta_2 * v(VT_alpha_4) *
    alpha_1 + beta_3 * beta_2 * alpha_1;
  const ML_CalcType E_1      = x_1 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_2      = x_2 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_3      = x_3 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType E_4      = x_4 / (x_1 + x_2 + x_3 + x_4);
  const ML_CalcType J_NaK_Na = 3.0 * (E_1 * alpha_3 - E_2 * beta_3);
  const ML_CalcType J_NaK_K  = 2.0 * (E_4 * v(VT_beta_1) - E_3 * alpha_1);
  double P_NaK               = 15.4509;
  I_NaK = v(VT_INaK_Multiplier)
    * P_NaK * (z_Na * J_NaK_Na + z_K * J_NaK_K);

  /// Calculate I_CaCl
  I_CaCl_junc = (0.2843 * (V_m - E_Cl)) / (1. + (0.1/Ca_ss));
  I_CaCl_sl = (0 * 0.2843 * (V_m - E_Cl)) / (1. + (0.1/Ca_i));
  I_CaCl = I_CaCl_junc + I_CaCl_sl;

  /// calculate I_Kb
  const ML_CalcType x_Kb = 1.0 / (1.0 + exp(-(V_m - 10.8968) / (23.9871)));
  double G_Kb            = 0.0189;
  I_Kb = v(VT_IKb_Multiplier) * G_Kb * x_Kb * (V_m - E_K);

  /// calculate I_Nab
  double P_Nab = 1.9239e-9;
  I_Nab = v(VT_INab_Multiplier) * P_Nab * z_Na * z_Na * VFFoverRT *
    ((Na_i * exp_z_Na - v(VT_Na_o)) / (exp_z_Na - 1.0));

  /// calculate I_Cab
  double P_Cab = 5.9194e-8;
  I_Cab = v(VT_ICab_Multiplier) * P_Cab * z_Ca * z_Ca * VFFoverRT *
    ((gamma_Ca_i * Ca_i * exp_z_Ca - gamma_Ca_o * v(VT_Ca_o)) / (exp_z_Ca - 1.0));

  /// calculate I_Clb
  double P_Clb = 1.98e-3;
  I_Clb = P_Clb * (V_m - E_Cl);

  /// calculate I_pCa
  double G_pCa = 0.0005;
  I_pCa = v(VT_IpCa_Multiplier) * ((G_pCa * Ca_i) / (0.0005 + Ca_i));

  /// stretch activated channel
 #ifdef ISAC
  I_SAC = v(VT_G_sac) * (V_m - v(VT_E_sac)) / (1. + v(VT_K_sac) * exp(-v(VT_alpha) * (stretch - 1.)));
  I_SAC *= v(VT_ISAC_SWITCH);
 #endif  // ifdef ISAC

  /// I_tot
   cur_I_tot =
    I_Na + I_Na_late + I_to + I_CaL + I_CaNa + I_CaK + I_Kr + I_Ks + I_K1 + I_NaCa_i + I_NaCa_ss + I_NaK +
    I_Nab + I_Cab + I_Kb + I_pCa + I_Clb + I_CaCl 
  #ifdef ISAC
    + I_SAC
  #endif  // ifdef ISAC
    + i_external;
  I_tot = -cur_I_tot;

  /// diffusion fluxes J_diff_Na, J_diff_Ca, J_diff_K
  double tau_diff_Na          = 2.0;
  double tau_diff_Ca          = 0.2;
  double tau_diff_K           = 2.0;
  J_diff_Na = (Na_ss - Na_i) / tau_diff_Na;
  J_diff_Ca = (Ca_ss - Ca_i) / tau_diff_Ca;
  J_diff_K  = (K_ss - K_i) / tau_diff_K;

  /// SR Calcuim release flux J_rel
  double K_inf_rel = 1.7;
  ML_CalcType J_rel_NP_inf = (v(VT_alpha_rel) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
  ML_CalcType J_rel_CaMK_inf = (v(VT_alpha_rel_CaMK) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
  if (v(VT_celltype) == 2.0) {
    J_rel_NP_inf   *= 1.7;
    J_rel_CaMK_inf *= 1.7;
  }
  const ML_CalcType tau_rel_NP_b = v(VT_beta_tau) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_NP   = tau_rel_NP_b < 0.001 ? 0.001 : tau_rel_NP_b;
  J_rel_NP = J_rel_NP_inf - (J_rel_NP_inf - J_rel_NP) * exp(-tinc / tau_rel_NP);
  
  const ML_CalcType tau_rel_CaMK_b = v(VT_beta_tau_CaMK) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_CaMK   = tau_rel_CaMK_b < 0.001 ? 0.001 : tau_rel_CaMK_b;
  J_rel_CaMK = J_rel_CaMK_inf - (J_rel_CaMK_inf - J_rel_CaMK) * exp(-tinc / tau_rel_CaMK);
  
  const ML_CalcType phi_rel_CaMK = phi_INa_CaMK;
  J_rel        = 1.5378 * (1.0 - phi_rel_CaMK) * J_rel_NP + phi_rel_CaMK * J_rel_CaMK;

  /// J_up
  ML_CalcType J_up_NP   = (0.005425 * Ca_i) / (0.00092 + Ca_i);
  ML_CalcType J_up_CaMK = (2.75 * 0.005425 * Ca_i) / (0.00092 - 0.00017 + Ca_i);
  if (v(VT_celltype) == 1.0) {
    J_up_NP   *= 1.3;
    J_up_CaMK *= 1.3;
  }
  const ML_CalcType phi_up_CaMK = (1.0 /(1.0 + (KmCaMK / CaMK_active)));
  J_leak      = (0.0048825 * Ca_nsr) / 15.0;

  J_up = ((1.0 - phi_up_CaMK) * J_up_NP + (phi_up_CaMK * J_up_CaMK) - J_leak);


  /// J_tr
  double tau_tr          = 60.0;
  J_tr = (Ca_nsr - Ca_jsr) / tau_tr;

  /// Concentrations
  const ML_CalcType cur_Nai = I_Na + I_Na_late + 3.0 * I_NaCa_i + I_CaNa_i + 3.0 * I_NaK + I_Nab;
  const ML_CalcType dNa_i = -(cur_Nai) * util_1 + J_diff_Na * util_2;
  Na_i += tinc * dNa_i;
  const ML_CalcType dNa_ss = -(I_CaNa_ss + 3.0 * I_NaCa_ss) * util_3 - J_diff_Na;
  Na_ss += tinc * dNa_ss;

  //const ML_CalcType cur_Ki = I_to + I_Kr + I_Ks + I_K1 + I_Kb - 2.0 * I_NaK + I_CaK_i;
  const ML_CalcType cur_Ki = I_to + I_Kr + I_Ks + I_K1 + I_Kb + i_external - 2.0 * I_NaK + I_CaK_i;
  const ML_CalcType dK_i = -(cur_Ki) * util_1 + J_diff_K * util_2;
  K_i += tinc * dK_i;
  const ML_CalcType dK_ss = -(I_CaK_ss) * util_3 - J_diff_K;
  K_ss += tinc * dK_ss;

  /// Ca buffer constants
  double cmdnmax = 0.05;
  if (v(VT_celltype) == 1) {
    cmdnmax *= 1.3;
  }
  double kmcmdn  = 0.00238;
  double trpnmax = 0.07;
  double kmtrpn  = 0.0005;
  double BSRmax  = 0.047;
  double KmBSR   = 0.00087;
  double BSLmax  = 1.124;
  double KmBSL   = 0.0087;
  double csqnmax = 10.0;
  double kmcsqn  = 0.8;
  #ifdef TRPN
  const ML_CalcType Ca_T50 = (v(VT_Ca_T50) + v(VT_beta1) * (lambda_m - 1.0));
  const ML_CalcType dCaTRPN = v(VT_k_TRPN)*(pow(1000*Ca_i/Ca_T50, v(VT_n_TRPN))*(1.-Ca_TRPN)-Ca_TRPN);// needs muM with Land parameters
  Ca_TRPN += tinc * dCaTRPN;
  const ML_CalcType I_Trpn = trpnmax * dCaTRPN;
  #endif  // ifdef TRPN
  
  const ML_CalcType beta_Cai = 1.0 / (1.0 + ((cmdnmax * kmcmdn) / pow((kmcmdn + Ca_i), 2.0))
  #ifndef TRPN
                             + ((trpnmax * kmtrpn) / pow((kmtrpn + Ca_i), 2.0))
  #endif  // ifndef TRPN
                                      );
  const ML_CalcType cur_Cai = I_CaL_i + I_pCa + I_Cab - 2.0 * I_NaCa_i;
  

  //const ML_CalcType dCa_i = beta_Cai * (-(cur_Cai) * 0.5 * util_1 - J_up * (v(VT_v_nsr) / v(VT_v_myo)) + J_diff_Ca * util_2
  const ML_CalcType dCa_i = beta_Cai * (-(cur_Cai) * 0.5 * util_1 - J_up * (v(VT_v_nsr) / v(VT_v_myo)) + J_diff_Ca * util_2
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
    (-(I_CaL_ss - 2.0 * I_NaCa_ss) * 0.5 * util_3 + J_rel * (v(VT_v_jsr) / v(VT_v_ss)) - J_diff_Ca);
  Ca_ss += tinc * dCa_ss;
  const ML_CalcType dCa_nsr = J_up - ((J_tr * v(VT_v_jsr)) / v(VT_v_nsr));
  Ca_nsr += tinc * dCa_nsr;
  const ML_CalcType beta_Cajsr = 1.0 / (1.0 + ((csqnmax * kmcsqn) / pow((kmcsqn + Ca_jsr), 2.0)));
  const ML_CalcType dCa_jsr    = beta_Cajsr * (J_tr - J_rel);
  Ca_jsr += tinc * dCa_jsr;


  return 0.001 * tinc * I_tot;
}  // Tomek::Calc

void Tomek::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr << tArg << ' ' << V << ' ' <<
    m << ' '  << A_h << ' '  << B_h << ' '  << h << ' ' << A_j << ' '  << B_j << ' '<< j << ' ' << h_p << ' ' << j_p << ' ' << m_L << ' ' <<
    h_L << ' ' << h_L_CaMK << ' ' << a << ' ' <<
    i_fast << ' ' << i_slow << ' ' << a_CaMK << ' ' << i_CaMK_fast << ' ' << i_CaMK_slow << ' ' << d << ' ' << f_fast <<
    ' ' << f_slow << ' ' << f_Ca_fast << ' ' << f_Ca_slow << ' ' <<
    j_Ca << ' ' << f_CaMK_fast << ' ' << f_Ca_CaMK_fast << ' ' << n_ss << ' ' << n_i << ' ' << C_0 << ' ' << C_1 << ' ' << C_2 << ' ' << O << ' ' << I << ' ' <<
    x_s1 << ' ' << x_s2 <<' ' <<
    Na_i << ' ' << Na_ss << ' ' << K_i << ' ' << K_ss << ' ' << Ca_i << ' ' << Ca_ss << ' ' << Ca_nsr << ' ' <<
    Ca_jsr << ' ' << Cl_i << ' ' << CaMK_trap << ' ' << J_rel_NP << ' ' << J_rel_CaMK << ' '
    #ifdef TRPN
    << 0.07*Ca_TRPN << ' ' // export in mM
    #endif  // ifdef TRPN
  ;
}

void Tomek::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);

  tempstr << I_Na << ' ' << I_Na_late << ' ' << I_to << ' ' << I_CaL_i << ' ' << I_CaL_ss << ' ' << I_CaNa_i << ' ' << I_CaNa_ss << ' ' << I_CaK_i << ' ' << I_CaK_ss << ' ' << I_CaL << ' ' << I_CaNa << ' ' << I_CaK << ' ' <<
    I_Kr << ' ' << I_Ks << ' ' << I_K1 << ' ' << I_NaCa_i << ' ' << I_NaCa_ss << ' ' << I_NaK << ' ' << I_CaCl << ' ' << I_Nab << ' ' <<
    I_Cab << ' ' << I_Kb << ' ' << I_Clb << ' '<< I_pCa << ' '<< J_diff_Na << ' '<< J_diff_Ca << ' '<< J_diff_K << ' '<< J_leak << ' '<< J_rel << ' '<< J_tr << ' '<< J_up << ' '<< cur_I_tot << ' '
  #ifdef ISAC
    << I_SAC << ' '
  #endif // ifdef ISAC
  ;
}  // Tomek::LongPrint

void Tomek::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {
    "m", "A_h", "B_h",          "h",  "A_j", "B_j",                "j",                            "h_p",
    "j_p",
    "m_L",         "h_L",
    "h_L_CaMK",    "a",                  "i_fast",                  "i_slow",                       "a_CaMK",
    "i_CaMK_fast", "i_CaMK_slow",        "d",                       "f_fast",                       "f_slow",
    "f_Ca_fast",
    "f_Ca_slow",
    "j_Ca",        "f_CaMK_fast",        "f_Ca_CaMK_fast",          "n_ss",          "n_i", "C_0", "C_1", "C_2", "O", "I",
    "x_s1",
    "x_s2",               "Na_i",                    "Na_ss",                        "K_i",
    "K_ss",
    "Ca_i",        "Ca_ss",
    "Ca_nsr",      "Ca_jsr", "Cl_i",             "CaMK_trap",               "J_rel_NP",
    "J_rel_CaMK"
   #ifdef TRPN
    ,              "Ca_TRPN"
   #endif  // ifdef TRPN
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}  // Tomek::GetParameterNames

void Tomek::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {
    "I_Na", "I_Na_late", "I_to", "I_CaL_i", "I_CaL_ss", "I_CaNa_i", "I_CaNa_ss", "I_CaK_i", "I_CaK_ss",      "I_CaL",      "I_CaNa",      "I_CaK",      "I_Kr",      "I_Ks",
    "I_K1",
    "I_NaCa_i",
    "I_NaCa_ss",
    "I_NaK",     "I_CaCl",     "I_Nab",     "I_Cab",     "I_Kb",     "I_Clb",       "I_pCa", "J_diff_Na", "J_diff_Ca", "J_diff_K", "J_leak", "J_rel", "J_tr", "J_up", "cur_I_tot"
#ifdef ISAC
    , "I_SAC"
#endif // ifdef ISAC
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
