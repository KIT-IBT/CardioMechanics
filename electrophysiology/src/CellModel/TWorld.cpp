/*
 * File: TWorld.cpp
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


#include <TWorld.h>
#include <algorithm>

TWorld::TWorld(TWorldParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_TWorldParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

TWorld::~TWorld() {}

#ifdef HETERO

inline bool TWorld::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool TWorld::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int TWorld::GetSize(void) {
  return sizeof(TWorld)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(TWorldParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char TWorld::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void TWorld::Init() {
#if KADEBUG
  cerr << "#initializing Class: TWorld ... " << endl;
        #endif  // if KADEBUG
    
    // Concentrations
    Na_o = v(VT_Na_o);
    Na_dyad = v(VT_Init_Na_dyad);
    Na_sl = v(VT_Init_Na_sl);
    Na_myo = v(VT_Init_Na_myo);
    K_myo = v(VT_Init_K_myo);
    Cl_myo = v(VT_Init_Cl_myo);
    Ca_dyad = v(VT_Init_Ca_dyad);
    Ca_sl = v(VT_Init_Ca_sl);
    Ca_i = v(VT_Init_Ca_i);
    Ca_SR = v(VT_Init_Ca_SR);
    
    // CaMK and Ca Signalling
    CaMK_trap = v(VT_Init_CaMK_trap);
    CaMK_f_ICaL = v(VT_Init_CaMK_f_ICaL);
    CaMK_f_RyR = v(VT_Init_CaMK_f_RyR);
    CaMK_f_PLB = v(VT_Init_CaMK_f_PLB);
    casig_serca_trap = v(VT_Init_casig_serca_trap);
    
    // Sodium current (INa, INaL)
    m = v(VT_Init_m);
    h = v(VT_Init_h);
    j = v(VT_Init_j);
    h_p = v(VT_Init_h_p);
    j_p = v(VT_Init_j_p);
    m_PKA = v(VT_Init_m_PKA);
    h_PKA = v(VT_Init_h_PKA);
    j_PKA = v(VT_Init_j_PKA);
    h_both = v(VT_Init_h_both);
    j_both = v(VT_Init_j_both);
    m_L = v(VT_Init_m_L);
    h_L = v(VT_Init_h_L);
    h_L_p = v(VT_Init_h_L_p);
    
    // L-type calcium current (I_CaL, I_CaNa, I_CaK)
    d = v(VT_Init_d);
    f_fast = v(VT_Init_f_fast);
    f_slow = v(VT_Init_f_slow);
    f_Ca_fast = v(VT_Init_f_Ca_fast);
    f_Ca_slow = v(VT_Init_f_Ca_slow);
    j_Ca = v(VT_Init_j_Ca);
    f_p_fast = v(VT_Init_f_p_fast);
    f_Ca_p_fast = v(VT_Init_f_Ca_p_fast);
    d_PKA = v(VT_Init_d_PKA);
    f_PKA_fast = v(VT_Init_f_PKA_fast);
    f_PKA_slow = v(VT_Init_f_PKA_slow);
    f_Ca_PKA_fast = v(VT_Init_f_Ca_PKA_fast);
    f_Ca_PKA_slow = v(VT_Init_f_Ca_PKA_slow);
    f_both_fast = v(VT_Init_f_both_fast);
    f_Ca_both_fast = v(VT_Init_f_Ca_both_fast);
    n_Ca_dyad = v(VT_Init_n_Ca_dyad);
    n_Ca_sl = v(VT_Init_n_Ca_sl);
    I_CaL_pureCDI_dyad = v(VT_Init_I_CaL_pureCDI_dyad);
    I_CaL_pureCDI_sl = v(VT_Init_I_CaL_pureCDI_sl);
    
    // Transient outward current (Ito)
    a_slow = v(VT_Init_a_slow);
    a_fast = v(VT_Init_a_fast);
    i_slow = v(VT_Init_i_slow);
    i_fast = v(VT_Init_i_fast);
    a_p_slow = v(VT_Init_a_p_slow);
    a_p_fast = v(VT_Init_a_p_fast);
    i_p_slow = v(VT_Init_i_p_slow);
    i_p_fast = v(VT_Init_i_p_fast);
    
    // Rapid delayed rectifier current (IKr)
    C_0 = v(VT_Init_C_0);
    C_1 = v(VT_Init_C_1);
    C_2 = v(VT_Init_C_2);
    O = v(VT_Init_O);
    I = v(VT_Init_I);
    
    // Slow delayed rectifier current (IKs)
    xs_dyad = v(VT_Init_xs_dyad);
    xs_sl = v(VT_Init_xs_sl);
    
    // Calcium release from SR (Jrel, Jleak)
    J_rel_ICaLdep_act = v(VT_Init_J_rel_ICaLdep_act);
    J_rel_ICaLdep_f1 = v(VT_Init_J_rel_ICaLdep_f1);
    J_rel_ICaLdep_f2 = v(VT_Init_J_rel_ICaLdep_f2);
    ryr_R = v(VT_Init_ryr_R);
    ryr_O = v(VT_Init_ryr_O);
    ryr_I = v(VT_Init_ryr_I);
    ryr_CaRI = v(VT_Init_ryr_CaRI);
    ryr_R_p = v(VT_Init_ryr_R_p);
    ryr_O_p = v(VT_Init_ryr_O_p);
    ryr_I_p = v(VT_Init_ryr_I_p);
    ryr_CaRI_p = v(VT_Init_ryr_CaRI_p);
    
    // Buffering
    Buffer_NaBj = v(VT_Init_Buffer_NaBj);
    Buffer_NaBsl = v(VT_Init_Buffer_NaBsl);
    Buffer_TnClow = v(VT_Init_Buffer_TnClow);
    Buffer_TnCHc = v(VT_Init_Buffer_TnCHc);
    Buffer_TnCHm = v(VT_Init_Buffer_TnCHm);
    Buffer_CaM = v(VT_Init_Buffer_CaM);
    Buffer_Myosin_ca = v(VT_Init_Buffer_Myosin_ca);
    Buffer_Myosin_mg = v(VT_Init_Buffer_Myosin_mg);
    Buffer_SRB = v(VT_Init_Buffer_SRB);
    Buffer_SLLj = v(VT_Init_Buffer_SLLj);
    Buffer_SLLsl = v(VT_Init_Buffer_SLLsl);
    Buffer_SLHj = v(VT_Init_Buffer_SLHj);
    Buffer_SLHsl = v(VT_Init_Buffer_SLHsl);
    Buffer_Csqn = v(VT_Init_Buffer_Csqn);
    
    // Land-Niederer model of contraction
    XS        = (v(VT_XS_init));
    XW        = (v(VT_XW_init));
    CaTRPN    = (v(VT_TRPN_init));
    TmBlocked = (v(VT_TmBlocked_init));
    ZETAS     = (v(VT_ZETAS_init));
    ZETAW     = (v(VT_ZETAW_init));
    Cd        = (v(VT_Cd_init));
    
}  // TWorld::Init

ML_CalcType TWorld::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
  tinc *= 1000.0;  // second to millisecond conversion
  ML_CalcType V_m = V * 1000.0;
  double VEL = 0.0; //velocity/1000.0;  // 1/s to 1/ms

  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_m)+.5);  // array position

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Model parameters
  ////////////////////////////////////////////////////////////////////////////////////////
  // Constants
  double R = 8314.0; // [J/kmol*K]
  double F = 96485.0; // [C/mol]
  double T = 310.0; // [K]
  double frt = F / (R * T);
  double vfrt = V_m  * frt;
  double vffrt = V_m * F * frt;
  double Cmem = 1.3810e-10; // [F] membrane capacitance
  double Qpow = (T - 310.0) / 10.0;
    
  // Cell Tometry
  double cellLength = 100.0; // cell length [um]
  double cellRadius = 10.25; // cell radius [um]
  double Vcell = M_PI * pow(cellRadius, 2) * cellLength * 1e-15; // [L]
  double Vmyo = 0.65 * Vcell;
  double Vsr = 0.035 * Vcell;
  double Vsl = 0.02 * Vcell;
  double Vdyad = 0.0539 * 0.01 * Vcell;
    
  // Fractional currents in compartments
  double Fdyad = 0.11;
  double Fsl = 1 - Fdyad;
    
  // Fixed ion concentrations
  double Cl_o = 150.0; // Extracellular Cl  [mM]
  double Mg_myo = 0.5; // Intracellular Mg  [mM]
    
  // Nerst Potentials
  const ML_CalcType E_Na_dyad = (1.0 / frt) * log(Na_o / Na_dyad);
  const ML_CalcType E_Na_sl = (1.0 / frt) * log(Na_o / Na_sl);
  const ML_CalcType E_K = (1.0 / frt) * log(v(VT_K_o) / K_myo);
  const ML_CalcType E_Ca_dyad = (1.0 / frt / 2.0) * log(v(VT_Ca_o) / Ca_dyad);
  const ML_CalcType E_Ca_sl = (1.0 / frt / 2.0) * log(v(VT_Ca_o) / Ca_sl);
  const ML_CalcType E_Cl = (1.0 / frt) * log(Cl_myo / Cl_o);

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Buffering parameters
  ////////////////////////////////////////////////////////////////////////////////////////
  double Bmax_Naj = 7.561;
  double Bmax_Nasl = 1.65;
  double koff_na = 1e-3;
  double kon_na = 0.1e-3;
  double Bmax_TnClow = 70.0e-3;
  double Bmax_TnChigh = 140.0e-3;
  double koff_tnchca = 0.032e-3;
  double kon_tnchca = 2.37;
  double koff_tnchmg = 3.33e-3;
  double kon_tnchmg = 3.0e-3;
  double Bmax_CaM = 24.0e-3 * v(VT_CMDN_Multiplier);
  double koff_cam = 238.0e-3;
  double kon_cam = 34.0;
  double Bmax_myosin = 140.0e-3;
  double koff_myoca = 0.46e-3;
  double kon_myoca = 13.8;
  double koff_myomg = 0.057e-3;
  double kon_myomg = 0.0157;
  double Bmax_SR = 17.85854e-3;
  double koff_sr = 60.0e-3;
  double kon_sr = 100.0;
  double Bmax_SLlowsl = 33.923e-3 * Vmyo / Vsl;
  double Bmax_SLlowj = 4.89983e-4 * Vmyo / Vdyad;
  double koff_sll = 1300.0e-3;
  double kon_sll = 100.0;
  double Bmax_SLhighsl = 12.15423e-3 * Vmyo / Vsl;
  double Bmax_SLhighj = 1.75755e-4 * Vmyo / Vdyad;
  double koff_slh = 30.0e-3;
  double kon_slh = 100.0;
  double Bmax_Csqn = 136.55214e-3 * Vmyo / Vsr;
  double koff_csqn = 65.0;
  double kon_csqn = 100.0;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        CaMK and Ca signalling
  ////////////////////////////////////////////////////////////////////////////////////////
  double PP1_tot = v(VT_PP1_tot); //0.13698
  double CaMK0  = 2.0 * 0.05; // Equilibrium fraction of active CaMKII binding sites
  double Km_CaMK_Ca = 5.0 * 0.0015; //[mmol/L] CaMKII affinity for Ca2+/CaM activation %Adjusted because of smaller cleft space
    
  const ML_CalcType CaMK_bound = CaMK0 * (1.0 - CaMK_trap) / (1 + Km_CaMK_Ca / Ca_dyad);
  const ML_CalcType CaMK_active = CaMK_bound + CaMK_trap; // Fraction of active CaMKII
    
  double alpha = 0.05;
  double beta  = 6.8e-4;
  const ML_CalcType dCaMK_trap = alpha * CaMK_bound * CaMK_active - beta * CaMK_trap * (0.1 + 0.9 * PP1_tot / 0.1371); //dCaMK_Trap/dt
  CaMK_trap += tinc * dCaMK_trap;
    
  double tau_plb = 100000.0; // [ms] Time constant of CaMKII PLB phosphorylation
  double tau_ryr = 10000.0; // [ms] Time constant of CaMKII RyR phosphorylation
  double tau_cal = tau_ryr; // Time constant of CaMKII ICaL phosphorylation
    
  double K_Phos_CaMK = 0.35;  // Affinity of PLB, ICaL etc for CaMKII
  const ML_CalcType CaMK_Phos_ss_ICaL =  CaMK_active / (CaMK_active + K_Phos_CaMK);
  const ML_CalcType CaMK_Phos_ss_RyR =  CaMK_active / (CaMK_active + 1.0);
  const ML_CalcType CaMK_Phos_ss_PLB =  CaMK_active / (CaMK_active + 10.0);

  const ML_CalcType dCaMK_f_ICaL = (CaMK_Phos_ss_ICaL - CaMK_f_ICaL) / tau_cal;
  CaMK_f_ICaL += tinc * dCaMK_f_ICaL;
  const ML_CalcType dCaMK_f_RyR = (CaMK_Phos_ss_RyR - CaMK_f_RyR)  / tau_ryr;
  CaMK_f_RyR += tinc * dCaMK_f_RyR;
  const ML_CalcType dCaMK_f_PLB = (CaMK_Phos_ss_PLB - CaMK_f_PLB)  / tau_plb;
  CaMK_f_PLB += tinc * dCaMK_f_PLB;
    
  double alpha_serca = 0.05;
  const ML_CalcType bound_serca = CaMK0 * (1.0 - casig_serca_trap) / (1.0 + Km_CaMK_Ca / Ca_dyad);
  const ML_CalcType casig_SERCA_act = bound_serca + casig_serca_trap; // Fraction of active CaMKII
  const ML_CalcType dcasig_serca_trap = alpha_serca * bound_serca * casig_SERCA_act - beta * casig_serca_trap * (0.1 + 0.9 * PP1_tot / 0.1371); //dCaMK_Trap/dt
  casig_serca_trap += tinc * dcasig_serca_trap;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        PKA phosphorylation
  ////////////////////////////////////////////////////////////////////////////////////////
  double fINa_PKA = v(VT_fINa_PKA);
  double fICaL_PKA = v(VT_fICaL_PKA);
  double fINaK_PKA = v(VT_fINaK_PKA);
  double fIKs_PKA = v(VT_fIKs_PKA);
  double fPLB_PKA = v(VT_fIfPLB_PKA);
  double fTnI_PKA = v(VT_fTnI_PKA);
  double fMyBPC_PKA = v(VT_fMyBPC_PKA);
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium current (INa, INaL)
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////// calulate I_Na //////////
  // m gate
  const ML_CalcType m_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86)/9.03))*(1.0 + exp(-(V_m + 56.86)/9.03)));
  const ML_CalcType tau_m = std::min(0.02, 0.1292 * exp(-(((V_m + 45.79)/15.54) * ((V_m + 45.79)/15.54))) + 0.06487 * exp(-(((V_m - 4.823)/51.12) * ((V_m - 4.823)/51.12)))); // min set based on discussion with Tomek himself 
  m = m_inf - (m_inf - m) * exp(-tinc / tau_m);
    
  // h and j gate
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

  // gating CaMK-P
  const ML_CalcType h_p_inf = (1.0/((1.0+exp(((V_m+77.55)/7.43)))*(1.0+exp(((V_m+77.55)/7.43)))));
  h_p = h_p_inf - (h_p_inf - h_p) * exp(-tinc / tau_h);
  const ML_CalcType tau_j_p = 1.46 * tau_j;
  j_p = j_inf - (j_inf - j_p) * exp(-tinc / tau_j_p);
    
  // gating PKA
  const ML_CalcType m_PKA_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86 + 5.0)/9.03))*(1.0 + exp(-(V_m + 56.86 + 5.0)/9.03)));
  m_PKA = m_PKA_inf - (m_PKA_inf - m_PKA) * exp(-tinc / tau_m);
      
  const ML_CalcType h_PKA_inf = (1.0/((1.0+exp(((V_m+76.55)/7.43)))*(1.0+exp(((V_m+76.55)/7.43)))));
  h_PKA = h_PKA_inf - (h_PKA_inf - h_PKA) * exp(-tinc / tau_h);
    
  const ML_CalcType j_PKA_inf = h_PKA_inf;
  j_PKA = j_PKA_inf - (j_PKA_inf - j_PKA) * exp(-tinc / tau_j_p);
      
  // Both Phosphorylated
  const ML_CalcType h_both_inf = (1.0/((1.0+exp(((V_m+82.55)/7.43)))*(1.0+exp(((V_m+82.55)/7.43)))));
  h_both = h_both_inf - (h_both_inf - h_both) * exp(-tinc / tau_h);
  const ML_CalcType j_both_inf = h_both_inf;
  j_both = j_both_inf - (j_both_inf - j_both) * exp(-tinc / tau_j_p);
      
  // Putting together the channels behavior and fraction
  double G_Na = 22.08788 * v(VT_INaF_Multiplier);
  double G_Na_PKA = G_Na * 1.25;
      
  const ML_CalcType phi_INa_CaMK = CaMK_f_RyR;
  const ML_CalcType phi_INa_PKA = fINa_PKA;
  const ML_CalcType phi_INa_Both = phi_INa_CaMK * phi_INa_PKA;
  const ML_CalcType phi_INa_CaMKonly = phi_INa_CaMK - phi_INa_Both;
  const ML_CalcType phi_INa_PKAonly = phi_INa_PKA - phi_INa_Both;

  I_Na_Base_NP = G_Na * m*m*m * h * j; //Non-Phosphorylated
  I_Na_Base_CaMK = G_Na * m*m*m * h_p * j_p;
  I_Na_Base_PKA = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_PKA * j_PKA;
  I_Na_Base_Both = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_both * j_both;
      
  // 4 population
  I_Na_Base = ((1-phi_INa_CaMKonly-phi_INa_PKAonly-phi_INa_Both)*I_Na_Base_NP + phi_INa_CaMKonly*I_Na_Base_CaMK + phi_INa_PKAonly*I_Na_Base_PKA + phi_INa_Both*I_Na_Base_Both);
  I_NaFast_dyad = Fdyad * I_Na_Base*(V_m - E_Na_dyad);
  I_NaFast_sl = (1 - Fdyad) * I_Na_Base * (V_m - E_Na_sl);

  ///////////// calulate I_NaL //////////
  // m gate
  const ML_CalcType m_L_inf = 1.0 / (1.0 + exp(-(V_m + 42.85) / 5.264));
  const ML_CalcType tau_m_L = tau_m;
  m_L = m_L_inf - (m_L_inf - m_L) * exp(-tinc / tau_m_L);
    
  // h gate
  const ML_CalcType h_L_inf = 1.0 / (1.0 + exp((V_m + 87.61) / 7.488));
  double tau_h_L = 145.0;
  h_L = h_L_inf - (h_L_inf - h_L) * exp(-tinc / tau_h_L);
    
  // gating CaMK-P
  const ML_CalcType h_L_p_inf = 1.0 / (1.0 + exp((V_m + 93.81) / 7.488));
  const ML_CalcType tau_h_L_p = 3.0 * tau_h_L;
  h_L_p = h_L_p_inf - (h_L_p_inf - h_L_p) * exp(-tinc / tau_h_L_p);
    
  // Putting together the channels behavior and fraction
  const ML_CalcType phi_INaL_CaMK = CaMK_f_RyR;
  double G_Na_L = 0.04229 * v(VT_INaL_Multiplier) * (1.0 + phi_INaL_CaMK);
  if (v(VT_celltype) == 1.0) { // epi
      G_Na_L = G_Na_L * 0.7;
    }
  I_NaL_dyad = Fdyad  * G_Na_L * (V_m - E_Na_dyad) * m_L * ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_p);
  I_NaL_sl = Fsl * G_Na_L * (V_m - E_Na_sl) * m_L * ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_p);

  ///////////// Combinations I_Na and I_NaL //////////
  I_Na_dyad = I_NaFast_dyad + I_NaL_dyad;
  I_Na_sl = I_NaFast_sl + I_NaL_sl ;
  I_NaFast = I_NaFast_dyad+I_NaFast_sl;
  I_NaL = I_NaL_dyad+I_NaL_sl;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        L-type calcium current (I_CaL, I_CaNa, I_CaK)
  ////////////////////////////////////////////////////////////////////////////////////////
  // d gate
  const ML_CalcType d_inf = min(1.0763*exp((-1.007*exp((-0.0829*(V_m+3.62483))))), 1.0);
  const ML_CalcType tau_d = 1.5 + 1.0 / (exp(-0.05 * (V_m + 6.0)) + exp(0.09 * (V_m + 14.0)));
  d = d_inf - (d_inf - d) * exp(-tinc / tau_d);
    
  // f gate
  const ML_CalcType f_inf = 1.0 / (1.0 + exp((V_m + 19.58) / 3.696));
  const ML_CalcType tau_f_fast = 6.17111 + 1.0 / (0.00126 * exp(-(V_m + 26.63596) / 9.69961) + 0.00126 * exp((V_m + 26.63596) / 9.69961));
  const ML_CalcType tau_f_slow = 2719.22489 + 1.0 / (7.19411e-05 * exp(-(V_m + 5.74631) / 10.87690) + 7.19411e-05 * exp((V_m + 5.74631) / 16.31535));
    
  const ML_CalcType A_f_fast = 0.52477;
  const ML_CalcType A_f_slow = 1.0 - A_f_fast;
  f_fast = f_inf - (f_inf - f_fast) * exp(-tinc / tau_f_fast);
  f_slow = f_inf - (f_inf - f_slow) * exp(-tinc / tau_f_slow);
  const ML_CalcType f = A_f_fast * f_fast + A_f_slow * f_slow;
    
  const ML_CalcType f_Ca_inf      = f_inf;
  const ML_CalcType tau_f_Ca_fast = 13.50673 + 1.0 / (0.15420 * exp(-(V_m - 1.31611) / 11.33960) + 0.15420 * exp((V_m - 1.31611) / 11.33960));
  const ML_CalcType tau_f_Ca_slow = 177.95813 + 1.0 / (4.73955e-04 * exp((-V_m+0.79049) / 0.81777) + 4.73955e-04 * exp((V_m+2.40474) /1.90812));
  const ML_CalcType A_f_Ca_fast = 0.3 + 0.6 / (1.0 + exp((V_m - 9.24247) / 27.96201));
  const ML_CalcType A_f_Ca_slow = 1.0 - A_f_Ca_fast;
  f_Ca_fast = f_Ca_inf - (f_Ca_inf - f_Ca_fast) * exp(-tinc / tau_f_Ca_fast);
  f_Ca_slow = f_Ca_inf - (f_Ca_inf - f_Ca_slow) * exp(-tinc / tau_f_Ca_slow);
  const ML_CalcType f_Ca = A_f_Ca_fast * f_Ca_fast + A_f_Ca_slow * f_Ca_slow;
    
  // j gate
  double tau_j_Ca = 66.0;
  const ML_CalcType j_Ca_inf = 1. / (1 + exp((V_m + 17.66945) / (3.21501)));
  j_Ca = j_Ca_inf - (j_Ca_inf - j_Ca) * exp(-tinc / tau_j_Ca);
    
  // gating CaMK-P
  const ML_CalcType tau_f_p_fast = 2.5 * tau_f_fast;
  const ML_CalcType f_p_inf = f_inf;
  f_p_fast = f_p_inf - (f_p_inf - f_p_fast) * exp(-tinc / tau_f_p_fast);
  const ML_CalcType f_p = A_f_fast * f_p_fast + A_f_slow * f_slow;
    
  const ML_CalcType tau_f_Ca_p_fast = 2.5 * tau_f_Ca_fast;
  const ML_CalcType f_Ca_p_inf = f_inf;
  f_Ca_p_fast = f_Ca_p_inf - (f_Ca_p_inf - f_Ca_p_fast) * exp(-tinc / tau_f_Ca_p_fast);
  const ML_CalcType f_Ca_p = A_f_Ca_fast * f_Ca_p_fast + A_f_Ca_slow * f_Ca_slow;
    
  // gating PKA
  const ML_CalcType d_PKA_inf = min(1.0323 * exp(-1.0553 * exp(-0.0810 * (V_m + 12.62483))), 1.0);
  d_PKA = d_PKA_inf - (d_PKA_inf - d_PKA) * exp(-tinc / tau_d);
    
  const ML_CalcType f_PKA_inf = 1.0 / (1.0 + exp((V_m + 19.58 + 6) / 3.696));
  f_PKA_fast = f_PKA_inf - (f_PKA_inf - f_PKA_fast) * exp(-tinc / tau_f_fast);
  f_PKA_slow = f_PKA_inf - (f_PKA_inf - f_PKA_slow) * exp(-tinc / tau_f_slow);
  const ML_CalcType f_PKA = A_f_fast * f_PKA_fast + A_f_slow * f_PKA_slow;
    
  const ML_CalcType f_Ca_PKA_inf = f_PKA_inf;
  f_Ca_PKA_fast = f_Ca_PKA_inf - (f_Ca_PKA_inf - f_Ca_PKA_fast) * exp(-tinc / tau_f_Ca_fast);
  f_Ca_PKA_slow = f_Ca_PKA_inf - (f_Ca_PKA_inf - f_Ca_PKA_slow) * exp(-tinc / tau_f_Ca_slow);
  const ML_CalcType f_Ca_PKA = A_f_Ca_fast * f_Ca_PKA_fast + A_f_Ca_slow * f_Ca_PKA_slow;
    
  // Both Phosphorylated
  const ML_CalcType f_both_inf = f_PKA_inf;
  f_both_fast = f_both_inf - (f_both_inf - f_both_fast) * exp(-tinc / tau_f_p_fast);
  const ML_CalcType f_both = A_f_fast * f_both_fast + A_f_slow * f_PKA_slow;
    
  const ML_CalcType f_Ca_both_inf = f_Ca_PKA_inf;
  f_Ca_both_fast = f_Ca_both_inf - (f_Ca_both_inf - f_Ca_both_fast) * exp(-tinc / tau_f_Ca_p_fast);
  const ML_CalcType f_Ca_both = A_f_Ca_fast * f_Ca_both_fast + A_f_Ca_slow * f_Ca_PKA_slow;
    
  // nca
  double Kmn = 0.00222;
  double k2n = 957.85903;
  const ML_CalcType k_m2_n = j_Ca * 0.84191;
    
  const ML_CalcType alpha_n_Ca_dyad   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_dyad)), 3.80763));
  const ML_CalcType dn_Ca_dyad = alpha_n_Ca_dyad * k2n - n_Ca_dyad * k_m2_n;
  n_Ca_dyad += tinc * dn_Ca_dyad;
    
  const ML_CalcType alpha_n_Ca_sl   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_sl)), 3.80763));
  const ML_CalcType dn_Ca_sl = alpha_n_Ca_sl * k2n - n_Ca_sl * k_m2_n;
  n_Ca_sl += tinc * dn_Ca_sl;
    
  // driving force
  const ML_CalcType I_o = (0.5 * (Na_o + v(VT_K_o) + v(VT_Cl_o) + (4.0 * v(VT_Ca_o))) / 1000.0);
  const ML_CalcType I_dyad = (0.5 * (Na_sl + K_myo + Cl_myo + (4.0 * Ca_sl))) / 1000.0;
  const ML_CalcType I_sl = (0.5 * (Na_myo + K_myo + Cl_myo + (4.0 * Ca_i))) / 1000.0;
    
  double dielConstant = 74.0; // water at 37°
  double temp = 310.0; // body temp in kelvins.
  double constA = 1.82e6 * pow((dielConstant * temp), -1.5);
    
  const ML_CalcType gamma_Ca_o = pow(10.0, -constA * 4.0 * (sqrt(I_o)/(1.0 + sqrt(I_o)) - 0.3 * I_o));
  const ML_CalcType gamma_Ca_dyad = pow(10.0, -constA * 4.0 * (sqrt(I_dyad)/(1.0 + sqrt(I_dyad)) - 0.3 * I_dyad));
  const ML_CalcType gamma_Ca_sl = pow(10.0, -constA * 4.0 * (sqrt(I_sl)/(1.0 + sqrt(I_sl)) - 0.3 * I_sl));
  const ML_CalcType gamma_Na_o = pow(10.0, -constA * 1.0 * (sqrt(I_o)/(1.0 + sqrt(I_o)) - 0.3 * I_o));
  const ML_CalcType gamma_Na_dyad = pow(10.0, -constA * 1.0 * (sqrt(I_dyad)/(1.0 + sqrt(I_dyad)) - 0.3 * I_dyad));
  const ML_CalcType gamma_Na_sl = pow(10.0, -constA * 1.0 * (sqrt(I_sl)/(1.0 + sqrt(I_sl)) - 0.3 * I_sl));
  const ML_CalcType gamma_K_o = pow(10.0, -constA * 1.0 * (sqrt(I_o)/(1.0 + sqrt(I_o)) - 0.3 * I_o));
  const ML_CalcType gamma_K_dyad = pow(10.0, -constA * 1.0 * (sqrt(I_dyad)/(1.0 + sqrt(I_dyad)) - 0.3 * I_dyad));
  const ML_CalcType gamma_K_sl = pow(10.0, -constA * 1.0 * (sqrt(I_sl)/(1.0 + sqrt(I_sl)) - 0.3 * I_sl));

  const ML_CalcType Psi_CaL_dyad = 4.0 * vffrt * (gamma_Ca_dyad * Ca_sl * exp(2.0 * vfrt) - gamma_Ca_o * v(VT_Ca_o)) / (exp(2.0 * vfrt) - 1.0);
  const ML_CalcType Psi_CaNa_dyad = 1.0 * vffrt * (gamma_Na_dyad * Na_sl * exp(1.0 * vfrt) - gamma_Na_o * Na_o) / (exp(1.0 * vfrt) - 1.0);
  const ML_CalcType Psi_CaK_dyad = 1.0 * vffrt * (gamma_K_dyad * K_myo * exp(1.0 * vfrt) - gamma_K_o * v(VT_K_o)) / (exp(1.0 * vfrt) - 1.0);
  const ML_CalcType Psi_CaL_sl = 4.0 * vffrt * (gamma_Ca_sl * Ca_i * exp(2.0 * vfrt) - gamma_Ca_o * v(VT_Ca_o)) / (exp(2.0 * vfrt) - 1.0);
  const ML_CalcType Psi_CaNa_sl = 1.0 * vffrt * (gamma_Na_sl * Na_myo * exp(1.0 * vfrt) - gamma_Na_o * Na_o) / (exp(1.0 * vfrt) - 1.0);
  const ML_CalcType Psi_CaK_sl = 1.0 * vffrt * (gamma_K_sl * K_myo * exp(1.0 * vfrt) - gamma_K_o * v(VT_K_o)) / (exp(1.0 * vfrt) - 1.0);
    
  // Calculating "pure" CDI
  double rateRecovery = 0.02313;
    const ML_CalcType sigmoidTransition_dyad = 1 - 1/(1 + (1.86532 * Ca_dyad/0.032));
  const ML_CalcType tauTransition_dyad = 1.09670 + (1.0 - sigmoidTransition_dyad) * 141.42990;
  const ML_CalcType dI_CaL_pureCDI_dyad = -I_CaL_pureCDI_dyad * sigmoidTransition_dyad / tauTransition_dyad + (1.0 - I_CaL_pureCDI_dyad) * rateRecovery;
  I_CaL_pureCDI_dyad += tinc * dI_CaL_pureCDI_dyad;
    const ML_CalcType sigmoidTransition_sl = 1 - 1/(1 + (1.86532 * Ca_sl/0.032));
  const ML_CalcType tauTransition_sl = 1.09670 + (1.0 - sigmoidTransition_sl) * 141.42990;
  const ML_CalcType dI_CaL_pureCDI_sl = -I_CaL_pureCDI_sl * sigmoidTransition_sl / tauTransition_sl + (1.0 - I_CaL_pureCDI_sl) * rateRecovery;
  I_CaL_pureCDI_sl += tinc * dI_CaL_pureCDI_sl;
    
  // Channel conductances
  double P_Ca = 1.5768e-04 * v(VT_ICaL_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      P_Ca = P_Ca * 1.025;
  } else if (v(VT_celltype) == 2.0) { // mid
      P_Ca = P_Ca * 1.1;
  }
  double P_Ca_p = 1.1 * P_Ca;
  double P_Ca_PKA = P_Ca * 1.9; // BetaAdrenergic;  (%Gong--> 1.9)
  if (v(VT_celltype) == 1.0) { // epi
      P_Ca_PKA = P_Ca_PKA * 1.025;
  } else if (v(VT_celltype) == 2.0) { // mid
      P_Ca_PKA = P_Ca_PKA * 1.1;
  }
  double P_CaNa = 1.1737 / 1.8969 * 0.00125 * P_Ca;
  double P_CaK = 1.1737 / 1.8969 * 3.574e-4 * P_Ca;
  double P_CaNa_p = 1.1737 / 1.8969 * 0.00125 * P_Ca_p;
  double P_CaK_p = 1.1737 / 1.8969 * 3.574e-4 * P_Ca_p;
  double P_CaNa_PKA = 1.1737 / 1.8969 * 0.00125 * P_Ca_PKA;
  double P_CaK_PKA = 1.1737 / 1.8969 * 3.574e-4 * P_Ca_PKA;
      
  // Putting together the channels behavior and fraction
  I_CaL_dyad_NP = P_Ca * Psi_CaL_dyad * d * (f * (1.0 - n_Ca_dyad) + j_Ca * f_Ca * n_Ca_dyad);
  I_CaL_dyad_CaMK = P_Ca_p * Psi_CaL_dyad * d * (f_p * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_p * n_Ca_dyad);
  I_CaL_dyad_PKA = P_Ca_PKA * Psi_CaL_dyad * d_PKA * (f_PKA * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_PKA * n_Ca_dyad);
  I_CaL_dyad_Both = P_Ca_PKA * Psi_CaL_dyad * d_PKA * (f_both * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_both * n_Ca_dyad);
  I_CaL_sl_NP = P_Ca * Psi_CaL_sl * d * (f * (1.0 - n_Ca_sl) + j_Ca * f_Ca * n_Ca_sl);
  I_CaL_sl_CaMK = P_Ca_p * Psi_CaL_sl * d * (f_p * (1.0 - n_Ca_sl) + j_Ca * f_Ca_p * n_Ca_sl);
  I_CaL_sl_PKA = P_Ca_PKA * Psi_CaL_sl * d_PKA * (f_PKA * (1.0 - n_Ca_sl) + j_Ca * f_Ca_PKA * n_Ca_sl);
  I_CaL_sl_Both = P_Ca_PKA * Psi_CaL_sl * d_PKA * (f_both * (1.0 - n_Ca_sl) + j_Ca * f_Ca_both * n_Ca_sl);
    
  I_CaNa_dyad_NP = P_CaNa * Psi_CaNa_dyad * d * (f * (1.0 - n_Ca_dyad) + j_Ca * f_Ca * n_Ca_dyad);
  I_CaNa_dyad_CaMK = P_CaNa_p * Psi_CaNa_dyad * d * (f_p * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_p * n_Ca_dyad);
  I_CaNa_dyad_PKA = P_CaNa_PKA * Psi_CaNa_dyad * d_PKA * (f_PKA * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_PKA * n_Ca_dyad);
  I_CaNa_dyad_Both = P_CaNa_PKA * Psi_CaNa_dyad * d_PKA * (f_both * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_both * n_Ca_dyad);
  I_CaNa_sl_NP = P_CaNa * Psi_CaNa_sl * d * (f * (1.0 - n_Ca_sl) + j_Ca * f_Ca * n_Ca_sl);
  I_CaNa_sl_CaMK = P_CaNa_p * Psi_CaNa_sl * d * (f_p * (1.0 - n_Ca_sl) + j_Ca * f_Ca_p * n_Ca_sl);
  I_CaNa_sl_PKA = P_CaNa_PKA * Psi_CaNa_sl * d_PKA * (f_PKA * (1.0 - n_Ca_sl) + j_Ca * f_Ca_PKA * n_Ca_sl);
  I_CaNa_sl_Both = P_CaNa_PKA * Psi_CaNa_sl * d_PKA * (f_both * (1.0 - n_Ca_sl) + j_Ca * f_Ca_both * n_Ca_sl);
    
  I_CaK_dyad_NP = P_CaK * Psi_CaK_dyad * d * (f * (1.0 - n_Ca_dyad) + j_Ca * f_Ca * n_Ca_dyad);
  I_CaK_dyad_CaMK = P_CaK_p * Psi_CaK_dyad * d * (f_p * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_p * n_Ca_dyad);
  I_CaK_dyad_PKA = P_CaK_PKA * Psi_CaK_dyad * d_PKA * (f_PKA * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_PKA * n_Ca_dyad);
  I_CaK_dyad_Both = P_CaK_PKA * Psi_CaK_dyad * d_PKA * (f_both * (1.0 - n_Ca_dyad) + j_Ca * f_Ca_both * n_Ca_dyad);
  I_CaK_sl_NP = P_CaK * Psi_CaK_sl * d * (f * (1.0 - n_Ca_sl) + j_Ca * f_Ca * n_Ca_sl);
  I_CaK_sl_CaMK = P_CaK_p * Psi_CaK_sl * d * (f_p * (1.0 - n_Ca_sl) + j_Ca * f_Ca_p * n_Ca_sl);
  I_CaK_sl_PKA = P_CaK_PKA * Psi_CaK_sl * d_PKA * (f_PKA * (1.0 - n_Ca_sl) + j_Ca * f_Ca_PKA * n_Ca_sl);
  I_CaK_sl_Both = P_CaK_PKA * Psi_CaK_sl * d_PKA * (f_both * (1.0 - n_Ca_sl) + j_Ca * f_Ca_both * n_Ca_sl);
    
  //4 population combination
  const ML_CalcType phi_ICaL_CaMK = CaMK_f_ICaL;
  const ML_CalcType phi_ICaL_PKA = fICaL_PKA;
  const ML_CalcType phi_ICaL_Both = phi_ICaL_CaMK * phi_ICaL_PKA;
  const ML_CalcType phi_ICaL_CaMKonly = phi_ICaL_CaMK - phi_ICaL_Both;
  const ML_CalcType phi_ICaL_PKAonly = phi_ICaL_PKA - phi_ICaL_Both;
    
  I_CaL_dyad = ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaL_dyad_NP + phi_ICaL_CaMKonly * I_CaL_dyad_CaMK + phi_ICaL_PKAonly * I_CaL_dyad_PKA + phi_ICaL_Both * I_CaL_dyad_Both) * I_CaL_pureCDI_dyad;
  I_CaNa_dyad= ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaNa_dyad_NP + phi_ICaL_CaMKonly * I_CaNa_dyad_CaMK + phi_ICaL_PKAonly * I_CaNa_dyad_PKA + phi_ICaL_Both * I_CaNa_dyad_Both) * I_CaL_pureCDI_dyad;
  I_CaK_dyad = ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaK_dyad_NP + phi_ICaL_CaMKonly * I_CaK_dyad_CaMK + phi_ICaL_PKAonly * I_CaK_dyad_PKA + phi_ICaL_Both * I_CaK_dyad_Both) * I_CaL_pureCDI_dyad;
  I_CaL_sl = ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaL_sl_NP + phi_ICaL_CaMKonly * I_CaL_sl_CaMK + phi_ICaL_PKAonly * I_CaL_sl_PKA + phi_ICaL_Both * I_CaL_sl_Both) * I_CaL_pureCDI_sl;
  I_CaNa_sl = ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaNa_sl_NP + phi_ICaL_CaMKonly * I_CaNa_sl_CaMK + phi_ICaL_PKAonly * I_CaNa_sl_PKA + phi_ICaL_Both * I_CaNa_sl_Both) * I_CaL_pureCDI_sl;
  I_CaK_sl = ((1.0 - phi_ICaL_CaMKonly - phi_ICaL_PKAonly - phi_ICaL_Both) * I_CaK_sl_NP + phi_ICaL_CaMKonly * I_CaK_sl_CaMK + phi_ICaL_PKAonly * I_CaK_sl_PKA + phi_ICaL_Both * I_CaK_sl_Both) * I_CaL_pureCDI_sl;
    
  // Weigh I_CaL in sl and dyad
  I_CaL_dyad = I_CaL_dyad * v(VT_ICaL_fractionSS);
  I_CaNa_dyad = I_CaNa_dyad * v(VT_ICaL_fractionSS);
  I_CaK_dyad = I_CaK_dyad * v(VT_ICaL_fractionSS);
  I_CaL_sl = I_CaL_sl * (1 - v(VT_ICaL_fractionSS));
  I_CaNa_sl = I_CaNa_sl * (1 - v(VT_ICaL_fractionSS));
  I_CaK_sl = I_CaK_sl * (1 - v(VT_ICaL_fractionSS));

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Transient outward current (Ito)
  ////////////////////////////////////////////////////////////////////////////////////////
  // activation (a gate)
  const ML_CalcType a_inf = 1.0 / (1.0 + exp(-(V_m - 19.0) / 13.0));; //Info: changed x -> a and y -> i to better match Tomek and OHara syntax
  const ML_CalcType tau_a_slow = 9.0 / (1.0 + exp((V_m + 3.0) / 15.0)) + 0.5;
  a_slow = a_inf - (a_inf - a_slow) * exp(-tinc / tau_a_slow);
    
  const ML_CalcType tau_a_fast = 8.5 * exp(-pow(((V_m + 45.0) / 50.0), 2)) + 0.5;
  a_fast = a_inf - (a_inf - a_fast) * exp(-tinc / tau_a_fast);
    
  // inactivation (i gate)
  const ML_CalcType i_inf = 1.0 / (1.0 + exp((V_m + 19.5) / 5.0));
  const ML_CalcType tau_i_slow = 800.0 / (1.0 +exp((V_m + 60.0) / 10.0)) + 30.0;
  i_slow = i_inf - (i_inf - i_slow) * exp(-tinc / tau_i_slow);
    
  const ML_CalcType tau_i_fast = 85.0 * exp(-pow((V_m + 40.0), 2.0) / 220) + 7.0;
  i_fast = i_inf - (i_inf - i_fast) * exp(-tinc / tau_i_fast);
    
  // gating CaMK-P
  const ML_CalcType a_p_inf = 1.0 / (1.0 + exp(-(V_m - 29.0) / 13.0));
  a_p_slow = a_p_inf - (a_p_inf - a_p_slow) * exp(-tinc / tau_a_slow);
  a_p_fast = a_p_inf - (a_p_inf - a_p_fast) * exp(-tinc / tau_a_fast);
    
  const ML_CalcType delta_p_develop = 1.354 + 1.0e-04 / (exp((V_m - 167.4) / 15.89) + exp(-(V_m - 12.23) / 0.2154));
  const ML_CalcType delta_p_recover = 1.0 - 0.5 / (1.0 + exp((V_m + 70.0) / 20.0));
  const ML_CalcType tau_i_p_slow = tau_i_slow * delta_p_develop * delta_p_recover;
  const ML_CalcType tau_i_p_fast = tau_i_fast * delta_p_develop * delta_p_recover;
  i_p_slow = i_inf - (i_inf - i_p_slow) * exp(-tinc / tau_i_p_slow);
  i_p_fast = i_inf - (i_inf - i_p_fast) * exp(-tinc / tau_i_p_fast);
    
  // Putting together the channels behavior and fraction
  double G_to_slow = v(VT_Ito_slow_Multiplier);
 double G_to_fast = v(VT_Ito_fast_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_to_slow = 0.02036 * G_to_slow;
      G_to_fast = 0.29856 * G_to_fast;
  } else if (v(VT_celltype) == 2.0) { // mid
      G_to_slow = 0.04632 * G_to_slow;
      G_to_fast = 0.14928 * G_to_fast;
  } else if (v(VT_celltype) == 0.0) { //endo
      G_to_slow = 0.07210 * G_to_slow;
      G_to_fast = 0.01276 * G_to_fast;
  }
  const ML_CalcType phi_Ito_CaMK = CaMK_f_RyR;
  I_to_slow = v(VT_Ito_Multiplier) * G_to_slow * (V_m - E_K) * ((1.0 - phi_Ito_CaMK) * a_slow * i_slow + phi_Ito_CaMK * a_p_slow * i_p_slow);
  I_to_fast = v(VT_Ito_Multiplier) * G_to_fast * (V_m - E_K) * ((1.0 - phi_Ito_CaMK) * a_fast * i_fast + phi_Ito_CaMK * a_p_fast * i_p_fast);
  I_to = I_to_slow + I_to_fast;

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Rapid delayed rectifier current (IKr)
  ////////////////////////////////////////////////////////////////////////////////////////
  const ML_CalcType alpha_Kr = 0.1161 * exp(0.299 * vfrt);
  const ML_CalcType beta_Kr = 0.2442 * exp(-1.604 * vfrt);
  double alpha_Kr_1 = 1.25 * 0.1235;
  double beta_Kr_1 = 0.1911;
  const ML_CalcType alpha_Kr_2 = 0.0578 * exp(0.971 * vfrt);
  const ML_CalcType beta_Kr_2 = 0.349e-3 * exp(-1.062 * vfrt);
  const ML_CalcType alpha_Kr_i = 0.2533 * exp(0.5953 * vfrt);
  const ML_CalcType beta_Kr_i = 0.04568 * exp(-0.8209 * vfrt);
  const ML_CalcType alpha_C2_to_I = 0.52e-4 * exp(1.525 * vfrt);
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
    
  // Putting together the channels behavior and fraction
  double G_Kr = 0.043 * sqrt(v(VT_K_o) / 5.0) * v(VT_IKr_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_Kr = G_Kr * 1.25;
  } else if (v(VT_celltype) == 2.0) { // mid
      G_Kr = G_Kr * 0.7;
  }
  I_Kr = G_Kr * O * (V_m - E_K);
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Slow delayed rectifier current (IKs)
  ////////////////////////////////////////////////////////////////////////////////////////
  // gating
  const ML_CalcType k_PKA = fIKs_PKA;
  const ML_CalcType V_h_0 = -1.0 - 10.0 * k_PKA;
  const ML_CalcType V_h_max = -12.0 - 9.0 * k_PKA;
  const ML_CalcType tau_0 = 26.0 + 9.0 * k_PKA;
  const ML_CalcType tau_max = 40.0 + 4.0 * k_PKA;
    
  const ML_CalcType V_dyad = V_h_0 + (V_h_max - V_h_0) / (1.0 + pow((350e-06 / Ca_dyad), 4.0)); // Regulated by PKA
  const ML_CalcType xs_dyad_inf = 1.0 / (1.0 + exp(-(V_m - V_dyad) / 25.0));
  const ML_CalcType V_tau_dyad = tau_0 + (tau_max - tau_0) / (1.0 + pow((150e-06/Ca_dyad), 3.0)); // Regulated by PKA
  const ML_CalcType tau_xs_dyad = 2.0 * (50.0 + (50.0 + 350.0 * exp(-(pow((V_m + 30.0), 2.0)) / 4000.0)) * 1.0 / (1.0 + exp(-(V_m + V_tau_dyad) / 10.0)));
  xs_dyad = xs_dyad_inf - (xs_dyad_inf - xs_dyad) * exp(-tinc / tau_xs_dyad);
    
  const ML_CalcType V_sl = V_h_0 + (V_h_max - V_h_0) / (1.0 + pow((350e-06 / Ca_sl), 4.0)); // Regulated by PKA
  const ML_CalcType xs_sl_inf = 1.0 / (1.0 + exp(-(V_m - V_sl) / 25.0));
  const ML_CalcType V_tau_sl = tau_0 + (tau_max - tau_0) / (1.0 + pow((150e-06/Ca_sl), 3.0)); // Regulated by PKA
  const ML_CalcType tau_xs_sl = 2.0 * (50.0 + (50.0 + 350.0 * exp(-(pow((V_m + 30.0), 2.0)) / 4000.0)) * 1.0 / (1.0 + exp(-(V_m + V_tau_sl) / 10.0)));
  xs_sl = xs_sl_inf - (xs_sl_inf - xs_sl) * exp(-tinc / tau_xs_sl);
    
  // conductances
  double G_Ks_factor_SA = 2.97002 * v(VT_IKs_Multiplier);
  
  if (v(VT_celltype) == 1.0) { // epi
      G_Ks_factor_SA = 1.4 * G_Ks_factor_SA;
  }
  else if (v(VT_celltype) == 2.0) { // mid
      G_Ks_factor_SA = 0.5 * G_Ks_factor_SA;
  }
  double G_Ks_factor = 0.01;
  const ML_CalcType G_Ks_0 = G_Ks_factor * (0.2 + 0.2 * k_PKA);
  const ML_CalcType G_Ks_max = G_Ks_factor * (0.8 + 7.0 * k_PKA);
  const ML_CalcType G_Ks_dyad = G_Ks_0 + (G_Ks_max - G_Ks_0) / (1.0 + pow((150e-06 / Ca_dyad), 1.3)); // Regulated by PKA
  const ML_CalcType G_Ks_sl = G_Ks_0 + (G_Ks_max - G_Ks_0) / (1.0 + pow((150e-06 / Ca_sl), 1.3)); // Regulated by PKA
    
  // Nernst potential
  double pNaK = 0.01833;
  const ML_CalcType E_Ks = (1.0/frt) * log((v(VT_K_o) + pNaK * Na_o) / (K_myo + pNaK * Na_myo));
    
  // Putting together the channels behavior and fraction
  I_Ks_dyad = Fdyad * G_Ks_factor_SA * G_Ks_dyad * xs_dyad * xs_dyad * (V_m - E_Ks);
  I_Ks_sl = Fsl * G_Ks_factor_SA * G_Ks_sl * xs_sl * xs_sl * (V_m - E_Ks);
  I_Ks = I_Ks_dyad + I_Ks_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Inward rectifier current (IK1)
  ////////////////////////////////////////////////////////////////////////////////////////
  // gating
  const ML_CalcType a_K1 = 4.094 / (1.0 + exp(0.1217 * (V_m - E_K - 49.934)));
  const ML_CalcType b_K1 = (15.72 * exp(0.0674 * (V_m - E_K - 3.257)) + exp(0.0618 * (V_m - E_K - 594.31))) / (1.0 + exp(-0.1629 * (V_m - E_K + 14.207)));
  const ML_CalcType K1_SS = (a_K1) / (a_K1 + b_K1);
    
  // Putting together the channels behavior and fraction
  double G_K1 = 0.6992 * v(VT_IK1_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_K1 = G_K1 * 1.1;
  } else if (v(VT_celltype) == 2.0) { // mid
      G_K1 = G_K1 * 1.3;
  }
  I_K1 = G_K1 * sqrt(v(VT_K_o)/5.0) * K1_SS * (V_m - E_K);
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium-calcium exchanger (INaCa)
  ////////////////////////////////////////////////////////////////////////////////////////
  double z_Ca = 2.0;
  double kna1 = 11.9712;
  double kna2 = 2.76;
  double kna3 = 88.767;
  double kasymm = 19.4258;
  double wna  = 3.2978e+04;
  double wca = 5.1756e+04;
  double wnaca  = 2.7763e+03;
  double kcaon = 3.4164e+06;
  double kcaoff = 3.8532e+03;
  double qna = 0.6718;
  double qca  = 0.0955;
  const ML_CalcType h_Ca = exp((qca * (V_m - 8.3117) * F) / (R * T));
  const ML_CalcType h_Na = exp((qna * (V_m - 8.3117) * F) / (R * T));
    
  // calculate I_NaCa_sl
  const ML_CalcType h_1_sl = 1.0 + Na_myo/kna3 * (1.0 + h_Na);
  const ML_CalcType h_2_sl = (Na_myo * h_Na) / (kna3 * h_1_sl);
  const ML_CalcType h_3_sl = 1.0 / h_1_sl;
  const ML_CalcType h_4_sl = 1.0 + (Na_myo * (1.0 + Na_myo / kna2)) / kna1;
  const ML_CalcType h_5_sl = (Na_myo * Na_myo) / (h_4_sl * kna1 * kna2);
  const ML_CalcType h_6_sl = 1.0 / h_4_sl;
  const ML_CalcType h_7_sl = 1.0 + (Na_o * (1.0 + 1.0 / h_Na)) / kna3;
  const ML_CalcType h_8_sl = Na_o / (kna3 * h_Na * h_7_sl);
  const ML_CalcType h_9_sl = 1.0 / h_7_sl;
  const ML_CalcType h_10_sl = kasymm + 1.0 + Na_o / kna1 * (1.0 + Na_o / kna2);
  const ML_CalcType h_11_sl = Na_o * Na_o / (h_10_sl * kna1 * kna2);
  const ML_CalcType h_12_sl = 1.0 / h_10_sl;
    
  const ML_CalcType k_1_sl = h_12_sl * v(VT_Ca_o) * kcaon;
  const ML_CalcType k_2_sl = kcaoff;
  const ML_CalcType k_3_d_sl = h_9_sl * wca;
  const ML_CalcType k_3_dd_sl = h_8_sl * wnaca;
  const ML_CalcType k_3_sl = k_3_d_sl + k_3_dd_sl;
  const ML_CalcType k_4_d_sl = (h_3_sl * wca) / h_Ca;
  const ML_CalcType k_4_dd_sl = h_2_sl * wnaca;
  const ML_CalcType k_4_sl = k_4_d_sl + k_4_dd_sl;
  const ML_CalcType k_5_sl = kcaoff;
  const ML_CalcType k_6_sl = h_6_sl * Ca_i * kcaon;
  const ML_CalcType k_7_sl = h_5_sl * h_2_sl * wna;
  const ML_CalcType k_8_sl = h_8_sl * h_11_sl * wna;
 
  const ML_CalcType x_1_sl = k_2_sl * k_4_sl * (k_7_sl + k_6_sl) + k_5_sl * k_7_sl * (k_2_sl + k_3_sl);
  const ML_CalcType x_2_sl = k_1_sl * k_7_sl * (k_4_sl + k_5_sl) + k_4_sl * k_6_sl * (k_1_sl + k_8_sl);
  const ML_CalcType x_3_sl = k_1_sl * k_3_sl * (k_7_sl + k_6_sl) + k_8_sl * k_6_sl * (k_2_sl + k_3_sl);
  const ML_CalcType x_4_sl = k_2_sl * k_8_sl * (k_4_sl + k_5_sl) + k_3_sl * k_5_sl * (k_1_sl + k_8_sl);
    
  const ML_CalcType E_1_sl = x_1_sl / (x_1_sl + x_2_sl + x_3_sl + x_4_sl);
  const ML_CalcType E_2_sl = x_2_sl / (x_1_sl + x_2_sl + x_3_sl + x_4_sl);
  const ML_CalcType E_3_sl = x_3_sl / (x_1_sl + x_2_sl + x_3_sl + x_4_sl);
  const ML_CalcType E_4_sl = x_4_sl / (x_1_sl + x_2_sl + x_3_sl + x_4_sl);
    
  double KmCaAct = 150.0e-6;
  const ML_CalcType allo_sl = 1.0 / (1.0 + ((KmCaAct / Ca_sl) * (KmCaAct / Ca_sl)));
  double z_Na = 1.0;
  const ML_CalcType J_NaCa_Na_sl = 3.0 * (E_4_sl * k_7_sl - E_1_sl * k_8_sl) + E_3_sl * k_4_dd_sl - E_2_sl * k_3_dd_sl;
  const ML_CalcType J_NaCa_Ca_sl = E_2_sl * k_2_sl - E_1_sl * k_1_sl;
  double G_NaCa = 0.00179 * v(VT_INaCa_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_NaCa = G_NaCa * 1.1;
    } else if (v(VT_celltype) == 2.0) { // mid
        G_NaCa = G_NaCa * 1.4;
    }
  I_NaCa_sl = G_NaCa * allo_sl * (z_Na * J_NaCa_Na_sl + z_Ca * J_NaCa_Ca_sl) * (1-v(VT_INaCa_fractionSS));
    
  // calculate I_NaCa_dyad
  const ML_CalcType h_1_dyad = 1.0 + (Na_sl * (1.0 + h_Na)) / kna3;
  const ML_CalcType h_2_dyad = (Na_sl * h_Na) / (kna3 * h_1_dyad);
  const ML_CalcType h_3_dyad = 1.0 / h_1_dyad;
  const ML_CalcType h_4_dyad = 1.0 + (Na_sl * (1.0 + Na_sl / kna2)) / kna1;
  const ML_CalcType h_5_dyad = (Na_sl * Na_sl) / (h_4_dyad * kna1 * kna2);
  const ML_CalcType h_6_dyad = 1.0 / h_4_dyad;
  const ML_CalcType h_7_dyad = 1.0 + Na_o / kna3 * (1.0 + 1.0 / h_Na);
  const ML_CalcType h_8_dyad = Na_o / (kna3 * h_Na * h_7_dyad);
  const ML_CalcType h_9_dyad = 1.0 / h_7_dyad;
  const ML_CalcType h_10_dyad = kasymm + 1.0 + Na_o / kna1 * (1 + Na_o / kna2);
  const ML_CalcType h_11_dyad = Na_o * Na_o / (h_10_dyad * kna1 * kna2);
  const ML_CalcType h_12_dyad = 1.0 / h_10_dyad;
    
  const ML_CalcType k_1_dyad = h_12_dyad * v(VT_Ca_o) * kcaon;
  const ML_CalcType k_2_dyad = kcaoff;
  const ML_CalcType k_3_d_dyad = h_9_dyad * wca;
  const ML_CalcType k_3_dd_dyad = h_8_dyad * wnaca;
  const ML_CalcType k_3_dyad = k_3_d_dyad + k_3_dd_dyad;
  const ML_CalcType k_4_d_dyad = (h_3_dyad * wca) / h_Ca;
  const ML_CalcType k_4_dd_dyad = h_2_dyad * wnaca;
  const ML_CalcType k_4_dyad = k_4_d_dyad + k_4_dd_dyad;
  const ML_CalcType k_5_dyad = kcaoff;
  const ML_CalcType k_6_dyad = h_6_dyad * Ca_sl * kcaon;
  const ML_CalcType k_7_dyad = h_5_dyad * h_2_dyad * wna;
  const ML_CalcType k_8_dyad = h_8_dyad * h_11_dyad * wna;
    
  const ML_CalcType x_1_dyad = k_2_dyad * k_4_dyad * (k_7_dyad + k_6_dyad) + k_5_dyad * k_7_dyad * (k_2_dyad + k_3_dyad);
  const ML_CalcType x_2_dyad = k_1_dyad * k_7_dyad * (k_4_dyad + k_5_dyad) + k_4_dyad * k_6_dyad * (k_1_dyad + k_8_dyad);
  const ML_CalcType x_3_dyad = k_1_dyad * k_3_dyad * (k_7_dyad + k_6_dyad) + k_8_dyad * k_6_dyad * (k_2_dyad + k_3_dyad);
  const ML_CalcType x_4_dyad = k_2_dyad * k_8_dyad * (k_4_dyad + k_5_dyad) + k_3_dyad * k_5_dyad * (k_1_dyad + k_8_dyad);

  const ML_CalcType E_1_dyad = x_1_dyad / (x_1_dyad + x_2_dyad + x_3_dyad + x_4_dyad);
  const ML_CalcType E_2_dyad = x_2_dyad / (x_1_dyad + x_2_dyad + x_3_dyad + x_4_dyad);
  const ML_CalcType E_3_dyad = x_3_dyad / (x_1_dyad + x_2_dyad + x_3_dyad + x_4_dyad);
  const ML_CalcType E_4_dyad = x_4_dyad / (x_1_dyad + x_2_dyad + x_3_dyad + x_4_dyad);
    
  const ML_CalcType allo_dyad = 1.0 / (1.0 + (KmCaAct / Ca_dyad) * (KmCaAct / Ca_dyad));
  const ML_CalcType J_NaCa_Na_dyad = 3.0 * (E_4_dyad * k_7_dyad - E_1_dyad * k_8_dyad) + E_3_dyad * k_4_dd_dyad - E_2_dyad * k_3_dd_dyad;
  const ML_CalcType J_NaCa_Ca_dyad = E_2_dyad * k_2_dyad - E_1_dyad * k_1_dyad;
    
  I_NaCa_dyad = v(VT_INaCa_fractionSS) * G_NaCa * allo_dyad * (z_Na * J_NaCa_Na_dyad + z_Ca * J_NaCa_Ca_dyad);

  I_NaCa = I_NaCa_dyad + I_NaCa_sl;

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium-potassium pump (INaK)
  ////////////////////////////////////////////////////////////////////////////////////////
  double I_bar_NaK = 2.10774 * v(VT_INaK_Multiplier);
  double KmNaip = 11;
  double KmNaip_PKA = 8.4615;
  double KmKo = 1.5;
    
  const ML_CalcType f_NaK = 0.75 + (0.00375 - ((140 - Na_o) / 50) * 0.001) * V_m; //Varying the slope mainly based on https://rupress.org/jgp/article-pdf/94/3/539/1814046/539.pdf
    
  I_NaK_dyad_noPKA = Fdyad * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip / Na_dyad), 4.0)) /(v(VT_K_o) + KmKo);
  I_NaK_dyad_PKA = Fdyad * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip_PKA / Na_dyad),4.0)) /(v(VT_K_o) + KmKo);
  I_NaK_sl_noPKA = Fsl * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip / Na_sl), 4.0)) / (v(VT_K_o) + KmKo);
  I_NaK_sl_PKA = Fsl * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip_PKA / Na_sl), 4.0)) / (v(VT_K_o) + KmKo);
    
  I_NaK_dyad = (1.0 - fINaK_PKA) * I_NaK_dyad_noPKA + fINaK_PKA * I_NaK_dyad_PKA;
  I_NaK_sl = (1.0 - fINaK_PKA) * I_NaK_sl_noPKA + fINaK_PKA * I_NaK_sl_PKA;
  I_NaK = I_NaK_dyad + I_NaK_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Chloride currents (ICaCl, IClb)
  ////////////////////////////////////////////////////////////////////////////////////////
  // calculate I_CaCl
  double G_CaCl = 0.01615 * v(VT_ICaCl_Multiplier);
  double Kd_CaCl = 100e-03;
    
  // I_CaCl_dyad = 0.5 * Fdyad * G_CaCl / (1.0 + Kd_CaCl / Ca_dyad)*(V_m - E_Cl);
  I_CaCl_dyad = 0.5 * G_CaCl / (1.0 + Kd_CaCl / Ca_dyad)*(V_m - E_Cl); //removed Fdyad based on suggested updates from Tomek himself
  // I_CaCl_sl = 0.5 * Fsl * G_CaCl / (1.0 + Kd_CaCl / Ca_sl)*(V_m - E_Cl);
  I_CaCl_sl = 0.5 * G_CaCl / (1.0 + Kd_CaCl / Ca_sl)*(V_m - E_Cl); //remove Fsl based on suggested updates from Tomek himself
  I_CaCl = I_CaCl_dyad + I_CaCl_sl;
    
  // calculate I_Clb
  double G_Clb = 0.00241 * v(VT_IClb_Multiplier);
  I_Clb = G_Clb * (V_m - E_Cl);
  
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Background currents (INab, ICab, IKb)
  ////////////////////////////////////////////////////////////////////////////////////////
  // calculate I_Nab
  double G_Nab = 2.0 * 0.297e-03 * v(VT_INab_Multiplier);
  I_Nab_dyad = Fdyad * G_Nab * (V_m - E_Na_dyad);
  I_Nab_sl = Fsl * G_Nab * (V_m - E_Na_sl);
  I_Nab = I_Nab_dyad + I_Nab_sl;
    
  // calculate I_Kb
  double G_Kb = 0.010879 * v(VT_IKb_Multiplier);
  const ML_CalcType x_Kb = 1.0 / (1.0 + exp(-(V_m - 10.8968) / (23.9871)));
  I_Kb = G_Kb * x_Kb * (V_m - E_K);
    
  // calculate I_Cab
  double G_Cab = 5.15575e-04 * v(VT_ICab_Multiplier);
  I_Cab_dyad = Fdyad * G_Cab * (V_m - E_Ca_dyad);
  I_Cab_sl = Fsl * G_Cab * (V_m - E_Ca_sl);
  I_Cab = I_Cab_dyad + I_Cab_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Calcium release from SR (Jrel, Jleak)
  ////////////////////////////////////////////////////////////////////////////////////////
  double ks = 26.6 * v(VT_Jrel_Multiplier);
  double koCa = 23.87221;
  double kom = 0.16219;
  double kiCa = 0.39871;
  double kim = 0.04311;
  double ec50SR = 0.75385;
  double steepnessCaSR = 5.09473;
  double caExpFactor = 0.68655 ;
  double caTransFactor = 0.94428 ;
  double caExpFactor2 = 2.06273 ;
  double caTransFactor2 = 0.52967 ;

  ///////////// calculate directly I_CaL-coupled RyRy //////////
  double directRelMidpoint = 0.95271;
  double bt = 12.47670;
  double a_rel=1.25 * bt;
  const ML_CalcType I_Ca_dyad_positive = abs(I_CaL_dyad);
  const ML_CalcType I_Ca_dyad_sigmoided = 1.0 - 1.0 / (1.0 + pow((I_Ca_dyad_positive / 0.45), 4.5));
  const ML_CalcType J_rel_inf = a_rel * (I_Ca_dyad_sigmoided) / (1.0 + pow((directRelMidpoint / Ca_SR), 7.72672));
  const ML_CalcType tau_rel = max(bt / (1.0 + 0.0123 / Ca_SR), 0.001);
  const ML_CalcType dJ_rel_ICaLdep_act = (J_rel_inf - J_rel_ICaLdep_act) / tau_rel;
  J_rel_ICaLdep_act += tinc * dJ_rel_ICaLdep_act;
    
  double tauInact = 64.11202;
  const ML_CalcType Jrel_inact_inf = 1.0 / (1.0 + (I_Ca_dyad_sigmoided / 1e-3));
  const ML_CalcType dJ_rel_ICaLdep_f1 = (Jrel_inact_inf - J_rel_ICaLdep_f1)/tauInact;
  J_rel_ICaLdep_f1 += tinc * dJ_rel_ICaLdep_f1;
    
  double tauInact2 = 119.48978;
  const ML_CalcType Jrel_inact_inf2 = 1.0 / (1.0 + (I_Ca_dyad_sigmoided / 0.6e-3));
  const ML_CalcType dJ_rel_ICaLdep_f2 = (Jrel_inact_inf2 - J_rel_ICaLdep_f2)/tauInact2;
  J_rel_ICaLdep_f2 += tinc * dJ_rel_ICaLdep_f2;
    
  J_rel_ICaLdep = (0.00174 * J_rel_ICaLdep_act * J_rel_ICaLdep_f1 * J_rel_ICaLdep_f2);
    
  ///////////// Main Ca-sensitive RyRs //////////
  double MaxSR = 15;
  double MinSR = 1;
  const ML_CalcType kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow((ec50SR / Ca_SR), steepnessCaSR));
  const ML_CalcType koSRCa = koCa / kCaSR;
  const ML_CalcType kiSRCa = kiCa * kCaSR;
    
  double ecCaI = 0.001;
  double steepnessCaI = 5.93447; // change compared to 5.7.4.5.4 was reflected here
  double minCaI = 0.93249;
  double maxCaI = 30.13294;
    
  double baseRateCaI = 3.02320e-04;
  const ML_CalcType sigmoidBaseCaI = minCaI + (maxCaI - minCaI) / (1.0 + pow((ecCaI / Ca_dyad), steepnessCaI));
  const ML_CalcType RI_to_CI = baseRateCaI * sigmoidBaseCaI;
  double CI_to_RI = 0.00248;
    
  // not phosphorylated by CaMKII
  const ML_CalcType RIcleft = 1.0 - ryr_R - ryr_O - ryr_I - ryr_CaRI;
  const ML_CalcType dryr_R = (kim * RIcleft - kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_R) - (caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * ryr_R - kom * ryr_O);
  ryr_R += tinc * dryr_R;
  const ML_CalcType dryr_O = (caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * ryr_R - kom * ryr_O) - (kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_O - kim * ryr_I);
  ryr_O += tinc * dryr_O;
  const ML_CalcType dryr_I = (kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_O - kim * ryr_I) - (kom * ryr_I - caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * RIcleft);
  ryr_I += tinc * dryr_I;
  const ML_CalcType dryr_CaRI =  RI_to_CI * RIcleft - CI_to_RI * ryr_CaRI; // shift of 29 state numbers
  ryr_CaRI += tinc * dryr_CaRI;
  
  J_SR_Carel_NP = ks * ryr_O * (Ca_SR - Ca_dyad) + J_rel_ICaLdep;
    
  // And also a version of phosphorylated
  caTransFactor2 = caTransFactor2 * 1.5; // *1.69
  const ML_CalcType RIcleftP = 1 - ryr_R_p - ryr_O_p - ryr_I_p - ryr_CaRI_p;
  const ML_CalcType dryr_R_p = (kim * RIcleftP - kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_R_p) - (caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * ryr_R_p - kom * ryr_O_p);
  ryr_R_p += tinc * dryr_R_p;
  const ML_CalcType dryr_O_p = (caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * ryr_R_p - kom * ryr_O_p) - (kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_O_p - kim * ryr_I_p);
  ryr_O_p += tinc * dryr_O_p;
  const ML_CalcType dryr_I_p = (kiSRCa * (caTransFactor * pow(Ca_dyad, caExpFactor)) * ryr_O_p - kim * ryr_I_p) - (kom * ryr_I_p - caTransFactor2 * koSRCa * pow(Ca_dyad, caExpFactor2) * RIcleftP);
  ryr_I_p += tinc * dryr_I_p;
  const ML_CalcType dryr_CaRI_p =  RI_to_CI * RIcleftP - CI_to_RI * ryr_CaRI_p;
  ryr_CaRI_p += tinc * dryr_CaRI_p;
  
  J_SR_Carel_CaMK = ks * ryr_O_p * (Ca_SR - Ca_dyad) + J_rel_ICaLdep;
    
  // Total release
  J_SR_Carel = J_SR_Carel_CaMK * CaMK_f_RyR + J_SR_Carel_NP * (1 - CaMK_f_RyR);
    
  ///////////// Additional leak //////////
  const ML_CalcType nonlinearModifier = 0.2144 * exp(1.83 * Ca_SR); // Leak should be nonlinear with load
  const ML_CalcType CaMKIILeakMultiplier = 1.0 + 2.0 * CaMK_f_RyR; // And it should be promoted by CaMKII
  J_SR_leak = 1.59306e-06 * (Ca_SR - Ca_dyad) * nonlinearModifier * CaMKIILeakMultiplier;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Calcium reuptake to the SR (Jup)
  ////////////////////////////////////////////////////////////////////////////////////////
  double Q10SRCaP = 2.6; // [none]
  double Vmax_SRCaP = 0.00543 * v(VT_Jup_Multiplier); // [mM/msec] (286 umol/L cytosol/sec)
  double Kmr = 2.31442; // [mM]L cytosol
  double hillSRCaP =  1.02809; // [mM]
  double Kmf = 0.30672e-03; // [mM] default
    
  if (v(VT_celltype) == 1.0) { // epi
      Vmax_SRCaP = Vmax_SRCaP * 1.2; // based loosely on https://www.ahajournals.org/doi/10.1161/01.RES.62468.25308.27?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed - there may be other datasets, and other values. % the 20% is just an initial guesstimate.
  }
  
  // PLB phosphorylation effect on affinity
  const ML_CalcType phosphorylationTotal = CaMK_f_PLB + fPLB_PKA - CaMK_f_PLB * fPLB_PKA; // we assume the same effect, just making sure we don't count it twice.
  double Kmf_Phospho = Kmf * 0.5; // Similar percentage effect as in Heijman 2011 % CHANGED JAKUB
    
  // Direct Ca-based acceleration
  double Km_SERCA_Ca = 0.4; // 0.03 in Heijman 2011; affinity for direct Vmax modulation by CaMKII
  double Max_Vmax_SERCA_Ca = 1.11142;
  double Vmax_mult = 1.0 + Max_Vmax_SERCA_Ca / (1.0 + pow((Km_SERCA_Ca / casig_SERCA_act), 2));
    
  J_up_NP = pow(Q10SRCaP, Qpow) * Vmax_SRCaP * Vmax_mult * (pow((Ca_i / Kmf), hillSRCaP) - pow((Ca_SR / Kmr), hillSRCaP)) / (1.0 + pow((Ca_i / Kmf), hillSRCaP) + pow((Ca_SR / Kmr), hillSRCaP));
  J_up_CaMK = pow(Q10SRCaP, Qpow) * Vmax_SRCaP * Vmax_mult * (pow((Ca_i / Kmf_Phospho), hillSRCaP) - pow((Ca_SR / Kmr), hillSRCaP)) / (1.0 + pow((Ca_i / Kmf_Phospho), hillSRCaP) + pow((Ca_SR / Kmr), hillSRCaP));
  J_up = J_up_NP * (1.0 - phosphorylationTotal) + J_up_CaMK * phosphorylationTotal;
    

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sarcolemmal calcium pump (pCa)
  ////////////////////////////////////////////////////////////////////////////////////////
  double IbarSLCaP = 0.02064  * v(VT_IpCa_Multiplier);
  double KmPCa = 0.5e-3; // [mM]
  double Q10SLCaP = 2.35; // [none]
  I_pCa_dyad = Fdyad * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Ca_dyad, 1.6) / (pow(KmPCa, 1.6) + pow(Ca_dyad, 1.6));
  I_pCa_sl = Fsl * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Ca_sl, 1.6) / (pow(KmPCa, 1.6) + pow(Ca_sl, 1.6));
  I_pCa = I_pCa_dyad + I_pCa_sl;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Land-Niederer model of contraction
  ////////////////////////////////////////////////////////////////////////////////////////
  double fracTnIpo = 0.0031; // Derived quantity (TnI_PKAp(baseline)/TnItot)
  double fPKA_TnI = (1.45 - 0.45 * (1.0 - fTnI_PKA)/(1.0 - fracTnIpo)); // multiplier for Ca unbinding from troponin.
  double PKAForceMultiplier = 1.0 + fMyBPC_PKA * 0.26;
  double PKAXBacceleration = 1.0 + fMyBPC_PKA/2;
    
  double lambda_m = std::min(1.2, stretch);
  double overlap  = 1.0 + v(VT_beta_0) * (lambda_m + std::min(0.87, lambda_m) - 1.87);
  double h        = std::max(0.0, overlap);
    
  // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
  double XU    = (1.0 - TmBlocked) - XW - XS;
  double xb_ws = v(VT_k_ws) * XW * PKAXBacceleration;
  double xb_uw = v(VT_k_uw) * XU * PKAXBacceleration;
  double xb_wu = v(VT_k_wu) * XW;
  double xb_su = v(VT_k_su) * XS;
    
  double zs_pos      = (ZETAS > 0.0) * ZETAS;
  double zs_neg      = (ZETAS < -1.0) * (-ZETAS - 1.0);
  double zs_         = std::max(zs_pos, zs_neg);  // should be zero if ZETAS + 1 is in [0,1] interval
  double gamma_rate  = (v(VT_gamma)) * zs_;
  double xb_su_gamma = gamma_rate * XS;
  double diff_XS     = xb_ws - xb_su - xb_su_gamma;
  XS += tinc * diff_XS;

  double gr_w_        = ((ZETAW < 0.0) ? -ZETAW : ZETAW); // use absolute value of ZetaW
  double gamma_rate_w = (v(VT_gamma_wu)) * gr_w_;  // weak xbs don't like being strained
  double xb_wu_gamma  = gamma_rate_w * XW;
  double diff_XW      = xb_uw - xb_wu - xb_ws - xb_wu_gamma;
  XW += tinc * diff_XW;

  double ca50_     = (v(VT_ca50) * fPKA_TnI) + v(VT_beta_1) * min(0.2,(lambda_m - 1.0));
  double diff_TRPN = v(VT_koff) * (pow((1000*Ca_i/ca50_), v(VT_TRPN_n)) * (1.0 - CaTRPN) - CaTRPN);
  CaTRPN += tinc * diff_TRPN;

  double trpn_np_       = pow(CaTRPN, -v(VT_nperm)/2.0);
  double trpn_np        = std::min(100.0, trpn_np_);
  double diff_TmBlocked = v(VT_ktm_block) * trpn_np * XU - v(VT_ktm_unblock) * pow(CaTRPN, v(VT_nperm)/2.0) * TmBlocked;
  TmBlocked += tinc * diff_TmBlocked;

  // Velocity dependence -- assumes distortion resets on W->S
  double diff_ZETAS = v(VT_A) * VEL - v(VT_cds) * ZETAS;  // - gamma_rate * ZETAS;
  double diff_ZETAW = v(VT_A) * VEL - v(VT_cdw) * ZETAW;  // - gamma_rate_w * ZETAW;
  ZETAS += tinc * diff_ZETAS;
  ZETAW += tinc * diff_ZETAW;

  // Active and Total Force
  Ta = h * PKAForceMultiplier * (v(VT_Tref) / v(VT_dr)) * ((ZETAS + 1.0) * XS + ZETAW * XW);

  // Minimal implementation of the passive cell model
  // Similar to a standard linear solid model. It is used for the viscoelastic response.
  double C_s     = (stretch - 1.0) - Cd;
  double eta     = (C_s > 0.0) ? v(VT_eta_l) : v(VT_eta_s);
  double diff_Cd = v(VT_k) * C_s / eta;
  Cd += tinc * diff_Cd;

  double F_d = v(VT_a) * v(VT_k) * C_s;

  // Total Tension
  Tension = Ta + F_d;
    

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Buffering
  ////////////////////////////////////////////////////////////////////////////////////////
  const ML_CalcType dBuffer_NaBj = kon_na * Na_dyad * (Bmax_Naj - Buffer_NaBj) - koff_na * Buffer_NaBj;
  Buffer_NaBj += tinc * dBuffer_NaBj;
  const ML_CalcType dBuffer_NaBsl = kon_na * Na_sl * (Bmax_Nasl - Buffer_NaBsl) - koff_na * Buffer_NaBsl;
  Buffer_NaBsl += tinc * dBuffer_NaBsl;
  const ML_CalcType dBuffer_TnClow = Bmax_TnClow * diff_TRPN;
  Buffer_TnClow += tinc * dBuffer_TnClow;
  const ML_CalcType dBuffer_TnCHc = kon_tnchca * Ca_i * (Bmax_TnChigh - Buffer_TnCHc - Buffer_TnCHm) - koff_tnchca * Buffer_TnCHc;
  Buffer_TnCHc += tinc * dBuffer_TnCHc;
  const ML_CalcType dBuffer_TnCHm = kon_tnchmg * Mg_myo * (Bmax_TnChigh - Buffer_TnCHc - Buffer_TnCHm) - koff_tnchmg * Buffer_TnCHm;
  Buffer_TnCHm += tinc * dBuffer_TnCHm;
  const ML_CalcType dBuffer_CaM = kon_cam * Ca_i * (Bmax_CaM - Buffer_CaM) - koff_cam * Buffer_CaM;
  Buffer_CaM += tinc * dBuffer_CaM;
  const ML_CalcType dBuffer_Myosin_ca = kon_myoca * Ca_i * (Bmax_myosin - Buffer_Myosin_ca - Buffer_Myosin_mg) - koff_myoca * Buffer_Myosin_ca;
  Buffer_Myosin_ca += tinc * dBuffer_Myosin_ca;
  const ML_CalcType dBuffer_Myosin_mg = kon_myomg * Mg_myo * (Bmax_myosin - Buffer_Myosin_ca - Buffer_Myosin_mg) - koff_myomg * Buffer_Myosin_mg;
  Buffer_Myosin_mg += tinc * dBuffer_Myosin_mg;
  const ML_CalcType dBuffer_SRB = kon_sr * Ca_i * (Bmax_SR - Buffer_SRB) - koff_sr * Buffer_SRB;
  Buffer_SRB += tinc * dBuffer_SRB;
  const ML_CalcType dBuffer_SLLj = kon_sll * Ca_dyad * (Bmax_SLlowj - Buffer_SLLj) - koff_sll * Buffer_SLLj;
  Buffer_SLLj += tinc * dBuffer_SLLj;
  const ML_CalcType dBuffer_SLLsl = kon_sll * Ca_sl * (Bmax_SLlowsl - Buffer_SLLsl) - koff_sll * Buffer_SLLsl;
  Buffer_SLLsl += tinc * dBuffer_SLLsl;
  const ML_CalcType dBuffer_SLHj = kon_slh * Ca_dyad * (Bmax_SLhighj - Buffer_SLHj) - koff_slh * Buffer_SLHj;
  Buffer_SLHj += tinc * dBuffer_SLHj;
  const ML_CalcType dBuffer_SLHsl = kon_slh * Ca_sl * (Bmax_SLhighsl - Buffer_SLHsl) - koff_slh * Buffer_SLHsl;
  Buffer_SLHsl += tinc * dBuffer_SLHsl;
  const ML_CalcType dBuffer_Csqn = kon_csqn * Ca_SR * (Bmax_Csqn - Buffer_Csqn) - koff_csqn * Buffer_Csqn;
  Buffer_Csqn += tinc * dBuffer_Csqn;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Diffusion
  ////////////////////////////////////////////////////////////////////////////////////////
  double J_Ca_dyad_sl = 1.0 / 3.06685e12;
  double J_Ca_sl_myo = 1.0 / 0.74556e11;
  double J_Na_dyad_sl = 1.0 / (1.6382e12 / 3.0 * 100.0);
  double J_Na_sl_myo = 1.0 / (1.8308e10 / 3.0 * 100.0);
    
  const ML_CalcType J_CaBuffer_myo = dBuffer_TnClow + dBuffer_TnCHc + dBuffer_CaM + dBuffer_Myosin_ca + dBuffer_SRB;
  const ML_CalcType J_CaBuffer_dyad = dBuffer_SLLj + dBuffer_SLHj;
  const ML_CalcType J_CaBuffer_sl = dBuffer_SLLsl + dBuffer_SLHsl;

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Total ion currents and concentration changes
  ////////////////////////////////////////////////////////////////////////////////////////
  // Sodium Concentration
  I_Na_tot_dyad = I_Na_dyad + I_Nab_dyad + 3 * I_NaCa_dyad + 3 * I_NaK_dyad + I_CaNa_dyad;
  I_Na_tot_sl = I_Na_sl + I_Nab_sl + 3 * I_NaCa_sl + 3 * I_NaK_sl + I_CaNa_sl;
  I_Na_tot = I_Na_tot_dyad + I_Na_tot_sl;
  const ML_CalcType dNa_dyad = -I_Na_tot_dyad * Cmem / (Vdyad * F) + J_Na_dyad_sl / Vdyad * (Na_sl - Na_dyad) - dBuffer_NaBj;
  Na_dyad += tinc * dNa_dyad;
  const ML_CalcType dNa_sl = -I_Na_tot_sl * Cmem / (Vsl * F) + J_Na_dyad_sl / Vsl * (Na_dyad - Na_sl) + J_Na_sl_myo / Vsl * (Na_myo - Na_sl) - dBuffer_NaBsl;
  Na_sl += tinc * dNa_sl;
  const ML_CalcType dNa_myo = J_Na_sl_myo / Vmyo * (Na_sl - Na_myo);
  Na_myo += tinc * dNa_myo;
    
  // Potassium Concentration
  I_K_tot = I_to + I_Kr + I_Ks + I_K1 - (2 * I_NaK) + (I_CaK_dyad + I_CaK_sl) + I_Kb + i_external;
  const ML_CalcType dK_myo = -I_K_tot * Cmem / (Vmyo * F);
  K_myo += tinc * dK_myo;
    
  // Cloride Concentration
  I_Cl_tot = I_CaCl + I_Clb;
  const ML_CalcType dCl_myo = -I_Cl_tot * Cmem / (-1.0 * Vmyo * F);
  Cl_myo += tinc * dCl_myo;
    
  // Calcium Concentration
  I_Ca_tot_dyad = I_CaL_dyad + I_Cab_dyad + I_pCa_dyad - (2 * I_NaCa_dyad);
  I_Ca_tot_sl = I_CaL_sl + I_Cab_sl + I_pCa_sl - (2 * I_NaCa_sl);
  I_Ca_tot = I_Ca_tot_dyad + I_Ca_tot_sl;
  const ML_CalcType dCa_dyad = -I_Ca_tot_dyad * Cmem / (Vdyad * 2.0 * F) + J_Ca_dyad_sl / Vdyad * (Ca_sl - Ca_dyad) - J_CaBuffer_dyad + J_SR_Carel * Vsr / Vdyad + J_SR_leak * Vmyo / Vdyad;
  Ca_dyad += tinc * dCa_dyad;
  const ML_CalcType dCa_sl = -I_Ca_tot_sl * Cmem / (Vsl * 2.0 * F) + J_Ca_dyad_sl / Vsl * (Ca_dyad - Ca_sl) + J_Ca_sl_myo / Vsl * (Ca_i - Ca_sl) - J_CaBuffer_sl;
  Ca_sl += tinc * dCa_sl;
  const ML_CalcType dCa_i = -J_up * Vsr / Vmyo - J_CaBuffer_myo + J_Ca_sl_myo / Vmyo * (Ca_sl - Ca_i);
  Ca_i += tinc * dCa_i;
  const ML_CalcType dCa_SR = J_up - (J_SR_leak * Vmyo / Vsr + J_SR_Carel) - dBuffer_Csqn;
  Ca_SR += tinc * dCa_SR;

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Change Membrane Potential (I_tot)
  ////////////////////////////////////////////////////////////////////////////////////////
  I_tot = -(I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot); // + i_external;

  return 0.001 * tinc * I_tot;
}  // TWorld::Calc

void TWorld::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr << tArg << ' ' << V << ' ' <<
    Na_dyad << ' ' << Na_sl << ' ' << Na_myo << ' ' << K_myo << ' ' << Cl_myo << ' ' << Ca_dyad << ' ' << Ca_sl << ' ' << Ca_i << ' ' << Ca_SR << ' ' << Ta << ' ' << Tension << ' '
    ;
}

void TWorld::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
    
  tempstr << CaMK_trap << ' ' << CaMK_f_ICaL << ' ' << CaMK_f_RyR << ' ' << CaMK_f_PLB << ' ' << casig_serca_trap << ' ' << m << ' ' << h << ' ' << j << ' ' << h_p << ' ' << j_p << ' ' << m_PKA << ' ' << h_PKA << ' ' << j_PKA << ' ' << h_both << ' ' << j_both << ' ' << m_L << ' ' << h_L << ' ' << h_L_p << ' ' << d << ' ' << f_fast << ' ' << f_slow << ' ' << f_Ca_fast << ' ' << f_Ca_slow << ' ' << j_Ca << ' ' << f_p_fast << ' ' << f_Ca_p_fast << ' ' << d_PKA << ' ' << f_PKA_fast << ' ' << f_PKA_slow << ' ' << f_Ca_PKA_fast << ' ' << f_Ca_PKA_slow << ' ' << f_both_fast << ' ' << f_Ca_both_fast << ' ' << n_Ca_dyad << ' ' << n_Ca_sl << ' ' << I_CaL_pureCDI_dyad << ' ' << I_CaL_pureCDI_sl << ' ' << a_slow << ' ' << a_fast << ' ' << i_slow << ' ' << i_fast << ' ' << a_p_slow << ' ' << a_p_fast << ' ' << i_p_slow << ' ' << i_p_fast << ' ' << C_0 << ' ' << C_1 << ' ' << C_2 << ' ' << O << ' ' << I << ' ' << xs_dyad << ' ' << xs_sl << ' ' << J_rel_ICaLdep_act << ' ' << J_rel_ICaLdep_f1 << ' ' << J_rel_ICaLdep_f2 << ' ' << ryr_R << ' ' << ryr_O << ' ' << ryr_I << ' ' << ryr_CaRI << ' ' << ryr_R_p << ' ' << ryr_O_p << ' ' << ryr_I_p << ' ' << ryr_CaRI_p << ' ' << Buffer_NaBj << ' ' << Buffer_NaBsl << ' ' << Buffer_TnClow << ' ' << Buffer_TnCHc << ' ' << Buffer_TnCHm << ' ' << Buffer_CaM << ' ' << Buffer_Myosin_ca << ' ' << Buffer_Myosin_mg << ' ' << Buffer_SRB << ' ' << Buffer_SLLj << ' ' << Buffer_SLLsl << ' ' << Buffer_SLHj << ' ' << Buffer_SLHsl << ' ' << Buffer_Csqn << ' ' << XS << ' ' << XW << ' ' << CaTRPN << ' ' << TmBlocked << ' ' << ZETAS << ' ' << ZETAW << ' ' << Cd << ' ';
}  // TWorld::LongPrint

void TWorld::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {
    "Na_dyad", "Na_sl", "Na_myo", "K_myo", "Cl_myo", "Ca_dyad", "Ca_sl", "Ca_i", "Ca_SR", "Ta", "Tension"
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}  // TWorld::GetParameterNames

void TWorld::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {
    "CaMK_trap","CaMK_f_ICaL","CaMK_f_RyR","CaMK_f_PLB","casig_serca_trap",
    "m","h","j","h_p","j_p","m_PKA","h_PKA","j_PKA","h_both","j_both","m_L","h_L","h_L_p",
    "d","f_fast","f_slow","f_Ca_fast","f_Ca_slow","j_Ca","f_p_fast","f_Ca_p_fast","d_PKA","f_PKA_fast","f_PKA_slow","f_Ca_PKA_fast","f_Ca_PKA_slow","f_both_fast","f_Ca_both_fast","n_Ca_dyad","n_Ca_sl","I_CaL_pureCDI_dyad","I_CaL_pureCDI_sl",
    "a_slow","a_fast","i_slow","i_fast","a_p_slow","a_p_fast","i_p_slow","i_p_fast",
    "C_0","C_1","C_2","O","I",
    "xs_dyad","xs_sl",
    "J_rel_ICaLdep_act","J_rel_ICaLdep_f1","J_rel_ICaLdep_f2","ryr_R","ryr_O","ryr_I","ryr_CaRI","ryr_R_p","ryr_O_p","ryr_I_p","ryr_CaRI_p",
    "Buffer_NaBj","Buffer_NaBsl","Buffer_TnClow","Buffer_TnCHc","Buffer_TnCHm","Buffer_CaM","Buffer_Myosin_ca","Buffer_Myosin_mg","Buffer_SRB","Buffer_SLLj","Buffer_SLLsl","Buffer_SLHj","Buffer_SLHsl","Buffer_Csqn",
    "XS","XW","CaTRPN","TmBlocked","ZETAS","ZETAW","Cd"
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
