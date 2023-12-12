/*
 * File: FabbriParameters.cpp
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


#include <FabbriParameters.h>

FabbriParameters::FabbriParameters(const char *initFile, ML_CalcType tinc) {
  // Konstruktor
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

FabbriParameters::~FabbriParameters() {
  // Destruktor
}

void FabbriParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the Fabbri parameter from " << initFile << " ...\n";
#endif  // if KADEBUG


  // Initialization of the Parameters ...
  P[VT_Vari_Nai].name = "Vari_Nai";
  P[VT_Vari_Ki].name = "Vari_Ki";
  P[VT_R].name = "R";
  P[VT_T].name = "T";
  P[VT_F].name = "F";
  P[VT_C_m].name = "C_m";
  P[VT_ACh].name = "ACh";
  P[VT_Iso_1_uM].name = "Iso_1_uM";
  P[VT_Nao].name = "Nao";
  P[VT_Init_Ki].name = "Init_Ki";
  P[VT_Ko].name = "Ko";
  P[VT_Cao].name = "Cao";
  P[VT_g_fNa].name = "g_fNa";
  P[VT_g_fK].name = "g_fK";
  P[VT_Km_Kp].name = "Km_Kp";
  P[VT_Km_Nap].name = "Km_Nap";
  P[VT_i_NaK_max].name = "i_NaK_max";
  P[VT_K_NaCa].name = "K_NaCa";
  P[VT_Qci].name = "Qci";
  P[VT_Qn].name = "Qn";
  P[VT_Qco].name = "Qco";
  P[VT_K3ni].name = "K3ni";
  P[VT_Kci].name = "Kci";
  P[VT_K1ni].name = "K1ni";
  P[VT_K2ni].name = "K2ni";
  P[VT_Kcni].name = "Kcni";
  P[VT_K3no].name = "K3no";
  P[VT_K1no].name = "K1no";
  P[VT_K2no].name = "K2no";
  P[VT_Kco].name = "Kco";
  P[VT_g_na].name = "g_na";
  P[VT_P_CaL].name = "P_CaL";
  P[VT_k_dl].name = "k_dl";
  P[VT_V_dl].name = "V_dl";
  P[VT_k_fl].name = "k_fl";
  P[VT_V_fl].name = "V_fl";
  P[VT_alpha_fCa].name = "alpha_fCa";
  P[VT_Km_fCa].name = "Km_fCa";
  P[VT_P_CaT].name = "P_CaT";
  P[VT_ks].name = "ks";
  P[VT_MaxSR].name = "MaxSR";
  P[VT_MinSR].name = "MinSR";
  P[VT_EC50_SR].name = "EC50_SR";
  P[VT_EC50_SK].name = "EC50_SK";
  P[VT_n_SK].name = "n_SK";
  P[VT_HSR].name = "HSR";
  P[VT_koCa].name = "koCa";
  P[VT_kiCa].name = "kiCa";
  P[VT_kim].name = "kim";
  P[VT_kom].name = "kom";
  P[VT_tau_dif_Ca].name = "tau_dif_Ca";
  P[VT_tau_tr].name = "tau_tr";
  P[VT_K_up].name = "K_up";
  P[VT_slope_up].name = "slope_up";
  P[VT_TC_tot].name = "TC_tot";
  P[VT_TMC_tot].name = "TMC_tot";
  P[VT_CM_tot].name = "CM_tot";
  P[VT_CQ_tot].name = "CQ_tot";
  P[VT_kf_TC].name = "kf_TC";
  P[VT_kf_TMM].name = "kf_TMM";
  P[VT_kf_TMC].name = "kf_TMC";
  P[VT_kf_CM].name = "kf_CM";
  P[VT_kf_CQ].name = "kf_CQ";
  P[VT_kb_TC].name = "kb_TC";
  P[VT_kb_TMC].name = "kb_TMC";
  P[VT_kb_TMM].name = "kb_TMM";
  P[VT_kb_CM].name = "kb_CM";
  P[VT_kb_CQ].name = "kb_CQ";
  P[VT_Init_Mgi].name = "Init_Mgi";
  P[VT_V_jsr_part].name = "V_jsr_part";
  P[VT_V_i_part].name = "V_i_part";
  P[VT_V_nsr_part].name = "V_nsr_part";
  P[VT_R_cell].name = "R_cell";
  P[VT_L_cell].name = "L_cell";
  P[VT_L_sub].name = "L_sub";
  P[VT_g_Kur].name = "g_Kur";
  P[VT_g_to].name = "g_to";
  P[VT_g_Kr].name = "g_Kr";
  P[VT_g_Ks].name = "g_Ks";
  P[VT_g_KACh].name = "g_KACh";
  P[VT_E_K].name = "E_K";
  P[VT_k34].name = "k34";
  P[VT_alpha_a].name = "alpha_a";
  P[VT_RTd2F].name = "RTd2F";
  P[VT_Init_Vm].name = "Init_Vm";
  P[VT_Init_Ca_sub].name = "Init_Ca_sub";
  P[VT_Init_Nai].name = "Init_Nai";
  P[VT_Init_y].name = "Init_y";
  P[VT_Init_m].name = "Init_m";
  P[VT_Init_h].name = "Init_h";
  P[VT_Init_dL].name = "Init_dL";
  P[VT_Init_fL].name = "Init_fL";
  P[VT_Init_fCa].name = "Init_fCa";
  P[VT_Init_dT].name = "Init_dT";
  P[VT_Init_fT].name = "Init_fT";
  P[VT_Init_R_Ca_rel].name = "Init_R_Ca_rel";
  P[VT_Init_O_Ca_rel].name = "Init_O_Ca_rel";
  P[VT_Init_I_Ca_rel].name = "Init_I_Ca_rel";
  P[VT_Init_RI_Ca_rel].name = "Init_RI_Ca_rel";
  P[VT_Init_Ca_jsr].name = "Init_Ca_jsr";
  P[VT_Init_Ca_nsr].name = "Init_Ca_nsr";
  P[VT_Init_Cai].name = "Init_Cai";
  P[VT_Init_fTMM].name = "Init_fTMM";
  P[VT_Init_fCMi].name = "Init_fCMi";
  P[VT_Init_fCMs].name = "Init_fCMs";
  P[VT_Init_fTC].name = "Init_fTC";
  P[VT_Init_fTMC].name = "Init_fTMC";
  P[VT_Init_fCQ].name = "Init_fCQ";
  P[VT_Init_r_Kur].name = "Init_r_Kur";
  P[VT_Init_s].name = "Init_s";
  P[VT_Init_q].name = "Init_q";
  P[VT_Init_r_to].name = "Init_r_to";
  P[VT_Init_paS].name = "Init_paS";
  P[VT_Init_paF].name = "Init_paF";
  P[VT_Init_piy].name = "Init_piy";
  P[VT_Init_n].name = "Init_n";
  P[VT_Init_a].name = "Init_a";
  P[VT_Init_x].name = "Init_x";
  P[VT_g_na_l].name = "g_na_l";
  P[VT_g_na_mut].name = "g_na_mut";
  P[VT_g_SK].name = "g_SK";

  P[VT_V_cell].name = "V_cell";
  P[VT_V_sub].name =  "V_sub";
  P[VT_V_jsr].name = "V_jsr";
  P[VT_V_i].name = "V_i";
  P[VT_V_nsr].name = "V_nsr";

  P[VT_V_cell].readFromFile = false;
  P[VT_V_sub].readFromFile = false;
  P[VT_V_jsr].readFromFile = false;
  P[VT_V_i].readFromFile = false;
  P[VT_V_nsr].readFromFile = false;

  P[VT_Iso_shift_If].name = "Iso_shift_If";
  P[VT_Iso_INak_increase].name = "Iso_INak_increase";
  P[VT_Iso_ICaL_increase].name = "Iso_ICaL_increase";
  P[VT_Iso_shift_dL].name = "Iso_shift_dL";
  P[VT_Iso_slope_dL].name = "Iso_slope_dL";
  P[VT_Iso_IKs_increase].name = "Iso_IKs_increase";
  P[VT_Iso_shift_n].name = "Iso_shift_n";
  P[VT_Iso_b_up].name = "Iso_b_up";
  P[VT_ACh_if_shift].name = "ACh_if_shift";
  P[VT_ACh_iCaL_block].name = "ACh_iCal_block";

  P[VT_Iso_shift_If].readFromFile = false;
  P[VT_Iso_INak_increase].readFromFile = false;
  P[VT_Iso_ICaL_increase].readFromFile = false;
  P[VT_Iso_shift_dL].readFromFile = false;
  P[VT_Iso_slope_dL].readFromFile = false;
  P[VT_Iso_IKs_increase].readFromFile = false;
  P[VT_Iso_shift_n].readFromFile = false;
  P[VT_Iso_b_up].readFromFile = false;
  P[VT_ACh_if_shift].readFromFile = false;
  P[VT_ACh_iCaL_block].readFromFile = false;

  P[VT_K_Iso_shift].name = "K_iso_shift";
  P[VT_K_Iso_increase].name = "K_iso_increase";
  P[VT_K_Iso_slope_dL].name = "K_iso_slope_dL";

  P[VT_y_inf_b].name = "y_inf_b";
  P[VT_y_inf_a].name = "y_inf_a";
  P[VT_y_inf_b].readFromFile = false;
  P[VT_y_inf_a].readFromFile = false;

  P[VT_V_0p5_y_inf].name = "V_0p5_y_inf";
  P[VT_tau_y_a_shift].name = "tau_y_a_shift";
  P[VT_tau_y_b_shift].name = "tau_y_b_shift";
  P[VT_V_0p5_m_inf].name = "V_0p5_m_inf";
  P[VT_k_m].name = "k_m";
  P[VT_tau_m_shift].name = "tau_m_shift";
  P[VT_V_0p5_h_inf].name = "V_0p5_h_inf";
  P[VT_k_h].name = "k_h";
  P[VT_tau_h_gain].name = "tau_h_gain";
  P[VT_V_0p5_n_inf].name = "V_0p5_n_inf";
  P[VT_k_n].name = "k_n";
  P[VT_R231C].name = "R231C";

  ParameterLoader EPL(initFile, EMT_FabbriEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile) {
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
    }
    // End Initialization of the Parameters ...
    else {
#if KADEBUG
#endif
    }
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ...\n";
#endif  // if KADEBUG
}  // FabbriParameters::Init

void FabbriParameters::Calculate() {
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
#if KADEBUG
  cerr << "FabbriParameters - Calculate ..." << endl;
#endif // if KADEBUG
  const double pi = 3.14159265359;


  P[VT_V_cell].value = 3.20150282549375e-06;  // 1e-9*pi*P[VT_R_cell].value*P[VT_R_cell].value*P[VT_L_cell].value;
  P[VT_V_sub].value =  3.27517317322043e-08;  // 1e-9*2*pi*P[VT_L_sub].value*(P[VT_R_cell].value-(P[VT_L_sub].value/2))*P[VT_L_cell].value;
  P[VT_V_jsr].value = 3.84180339059250e-09;   // P[VT_V_jsr_part].value*P[VT_V_cell].value;
  P[VT_V_i].value = 1.43993956799492e-06;     // P[VT_V_i_part].value*P[VT_V_cell].value-P[VT_V_sub].value;
  P[VT_V_nsr].value = 3.71374327757275e-08;   // P[VT_V_nsr_part].value*P[VT_V_cell].value;

  // Weighting the ISO Influence
  if ((P[VT_Iso_1_uM].value == 0) || (P[VT_K_Iso_increase].value == 0)) {
    P[VT_Iso_INak_increase].value = 1;
    P[VT_Iso_ICaL_increase].value = 1;
    P[VT_Iso_IKs_increase].value = 1;
  } else {
    P[VT_Iso_INak_increase].value = 1 + 0.2*P[VT_K_Iso_increase].value*P[VT_Iso_1_uM].value;
    P[VT_Iso_ICaL_increase].value = 1 + 0.23*P[VT_K_Iso_increase].value*P[VT_Iso_1_uM].value;
    P[VT_Iso_IKs_increase].value = 1 + 0.2*P[VT_K_Iso_increase].value*P[VT_Iso_1_uM].value;
  }
  P[VT_Iso_shift_If].value = 7.5*P[VT_K_Iso_shift].value*P[VT_Iso_1_uM].value;
  P[VT_Iso_shift_dL].value = -8*P[VT_K_Iso_shift].value*P[VT_Iso_1_uM].value; // -8*P[VT_K_Iso_shift].value*P[VT_Iso_1_uM].value;
  P[VT_Iso_slope_dL].value = -27*P[VT_K_Iso_slope_dL].value*P[VT_Iso_1_uM].value; // -27*P[VT_K_Iso_slope_dL].value*P[VT_Iso_1_uM].value;
  P[VT_Iso_shift_n].value = -14*P[VT_K_Iso_shift].value*P[VT_Iso_1_uM].value;
  P[VT_Iso_b_up].value = 0.25*P[VT_K_Iso_shift].value*P[VT_Iso_1_uM].value;

  // Weighting ACh - Influence
  if (P[VT_ACh].value > 0) {
    P[VT_ACh_iCaL_block].value = 0.31*P[VT_ACh].value/(P[VT_ACh].value+0.00009);
    P[VT_ACh_if_shift].value = (-1-(9.898*pow(P[VT_ACh].value, 0.618)/(pow(P[VT_ACh].value, 0.618)+0.00122423)));
  } else {
    P[VT_ACh_iCaL_block].value = 0;
    P[VT_ACh_if_shift].value = 0;
  }

  // calculate y_inf Parameter
  double slope_y_inf = -0.038877;
  P[VT_y_inf_b].value = -1/(2*slope_y_inf);
  P[VT_y_inf_a].value = 1/(2*exp(-P[VT_V_0p5_y_inf].value/P[VT_y_inf_b].value));

  // P[VT_y_inf_b].value= 12.861;
  // P[VT_y_inf_a].value= 0.0002501;
} // FabbriParameters::Calculate

void FabbriParameters::InitTable(ML_CalcType tinc) {
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);

    // Berechnung von RT/F

    const double RTdF = (P[VT_R].value*P[VT_T].value)/P[VT_F].value;


    // Gates
    //
    // y von I_f:
    // y_inf
    double y_inf_shift = P[VT_V_0p5_y_inf].value + 97.75063;
    if (V < -80+y_inf_shift+P[VT_Iso_shift_If].value+P[VT_ACh_if_shift].value) {
      m_y[Vi] = 0.01329+
        (0.99921/(1+exp(((V+97.134-y_inf_shift-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value)/8.1752))));
    } else {
      m_y[Vi] = P[VT_y_inf_a].value*exp(-(V-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value)/P[VT_y_inf_b].value);
    }
    double tau_y =
      (1/
       (((0.36*(V-P[VT_tau_y_a_shift].value-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value))/
         (exp(0.066*(V-P[VT_tau_y_a_shift].value-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value))-1))+
        ((0.1*(V-P[VT_tau_y_b_shift].value-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value))/
         (1-(exp((-0.2)*(V-P[VT_tau_y_b_shift].value-P[VT_Iso_shift_If].value-P[VT_ACh_if_shift].value)))))))-0.054;
    exptau_y[Vi] = exp(-tinc/tau_y);

    // do von I_NaCa
    doo[Vi] = 1+ (P[VT_Cao].value/P[VT_Kco].value)*(1+exp((P[VT_Qco].value*V)/RTdF))+(P[VT_Nao].value/P[VT_K1no].value)*
      (1+(P[VT_Nao].value/P[VT_K2no].value)*(1+(P[VT_Nao].value/P[VT_K3no].value)));

    // k41 von I_NaCa
    k41[Vi]  = exp(((-1*P[VT_Qn].value)*V)/(2*RTdF));

    //
    // k21 von I_NaCa
    k21[Vi] = ((P[VT_Cao].value/P[VT_Kco].value)*exp((P[VT_Qco].value*V)/RTdF))/doo[Vi];

    //

    // k23 von I_NaCa
    k23[Vi] =
      ((((P[VT_Nao].value/P[VT_K1no].value)*P[VT_Nao].value)/P[VT_K2no].value)*(1+(P[VT_Nao].value/P[VT_K3no].value))*
       exp(((-P[VT_Qn].value)*V)/(2*RTdF)))/doo[Vi];

    //
    // k32 von I_NaCa
    k32[Vi] = exp((P[VT_Qn].value*V)/(2*RTdF));

    //

    // m von I_Na:
    const double E0_m = V+41;

    double b = 8000*exp(-0.056*(V+66));
    double a = (200*E0_m)/(1-exp(-0.1*E0_m));
    m_m[Vi] = 1/(1+exp(-(V+42.0504)/8.3106));
    double tau_m = 1/(a+b);
    exptau_m[Vi] = exp(-tinc/tau_m);

    // m_mut von I_Na
    const double E0_m_mut = V-P[VT_tau_m_shift].value;

    a = (200*E0_m_mut)/(1-exp(-0.1*E0_m_mut));
    b = 8000*exp(-0.056*(E0_m_mut+25));
    m_m_mut[Vi] = 1/(1+exp(-(V-P[VT_V_0p5_m_inf].value)/P[VT_k_m].value));
    double tau_m_mut = 1/(a+b);
    exptau_m_mut[Vi] = exp(-tinc/tau_m_mut);

    // h von I_Na:
    a = 20*exp(-0.125*(V+75));
    b = 2000/(320*exp((-0.1)*(V+75))+1);
    m_h[Vi] = 1/(1+exp((V+69.804)/4.4565));
    double tau_h = 1/(a+b);
    exptau_h[Vi] = exp(-tinc/tau_h);

    // h_mut von I_Na
    m_h_mut[Vi] = 1/(1+exp((V-P[VT_V_0p5_h_inf].value)/P[VT_k_h].value));
    double tau_h_mut = P[VT_tau_h_gain].value/(a+b);
    exptau_h_mut[Vi] = exp(-tinc/tau_h_mut);


    // dL von I_CaL
    a = ((-0.02839*(V+41.8-P[VT_Iso_shift_dL].value))/(exp((-(V+41.8-P[VT_Iso_shift_dL].value))/2.5)-1))-
      ((0.0849*(V+6.8-P[VT_Iso_shift_dL].value))/(exp((-(V+6.8-P[VT_Iso_shift_dL].value))/4.8)-1));
    b = (0.01143*(V+1.8-P[VT_Iso_shift_dL].value))/(exp((V+1.8-P[VT_Iso_shift_dL].value)/2.5)-1);
    double tau_dL = 0.001/(a+b);
    m_dL[Vi] = 1/
      (1+exp(-(V-P[VT_V_dl].value-P[VT_Iso_shift_dL].value)/(P[VT_k_dl].value*(1+P[VT_Iso_slope_dL].value/100))));
    exptau_dL[Vi] = exp(-tinc/tau_dL);

    //
    // fL con I_CaL:
    m_fL[Vi] = 1/(1+exp((V-P[VT_V_fl].value)/P[VT_k_fl].value));
    double tau_fL = 0.001*(44.3+230*exp(-(pow(((V+36)/10), 2))));
    exptau_fL[Vi] = exp(-tinc/tau_fL);

    //
    // dT von I_CaT:
    m_dT[Vi] = 1/(1+exp(-(V+38.3)/5.5));
    double tau_dT = 0.001/(1.068*exp((V+38.3)/30)+1.068*exp(-(V+38.3)/30));
    exptau_dT[Vi] = exp(-tinc/tau_dT);

    //
    // fT von I_CaT:
    m_fT[Vi] = 1/(1+exp((V+58.7)/3.8));
    double tau_fT =  (1/(16.67*exp(-(V+75)/83.3)+16.67*exp((V+75)/15.38)));
    exptau_fT[Vi] = exp(-tinc/tau_fT);

    //
    // r_Kur von I_Kur:
    m_r_Kur[Vi] = 1/(1+exp((V+6)/ -8.6));
    double tau_r_Kur = (0.009/(1+exp((V+5)/12)))+0.0005;
    exptau_r_Kur[Vi] = exp(-tinc/tau_r_Kur);

    //
    // s von I_Kur:
    m_s[Vi] = 1/(1+exp((V+7.5)/10));
    double tau_s = (0.59/(1+exp((V+60)/10)))+3.05;
    exptau_s[Vi] = exp(-tinc/tau_s);

    //
    // q von I_to:
    m_q[Vi] = 1/(1+exp((V+49)/13));
    double tau_q = 0.001*0.6*(65.17/(0.57*exp(-0.08*(V+44))+0.065*exp(0.1*(V+45.93)))+10.1);
    exptau_q[Vi] = exp(-tinc/tau_q);

    //
    // r_to von I_to:
    m_r_to[Vi] = 1/(1+exp(-(V-19.3)/15));
    double tau_r_to = 0.001*0.66*1.4*(15.59/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84)))+2.98);
    exptau_r_to[Vi] = exp(-tinc/tau_r_to);

    //
    // paF von I_Kr:
    m_pa[Vi] = 1/(1+exp(-(V+10.0144)/7.6607));
    m_pi[Vi] = 1/(1+exp((V+28.6)/17.1));
    double tau_paS = 0.84655354/(4.2*exp(V/17)+0.15*exp(-V/21.6));
    exptau_paS[Vi] = exp(-tinc/tau_paS);
    double tau_paF = 1/(30*exp(V/10)+exp(-V/12));
    exptau_paF[Vi] = exp(-tinc/tau_paF);
    double tau_pi = 1/(100*exp(-V/54.645)+656*exp(V/106.157));
    exptau_pi[Vi] = exp(-tinc/tau_pi);

    //
    // n von I_Ks:
    a = 28/(1+exp(-(V-40-P[VT_Iso_shift_n].value)/3));
    b = exp(-(V-5-P[VT_Iso_shift_n].value)/25);

    // m_n[Vi]=sqrt(1/(1+exp(-(V+0.6383)/10.7071)));
    m_n[Vi] = sqrt(1/(1+exp(-(V-P[VT_V_0p5_n_inf].value-P[VT_Iso_shift_n].value)/P[VT_k_n].value)));
    double tau_n = 1/(a+b);
    exptau_n[Vi] = exp(-tinc/tau_n);

    //

    // a von I_KACh:
    a = (3.5988-0.025641)/(1+(0.0000012155/pow(P[VT_ACh].value, 1.6951)))+0.025641;
    b = 10*exp(0.0133*(V+40));
    m_a[Vi] = a/(a+b);
    double tau_a = 1/(a+b);
    exptau_a[Vi] = exp(-tinc/tau_a);
  }
}  // FabbriParameters::InitTable

void FabbriParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"FabbriParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}
