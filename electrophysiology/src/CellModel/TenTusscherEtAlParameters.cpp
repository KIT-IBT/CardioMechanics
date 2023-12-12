/*
 * File: TenTusscherEtAlParameters.cpp
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



#include <TenTusscherEtAlParameters.h>
TenTusscherEtAlParameters::TenTusscherEtAlParameters(const char *initFile, ML_CalcType tinc) {
  // Konstruktor
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

TenTusscherEtAlParameters::~TenTusscherEtAlParameters() {
  // Destruktor
}

void TenTusscherEtAlParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"TenTusscherEtAlParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout <<     "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void TenTusscherEtAlParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "#Loading the TenTusscherEtAl parameter from " << initFile << " ...\n";
#endif  // if KADEBUG


  // Initialization of the Parameters ...
  P[VT_R].name            =   "R";
  P[VT_Tx].name           =  "Tx";
  P[VT_F].name            =   "F";
  P[VT_m_init].name       =      "init_m";
  P[VT_h_init].name       =      "init_h";
  P[VT_j_init].name       =      "init_j";
  P[VT_xr1_init].name     =    "init_Xr1";
  P[VT_xr2_init].name     =    "init_Xr2";
  P[VT_xs_init].name      =     "init_Xs";
  P[VT_r_init].name       =      "init_r";
  P[VT_s_init].name       =      "init_s";
  P[VT_d_init].name       =      "init_d";
  P[VT_f_init].name       =      "init_f";
  P[VT_fCa_init].name     =    "init_fCa";
  P[VT_g_init].name       =      "init_g";
  P[VT_K_o].name          = "K_o";
  P[VT_Ca_o].name         =        "Ca_o";
  P[VT_Na_o].name         =        "Na_o";
  P[VT_Vc].name           =  "Vc";
  P[VT_Vsr].name          = "Vsr";
  P[VT_Bufc].name         =        "Bufc";
  P[VT_Kbufc].name        =       "Kbufc";
  P[VT_Bufsr].name        =       "Bufsr";
  P[VT_Kbufsr].name       =      "Kbufsr";
  P[VT_taufca].name       =      "taufca";
  P[VT_taug].name         =        "taug";
  P[VT_Vmaxup].name       =      "Vmaxup";
  P[VT_Kup].name          = "Kup";
  P[VT_C].name            =   "C";
  P[VT_g_Kr].name         =        "g_Kr";
  P[VT_pKNa].name         =        "pKNa";
  P[VT_g_Ks].name         =        "g_Ks";
  P[VT_g_K1].name         =        "g_K1";
  P[VT_g_to].name         =        "g_to";
  P[VT_g_Na].name         =        "g_Na";
  P[VT_g_bNa].name        =       "g_bNa";
  P[VT_KmK].name          = "KmK";
  P[VT_KmNa].name         =        "KmNa";
  P[VT_knak].name         =        "knak";
  P[VT_g_CaL].name        =       "g_CaL";
  P[VT_g_bCa].name        =       "g_bCa";
  P[VT_kNaCa].name        =       "knaca";
  P[VT_KmNai].name        =       "KmNai";
  P[VT_KmCa].name         =        "KmCa";
  P[VT_ksat].name         =        "ksat";
  P[VT_n].name            =   "n";
  P[VT_g_pCa].name        =       "g_pCa";
  P[VT_KpCa].name         =        "KpCa";
  P[VT_g_pK].name         =        "g_pK";
  P[VT_V_init].name       =      "init_V_m";
  P[VT_Cai_init].name     =    "init_Ca_i";
  P[VT_CaSR_init].name    =   "init_Ca_SR";
  P[VT_Nai_init].name     =    "init_Na_i";
  P[VT_Ki_init].name      =     "init_K_i";
  P[VT_s_inf_vHalf].name  = "s_inf_vHalf";
  P[VT_tau_s_f1].name     =    "tau_s_f1";
  P[VT_tau_s_slope1].name =        "tau_s_slope1";
  P[VT_tau_s_vHalf1].name =        "tau_s_vHalf1";
  P[VT_tau_s_f2].name     =    "tau_s_f2";
  P[VT_tau_s_f3].name     =    "tau_s_f3";
  P[VT_m_Xr1_1].name      =     "m_Xr1_1";
  P[VT_m_Xr1_2].name      =     "m_Xr1_2";
  P[VT_a_Xr1_1].name      =     "a_Xr1_1";
  P[VT_a_Xr1_2].name      =     "a_Xr1_2";
  P[VT_b_Xr1_1].name      =     "b_Xr1_1";
  P[VT_b_Xr1_2].name      =     "b_Xr1_2";
  P[VT_K_Q10Xr1].name     =    "K_Q10Xr1";
  P[VT_m_Xr2_1].name      =     "m_Xr2_1";
  P[VT_m_Xr2_2].name      =     "m_Xr2_2";
  P[VT_a_Xr2_1].name      =     "a_Xr2_1";
  P[VT_a_Xr2_2].name      =     "a_Xr2_2";
  P[VT_b_Xr2_1].name      =     "b_Xr2_1";
  P[VT_b_Xr2_2].name      =     "b_Xr2_2";
  P[VT_K_Q10Xr2].name     =    "K_Q10Xr2";
  P[VT_inverseVcF2].name  = "inverseVcF2";
  P[VT_inverseVcF2C].name =        "inverseVcF2C";
  P[VT_inverseVcFC].name  = "inverseVcFC";
  P[VT_VcdVsr].name       =      "VcdVsr";
  P[VT_Kupsquare].name    =   "Kupsquare";
  P[VT_BufcPKbufc].name   =  "BufcPKbufc";
  P[VT_Kbufcsquare].name  = "Kbufcsquare";
  P[VT_Kbufc2].name       =      "Kbufc2";
  P[VT_BufsrPKbufsr].name =        "BufsrPKbufsr";
  P[VT_Kbufsrsquare].name =        "Kbufsrsquare";
  P[VT_Kbufsr2].name      =     "Kbufsr2";
  P[VT_KopKNaNao].name    =   "KopKNaNao";
  P[VT_KmNai3].name       =      "KmNai3";
  P[VT_Nao3].name         =        "Nao3";
  P[VT_inverseRTONF].name =        "inverseRTONF";
  P[VT_RTONF].name        =       "RTONF";
  P[VT_Amp].name          =         "Amp";

  P[VT_RTONF].readFromFile        = false;
  P[VT_inverseRTONF].readFromFile = false;
  P[VT_inverseVcF2C].readFromFile = false;
  P[VT_inverseVcF2].readFromFile  = false;
  

  // inverseVcF2 is a unused value! (dw)
  P[VT_inverseVcFC].readFromFile  = false;
  P[VT_VcdVsr].readFromFile       = false;
  P[VT_Kupsquare].readFromFile    = false;
  P[VT_BufcPKbufc].readFromFile   = false;
  P[VT_Kbufcsquare].readFromFile  = false;
  P[VT_Kbufc2].readFromFile       = false;
  P[VT_BufsrPKbufsr].readFromFile = false;
  P[VT_Kbufsrsquare].readFromFile = false;
  P[VT_Kbufsr2].readFromFile      = false;
  P[VT_KopKNaNao].readFromFile    = false;
  P[VT_KmNai3].readFromFile       = false;
  P[VT_Nao3].readFromFile         = false;

  ParameterLoader EPL(initFile, EMT_TenTusscher);
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
}  // TenTusscherEtAlParameters::Init

void TenTusscherEtAlParameters::Calculate() {
#if KADEBUG
  cerr << "#TenTusscherEtAlParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_RTONF].value        = P[VT_R].value*(P[VT_Tx].value)/(P[VT_F].value);
  P[VT_inverseRTONF].value = 1/(P[VT_RTONF].value);
  P[VT_inverseVcF2C].value = (1/(2*(P[VT_Vc].value)*(P[VT_F].value)))*(P[VT_C].value);
  P[VT_inverseVcFC].value  = (1./((P[VT_Vc].value)*(P[VT_F].value)))*(P[VT_C].value);
  P[VT_VcdVsr].value       = P[VT_Vc].value/(P[VT_Vsr].value);
  P[VT_Kupsquare].value    = P[VT_Kup].value*(P[VT_Kup].value);
  P[VT_BufcPKbufc].value   = P[VT_Bufc].value+(P[VT_Kbufc].value);
  P[VT_Kbufcsquare].value  = P[VT_Kbufc].value*(P[VT_Kbufc].value);
  P[VT_Kbufc2].value       = 2*(P[VT_Kbufc].value);
  P[VT_BufsrPKbufsr].value = P[VT_Bufsr].value+(P[VT_Kbufsr].value);
  P[VT_Kbufsrsquare].value = P[VT_Kbufsr].value*(P[VT_Kbufsr].value);
  P[VT_Kbufsr2].value      = 2*(P[VT_Kbufsr].value);
  P[VT_KopKNaNao].value    = P[VT_K_o].value+(P[VT_pKNa].value)*(P[VT_Na_o].value);
  P[VT_KmNai3].value       = P[VT_KmNai].value*(P[VT_KmNai].value)*(P[VT_KmNai].value);
  P[VT_Nao3].value         = P[VT_Na_o].value*(P[VT_Na_o].value)*(P[VT_Na_o].value);
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
}

void TenTusscherEtAlParameters::InitTable(ML_CalcType tinc) {
  tinc *= 1000.;  // sec. -> ms

#if KADEBUG
  cerr << "#TenTusscherEtAlParameters - InitTable ..." << endl;
#endif  // if KADEBUG
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    // V in mV
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);

    // XP[VT_F].value/(R*T)=X/((R*T)/F)=X/RTONF=XP[VT_inverseRTONF].value
    const ML_CalcType rec_iNaK =
      (1./(1.+0.1245*exp(-0.1*V*(P[VT_inverseRTONF].value))+0.0353*exp(-V*(P[VT_inverseRTONF].value))));
    NaK_P1[Vi] = P[VT_knak].value*(P[VT_K_o].value/(P[VT_K_o].value+(P[VT_KmK].value)))*rec_iNaK;

    //  const T rec_iNaK=(1./(1.+0.1245*exp(-0.1*svoltP[VT_F].value/(R*T))+0.0353*exp(-svoltP[VT_F].value/(R*T))));
    rec_ipK[Vi] = 1./(1.+exp((25-V)/5.98));

    //  const T rec_ipK=1./(1.+exp((25-svolt)/5.98));
    const ML_CalcType a_m = 1./(1.+exp((-60.-V)/5.));
    const ML_CalcType b_m = 0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    exptau_m[Vi] = exp(-tinc/(a_m*b_m));
    m_inf[Vi]    = 1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
    if (V >= -40.) {
      const ML_CalcType AH_1 = 0.;
      const ML_CalcType BH_1 = (0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      exptau_h[Vi] = exp(-tinc * (AH_1+BH_1));
    } else {
      const ML_CalcType AH_2 = (0.057*exp(-(V+80.)/6.8));
      const ML_CalcType BH_2 = (2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
      exptau_h[Vi] = exp(-tinc * (AH_2+BH_2));
    }
    h_inf[Vi] = 1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
    if (V >= -40.) {
      const ML_CalcType AJ_1 = 0.;
      const ML_CalcType BJ_1 = (0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
      exptau_j[Vi] = exp(-tinc * (AJ_1+BJ_1));
    } else {
      const ML_CalcType AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                                 exp(-0.04391*V))*(V+37.78)/
                                (1.+exp(0.311*(V+79.23))));
      const ML_CalcType BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
      exptau_j[Vi] = exp(-tinc * (AJ_2+BJ_2));
    }
    j_inf[Vi]   = h_inf[Vi];
    Xr1_inf[Vi] = 1./(1.0+exp((P[VT_m_Xr1_1].value-V)/(P[VT_m_Xr1_2].value)));
    const ML_CalcType a_Xr1 = 450./(1.+exp((P[VT_a_Xr1_1].value-V)/(P[VT_a_Xr1_2].value)));
    const ML_CalcType b_Xr1 = 6./(1.+exp((V-(P[VT_b_Xr1_1].value))/(P[VT_b_Xr1_2].value)));
    exptau_Xr1[Vi] = exp(-tinc/(a_Xr1*b_Xr1*fabs(P[VT_K_Q10Xr1].value)));
    Xr2_inf[Vi]    = 1./(1.+exp((V-(P[VT_m_Xr2_1].value))/(P[VT_m_Xr2_2].value)));
    const ML_CalcType a_Xr2 = 3./(1.+exp((P[VT_a_Xr2_1].value-V)/(P[VT_a_Xr2_2].value)));
    const ML_CalcType b_Xr2 = 1.12/(1.+exp((V-(P[VT_b_Xr2_1].value))/(P[VT_b_Xr2_2].value)));
    exptau_Xr2[Vi] = exp(-tinc / (a_Xr2*b_Xr2*fabs(P[VT_K_Q10Xr2].value)));
    Xs_inf[Vi]     = 1./(1.+exp((-5.-V)/14.));
    const ML_CalcType a_Xs = 1100./(sqrt(1.+exp((-10.-V)/6)));
    const ML_CalcType b_Xs = 1./(1.+exp((V-60.)/20.));
    exptau_Xs[Vi] = exp(-tinc / (a_Xs*b_Xs));
    r_inf[Vi]     = 1./(1.+exp((20-V)/6.));
    s_inf[Vi]     = 1./(1.+exp((V+(P[VT_s_inf_vHalf].value))/5.));
    exptau_r[Vi]  = exp(-tinc / (9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8));
    exptau_s[Vi]  =
      exp(-tinc /
          ((P[VT_tau_s_f1].value)*
           exp(-(V+(P[VT_tau_s_vHalf1].value))*(V+(P[VT_tau_s_vHalf1].value))/(P[VT_tau_s_slope1].value))+
           (P[VT_tau_s_f2].value)/(1.+exp((V-20.)/5.))+(P[VT_tau_s_f3].value)));
    d_inf[Vi] = 1./(1.+exp((-5-V)/7.5));
    const ML_CalcType a_d = 1.4/(1.+exp((-35-V)/13))+0.25;
    const ML_CalcType b_d = 1.4/(1.+exp((V+5)/5));
    const ML_CalcType c_d = 1./(1.+exp((50-V)/20));
    exptau_d[Vi] = exp(-tinc / (a_d*b_d+c_d));
    f_inf[Vi]    = 1./(1.+exp((V+20)/7));
    exptau_f[Vi] = exp(-tinc/(1125*exp(-(V+27)*(V+27)/300)+80+165/(1.+exp((25-V)/10))));
    const ML_CalcType NaCaP1 = (P[VT_kNaCa].value)*(1./((P[VT_KmNai3].value)+(P[VT_Nao3].value)))*
      (1./((P[VT_KmCa].value)+(P[VT_Ca_o].value)))*
      (1./(1.+(P[VT_ksat].value)*exp(((P[VT_n].value)-1.)*V*(P[VT_inverseRTONF].value))));
    const ML_CalcType NaCaP2 = exp((P[VT_n].value)*V*(P[VT_inverseRTONF].value))*(P[VT_Ca_o].value);
    const ML_CalcType NaCaP3 = exp(((P[VT_n].value)-1)*V*(P[VT_inverseRTONF].value))*(P[VT_Nao3].value)*2.5;

    // INaCa=NaCaP1*(NaCaP2*Na_I^3-NaCaP3P[VT_C].valuea_i)=NaCa_P1[Vi]*Na_i^3-NaCa_P2[Vi]P[VT_C].valuea_i ,
    // NaCa_P1[Vi]=NaCaP1*NaCaP2, NaCa_P2[Vi]=NaCaP1*NaCaP3
    NaCa_P1[Vi] = NaCaP1*NaCaP2;
    NaCa_P2[Vi] = NaCaP1*NaCaP3;
    const ML_CalcType CaLP1 = 2*V*(P[VT_inverseRTONF].value);
    const ML_CalcType CaLP2 = 2*(P[VT_g_CaL].value)*CaLP1*(P[VT_F].value)/(exp(CaLP1)-1.);
    CaL_P1[Vi] = CaLP2*exp(CaLP1);
    CaL_P2[Vi] = CaLP2* -0.341*(P[VT_Ca_o].value);  // ICaL=d*f*fCa*(CaL_P1[Vi]P[VT_C].valuea_i+CaL_P2[Vi])
  }
}  // TenTusscherEtAlParameters::InitTable
