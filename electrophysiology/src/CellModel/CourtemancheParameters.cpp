/* -------------------------------------------------------

   CourtemancheParameters.cpp

   Ver. 1.1.0

   Created:       dw (27.02.2007)
   Last modified: Tobias Gerach (07.05.2023)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */


#include <CourtemancheParameters.h>
CourtemancheParameters::CourtemancheParameters(const char *initFile, ML_CalcType tinc) {
  // Konstruktor
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

CourtemancheParameters::~CourtemancheParameters() {
  // Destruktor
}

void CourtemancheParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the Courtemanche parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  // Initialization of the Parameters ...
  P[VT_R].name             =   "R";
  P[VT_Tx].name            =  "Tx";
  P[VT_F].name             =   "F";
  P[VT_C_m].name           = "C_m";
  P[VT_Vcell].name         =       "Vcell";
  P[VT_Vcell_outside].name =     "Vcell_outside";
  P[VT_part_i].name        =      "part_i";
  P[VT_part_up].name       =     "part_up";
  P[VT_part_rel].name      =    "part_rel";
  P[VT_Tfac].name          =        "Tfac";
  P[VT_K_o].name           = "K_o";
  P[VT_Na_o].name          =        "Na_o";
  P[VT_Ca_o].name          =        "Ca_o";
  P[VT_g_Na].name          =        "g_Na";
  P[VT_g_to].name          =        "g_to";
  P[VT_g_Kr].name          =        "g_Kr";
  P[VT_g_Ks].name          =        "g_Ks";
  P[VT_g_CaL].name         =       "g_CaL";
  P[VT_g_K1].name          =        "g_K1";
  P[VT_g_bNa].name         =       "g_bNa";
  P[VT_g_bCa].name         =       "g_bCa";
  P[VT_I_NaKmax].name      =    "I_NaKmax";
  P[VT_k_mNai].name        =      "k_mNai";
  P[VT_k_mKo].name         =       "k_mKo";
  P[VT_I_pCamax].name      =    "I_pCamax";
  P[VT_k_mpCa].name        =      "k_mpCa";
  P[VT_k_NaCa].name        =      "k_NaCa";
  P[VT_k_mNa].name         =       "k_mNa";
  P[VT_k_mCa].name         =       "k_mCa";
  P[VT_k_sat].name         =       "k_sat";
  P[VT_gamm].name          =        "gamm";
  P[VT_CMDNmax].name       =     "CMDNmax";
  P[VT_TRPNmax].name       =     "TRPNmax";
  P[VT_CSQNmax].name       =     "CSQNmax";
  P[VT_K_mcmdn].name       =     "K_mcmdn";
  P[VT_K_mtrpn].name       =     "K_mtrpn";
  P[VT_K_mcsqn].name       =     "K_mcsqn";
  P[VT_k_rel].name         =       "k_rel";
  P[VT_k_up].name          =        "k_up";
  P[VT_I_upmax].name       =     "I_upmax";
  P[VT_Ca_upmax].name      =    "Ca_upmax";
  P[VT_t_tr].name          =        "t_tr";
  P[VT_shiftm].name        =      "shiftm";
  P[VT_shifth].name        =      "shifth";
  P[VT_shiftj].name        =      "shiftj";
  P[VT_shiftoa].name       =     "shiftoa";
  P[VT_C_oa_1].name        =      "C_oa_1";
  P[VT_C_oa_2].name        =      "C_oa_2";
  P[VT_C_oa_3].name        =      "C_oa_3";
  P[VT_C_oa_4].name        =      "C_oa_4";
  P[VT_C_oa_5].name        =  "C_oa_5";
  P[VT_C_oa_6].name        =  "C_oa_6";
  P[VT_CaL_fT_1].name      =        "CaL_fT_1";
  P[VT_Init_Na_i].name     =   "Init_Na_i";
  P[VT_Init_K_i].name      =    "Init_K_i";
  P[VT_Init_m].name        =      "Init_m";
  P[VT_Init_oa].name       =     "Init_oa";
  P[VT_Init_ua].name       =     "Init_ua";
  P[VT_Init_d].name        =      "Init_d";
  P[VT_Init_f_Ca].name     =   "Init_f_Ca";
  P[VT_Init_Ca_i].name     =   "Init_Ca_i";
  P[VT_Init_Ca_up].name    =  "Init_Ca_up";
  P[VT_Init_Ca_rel].name   = "Init_Ca_rel";
  P[VT_Init_h].name        =      "Init_h";
  P[VT_Init_j].name        =      "Init_j";
  P[VT_Init_oi].name       =     "Init_oi";
  P[VT_Init_ui].name       =     "Init_ui";
  P[VT_Init_Xr].name       =     "Init_Xr";
  P[VT_Init_Xs].name       =     "Init_Xs";
  P[VT_Init_f].name        =      "Init_f";
  P[VT_Init_u].name        =      "Init_u";
  P[VT_Init_v].name        =      "Init_v";
  P[VT_Init_w].name        =      "Init_w";
  P[VT_Init_Vm].name       =     "Init_Vm";
  P[VT_ua_a1].name         =       "ua_a1";
  P[VT_ua_a2].name         =       "ua_a2";
  P[VT_ua_a3].name         =       "ua_a3";
  P[VT_ua_a4].name         =       "ua_a4";
  P[VT_ua_a5].name         =       "ua_a5";
  P[VT_ua_b1].name         =       "ua_b1";
  P[VT_ua_b2].name         =       "ua_b2";
  P[VT_ua_b3].name         =       "ua_b3";
  P[VT_ua_m1].name         =       "ua_m1";
  P[VT_ua_m2].name         =       "ua_m2";
  P[VT_ua_KQ10].name       =     "ua_KQ10";
  P[VT_ui_a1].name         =       "ui_a1";
  P[VT_ui_a2].name         =       "ui_a2";
  P[VT_ui_a3].name         =       "ui_a3";
  P[VT_ui_a4].name         =       "ui_a4";
  P[VT_ui_b1].name         =       "ui_b1";
  P[VT_ui_b2].name         =       "ui_b2";
  P[VT_ui_m1].name         =       "ui_m1";
  P[VT_ui_m2].name         =       "ui_m2";
  P[VT_ui_KQ10].name       =     "ui_KQ10";
  P[VT_g_Kur1].name        =      "g_Kur1";
  P[VT_g_Kur2].name        =      "g_Kur2";
  P[VT_g_Kur3].name        =      "g_Kur3";
  P[VT_g_Kur4].name        =      "g_Kur4";
  P[VT_g_K1_1].name        =      "g_K1_1";
  P[VT_g_K1_2].name        =      "g_K1_2";
  P[VT_Xr_a1].name         =       "Xr_a1";
  P[VT_Xr_a2].name         =       "Xr_a2";
  P[VT_Xr_a3].name         =       "Xr_a3";
  P[VT_Xr_b1].name         =       "Xr_b1";
  P[VT_Xr_b2].name         =       "Xr_b2";
  P[VT_Xr_KQ10].name       =     "Xr_KQ10";
  P[VT_Xr_m1].name         =       "Xr_m1";
  P[VT_Xr_m2].name         =       "Xr_m2";
  P[VT_g_Kr1].name         =       "g_Kr1";
  P[VT_g_Kr2].name         =       "g_Kr2";
  P[VT_Xs_a1].name         =       "Xs_a1";
  P[VT_Xs_a2].name         =       "Xs_a2";
  P[VT_Xs_a3].name         =       "Xs_a3";
  P[VT_Xs_b1].name         =       "Xs_b1";
  P[VT_Xs_b2].name         =       "Xs_b2";
  P[VT_Xs_KQ10].name       =     "Xs_KQ10";
  P[VT_Xs_m1].name         =       "Xs_m1";
  P[VT_Xs_m2].name         =       "Xs_m2";
  P[VT_m_a1].name          =        "m_a1";
  P[VT_m_a2].name          =        "m_a2";
  P[VT_m_b1].name          =        "m_b1";
  P[VT_m_b2].name          =        "m_b2";
  P[VT_h_a1].name          =        "h_a1";
  P[VT_h_a2].name          =        "h_a2";
  P[VT_h_b1].name          =        "h_b1";
  P[VT_h_b2].name          =        "h_b2";
  P[VT_h_b3].name          =        "h_b3";
  P[VT_h_b4].name          =        "h_b4";
  P[VT_h_b5].name          =        "h_b5";
  P[VT_h_b6].name          =        "h_b6";
  P[VT_j_a1].name          =        "j_a1";
  P[VT_j_a2].name          =        "j_a2";
  P[VT_j_a3].name          =        "j_a3";
  P[VT_j_a4].name          =        "j_a4";
  P[VT_j_a5].name          =        "j_a5";
  P[VT_j_b1].name          =        "j_b1";
  P[VT_j_b2].name          =        "j_b2";
  P[VT_j_b3].name          =        "j_b3";
  P[VT_j_b4].name          =        "j_b4";
  P[VT_j_b5].name          =        "j_b5";
  P[VT_j_b6].name          =        "j_b6";
  P[VT_tau_j_mult].name    =        "tau_j_mult";
  P[VT_tau_f_mult].name    =        "tau_f_mult";
  P[VT_E_rL].name          =        "E_rL";
#ifdef ISAC
  P[VT_ISAC_SWITCH].name   =     "ISAC_SWITCH";
  P[VT_alpha].name         =     "alpha";
  P[VT_G_sac].name             =     "G_sac";
  P[VT_K_sac].name             =     "K_sac";
  P[VT_E_sac].name         =     "E_sac";
#endif // ifdef ISAC

#ifdef TRPN
  P[VT_Init_Ca_TRPN].name = "Init_Ca_TRPN";
  P[VT_beta1].name        = "beta1";
  P[VT_Ca_T50].name       = "Ca_T50";
  P[VT_k_TRPN].name       = "k_TRPN";
  P[VT_n_TRPN].name       = "n_TRPN";
#endif  // if TRPN

  P[VT_RTdF].name     =        "RTdF";
  P[VT_FdRT].name     =        "FdRT";
  P[VT_RTd2F].name    =       "RTd2F";
  P[VT_Vi].name       =          "Vi";
  P[VT_Vup].name      =         "Vup";
  P[VT_Vrel].name     =        "Vrel";
  P[VT_VupdVi].name   =      "VupdVi";
  P[VT_VreldVi].name  =     "VreldVi";
  P[VT_VreldVup].name =    "VreldVup";
  P[VT_CmdFvi].name   =      "CmdFvi";
  P[VT_Cmd2F].name    =       "Cmd2F";
  P[VT_csqnkm].name   =      "csqnkm";
  P[VT_kmcsqnm2].name =    "kmcsqnm2";
  P[VT_kkmcsqn].name  =     "kkmcsqn";
  P[VT_cmdnkm].name   =      "cmdnkm";
  P[VT_kmcmdnm2].name =    "kmcmdnm2";
  P[VT_kkmcmdn].name  =     "kkmcmdn";
  P[VT_trpnkm].name   =      "trpnkm";
  P[VT_kmtrpnm2].name =    "kmtrpnm2";
  P[VT_kkmtrpn].name  =     "kkmtrpn";
  P[VT_kupleak].name  =     "kupleak";


  P[VT_RTdF].readFromFile     = false;
  P[VT_FdRT].readFromFile     = false;
  P[VT_RTd2F].readFromFile    = false;
  P[VT_Vi].readFromFile       = false;
  P[VT_Vup].readFromFile      = false;
  P[VT_Vrel].readFromFile     = false;
  P[VT_VupdVi].readFromFile   = false;
  P[VT_VreldVi].readFromFile  = false;
  P[VT_VreldVup].readFromFile = false;
  P[VT_CmdFvi].readFromFile   = false;
  P[VT_Cmd2F].readFromFile    = false;
  P[VT_csqnkm].readFromFile   = false;
  P[VT_kmcsqnm2].readFromFile = false;
  P[VT_kkmcsqn].readFromFile  = false;
  P[VT_cmdnkm].readFromFile   = false;
  P[VT_kmcmdnm2].readFromFile = false;
  P[VT_kkmcmdn].readFromFile  = false;
  P[VT_trpnkm].readFromFile   = false;
  P[VT_kmtrpnm2].readFromFile = false;
  P[VT_kkmtrpn].readFromFile  = false;
  P[VT_kupleak].readFromFile  = false;
  P[VT_dt_tr].readFromFile    = false;

  ParameterLoader EPL(initFile, EMT_CourtemancheEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  // End Initialization of the Parameters ...

  exp250 = exp(250.0);
  exp50  = exp(169.5);
  Calculate();
  InitTable(tinc);
}  // CourtemancheParameters::Init

void CourtemancheParameters::Calculate() {
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
#if KADEBUG
  cerr << "CourtemancheParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_RTdF].value     = P[VT_R].value*(P[VT_Tx].value)/(P[VT_F].value);
  P[VT_FdRT].value     = 1.0/(P[VT_RTdF].value);
  P[VT_RTd2F].value    = P[VT_RTdF].value*.5;
  P[VT_Vi].value       = P[VT_Vcell].value*(P[VT_part_i].value);
  P[VT_Vup].value      = P[VT_Vcell].value*(P[VT_part_up].value);
  P[VT_Vrel].value     = P[VT_Vcell].value*(P[VT_part_rel].value);
  P[VT_VupdVi].value   = P[VT_Vup].value/(P[VT_Vi].value);
  P[VT_VreldVi].value  = P[VT_Vrel].value/(P[VT_Vi].value);
  P[VT_VreldVup].value = P[VT_Vrel].value/(P[VT_Vup].value);
  P[VT_CmdFvi].value   = P[VT_C_m].value/(P[VT_F].value*(P[VT_Vi].value));
  P[VT_Cmd2F].value    = P[VT_C_m].value/(P[VT_F].value*2.0);
  P[VT_csqnkm].value   = P[VT_CSQNmax].value*P[VT_K_mcsqn].value;
  P[VT_kmcsqnm2].value = P[VT_K_mcsqn].value*2.0;
  P[VT_kkmcsqn].value  = P[VT_K_mcsqn].value*P[VT_K_mcsqn].value;
  P[VT_cmdnkm].value   = P[VT_CMDNmax].value*P[VT_K_mcmdn].value;
  P[VT_kmcmdnm2].value = P[VT_K_mcmdn].value*2.0;
  P[VT_kkmcmdn].value  = P[VT_K_mcmdn].value*P[VT_K_mcmdn].value;
  P[VT_trpnkm].value   = P[VT_TRPNmax].value*P[VT_K_mtrpn].value;
  P[VT_kmtrpnm2].value = P[VT_K_mtrpn].value*2.0;
  P[VT_kkmtrpn].value  = P[VT_K_mtrpn].value*P[VT_K_mtrpn].value;
  P[VT_kupleak].value  = P[VT_I_upmax].value/(P[VT_Ca_upmax].value);
  P[VT_dt_tr].value    = 1.0/(P[VT_t_tr].value);
}  // CourtemancheParameters::Calculate

void CourtemancheParameters::InitTable(ML_CalcType tinc) {
  tinc *= 1000;

  double sigma = (exp(P[VT_Na_o].value/67.3)-1.0)/7.0;

  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi   = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    double a =
      ((fabs(V+47.13-P[VT_shiftm].value) <
        1e-10 ? 3.2 : P[VT_m_a1].value*(V+47.13-P[VT_shiftm].value)/
        (1.0-exp(-(V+47.13-P[VT_shiftm].value)*P[VT_m_a2].value)) ));
    double b = P[VT_m_b1].value*exp(-(V-P[VT_shiftm].value)/P[VT_m_b2].value);
    exptau_m[Vi] = exp(-tinc * (a+b));
    m_m[Vi]      = a/(a+b);

    if (V >= -40.0+P[VT_shifth].value) {
      a = .0;
      b = 1.0/(P[VT_h_b5].value*(1.0+exp(-(V+10.66-P[VT_shifth].value)/P[VT_h_b6].value)));
    } else {
      a = P[VT_h_a1].value*exp(-(V+80.0-P[VT_shifth].value)/P[VT_h_a2].value);
      b = P[VT_h_b1].value*exp(P[VT_h_b2].value*(V-P[VT_shifth].value))+P[VT_h_b3].value*
        exp(P[VT_h_b4].value*(V-P[VT_shifth].value));
    }
    exptau_h[Vi] = exp(-tinc * (a+b));
    m_h[Vi]      = a/(a+b);

    if (V >= -40.0+P[VT_shiftj].value) {
      a = .0;
      b = P[VT_j_b4].value*exp(P[VT_j_b5].value*(V-P[VT_shiftj].value))/
        (1.0+exp(P[VT_j_b6].value*(V+32.0-P[VT_shiftj].value)));
    } else {
      a =
        (P[VT_j_a1].value*exp(P[VT_j_a2].value*(V-P[VT_shiftj].value))-(P[VT_j_a3].value)*
         exp(P[VT_j_a4].value*(V-P[VT_shiftj].value)))*(V+37.78-P[VT_shiftj].value)/
        (1.0+exp(P[VT_j_a5].value*(V+79.23-P[VT_shiftj].value)));
      b = P[VT_j_b1].value*exp(P[VT_j_b2].value*(V-P[VT_shiftj].value))/
        (1.0+exp(P[VT_j_b3].value*(V+40.14-P[VT_shiftj].value)));
    }
    exptau_j[Vi] = exp(-tinc * (a+b) / P[VT_tau_j_mult].value);
    m_j[Vi]      = a/(a+b);

    a = P[VT_C_oa_1].value/
      ((exp(-(V-P[VT_shiftoa].value+10.0)/(P[VT_C_oa_2].value)))+
       (exp(-(V-P[VT_shiftoa].value-30.0)/(P[VT_C_oa_3].value))));
    b             = P[VT_C_oa_1].value/(2.5+exp((V-P[VT_shiftoa].value+82.0)/(P[VT_C_oa_4].value)));
    m_oa[Vi]      = 1.0/(1.0+exp(-(V-P[VT_shiftoa].value+P[VT_C_oa_6].value)/(P[VT_C_oa_5].value)));
    exptau_oa[Vi] = exp(-tinc * ((a+b)*P[VT_Tfac].value));

    a             = 1.0/(18.53+exp((V+113.7)/10.95));
    b             = 1.0/(35.56+exp(-(V+1.26)/7.44));
    m_oi[Vi]      = 1.0/(1.0+exp((V+43.1)/5.3));
    exptau_oi[Vi] = exp(-tinc * ((a+b)*P[VT_Tfac].value));

    a = P[VT_ua_a1].value/
      ((exp((V+P[VT_ua_a2].value)/(P[VT_ua_a3].value)))+(exp((V+P[VT_ua_a4].value)/(P[VT_ua_a5].value))));
    b             = .65/(P[VT_ua_b1].value+exp((V+P[VT_ua_b2].value)/(P[VT_ua_b3].value)));
    m_ua[Vi]      = 1.0/(1.0+exp((V+P[VT_ua_m1].value)/(P[VT_ua_m2].value)));
    exptau_ua[Vi] = exp(-tinc * ((a+b)*P[VT_ua_KQ10].value));

    a             = P[VT_ui_a1].value/(P[VT_ui_a2].value+exp((V+P[VT_ui_a3].value)/(P[VT_ui_a4].value)));
    b             = 1.0/exp((V+P[VT_ui_b1].value)/(P[VT_ui_b2].value));
    m_ui[Vi]      = 1.0/(1.0+exp((V+P[VT_ui_m1].value)/(P[VT_ui_m2].value)));
    exptau_ui[Vi] = exp(-tinc * ((a+b)*P[VT_ui_KQ10].value));

    a =
      (fabs(V+P[VT_Xr_a2].value) <
       1e-10 ? 0.0015 : P[VT_Xr_a1].value*(V+P[VT_Xr_a2].value)/(1.0-exp((V+P[VT_Xr_a2].value)/(P[VT_Xr_a3].value))));
    b =
      (fabs(V-P[VT_Xr_b1].value) <
       1e-10 ? 3.7836118e-4 : 7.3898e-5*(V+P[VT_Xr_b1].value)/(exp((V+P[VT_Xr_b1].value)/(P[VT_Xr_b2].value))-1.0));
    m_Xr[Vi]      = 1.0/(1.0+exp((V+P[VT_Xr_m1].value)/(P[VT_Xr_m2].value)));
    exptau_Xr[Vi] = exp(-tinc * (P[VT_Xr_KQ10].value*(a+b)));

    a =
      (fabs(V+P[VT_Xs_a2].value) <
       1e-10 ? 0.00068 : P[VT_Xs_a1].value*(V+P[VT_Xs_a2].value)/(1.0-exp((V+P[VT_Xs_a2].value)/(P[VT_Xs_a3].value))));
    b =
      (fabs(V+P[VT_Xs_b1].value) <
       1e-10 ? 0.000315 : (3.5e-5)*(V+P[VT_Xs_b1].value)/(exp((V+P[VT_Xs_b1].value)/(P[VT_Xs_b2].value))-1.0));
    exptau_Xs[Vi] = exp(-tinc * (P[VT_Xs_KQ10].value*(a+b)));
    m_Xs[Vi]      = 1.0/sqrt(1.0+exp((V+P[VT_Xs_m1].value)/(P[VT_Xs_m2].value)));

    a            = 1.0/(1.0+exp(-(V+10.0)/6.24));
    exptau_d[Vi] =
      exp(-tinc * (fabs(V+10.0) < 1e-10 ? 1.0/(a*4.579) : 1.0/(a*(1.0-exp(-(V+10.0)/6.24))/(0.035*(V+10.0)))));
    m_d[Vi] = 1.0/(1.0+exp(-(V+10.0)/8));
    if (V > 0) {
      exptau_f[Vi] =
        exp(-tinc *
            (1.0/(P[VT_CaL_fT_1].value/(0.0197*exp(-0.0337*0.0337*(V+10.0)*(V+10.0))+0.02)*P[VT_tau_f_mult].value)));
    } else {
      exptau_f[Vi] = exp(-tinc * (1.0/(P[VT_CaL_fT_1].value/(0.0197*exp(-0.0337*0.0337*(V+10.0)*(V+10.0))+0.02))));
    }

    m_f[Vi] = exp(-(V+28.0)/6.9)/(1.0+exp(-(V+28.0)/6.9));

    exptau_w[Vi]  =
      exp(-tinc *
          (fabs(V-7.9) <
           1e-10 ? 1.0/(6.0*0.2/1.3) : 1.0/(6.0*(1.0-exp(-(V-7.9)/5.0))/(1.0+0.3*exp(-(V-7.9)/5.0))/(V-7.9))));
    m_w[Vi] = 1.0-1.0/(1.0+exp(-(V-40.0)/17.0));

    double VFdRT = V*P[VT_FdRT].value;
    expVm[Vi] = exp(-VFdRT)*P[VT_Na_o].value*P[VT_Na_o].value*P[VT_Na_o].value;
    g_Kur[Vi] = P[VT_g_Kur1].value+P[VT_g_Kur2].value/(1.0+exp((V+P[VT_g_Kur3].value)/(P[VT_g_Kur4].value)));
    CKr[Vi]   = P[VT_g_Kr].value/(1.0+exp((V+P[VT_g_Kr1].value)/(P[VT_g_Kr2].value)));
    CK1[Vi]   = P[VT_g_K1].value/(1.0+exp((V+P[VT_g_K1_1].value)*P[VT_g_K1_2].value));
    CNaK[Vi]  = P[VT_I_NaKmax].value/(1.0+.1245*exp(-.1*VFdRT)+.0365*sigma*exp(-VFdRT))*P[VT_K_o].value/
      (P[VT_K_o].value+P[VT_k_mKo].value);
    CNaCa[Vi] = P[VT_k_NaCa].value*exp(P[VT_gamm].value*VFdRT)/
      ((P[VT_k_mNa].value*P[VT_k_mNa].value*P[VT_k_mNa].value+P[VT_Na_o].value*P[VT_Na_o].value*P[VT_Na_o].value)*
       (P[VT_k_mCa].value+P[VT_Ca_o].value)*(1.0+P[VT_k_sat].value*exp((P[VT_gamm].value-1.0)*VFdRT)));
  }
}  // CourtemancheParameters::InitTable

void CourtemancheParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"CourtemancheParameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}
