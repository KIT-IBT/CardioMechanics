/*
 * File: TenTusscher2Parameters.cpp
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


#include <TenTusscher2Parameters.h>

TenTusscher2Parameters::TenTusscher2Parameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

TenTusscher2Parameters::~TenTusscher2Parameters() {}

void TenTusscher2Parameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"TenTusscher2Parameters:"<<endl;

  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}

void TenTusscher2Parameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "#Loading the TenTusscherEtAl2 parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  // Initialization of the Parameters ...
  P[VT_R].name                  =   "R";
  P[VT_Tx].name                 =  "Tx";
  P[VT_F].name                  =   "F";
  P[VT_m_init].name             =      "init_m";
  P[VT_h_init].name             =      "init_h";
  P[VT_j_init].name             =      "init_j";
  P[VT_xr1_init].name           =    "init_Xr1";
  P[VT_xr2_init].name           =    "init_Xr2";
  P[VT_xs_init].name            =     "init_Xs";
  P[VT_r_init].name             =      "init_r";
  P[VT_s_init].name             =      "init_s";
  P[VT_d_init].name             =      "init_d";
  P[VT_f_init].name             =      "init_f";
  P[VT_f2_init].name            =     "init_f2";
  P[VT_fCa_init].name           =    "init_fCa";
  P[VT_Rq_init].name            =     "init_Rq";
  P[VT_O_init].name             =      "init_O";
  P[VT_K_o].name                = "K_o";
  P[VT_Ca_o].name               =        "Ca_o";
  P[VT_Na_o].name               =        "Na_o";
  P[VT_Vc].name                 =  "Vc";
  P[VT_Vsr].name                = "Vsr";
  P[VT_Vss].name                = "Vss";
  P[VT_Vrel].name               =        "Vrel";
  P[VT_Vcell].name              =       "Vcell";
  P[VT_ks1].name                = "ks1";
  P[VT_ks2].name                = "ks2";
  P[VT_k3].name                 =  "k3";
  P[VT_k4].name                 =  "k4";
  P[VT_EC].name                 =  "EC";
  P[VT_max_sr].name             =      "max_sr";
  P[VT_min_sr].name             =      "min_sr";
  P[VT_Vxfer].name              =       "Vxfer";
  P[VT_Bufc].name               =        "Bufc";
  P[VT_Kbufc].name              =       "Kbufc";
  P[VT_Bufsr].name              =       "Bufsr";
  P[VT_Kbufsr].name             =      "Kbufsr";
  P[VT_Bufss].name              =       "Bufss";
  P[VT_Kbufss].name             =      "Kbufss";
  P[VT_taug].name               =        "taug";
  P[VT_Vmaxup].name             =      "Vmaxup";
  P[VT_Kup].name                = "Kup";
  P[VT_C].name                  =   "C";
  P[VT_g_Kr].name               =        "g_Kr";
  P[VT_pKNa].name               =        "pKNa";
  P[VT_g_Ks].name               =        "g_Ks";
  P[VT_g_K1].name               =        "g_K1";
  P[VT_g_to].name               =        "g_to";
  P[VT_g_Na].name               =        "g_Na";
  P[VT_g_bNa].name              =       "g_bNa";
  P[VT_KmK].name                = "KmK";
  P[VT_KmNa].name               =        "KmNa";
  P[VT_knak].name               =        "knak";
  P[VT_g_CaL].name              =       "g_CaL";
  P[VT_g_bCa].name              =       "g_bCa";
  P[VT_kNaCa].name              =       "kNaCa";
  P[VT_KmNai].name              =       "KmNai";
  P[VT_KmCa].name               =        "KmCa";
  P[VT_ksat].name               =        "ksat";
  P[VT_n].name                  =   "n";
  P[VT_g_pCa].name              =       "g_pCa";
  P[VT_KpCa].name               =        "KpCa";
  P[VT_g_pK].name               =        "g_pK";
  P[VT_V_init].name             =      "init_V_m";
  P[VT_Cai_init].name           =    "init_Cai";
  P[VT_CaSR_init].name          =   "init_CaSR";
  P[VT_CaSS_init].name          =   "init_CaSS";
  P[VT_Nai_init].name           =    "init_Nai";
  P[VT_Ki_init].name            =     "init_Ki";
  P[VT_lenght].name             =      "lenght";
  P[VT_radius].name             =      "radius";
  P[VT_pi].name                 =  "pi";
  P[VT_CaiMin].name             =      "CaiMin";
  P[VT_CaiMax].name             =      "CaiMax";
  P[VT_s_inf_vHalf].name        = "s_inf_vHalf";
  P[VT_tau_s_f1].name           =    "tau_s_f1";
  P[VT_tau_s_slope1].name       =        "tau_s_slope1";
  P[VT_tau_s_vHalf1].name       =        "tau_s_vHalf1";
  P[VT_tau_s_enable].name       =        "tau_s_enable";
  P[VT_tau_s_add].name          =   "tau_s_add";
  P[VT_m_Xr1_1].name            =     "m_Xr1_1";
  P[VT_m_Xr1_2].name            =     "m_Xr1_2";
  P[VT_a_Xr1_1].name            =     "a_Xr1_1";
  P[VT_a_Xr1_2].name            =     "a_Xr1_2";
  P[VT_b_Xr1_1].name            =     "b_Xr1_1";
  P[VT_b_Xr1_2].name            =     "b_Xr1_2";
  P[VT_K_Q10Xr1].name           =    "K_Q10Xr1";
  P[VT_m_Xr2_1].name            =     "m_Xr2_1";
  P[VT_m_Xr2_2].name            =     "m_Xr2_2";
  P[VT_a_Xr2_1].name            =     "a_Xr2_1";
  P[VT_a_Xr2_2].name            =     "a_Xr2_2";
  P[VT_b_Xr2_1].name            =     "b_Xr2_1";
  P[VT_b_Xr2_2].name            =     "b_Xr2_2";
  P[VT_K_Q10Xr2].name           =    "K_Q10Xr2";
  P[VT_m_Xs_1].name             =      "m_Xs_1";
  P[VT_m_Xs_2].name             =      "m_Xs_2";
  P[VT_a_Xs_1].name             =      "a_Xs_1";
  P[VT_a_Xs_2].name             =      "a_Xs_2";
  P[VT_b_Xs_1].name             =      "b_Xs_1";
  P[VT_b_Xs_2].name             =      "b_Xs_2";
  P[VT_tau_Xs_add].name         =  "tau_Xs_add";
  P[VT_K_Q10Xs].name            =     "K_Q10Xs";
  P[VT_vol].name                = "vol";
  P[VT_vi].name                 =  "vi";
  P[VT_inverseviF].name         =  "inverseviF";
  P[VT_inverseviF2].name        = "inverseviF2";
  P[VT_inversevssF2].name       =        "inversevssF2";
  P[VT_volforCall].name         =  "volforCall";
  P[VT_INVERSECAPACITANCE].name =  "INVERSECAPACITANCE";
  P[VT_VcdVsr].name             =      "VcdVsr";
  P[VT_Kupsquare].name          =   "Kupsquare";
  P[VT_BufcPKbufc].name         =  "BufcPKbufc";
  P[VT_Kbufcsquare].name        = "Kbufcsquare";
  P[VT_Kbufc2].name             =      "Kbufc2";
  P[VT_BufsrPKbufsr].name       =        "BufsrPKbufsr";
  P[VT_Kbufsrsquare].name       =        "Kbufsrsquare";
  P[VT_Kbufsr2].name            =     "Kbufsr2";
  P[VT_KmNai3].name             =      "KmNai3";
  P[VT_Nao3].name               =        "Nao3";
  P[VT_StepCai].name            =     "StepCai";
  P[VT_inverseRTONF].name       =        "inverseRTONF";
  P[VT_RTONF].name              =       "RTONF";
  P[VT_stim_duration].name      =       "stim_duration";
  P[VT_Amp].name                =    "Amp";
  P[VT_alpha3_div].name         =  "alpha3_div";
  P[VT_alpha5_div].name         =  "alpha5_div";
  P[VT_beta5_div].name          =   "beta5_div";
  P[VT_alpha6].name             =      "alpha6";
  P[VT_beta6].name              =       "beta6";
  P[VT_initMINALC3].name        = "initMINALC3";
  P[VT_initMINALC2].name        = "initMINALC2";
  P[VT_initMINALC1].name        = "initMINALC1";
  P[VT_initMINALO].name         =  "initMINALO";
  P[VT_initMINAUC3].name        = "initMINAUC3";
  P[VT_initMINAUC2].name        = "initMINAUC2";
  P[VT_initMINAUC1].name        = "initMINAUC1";
  P[VT_initMINAUO].name         =  "initMINAUO";
  P[VT_initMINAUIC3].name       =        "initMINAUIC3";
  P[VT_initMINAUIC2].name       =        "initMINAUIC2";
  P[VT_initMINAUIF].name        = "initMINAUIF";
  P[VT_initMINAUIM1].name       =        "initMINAUIM1";
  P[VT_initMINAUIM2].name       =        "initMINAUIM2";
#ifdef ACTIVATE_ISAC_CHANNEL
#ifdef SAC_SACHS
    P[VT_g_SAC].name              =  "g_SAC";
    P[VT_alphaSAC].name           =  "alphaSAC";
    P[VT_ESAC].name               =  "ESAC";
    P[VT_KSAC].name               =  "KSAC";
#endif //SAC_SACHS
#ifdef SAC_KUIJPERS
    P[VT_PsacNa].name       =      "PsacNa";
    P[VT_PsacK].name        =       "PsacK";
    P[VT_PsacCa].name       =      "PsacCa";
    P[VT_Gsac].name         =        "Gsac";
    P[VT_Ksac].name         =        "Ksac";
    P[VT_alphasac].name     =    "alphasac";
    P[VT_lambdasac].name    =   "lambdasac";
    P[VT_lambdasac].readFromFile    = false; 
    P[VT_gsac].name         =        "gsac";
    P[VT_gsac].readFromFile = false;
#endif  // ifdef SAC_KUIJPERS
#endif  // ifdef ACTIVATE_ISAC_CHANNEL
#ifdef ACTIVATE_IKATP_CHANNEL
  P[VT_nicholsarea].name = "nicholsarea";
  P[VT_gkatp].name       = "gkatp";
  P[VT_Mgi].name         = "Mgi";
  P[VT_atpi].name        = "atpi";
  P[VT_adpi].name        = "adpi";
  P[VT_Km_factor].name   = "Km_factor";
# ifdef ISCHEMIA
  P[VT_IschemiaStart].name       = "IschemiaStart";
  P[VT_IschemiaStage1].name      = "IschemiaStage1";
  P[VT_IschemiaStage2].name      = "IschemiaStage2";
  P[VT_ZoneFactor].name          = "ZoneFactor";
  P[VT_DiffusionFactor].name     = "DiffusionFactor";
  P[VT_RestoreIschemia].name     = "RestoreIschemia";
  P[VT_Ko_ZoneFactor_Begin].name = "Ko_ZoneFactorVariance_Begin";
  P[VT_Ko_ZoneFactor_End].name   = "Ko_ZoneFactorVariance_End";

  // P[VT_dVmNa_ZoneFactor_Begin].name = "dVmNa_ZoneFactorVariance_Begin";
  // P[VT_dVmNa_ZoneFactor_End].name = "dVmNa_ZoneFactorVariance_End";
  P[VT_fpH_ZoneFactor_Begin].name = "fpH_ZoneFactorVariance_Begin";
  P[VT_fpH_ZoneFactor_End].name   = "fpH_ZoneFactorVariance_End";
  P[VT_pO_ZoneFactor_Begin].name  = "pO_ZoneFactorVariance_Begin";
  P[VT_pO_ZoneFactor_End].name    = "pO_ZoneFactorVariance_End";
  P[VT_K_o_stage1].name           = "K_o_stage1";
  P[VT_K_o_stage2].name           = "K_o_stage2";
  P[VT_gCaL_stage1].name          = "gCaL_stage1";
  P[VT_gCaL_stage2].name          = "gCaL_stage2";
  P[VT_gNa_stage1].name           = "gNa_stage1";
  P[VT_gNa_stage2].name           = "gNa_stage2";
  P[VT_Mgi_stage1].name           = "Mgi_stage1";
  P[VT_Mgi_stage2].name           = "Mgi_stage2";
  P[VT_ATP_stage1].name           = "ATP_stage1";
  P[VT_ATP_stage2].name           = "ATP_stage2";
  P[VT_ADP_stage1].name           = "ADP_stage1";
  P[VT_ADP_stage2].name           = "ADP_stage2";
  P[VT_ADP_stage2].name           = "ADP_stage2";
  P[VT_dVmNa_stage0].name         = "dVmNa_stage0";
  P[VT_dVmNa_stage1].name         = "dVmNa_stage1";
  P[VT_dVmNa_stage2].name         = "dVmNa_stage2";
  P[VT_knak_1b].name              = "knak_1b";
  P[VT_kNaCa_1b].name             = "kNaCa_1b";
  P[VT_Vmaxup_1b].name            = "Vmaxup_1b";
  P[VT_Vrel_1b].name              = "Vrel_1b";
# endif  // ifdef ISCHEMIA
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

  P[VT_RTONF].readFromFile              = false;
  P[VT_inverseRTONF].readFromFile       = false;
  P[VT_vol].readFromFile                = false;
  P[VT_vi].readFromFile                 = false;
  P[VT_vi].readFromFile                 = false;
  P[VT_inverseviF].readFromFile         = false;
  P[VT_inverseviF2].readFromFile        = false;
  P[VT_inversevssF2].readFromFile       = false;
  P[VT_volforCall].readFromFile         = false;
  P[VT_INVERSECAPACITANCE].readFromFile = false;
  P[VT_VcdVsr].readFromFile             = false;
  P[VT_Kupsquare].readFromFile          = false;
  P[VT_BufcPKbufc].readFromFile         = false;
  P[VT_Kbufcsquare].readFromFile        = false;
  P[VT_Kbufc2].readFromFile             = false;
  P[VT_BufsrPKbufsr].readFromFile       = false;
  P[VT_Kbufsrsquare].readFromFile       = false;
  P[VT_Kbufsr2].readFromFile            = false;
  P[VT_KmNai3].readFromFile             = false;
  P[VT_Nao3].readFromFile               = false;
  P[VT_StepCai].readFromFile            = false;
  P[VT_pi].readFromFile                 = false;

#ifdef ACTIVATE_IKATP_CHANNEL

  // Ischemia
  P[VT_f_T].readFromFile        = false;
  P[VT_gammaconst].readFromFile = false;
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

  // #ifdef ISCHEMIA
  //    P[VT_Ko_ZoneFactor_Begin].readFromFile=false;
  //    P[VT_Ko_ZoneFactor_End].readFromFile=false;
  // P[VT_dVmNa_ZoneFactor_Begin].readFromFile=false;
  // P[VT_dVmNa_ZoneFactor_End].readFromFile=false;
  //    P[VT_fpH_ZoneFactor_Begin].readFromFile=false;
  //    P[VT_fpH_ZoneFactor_End].readFromFile=false;
  //    P[VT_pO_ZoneFactor_Begin].readFromFile=false;
  //    P[VT_pO_ZoneFactor_End].readFromFile=false;
  //    P[VT_K_o_stage1].readFromFile=false;
  //    P[VT_K_o_stage2].readFromFile=false;
  //    P[VT_gCaL_stage1].readFromFile=false;
  //    P[VT_gCaL_stage2].readFromFile=false;
  //    P[VT_gNa_stage1].readFromFile=false;
  //    P[VT_gNa_stage2].readFromFile=false;
  //    P[VT_Mgi_stage1].readFromFile=false;
  //    P[VT_Mgi_stage2].readFromFile=false;
  //    P[VT_ATP_stage1].readFromFile=false;
  //    P[VT_ATP_stage2].readFromFile=false;
  //    P[VT_ADP_stage1].readFromFile=false;
  //    P[VT_ADP_stage2].readFromFile=false;
  //    P[VT_dVmNa_stage0].readFromFile=false;
  //    P[VT_dVmNa_stage1].readFromFile=false;
  //    P[VT_dVmNa_stage2].readFromFile=false;
  // #endif
  ParameterLoader EPL(initFile, EMT_TenTusscher2);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  // End Initialization of the Parameters ...

  Calculate();
  InitTable(tinc);
}  // TenTusscher2Parameters::Init

void TenTusscher2Parameters::Calculate() {
#if KADEBUG
  cerr << "#TenTusscherEtAl2Parameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_pi].value                 = M_PI;
  P[VT_RTONF].value              = (P[VT_R].value)*(P[VT_Tx].value)/(P[VT_F].value);
  P[VT_inverseRTONF].value       = 1/(P[VT_RTONF].value);
  P[VT_vol].value                = (P[VT_pi].value)*(P[VT_lenght].value)*(P[VT_radius].value)*(P[VT_radius].value);
  P[VT_vi].value                 = 0.49*(P[VT_vol].value);
  P[VT_vi].value                 = 0.016404;
  P[VT_inverseviF].value         = 1/((P[VT_vi].value)*(P[VT_F].value));
  P[VT_inverseviF2].value        = 1./(2*(P[VT_vi].value)*(P[VT_F].value));
  P[VT_inversevssF2].value       = 1/(2*((P[VT_vi].value)/300.)*(P[VT_F].value));
  P[VT_volforCall].value         = P[VT_vol].value/1000000000;
  P[VT_INVERSECAPACITANCE].value = 1./(P[VT_C].value);
  P[VT_VcdVsr].value             = (P[VT_Vc].value)/(P[VT_Vsr].value);
  P[VT_Kupsquare].value          = (P[VT_Kup].value)*(P[VT_Kup].value);
  P[VT_BufcPKbufc].value         = (P[VT_Bufc].value)+(P[VT_Kbufc].value);
  P[VT_Kbufcsquare].value        = (P[VT_Kbufc].value)*(P[VT_Kbufc].value);
  P[VT_Kbufc2].value             = 2*(P[VT_Kbufc].value);
  P[VT_BufsrPKbufsr].value       = (P[VT_Bufsr].value)+(P[VT_Kbufsr].value);
  P[VT_Kbufsrsquare].value       = (P[VT_Kbufsr].value)*(P[VT_Kbufsr].value);
  P[VT_Kbufsr2].value            = 2*(P[VT_Kbufsr].value);
  P[VT_KmNai3].value             = (P[VT_KmNai].value)*(P[VT_KmNai].value)*(P[VT_KmNai].value);
  P[VT_Nao3].value               = (P[VT_Na_o].value)*(P[VT_Na_o].value)*(P[VT_Na_o].value);
  P[VT_StepCai].value            = (P[VT_CaiMax].value-P[VT_CaiMin].value)/RTDT;

#ifdef ACTIVATE_IKATP_CHANNEL

  // temperature effect of I_Katp channel
  /// (E16)
  P[VT_f_T].value = pow(1.3, ( (P[VT_Tx].value - 273.15) - 36) / 10); // All simulations correspond to a
                                                                      // temperature of 37Â¡C
  P[VT_gammaconst].value = (P[VT_gkatp].value / P[VT_nicholsarea].value) * P[VT_INVERSECAPACITANCE].value;
#endif  // ifdef ACTIVATE_IKATP_CHANNEL
//#ifdef SAC_KUIJPERS
   // P[VT_gsac].value = P[VT_Gsac].value /
   // (1.0 + P[VT_Ksac].value * exp(-P[VT_alphasac].value * (stretch - 1.0) ) );
// #endif  // ifdef SAC_KUIJPERS
  // #ifdef ISCHEMIA
  //    // the following values should be transferred to the *.ev file - if necessary ...
  //    P[VT_Ko_ZoneFactor_Begin].value=0;
  //    P[VT_Ko_ZoneFactor_End].value=1;
  // P[VT_dVmNa_ZoneFactor_Begin].value=0.5;
  // P[VT_dVmNa_ZoneFactor_End].value=1;
  //    P[VT_fpH_ZoneFactor_Begin].value=0.5;
  //    P[VT_fpH_ZoneFactor_End].value=1;
  //    P[VT_pO_ZoneFactor_Begin].value=0;
  //    P[VT_pO_ZoneFactor_End].value=0.1;
  //    P[VT_K_o_stage1].value=8.7;
  //    P[VT_K_o_stage2].value=12.5;
  //    P[VT_gCaL_stage1].value=P[VT_g_CaL].value * 0.875;
  //    P[VT_gCaL_stage2].value=P[VT_g_CaL].value * 0.75;
  //    P[VT_gNa_stage1].value=P[VT_g_Na].value * 0.875;
  //    P[VT_gNa_stage2].value=P[VT_g_Na].value * 0.75;
  //    P[VT_Mgi_stage1].value=3;
  //    P[VT_Mgi_stage2].value=6;
  //    P[VT_ATP_stage1].value=5.7;
  //    P[VT_ATP_stage2].value=4.6;
  //    P[VT_ADP_stage1].value=57;
  //    P[VT_ADP_stage2].value=99;

  //    P[VT_dVmNa_stage0].value=0;
  //    P[VT_dVmNa_stage1].value=1.7;
  //    P[VT_dVmNa_stage2].value=3.4;
  //                    //see shaw97a p 270 - added by dw
  // #endif

  // if (PrintParameterMode == PrintParameterModeOn) PrintParameters();
}  // TenTusscher2Parameters::Calculate

void TenTusscher2Parameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#TenTusscherEtAl2Parameters - InitTable ..." << endl;
#endif  // if KADEBUG
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi                     = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    const ML_CalcType rec_iNaK =
      (1./(1.+0.1245*exp(-0.1*V*(P[VT_inverseRTONF].value))+0.0353*exp(-V*(P[VT_inverseRTONF].value))));
    NaK_P1[Vi]  = rec_iNaK;
    rec_ipK[Vi] = P[VT_g_pK].value/(1.+exp((25.-V)/5.98));
    const double a_m        = 1./(1.+exp((-60.-V)/5.));
    const double b_m        = 0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    const ML_CalcType tau_m = a_m*b_m;
    m_inf[Vi] = 1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
    ML_CalcType tau_h, tau_j;
    if (V >= -40.) {
      const ML_CalcType AH_1 = 0.;
      const ML_CalcType BH_1 = (0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      tau_h = 1.0/(AH_1+BH_1);
    } else {
      const ML_CalcType AH_2 = (0.057*exp(-(V+80.)/6.8));
      const ML_CalcType BH_2 = (2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
      tau_h = 1.0/(AH_2+BH_2);
    }
    h_inf[Vi] = 1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
    if (V >= -40.) {
      const ML_CalcType AJ_1 = 0.;
      const ML_CalcType BJ_1 = (0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
      tau_j = 1.0/(AJ_1+BJ_1);
    } else {
      const ML_CalcType AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                                 exp(-0.04391*V))*(V+37.78)/
                                (1.+exp(0.311*(V+79.23))));
      const ML_CalcType BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
      tau_j = 1.0/(AJ_2+BJ_2);
    }
    j_inf[Vi]   = h_inf[Vi];
    Xr1_inf[Vi] = 1./(1.+exp((P[VT_m_Xr1_1].value-V)/(P[VT_m_Xr1_2].value)));
    const ML_CalcType a_Xr1   = 450./(1.+exp((P[VT_a_Xr1_1].value-V)/(P[VT_a_Xr1_2].value)));
    const ML_CalcType b_Xr1   = 6./(1.+exp((V-(P[VT_b_Xr1_1].value))/(P[VT_b_Xr1_2].value)));
    const ML_CalcType tau_Xr1 = a_Xr1*b_Xr1*fabs(P[VT_K_Q10Xr1].value);
    Xr2_inf[Vi] = 1./(1.+exp((V-(P[VT_m_Xr2_1].value))/(P[VT_m_Xr2_2].value)));
    const ML_CalcType a_Xr2   = 3./(1.+exp((P[VT_a_Xr2_1].value-V)/(P[VT_a_Xr2_2].value)));
    const ML_CalcType b_Xr2   = 1.12/(1.+exp((V-(P[VT_b_Xr2_1].value))/(P[VT_b_Xr2_2].value)));
    const ML_CalcType tau_Xr2 = a_Xr2*b_Xr2*fabs(P[VT_K_Q10Xr2].value);
    Xs_inf[Vi] = 1./(1.+exp((P[VT_m_Xs_1].value-V)/P[VT_m_Xs_2].value));
    const ML_CalcType a_Xs   = (1400./(sqrt(1.+exp((P[VT_a_Xs_1].value-V)/P[VT_a_Xs_2].value))));
    const ML_CalcType b_Xs   = (1./(1.+exp((V-P[VT_b_Xs_1].value)/P[VT_b_Xs_2].value)));
    const ML_CalcType tau_Xs = (a_Xs*b_Xs+P[VT_tau_Xs_add].value)*fabs(P[VT_K_Q10Xs].value);
    r_inf[Vi] = 1./(1.+exp((20-V)/6.));                            // 56
    s_inf[Vi] = 1./(1.+exp((V+(P[VT_s_inf_vHalf].value))/5.));     // =1./(1.+exp((V+20)/5.));  //58
    const ML_CalcType tau_r = 9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;  // 57
    const ML_CalcType tau_s = P[VT_tau_s_f1].value*
      exp(-(V+(P[VT_tau_s_vHalf1].value))*(V+(P[VT_tau_s_vHalf1].value))/P[VT_tau_s_slope1].value)+
      P[VT_tau_s_enable].value*5./(1.+exp((V-20.)/5.))+P[VT_tau_s_add].value;
    d_inf[Vi] = 1./(1.+exp((-8.-V)/7.5));  // 41 n7
    const ML_CalcType a_d   = 1.4/(1.+exp((-35-V)/13))+0.25;
    const ML_CalcType b_d   = 1.4/(1.+exp((V+5)/5));
    const ML_CalcType c_d   = 1./(1.+exp((50-V)/20));
    const ML_CalcType tau_d = a_d*b_d+c_d;
    f_inf[Vi] = 1./(1.+exp((V+20)/7));  // 46 n12
    const ML_CalcType V27square = (V+27)*(V+27);                   // new
    const ML_CalcType exsquare  = V27square/(15*15);               // new
    const ML_CalcType a_f       = 1102.5*exp(-exsquare);           // new n13
    const ML_CalcType b_f       = 200./(1.+exp((13-V)/10.));       // new n14
    const ML_CalcType g_f       = (180./(1.+exp((V+30)/10.)))+20.; // new n15
    const ML_CalcType tau_f     = a_f+b_f+g_f;                     // new n16
    f2_inf[Vi] = (.67/(1.+exp((V+35.)/7.)))+.33;                   // new n17
    const ML_CalcType a_f2   = 600.*exp(-(V+25)*(V+25)/170.);      // new n18
    const ML_CalcType b_f2   = 31./(1.+exp((25-V)/10.));           // new n19
    const ML_CalcType c_f2   = 16./(1.+exp((V+30)/10.));           // new n20
    const ML_CalcType tau_f2 = a_f2+b_f2+c_f2;                     // new n21
    NaCa_P1[Vi] =
      ((1./((P[VT_KmNai3].value)+(P[VT_Nao3].value)))*(1./((P[VT_KmCa].value)+(P[VT_Ca_o].value)))*
       (1./(1.+(P[VT_ksat].value)*exp(((P[VT_n].value)-1.)*V*(P[VT_inverseRTONF].value)))))*
      (exp((P[VT_n].value)*V*(P[VT_inverseRTONF].value))*(P[VT_Ca_o].value));  // (exp(n*vv/RTONF)*Nai*Nai*NaiP[VT_C].valueao);
    NaCa_P2[Vi] =
      ((1./((P[VT_KmNai3].value)+(P[VT_Nao3].value)))*(1./((P[VT_KmCa].value)+(P[VT_Ca_o].value)))*
       (1./(1.+(P[VT_ksat].value)*exp(((P[VT_n].value)-1.)*V*(P[VT_inverseRTONF].value)))))*
      (exp(((P[VT_n].value)-1)*V/(P[VT_RTONF].value))*(P[VT_Na_o].value)*(P[VT_Na_o].value)*(P[VT_Na_o].value)*2.5);  //
                                                                                                                      //
                                                                                                                      // *7.5);//*2.5);
    CaL_P1[Vi] = 4*(V-15.)*(P[VT_F].value)/(P[VT_RTONF].value)*(0.25*exp(2*(V-15.)/(P[VT_RTONF].value)))/
      (exp(2*(V-15.)/(P[VT_RTONF].value))-1.);
    CaL_P2[Vi] = 4*(V-15.)*((P[VT_F].value)/(P[VT_RTONF].value))*(1*(P[VT_Ca_o].value))/
      (exp(2*(V-15.)/(P[VT_RTONF].value))-1.);

    const ML_CalcType HT = tinc*1000;
    exptau_m[Vi]   =  exp(-HT/ tau_m);
    exptau_h[Vi]   =  exp(-HT/ tau_h);
    exptau_j[Vi]   =  exp(-HT/ tau_j);
    exptau_Xr1[Vi] =  exp(-HT/ tau_Xr1);
    exptau_Xr2[Vi] =  exp(-HT/ tau_Xr2);
    exptau_Xs[Vi]  =  exp(-HT/ tau_Xs);
    exptau_s[Vi]   =  exp(-HT/ tau_s);
    exptau_r[Vi]   =  exp(-HT/ tau_r);
    exptau_d[Vi]   =  exp(-HT/ tau_d);
    exptau_f[Vi]   =  exp(-HT/ tau_f);
    exptau_f2[Vi]  =  exp(-HT/ tau_f2);
    const double Cai = Vi * (P[VT_StepCai].value);
    ECA[Vi] = ((P[VT_RTONF].value)*0.5)*(log(((P[VT_Ca_o].value) / Cai)));

#ifdef ACTIVATE_IKATP_CHANNEL

    /* Stuff used by the ATP dependent potassium channel I_Katp */

    // intracellular Na+ ions
    /// (E15)
    KhNa[Vi] = 25.9 * exp(-0.35 * P[VT_inverseRTONF].value * V);
#endif  // ifdef ACTIVATE_IKATP_CHANNEL

#ifdef MARKOV_I_NA
    a11[Vi] = 3.802/(0.1027*exp(-V/17.0)+0.2*exp(-V/150.0));
    a12[Vi] = 3.802/(0.1027*exp(-V/15.0)+0.23*exp(-V/150.0));
    a13[Vi] = 3.802/(0.1027*exp(-V/12.0)+0.25*exp(-V/150.0));

    b11[Vi] = 0.1917*exp(-V/20.3);
    b12[Vi] = 0.2*exp(-(V-5.0)/20.3);
    b13[Vi] = 0.22*exp(-(V-10.0)/20.3);

    a3[Vi] = 3.7933e-07*exp(-V/7.7)/P[VT_alpha3_div].value;
    b3[Vi] = 0.0084+0.00002*V;

    a2[Vi] = 9.178*exp(V/29.68);
    b2[Vi] = (a13[Vi]*a2[Vi]*a3[Vi])/(b13[Vi]*b3[Vi]);

    a4[Vi] = a2[Vi]/100.0;
    b4[Vi] = a3[Vi];

    a5[Vi] = a2[Vi]/P[VT_alpha5_div].value;
    b5[Vi] = a3[Vi]/P[VT_beta5_div].value;
#endif  // ifdef MARKOV_I_NA
//#ifdef SAC_KUIJPERS
      //double VinverseRTONF = V*P[VT_inverseRTONF].value;
      //Csac1[Vi] = P[VT_gsac].value*VinverseRTONF*P[VT_F].value/(1.0 - exp(-VinverseRTONF) );
      //Csac2[Vi] = P[VT_gsac].value*4.0*VinverseRTONF*P[VT_F].value/(1.0 - exp(-2.0*VinverseRTONF) );
//#endif  // ifdef SAC_KUIJPERS
      
    
  }
}  // TenTusscher2Parameters::InitTable
