/*
 * File: MaleckarEtAlParameters.cpp
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


#include <MaleckarEtAlParameters.h>
MaleckarEtAlParameters::MaleckarEtAlParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

MaleckarEtAlParameters::~MaleckarEtAlParameters() {}

void MaleckarEtAlParameters::PrintParameters() {
  cout << "MaleckarEtAlParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void MaleckarEtAlParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the MaleckarEtAl parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  P[VT_R].name               = "R";
  P[VT_T].name               = "T";
  P[VT_F].name               = "F";
  P[VT_Cm].name              = "Cm";
  P[VT_stim_duration].name   = "stim_duration";
  P[VT_Amp].name             = "Amp";
  P[VT_P_Na].name            = "P_Na";
  P[VT_g_Ca_L].name          = "g_Ca_L";
  P[VT_E_Ca_app].name        = "E_Ca_app";
  P[VT_k_Ca].name            = "k_Ca";
  P[VT_g_t].name             = "g_t";
  P[VT_g_kur].name           = "g_kur";
  P[VT_g_K1].name            = "g_K1";
  P[VT_g_Ks].name            = "g_Ks";
  P[VT_g_Kr].name            = "g_Kr";
  P[VT_g_B_Na].name          = "g_B_Na";
  P[VT_g_B_Ca].name          = "g_B_Ca";
  P[VT_K_NaK_K].name         = "K_NaK_K";
  P[VT_i_NaK_max].name       = "i_NaK_max";
  P[VT_pow_K_NaK_Na_15].name = "pow_K_NaK_Na_15";
  P[VT_i_CaP_max].name       = "i_CaP_max";
  P[VT_k_CaP].name           = "k_CaP";
  P[VT_K_NaCa].name          = "K_NaCa";
  P[VT_d_NaCa].name          = "d_NaCa";
  P[VT_gamma_Na].name        = "gamma_Na";
  P[VT_ACh].name             = "ACh";
  P[VT_phi_Na_en].name       = "phi_Na_en";
  P[VT_Vol_i].name           = "Vol_i";
  P[VT_V_corrcell].name      = "V_corrcell";
  P[VT_Vol_d].name           = "Vol_d";
  P[VT_tau_di].name          = "tau_di";
  P[VT_Mg_i].name            = "Mg_i";
  P[VT_Vol_c].name           = "Vol_c";
  P[VT_tau_Na].name          = "tau_Na";
  P[VT_tau_K].name           = "tau_K";
  P[VT_tau_Ca].name          = "tau_Ca";
  P[VT_Na_b].name            = "Na_b";
  P[VT_Ca_b].name            = "Ca_b";
  P[VT_K_b].name             = "K_b";
  P[VT_I_up_max].name        = "I_up_max";
  P[VT_k_cyca].name          = "k_cyca";
  P[VT_k_srca].name          = "k_srca";
  P[VT_k_xcs].name           = "k_xcs";
  P[VT_alpha_rel].name       = "alpha_rel";
  P[VT_Vol_up].name          = "Vol_up";
  P[VT_Vol_rel].name         = "Vol_rel";
  P[VT_r_recov].name         = "r_recov";
  P[VT_tau_tr].name          = "tau_tr";
  P[VT_k_rel_i].name         = "k_rel_i";
  P[VT_k_rel_d].name         = "k_rel_d";
  P[VT_V_init].name          = "V_init";
  P[VT_Na_c_init].name       = "Na_c_init";
  P[VT_Na_i_init].name       = "Na_i_init";
  P[VT_m_init].name          = "m_init";
  P[VT_h1_init].name         = "h1_init";
  P[VT_h2_init].name         = "h2_init";
  P[VT_Ca_d_init].name       = "Ca_d_init";
  P[VT_d_L_init].name        = "d_L_init";
  P[VT_f_L1_init].name       = "f_L1_init";
  P[VT_f_L2_init].name       = "f_L2_init";
  P[VT_K_c_init].name        = "K_c_init";
  P[VT_K_i_init].name        = "K_i_init";
  P[VT_r_init].name          = "r_init";
  P[VT_s_init].name          = "s_init";
  P[VT_a_ur_init].name       = "a_ur_init";
  P[VT_i_ur_init].name       = "i_ur_init";
  P[VT_n_init].name          = "n_init";
  P[VT_pa_init].name         = "pa_init";
  P[VT_Ca_c_init].name       = "Ca_c_init";
  P[VT_Ca_i_init].name       = "Ca_i_init";
  P[VT_O_C_init].name        = "O_C_init";
  P[VT_O_TC_init].name       = "O_TC_init";
  P[VT_O_TMgC_init].name     = "O_TMgC_init";
  P[VT_O_TMgMg_init].name    = "O_TMgMg_init";
  P[VT_O_init].name          = "O_init";
  P[VT_Ca_rel_init].name     = "Ca_rel_init";
  P[VT_Ca_up_init].name      = "Ca_up_init";
  P[VT_O_Calse_init].name    = "O_Calse_init";
  P[VT_F1_init].name         = "F1_init";
  P[VT_F2_init].name         = "F2_init";
  P[VT_dT_yKur1].name        =        "dT_yKur1";
  P[VT_dT_yKur2].name        =        "dT_yKur2";
  P[VT_Init_yKur].name       =       "Init_yKur";

  P[VT_InitTableDone].name =       "InitTableDone";
  P[VT_tincforTab].name    =  "tincforTab";
  P[VT_RTdF].name          =        "RTdF";
  P[VT_FdRT].name          =        "FdRT";
  P[VT_FFdRT].name         =       "FFdRT";
  P[VT_RTd2F].name         =       "RTd2F";
  P[VT_Ci_di].name         =       "Ci_di";
  P[VT_Ci_tr].name         =       "Ci_tr";
  P[VT_dVol_iF].name       =     "dVol_iF";
  P[VT_dVol_cF].name       =     "dVol_cF";
  P[VT_dVol_dF].name       =     "dVol_dF";
  P[VT_dVol_upF].name      =    "dVol_upF";
  P[VT_dVol_relF].name     =   "dVol_relF";
  P[VT_k_xcs2ds].name      =    "k_xcs2ds";

  P[VT_InitTableDone].readFromFile = false;
  P[VT_tincforTab].readFromFile    = false;
  P[VT_RTdF].readFromFile          = false;
  P[VT_FdRT].readFromFile          = false;
  P[VT_FFdRT].readFromFile         = false;
  P[VT_RTd2F].readFromFile         = false;
  P[VT_Ci_di].readFromFile         = false;
  P[VT_Ci_tr].readFromFile         = false;
  P[VT_dVol_iF].readFromFile       = false;
  P[VT_dVol_cF].readFromFile       = false;
  P[VT_dVol_dF].readFromFile       = false;
  P[VT_dVol_upF].readFromFile      = false;
  P[VT_dVol_relF].readFromFile     = false;
  P[VT_k_xcs2ds].readFromFile      = false;

  ParameterLoader EPL(initFile, EMT_MaleckarEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // MaleckarEtAlParameters::Init

void MaleckarEtAlParameters::Calculate() {
#if KADEBUG
  cerr << "#MaleckarEtAlParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_RTdF].value      = P[VT_R].value*(P[VT_T].value)/(P[VT_F].value);
  P[VT_FdRT].value      = 1.0/(P[VT_RTdF].value);
  P[VT_FFdRT].value     = (P[VT_F].value)/(P[VT_RTdF].value);
  P[VT_RTd2F].value     = P[VT_RTdF].value*.5;
  P[VT_Ci_di].value     = 2.00000*P[VT_Vol_d].value*P[VT_F].value/P[VT_tau_di].value;
  P[VT_Ci_tr].value     = 2.00000*P[VT_Vol_rel].value*P[VT_F].value/P[VT_tau_tr].value;
  P[VT_dVol_iF].value   = 1.0/(P[VT_Vol_i].value*P[VT_F].value);
  P[VT_dVol_cF].value   = 1.0/(P[VT_Vol_c].value*P[VT_F].value);
  P[VT_dVol_dF].value   = 1.0/(2.0*P[VT_Vol_d].value*P[VT_F].value);
  P[VT_dVol_upF].value  = 1.0/(2.0*P[VT_Vol_up].value*P[VT_F].value);
  P[VT_dVol_relF].value = 1.0/(2.0*P[VT_Vol_rel].value*P[VT_F].value);
  P[VT_k_xcs2ds].value  = P[VT_k_xcs].value*P[VT_k_xcs].value/P[VT_k_srca].value;
}

void MaleckarEtAlParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#MaleckarEtAlParameters - InitTable ..." << endl;
#endif  // if KADEBUG
  ML_CalcType HT = tinc;
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    exptau_r[Vi]      =  exp(-HT/ (0.00350000*(exp((((-V*V)/30.0000)/30.0000)))+0.00150000));
    r_infinity[Vi]    = 1.00000/(1.00000+(exp(((V - 1.00000)/ -11.0000))));
    a_ur_infinity[Vi] = 1.00000/(1.00000+(exp((-(V+6.00000)/8.60000))));
    exptau_a_ur[Vi]   =  exp(-HT/ (0.00900000/(1.00000+(exp(((V+5.00000)/12.0000))))+0.000500000));
    i_ur_infinity[Vi] = 1.00000/(1.00000+(exp(((V+7.50000)/10.0000))));
    exptau_i_ur[Vi]   =  exp(-HT/ (0.590000/(1.00000+(exp(((V+60.0000)/10.0000))))+3.05000));
    m_infinity[Vi]    = 1.00000/(1.00000+(exp(((V+27.1200)/ -8.21000))));
    ML_CalcType m_factor = (V+25.5700)/28.8000;
    exptau_m[Vi]   =  exp(-HT/ (4.20000e-05*(exp((-m_factor*m_factor)))+2.40000e-05));
    h_infinity[Vi] = 1.00000/(1.00000+(exp(((V+63.6000)/5.30000))));
    ML_CalcType h_factor = 1.00000/(1.00000+(exp(((V+35.1000)/3.20000))));
    exptau_h1[Vi]    =  exp(-HT/ (0.0300000*h_factor+0.000300000));
    exptau_h2[Vi]    =  exp(-HT/ (0.120000*h_factor+0.00300000));
    d_L_infinity[Vi] = 1.00000/(1.00000+(exp(((V+9.00000)/ -5.80000))));
    ML_CalcType d_L_factor = (V+35.0000)/30.0000;
    exptau_d_L[Vi]   =  exp(-HT/ (0.00270000*(exp((-d_L_factor*d_L_factor)))+0.00200000));
    f_L_infinity[Vi] = 1.00000/(1.00000+(exp(((V+27.4000)/7.10000))));
    ML_CalcType f_L_factor = V+40.0000;
    exptau_f_L1[Vi] =  exp(-HT/ (0.161000*(exp((((-f_L_factor*f_L_factor)/14.4000)/14.4000)))+0.0100000));
    exptau_f_L2[Vi] =  exp(-HT/ (1.33230*(exp((((-f_L_factor*f_L_factor)/14.2000)/14.2000)))+0.0626000));
    ML_CalcType s_factor = (V+52.4500)/15.8827;
    s_infinity[Vi] = 1.00000/(1.00000+(exp(((V+40.5000)/11.5000))));
    exptau_s[Vi]   =  exp(-HT/ (0.0256350*(exp((-s_factor*s_factor)))+0.0141400));
    ML_CalcType n_factor = (V - 20.0000)/20.0000;
    n_infinity[Vi] = 1.00000/(1.00000+(exp(((V - 19.9000)/ -12.7000))));
    exptau_n[Vi]   = exp(-HT/ (0.700000+ 0.400000*(exp((-n_factor*n_factor)))));
    ML_CalcType pa_factor = (V+20.1376)/22.1996;
    p_a_infinity[Vi] = 1.00000/(1.00000+(exp(((V+15.0000)/ -6.00000))));
    exptau_pa[Vi]    = exp(-HT/ (0.0311800+ 0.217180*(exp((-pa_factor*pa_factor)))));
    pip[Vi]          = 1.00000/(1.00000+(exp(((V+55.0000)/24.0000))));
    KIKACH[Vi]       =
      (10.0000/
       (1.00000+
        (9.13652*
         (pow((float)1.00000,
              (float)0.477811)))/
        (pow((float)P[VT_ACh].value,
             (float)0.477811))))*(0.0517000+0.451600/(1.00000+(exp(((V+59.5300)/17.1800)))))*P[VT_Cm].value;
    Q_tot[Vi]    =  0.0500000*V;
    CexpINa[Vi]  = V*(P[VT_FFdRT].value)/((exp(((V*(P[VT_FdRT].value))))) - 1.00000);
    CI_NaCa1[Vi] = exp(P[VT_FdRT].value*V*P[VT_gamma_Na].value);
    CI_NaCa2[Vi] = exp( (P[VT_gamma_Na].value - 1.00000)*V*P[VT_FdRT].value);

    // IKur Blockade
    yKur_infinity[Vi] = 0.1 + 0.9/(1.0+exp((V+40.0)/5.0));
    exptau_yKur[Vi]   =
      exp(-HT/(P[VT_dT_yKur1].value + (P[VT_dT_yKur2].value-P[VT_dT_yKur1].value)/(1.0+exp((V+40.0)/5.0))));
  }
}  // MaleckarEtAlParameters::InitTable
