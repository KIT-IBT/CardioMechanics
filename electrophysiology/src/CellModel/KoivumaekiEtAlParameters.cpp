/* File: KoivumaekiEtAlParameters.cpp
        Institute of Biomedical Engineering, Karlsruhe Institute of Technology (KIT) */

#include <KoivumaekiEtAlParameters.h>
KoivumaekiEtAlParameters::KoivumaekiEtAlParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

KoivumaekiEtAlParameters::~KoivumaekiEtAlParameters() {}

void KoivumaekiEtAlParameters::PrintParameters() {
  cout << "KoivumaekiEtAlParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void KoivumaekiEtAlParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the KoivumaekiEtAl parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  P[VT_R].name               = "R";
  P[VT_T].name               = "T";
  P[VT_F].name               = "F";
  P[VT_Cm].name              = "Cm";
  P[VT_stim_duration].name   = "stim_duration";
  P[VT_Amp].name             = "Amp";
  P[VT_Na_o].name            = "Na_o";
  P[VT_Ca_o].name            = "Ca_o";
  P[VT_K_o].name             = "K_o";
  P[VT_V_ss].name            = "V_ss";
  P[VT_V_corrcell].name      = "V_corrcell";
  P[VT_r_junct].name         = "r_junct";
  P[VT_l_cell].name          = "l_cell";
  P[VT_cell_dillation].name  = "cell_dillation";
  P[VT_dx].name              = "dx";
  P[VT_totCMDN].name         = "totCMDN";
  P[VT_totSLlow].name        = "totSLlow";
  P[VT_totSLhigh].name       = "totSLhigh";
  P[VT_totFura2].name        = "totFura2";
  P[VT_totHTRPN].name        = "totHTRPN";
  P[VT_totLTRPN].name        = "totLTRPN";
  P[VT_KdCMDN].name          = "KdCMDN";
  P[VT_KdSLlow].name         = "KdSLlow";
  P[VT_KdSLhigh].name        = "KdSLhigh";
  P[VT_KdFura2].name         = "KdFura2";
  P[VT_kHTRPNon].name        = "kHTRPNon";
  P[VT_kHTRPNoff].name       = "kHTRPNoff";
  P[VT_kLTRPNon].name        = "kLTRPNon";
  P[VT_kLTRPNoff].name       = "kLTRPNoff";
  P[VT_totCSQN].name         = "totCSQN";
  P[VT_KdCSQN].name          = "KdCSQN";
  P[VT_totBNa].name          = "totBNa";
  P[VT_KdBNa].name           = "KdBNa";
  P[VT_PNa].name             = "PNa";
  P[VT_ECa_app].name         = "ECa_app";
  P[VT_gCaL].name            = "gCaL";
  P[VT_kCan].name            = "kCan";
  P[VT_kCa].name             = "kCa";
  P[VT_gt].name              = "gt";
  P[VT_gsus].name            = "gsus";
  P[VT_gKs].name             = "gKs";
  P[VT_gKr].name             = "gKr";
  P[VT_gK1].name             = "gK1";
  P[VT_gNab].name            = "gNab";
  P[VT_gCab].name            = "gCab";
  P[VT_INaKmax].name         = "INaKmax";
  P[VT_kNaKK].name           = "kNaKK";
  P[VT_kNaKNa].name          = "kNaKNa";
  P[VT_ICaPmax].name         = "ICaPmax";
  P[VT_kCaP].name            = "kCaP";
  P[VT_kNaCa].name           = "kNaCa";
  P[VT_gam].name             = "gam";
  P[VT_dNaCa].name           = "dNaCa";
  P[VT_gIf].name             = "gIf";
  P[VT_DCa].name             = "DCa";
  P[VT_DCaSR].name           = "DCaSR";
  P[VT_DCaBm].name           = "DCaBm";
  P[VT_DNa].name             = "DNa";
  P[VT_SERCAKmf].name        = "SERCAKmf";
  P[VT_SERCAKmr].name        = "SERCAKmr";
  P[VT_k4].name              = "k4";
  P[VT_cpumps].name          = "cpumps";
  P[VT_kSRleak].name         = "kSRleak";
  P[VT_tau_fca].name         =     "tau_fca";
  P[VT_RyRtauadapt].name     = "RyRtauadapt";
  P[VT_RyRtauactss].name     = "RyRtauactss";
  P[VT_RyRtauinactss].name   =       "RyRtauinactss";
  P[VT_RyRtauact].name       =   "RyRtauact";
  P[VT_RyRtauinact].name     = "RyRtauinact";
  P[VT_V_init].name          = "V_init";
  P[VT_m_init].name          = "m_init";
  P[VT_h1_init].name         = "h1_init";
  P[VT_h2_init].name         = "h2_init";
  P[VT_d_init].name          = "d_init";
  P[VT_f1_init].name         = "f1_init";
  P[VT_fca_init].name        = "fca_init";
  P[VT_r_init].name          = "r_init";
  P[VT_s_init].name          = "s_init";
  P[VT_sus_r_init].name      = "sus_r_init";
  P[VT_sus_s_init].name      = "sus_s_init";
  P[VT_n_init].name          = "n_init";
  P[VT_pa_init].name         = "pa_init";
  P[VT_y_init].name          = "y_init";
  P[VT_RyR_oss_init].name    = "RyR_oss_init";
  P[VT_RyR_css_init].name    = "RyR_css_init";
  P[VT_RyR_ass_init].name    = "RyR_ass_init";
  P[VT_RyR_o1_init].name     = "RyR_o1_init";
  P[VT_RyR_c1_init].name     = "RyR_c1_init";
  P[VT_RyR_a1_init].name     = "RyR_a1_init";
  P[VT_RyR_o2_init].name     = "RyR_o2_init";
  P[VT_RyR_c2_init].name     = "RyR_c2_init";
  P[VT_RyR_a2_init].name     = "RyR_a2_init";
  P[VT_RyR_o3_init].name     = "RyR_o3_init";
  P[VT_RyR_c3_init].name     = "RyR_c3_init";
  P[VT_RyR_a3_init].name     = "RyR_a3_init";
  P[VT_SERCA_Ca1_init].name  = "SERCA_Ca1_init";
  P[VT_SERCA_Ca2_init].name  = "SERCA_Ca2_init";
  P[VT_SERCA_Ca3_init].name  = "SERCA_Ca3_init";
  P[VT_SERCA_Cass_init].name = "SERCA_Cass_init";
  P[VT_Na_ss_init].name      = "Na_ss_init";
  P[VT_Na_i_init].name       = "Na_i_init";
  P[VT_K_i_init].name        = "K_i_init";
  P[VT_HTRPNCa1_init].name   = "HTRPNCa1_init";
  P[VT_HTRPNCa2_init].name   = "HTRPNCa2_init";
  P[VT_HTRPNCa3_init].name   = "HTRPNCa3_init";
  P[VT_HTRPNCa4_init].name   = "HTRPNCa4_init";
  P[VT_LTRPNCa1_init].name   = "LTRPNCa1_init";
  P[VT_LTRPNCa2_init].name   = "LTRPNCa2_init";
  P[VT_LTRPNCa3_init].name   = "LTRPNCa3_init";
  P[VT_LTRPNCa4_init].name   = "LTRPNCa4_init";
  P[VT_Ca_ss_init].name      = "Ca_ss_init";
  P[VT_Ca_i1_init].name      = "Ca_i1_init";
  P[VT_Ca_i2_init].name      = "Ca_i2_init";
  P[VT_Ca_i3_init].name      = "Ca_i3_init";
  P[VT_Ca_i4_init].name      = "Ca_i4_init";
  P[VT_Ca_SR1_init].name     = "Ca_SR1_init";
  P[VT_Ca_SR2_init].name     = "Ca_SR2_init";
  P[VT_Ca_SR3_init].name     = "Ca_SR3_init";
  P[VT_Ca_SR4_init].name     = "Ca_SR4_init";
  P[VT_RTdF].name            =        "RTdF";
  P[VT_FdRT].name            =        "FdRT";
  P[VT_RTd2F].name           =       "RTd2F";
  P[VT_Dlcell].name          =      "Dlcell";
  P[VT_Ddcell].name          =      "Ddcell";
  P[VT_Dvcell].name          =      "Dvcell";
  P[VT_rstart].name          =      "rstart";
  P[VT_rend].name            =        "rend";
  P[VT_Aj_nj].name           =       "Aj_nj";
  P[VT_xj_nj].name           =       "xj_nj";
  P[VT_xj_nj_Nai].name       =   "xj_nj_Nai";
  P[VT_V_nonjunct1].name     = "V_nonjunct1";
  P[VT_V_nonjunct2].name     = "V_nonjunct2";
  P[VT_V_nonjunct3].name     = "V_nonjunct3";
  P[VT_V_nonjunct4].name     = "V_nonjunct4";
  P[VT_V_nonjunct_Nai].name  =      "V_nonjunct_Nai";
  P[VT_V_cytosol].name       =   "V_cytosol";
  P[VT_V_SR1].name           =       "V_SR1";
  P[VT_V_SR2].name           =       "V_SR2";
  P[VT_V_SR3].name           =       "V_SR3";
  P[VT_V_SR4].name           =       "V_SR4";
  P[VT_exp_tau_fca].name     = "exp_tau_fca";
  P[VT_pow_kNaKNa].name      =  "pow_kNaKNa";
  P[VT_k1].name              = "k1";
  P[VT_k2].name              = "k2";
  P[VT_k3].name              = "k3";
  P[VT_powKo].name           = "powKo";
  P[VT_Na_o3].name           = "Na_o3";
  P[VT_CJj_nj].name          = "CJj_nj";
  P[VT_CJNa].name            = "CJNa";

  P[VT_RTdF].readFromFile           = false;
  P[VT_FdRT].readFromFile           = false;
  P[VT_RTd2F].readFromFile          = false;
  P[VT_Dlcell].readFromFile         = false;
  P[VT_Ddcell].readFromFile         = false;
  P[VT_Dvcell].readFromFile         = false;
  P[VT_rstart].readFromFile         = false;
  P[VT_rend].readFromFile           = false;
  P[VT_Aj_nj].readFromFile          = false;
  P[VT_xj_nj].readFromFile          = false;
  P[VT_xj_nj_Nai].readFromFile      = false;
  P[VT_V_nonjunct1].readFromFile    = false;
  P[VT_V_nonjunct2].readFromFile    = false;
  P[VT_V_nonjunct3].readFromFile    = false;
  P[VT_V_nonjunct4].readFromFile    = false;
  P[VT_V_nonjunct_Nai].readFromFile = false;
  P[VT_V_cytosol].readFromFile      = false;
  P[VT_V_SR1].readFromFile          = false;
  P[VT_V_SR2].readFromFile          = false;
  P[VT_V_SR3].readFromFile          = false;
  P[VT_V_SR4].readFromFile          = false;
  P[VT_exp_tau_fca].readFromFile    = false;
  P[VT_pow_kNaKNa].readFromFile     = false;
  P[VT_k1].readFromFile             = false;
  P[VT_k2].readFromFile             = false;
  P[VT_k3].readFromFile             = false;
  P[VT_powKo].readFromFile          = false;
  P[VT_Na_o3].readFromFile          = false;
  P[VT_CJj_nj].readFromFile         = false;
  P[VT_CJNa].readFromFile           = false;


  ParameterLoader EPL(initFile, EMT_KoivumaekiEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // KoivumaekiEtAlParameters::Init

void KoivumaekiEtAlParameters::Calculate() {
#if KADEBUG
  cerr << "#KoivumaekiEtAlParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_RTdF].value   = P[VT_R].value*(P[VT_T].value)/(P[VT_F].value);
  P[VT_FdRT].value   = 1.0/(P[VT_RTdF].value);
  P[VT_RTd2F].value  = P[VT_RTdF].value*.5;
  P[VT_Dlcell].value = P[VT_cell_dillation].value;
  P[VT_Ddcell].value = (P[VT_Dlcell].value - 1.0)*(21.0/12.0) + 1.0;
  P[VT_Dvcell].value = P[VT_Dlcell].value*P[VT_Ddcell].value*P[VT_Ddcell].value;

  P[VT_r_junct].value   = P[VT_r_junct].value*P[VT_Ddcell].value;
  P[VT_l_cell].value    = P[VT_l_cell].value*P[VT_Dlcell].value;
  P[VT_dx].value        = P[VT_dx].value*P[VT_Ddcell].value;
  P[VT_totCMDN].value   = P[VT_totCMDN].value/P[VT_Ddcell].value;
  P[VT_totSLlow].value  = P[VT_totSLlow].value/P[VT_Ddcell].value;
  P[VT_totSLhigh].value = P[VT_totSLhigh].value/P[VT_Ddcell].value;
  P[VT_totHTRPN].value  = P[VT_totHTRPN].value/((P[VT_Dvcell].value-1.0)/4.0+1.0);
  P[VT_totLTRPN].value  = P[VT_totLTRPN].value/((P[VT_Dvcell].value-1.0)/4.0+1.0);
  P[VT_totBNa].value    = P[VT_totBNa].value/P[VT_Ddcell].value;
  P[VT_cpumps].value    = P[VT_cpumps].value/P[VT_Dvcell].value;

  P[VT_rstart].value = 0.5*P[VT_dx].value;
  P[VT_rend].value   = P[VT_r_junct].value - P[VT_rstart].value;
  P[VT_Aj_nj].value  = M_PI*P[VT_r_junct].value*2.0*P[VT_l_cell].value*0.5; // Area between junct and nonjunct
  P[VT_xj_nj].value  = 0.01*P[VT_Ddcell].value + P[VT_rstart].value; // diffusion distance from center to center of
                                                                     // junct to first njunct
  P[VT_xj_nj_Nai].value = 0.01*P[VT_Ddcell].value + 2*P[VT_dx].value; // diffusion distance from center of junct to
                                                                      // center of njunct (between 2nd and 3rd njunct)
  P[VT_V_nonjunct1].value = (M_PI*(P[VT_dx].value)*(P[VT_dx].value)*P[VT_l_cell].value)*1e-6*0.5;
  P[VT_V_nonjunct2].value =
    (M_PI*(2.*P[VT_dx].value)*(2.*P[VT_dx].value)*P[VT_l_cell].value-M_PI*(1.*P[VT_dx].value)*(1.*P[VT_dx].value)*
     P[VT_l_cell].value)*1e-6*0.5;
  P[VT_V_nonjunct3].value =
    (M_PI*(3.*P[VT_dx].value)*(3.*P[VT_dx].value)*P[VT_l_cell].value-M_PI*(2.*P[VT_dx].value)*(2.*P[VT_dx].value)*
     P[VT_l_cell].value)*1e-6*0.5;
  P[VT_V_nonjunct4].value =
    (M_PI*(4.*P[VT_dx].value)*(4.*P[VT_dx].value)*P[VT_l_cell].value-M_PI*(3.*P[VT_dx].value)*(3.*P[VT_dx].value)*
     P[VT_l_cell].value)*1e-6*0.5;
  P[VT_V_nonjunct_Nai].value = P[VT_V_nonjunct1].value+P[VT_V_nonjunct2].value+P[VT_V_nonjunct3].value+
    P[VT_V_nonjunct4].value;
  P[VT_V_cytosol].value  = P[VT_V_nonjunct_Nai].value+P[VT_V_ss].value;
  P[VT_V_SR1].value      = 0.05*P[VT_V_nonjunct1].value/2.0*0.9/P[VT_Dvcell].value;
  P[VT_V_SR2].value      = 0.05*P[VT_V_nonjunct2].value/2.0*0.9/P[VT_Dvcell].value;
  P[VT_V_SR3].value      = 0.05*P[VT_V_nonjunct3].value/2.0*0.9/P[VT_Dvcell].value;
  P[VT_V_SR4].value      = 0.05*P[VT_V_nonjunct4].value/2.0*0.9/P[VT_Dvcell].value;
  P[VT_pow_kNaKNa].value = pow(P[VT_kNaKNa].value, 1.5);
  P[VT_k3].value         = P[VT_k4].value / (P[VT_SERCAKmr].value*P[VT_SERCAKmr].value);
  P[VT_k1].value         = 1000000.0 * P[VT_k4].value;
  P[VT_k2].value         = P[VT_k1].value * (P[VT_SERCAKmf].value*P[VT_SERCAKmf].value);
  P[VT_powKo].value      = pow(P[VT_K_o].value, 0.445700);
  P[VT_Na_o3].value      = P[VT_Na_o].value*P[VT_Na_o].value*P[VT_Na_o].value;
  P[VT_CJj_nj].value     = P[VT_DCa].value * P[VT_Aj_nj].value / P[VT_xj_nj].value*1e-6;
  P[VT_CJNa].value       = P[VT_DNa].value * P[VT_Aj_nj].value / P[VT_xj_nj_Nai].value* 1e-6;
}  // KoivumaekiEtAlParameters::Calculate

void KoivumaekiEtAlParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#KoivumaekiEtAlParameters - InitTable ..." << endl;
#endif  // if KADEBUG
  ML_CalcType HT = tinc;
  P[VT_exp_tau_fca].value = exp(-HT/P[VT_tau_fca].value);
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    exptau_r[Vi]       =  exp(-HT/ (0.00350000*(exp((((-V*V)/30.0000)/30.0000)))+0.00150000));
    r_infinity[Vi]     = 1.00000/(1.00000+(exp(((V - 1.00000)/ -11.0000))));
    sus_r_infinity[Vi] = 1.00000/(1.00000+(exp((-(V+6.00000)/8.60000))));
    exptau_sus_r[Vi]   =  exp(-HT/ (0.00900000/(1.00000+(exp(((V+5.00000)/12.0000))))+0.000500000));
    sus_s_infinity[Vi] = 1.00000/(1.00000+(exp(((V+7.50000)/10.0000))));
    exptau_sus_s[Vi]   =  exp(-HT/ (0.590000/(1.00000+(exp(((V+60.0000)/10.0000))))+3.05000));
    m_infinity[Vi]     = 1.00000/(1.00000+(exp(((V+27.1200)/ -8.21000))));
    ML_CalcType m_factor = (V+25.5700)/28.8000;
    exptau_m[Vi]   =  exp(-HT/ (4.20000e-05*(exp((-m_factor*m_factor)))+2.40000e-05));
    h_infinity[Vi] = 1.00000/(1.00000+(exp(((V+63.6000)/5.30000))));
    ML_CalcType h_factor = 1.00000/(1.00000+(exp(((V+35.1000)/3.20000))));
    exptau_h1[Vi]  =  exp(-HT/ (0.0300000*h_factor+0.000300000));
    exptau_h2[Vi]  =  exp(-HT/ (0.120000*h_factor+0.00300000));
    d_infinity[Vi] = 1.00000/(1.00000+(exp(((V+9.00000)/ -5.80000))));
    ML_CalcType d_L_factor = (V+35.0000)/30.0000;
    exptau_d[Vi]   =  exp(-HT/ (0.00270000*(exp((-d_L_factor*d_L_factor)))+0.00200000));
    f_infinity[Vi] = 1.00000/(1.00000+(exp(((V+27.4000)/7.10000))));
    ML_CalcType f_L_factor = (V+30.16047)/7.09396;
    exptau_f1[Vi] =
      exp(-HT/
          (0.98698*exp(-f_L_factor*f_L_factor) + 0.04275/(1.0+exp((V-51.61555)/ -80.61331)) + 0.03576/
           (1+exp((V+29.57272)/13.21758)) - 0.00821));
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
    exptau_y[Vi]     =  exp(-HT * (0.00332*exp(-V/16.54103)+23.71839*exp(V/16.54103)));
    y_infinity[Vi]   = 1.00000/(1.00000+(exp(((V+97.82874)/12.48025))));
    CexpINa[Vi]      = P[VT_Na_o].value * V * P[VT_F].value*P[VT_FdRT].value / (exp(V*P[VT_FdRT].value) - 1.0);
    CI_NaCa1[Vi]     = P[VT_Ca_o].value*exp(P[VT_FdRT].value*V*P[VT_gam].value);
    CI_NaCa2[Vi]     = P[VT_Na_o3].value*exp( (P[VT_gam].value - 1.00000)*V*P[VT_FdRT].value);
  }
}  // KoivumaekiEtAlParameters::InitTable
