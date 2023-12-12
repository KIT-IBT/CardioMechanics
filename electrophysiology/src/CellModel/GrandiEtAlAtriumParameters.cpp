/*
 * File: GrandiEtAlAtriumParameters.cpp
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


#include <GrandiEtAlAtriumParameters.h>
GrandiEtAlAtriumParameters::GrandiEtAlAtriumParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

GrandiEtAlAtriumParameters::~GrandiEtAlAtriumParameters() {}

void GrandiEtAlAtriumParameters::PrintParameters() {
  cout << "GrandiEtAlAtriumParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void GrandiEtAlAtriumParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the GrandiEtAlAtrium parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

  P[VT_stim_amplitude].name = "stim_amplitude";
  P[VT_stim_duration].name  = "stim_duration";
  P[VT_R].name              = "R";
  P[VT_Cmem].name           = "Cmem";
  P[VT_cellLength].name     = "cellLength";
  P[VT_junctionLength].name = "junctionLength";
  P[VT_distSLcyto].name     = "distSLcyto";
  P[VT_distJuncSL].name     = "distJuncSL";
  P[VT_DcaJuncSL].name      = "DcaJuncSL";
  P[VT_DcaSLcyto].name      = "DcaSLcyto";
  P[VT_DnaJuncSL].name      = "DnaJuncSL";
  P[VT_DnaSLcyto].name      = "DnaSLcyto";
  P[VT_J_ca_juncsl].name    = "J_ca_juncsl";
  P[VT_J_ca_slmyo].name     = "J_ca_slmyo";
  P[VT_J_na_juncsl].name    = "J_na_juncsl";
  P[VT_J_na_slmyo].name     = "J_na_slmyo";
  P[VT_Fjunc].name          = "Fjunc";
  P[VT_Fjunc_CaL].name      = "Fjunc_CaL";
  P[VT_K_o].name            = "K_o";
  P[VT_Na_o].name           = "Na_o";
  P[VT_Ca_o].name           = "Ca_o";
  P[VT_Mgi].name            = "Mgi";
  P[VT_AF].name             = "AF";
  P[VT_ISO].name            = "ISO";
  P[VT_RA].name             = "RA";
  P[VT_GNaL].name           = "GNaL";
  P[VT_GNa].name            = "GNa";
  P[VT_GK1].name            = "GK1";
  P[VT_GNaB].name           = "GNaB";
  P[VT_IbarNaK].name        = "IbarNaK";
  P[VT_KmNaip].name         = "KmNaip";
  P[VT_KmKo].name           = "KmKo";
  P[VT_Q10NaK].name         = "Q10NaK";
  P[VT_Q10KmNai].name       = "Q10KmNai";
  P[VT_pNaK].name           = "pNaK";
  P[VT_gkp].name            = "gkp";
  P[VT_GClCa].name          = "GClCa";
  P[VT_GClB].name           = "GClB";
  P[VT_KdClCa].name         = "KdClCa";
  P[VT_pNa].name            = "pNa";
  P[VT_pCa].name            = "pCa";
  P[VT_pK].name             = "pK";
  P[VT_Q10CaL].name         = "Q10CaL";
  P[VT_IbarNCX].name        = "IbarNCX";
  P[VT_KmCai].name          = "KmCai";
  P[VT_KmCao].name          = "KmCao";
  P[VT_KmNai].name          = "KmNai";
  P[VT_KmNao].name          = "KmNao";
  P[VT_ksat].name           = "ksat";
  P[VT_nu].name             = "nu";
  P[VT_Kdact].name          = "Kdact";
  P[VT_Q10NCX].name         = "Q10NCX";
  P[VT_IbarSLCaP].name      = "IbarSLCaP";
  P[VT_KmPCa].name          = "KmPCa";
  P[VT_GCaB].name           = "GCaB";
  P[VT_Q10SLCaP].name       = "Q10SLCaP";
  P[VT_Q10SRCaP].name       = "Q10SRCaP";
  P[VT_GtoFast].name        = "GtoFast";
  P[VT_Gkur].name           = "Gkur";
  P[VT_Vmax_SRCaP].name     = "Vmax_SRCaP";
  P[VT_Kmf].name            = "Kmf";
  P[VT_Kmr].name            = "Kmr";
  P[VT_hillSRCaP].name      = "hillSRCaP";
  P[VT_ks].name             = "ks";
  P[VT_koCa].name           = "koCa";
  P[VT_kom].name            = "kom";
  P[VT_kim].name            = "kim";
  P[VT_ec50SR].name         = "ec50SR";
  P[VT_Bmax_Naj].name       = "Bmax_Naj";
  P[VT_Bmax_Nasl].name      = "Bmax_Nasl";
  P[VT_koff_na].name        = "koff_na";
  P[VT_kon_na].name         = "kon_na";
  P[VT_Bmax_TnClow].name    = "Bmax_TnClow";
  P[VT_koff_tncl].name      = "koff_tncl";
  P[VT_kon_tncl].name       = "kon_tncl";
  P[VT_Bmax_TnChigh].name   = "Bmax_TnChigh";
  P[VT_koff_tnchca].name    = "koff_tnchca";
  P[VT_kon_tnchca].name     = "kon_tnchca";
  P[VT_koff_tnchmg].name    = "koff_tnchmg";
  P[VT_kon_tnchmg].name     = "kon_tnchmg";
  P[VT_Bmax_CaM].name       = "Bmax_CaM";
  P[VT_koff_cam].name       = "koff_cam";
  P[VT_kon_cam].name        = "kon_cam";
  P[VT_Bmax_myosin].name    = "Bmax_myosin";
  P[VT_koff_myoca].name     = "koff_myoca";
  P[VT_kon_myoca].name      = "kon_myoca";
  P[VT_koff_myomg].name     = "koff_myomg";
  P[VT_kon_myomg].name      = "kon_myomg";
  P[VT_Bmax_SR].name        = "Bmax_SR";
  P[VT_koff_sr].name        = "koff_sr";
  P[VT_kon_sr].name         = "kon_sr";
  P[VT_koff_sll].name       = "koff_sll";
  P[VT_kon_sll].name        = "kon_sll";
  P[VT_koff_slh].name       = "koff_slh";
  P[VT_kon_slh].name        = "kon_slh";
  P[VT_koff_csqn].name      = "koff_csqn";
  P[VT_kon_csqn].name       = "kon_csqn";
  P[VT_gks_junc].name       = "gks_junc";
  P[VT_gks_sl].name         = "gks_sl";
  P[VT_fcaCaMSL].name       = "fcaCaMSL";
  P[VT_fcaCaj].name         = "fcaCaj";
  P[VT_MaxSR].name          = "MaxSR";
  P[VT_MinSR].name          = "MinSR";
  P[VT_Frdy].name           = "Frdy";
  P[VT_cellRadius].name     = "cellRadius";
  P[VT_junctionRadius].name = "junctionRadius";
  P[VT_Fsl].name            = "Fsl";
  P[VT_Fsl_CaL].name        = "Fsl_CaL";
  P[VT_sigma].name          = "sigma";
  P[VT_gkr].name            = "gkr";
  P[VT_gkr0].name           = "Gkr0";
  P[VT_kiCa].name           = "kiCa";
  P[VT_Temp].name           = "Temp";
  P[VT_Vcell].name          = "Vcell";
  P[VT_SAjunc].name         = "SAjunc";
  P[VT_SAsl].name           = "SAsl";
  P[VT_FoRT].name           = "FoRT";
  P[VT_Qpow].name           = "Qpow";
  P[VT_Vmyo].name           = "Vmyo";
  P[VT_Vsr].name            = "Vsr";
  P[VT_Vsl].name            = "Vsl";
  P[VT_Vjunc].name          = "Vjunc";
  P[VT_Cli].name            = "Cli";
  P[VT_Clo].name            = "Clo";
  P[VT_Bmax_SLlowsl].name   = "Bmax_SLlowsl";
  P[VT_Bmax_SLlowj].name    = "Bmax_SLlowj";
  P[VT_Bmax_SLhighsl].name  = "Bmax_SLhighsl";
  P[VT_Bmax_SLhighj].name   = "Bmax_SLhighj";
  P[VT_Bmax_Csqn].name      = "Bmax_Csqn";
  P[VT_ecl].name            = "ecl";
  P[VT_V_init].name         = "V_init";
  P[VT_Na_j_init].name      = "Na_j_init";
  P[VT_Na_sl_init].name     = "Na_sl_init";
  P[VT_K_i_init].name       = "K_i_init";
  P[VT_Ca_j_init].name      = "Ca_j_init";
  P[VT_Ca_sl_init].name     = "Ca_sl_init";
  P[VT_m_init].name         = "m_init";
  P[VT_h_init].name         = "h_init";
  P[VT_j_init].name         = "j_init";
  P[VT_ml_init].name        = "ml_init";
  P[VT_hl_init].name        = "hl_init";
  P[VT_x_kr_init].name      = "x_kr_init";
  P[VT_x_ks_init].name      = "x_ks_init";
  P[VT_Na_i_init].name      = "Na_i_init";
  P[VT_x_kur_init].name     = "x_kur_init";
  P[VT_y_kur_init].name     = "y_kur_init";
  P[VT_x_to_f_init].name    = "x_to_f_init";
  P[VT_y_to_f_init].name    = "y_to_f_init";
  P[VT_d_init].name         = "d_init";
  P[VT_f_init].name         = "f_init";
  P[VT_f_Ca_Bj_init].name   = "f_Ca_Bj_init";
  P[VT_f_Ca_Bsl_init].name  = "f_Ca_Bsl_init";
  P[VT_Ry_Rr_init].name     = "Ry_Rr_init";
  P[VT_Ry_Ro_init].name     = "Ry_Ro_init";
  P[VT_Ry_Ri_init].name     = "Ry_Ri_init";
  P[VT_Ca_sr_init].name     = "Ca_sr_init";
  P[VT_Ca_i_init].name      = "Ca_i_init";
  P[VT_Na_Bj_init].name     = "Na_Bj_init";
  P[VT_Na_Bsl_init].name    = "Na_Bsl_init";
  P[VT_Tn_CL_init].name     = "Tn_CL_init";
  P[VT_Tn_CHc_init].name    = "Tn_CHc_init";
  P[VT_Tn_CHm_init].name    = "Tn_CHm_init";
  P[VT_CaM_init].name       = "CaM_init";
  P[VT_Myo_c_init].name     = "Myo_c_init";
  P[VT_Myo_m_init].name     = "Myo_m_init";
  P[VT_SRB_init].name       = "SRB_init";
  P[VT_SLL_j_init].name     = "SLL_j_init";
  P[VT_SLL_sl_init].name    = "SLL_sl_init";
  P[VT_SLH_j_init].name     = "SLH_j_init";
  P[VT_SLH_sl_init].name    = "SLH_sl_init";
  P[VT_Csqn_b_init].name    = "Csqn_b_init";
  P[VT_powncx].name         = "powncx";
  P[VT_powcal].name         = "powcal";
  P[VT_powsrcap].name       = "powsrcap";
  P[VT_powslcap].name       = "powslcap";
  P[VT_powKmPCa].name       = "powKmPCa";
  P[VT_Ikiconst].name       = "Ikiconst";

  P[VT_Fsl].readFromFile           = false;
  P[VT_Fsl_CaL].readFromFile       = false;
  P[VT_sigma].readFromFile         = false;
  P[VT_gkr].readFromFile           = false;
  P[VT_Vcell].readFromFile         = false;
  P[VT_SAjunc].readFromFile        = false;
  P[VT_SAsl].readFromFile          = false;
  P[VT_FoRT].readFromFile          = false;
  P[VT_Qpow].readFromFile          = false;
  P[VT_Vmyo].readFromFile          = false;
  P[VT_Vsr].readFromFile           = false;
  P[VT_Vsl].readFromFile           = false;
  P[VT_Vjunc].readFromFile         = false;
  P[VT_Bmax_SLlowsl].readFromFile  = false;
  P[VT_Bmax_SLlowj].readFromFile   = false;
  P[VT_Bmax_SLhighsl].readFromFile = false;
  P[VT_Bmax_SLhighj].readFromFile  = false;
  P[VT_Bmax_Csqn].readFromFile     = false;
  P[VT_ecl].readFromFile           = false;
  P[VT_powncx].readFromFile        = false;
  P[VT_powcal].readFromFile        = false;
  P[VT_powsrcap].readFromFile      = false;
  P[VT_powslcap].readFromFile      = false;
  P[VT_powKmPCa].readFromFile      = false;
  P[VT_Ikiconst].readFromFile      = false;

  ParameterLoader EPL(initFile, EMT_GrandiEtAlAtrium);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // GrandiEtAlAtriumParameters::Init

void GrandiEtAlAtriumParameters::Calculate() {
#if KADEBUG
  cerr << "#GrandiEtAlAtriumParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_Fsl].value     = 1.00000 - (P[VT_Fjunc].value);
  P[VT_Fsl_CaL].value = 1.00000 - (P[VT_Fjunc_CaL].value);
  P[VT_sigma].value   = ((exp(((P[VT_Na_o].value)/67.3000))) - 1.00000)/7.00000;
  P[VT_gkr].value     =  P[VT_gkr0].value* pow(((P[VT_K_o].value)/5.40000), 1.0 / 2);
  P[VT_Vcell].value   =   3.14159265358979*
    (pow((P[VT_cellRadius].value), 2.00000))*(P[VT_cellLength].value)*1.00000e-15;
  P[VT_SAjunc].value =  20150.0* 3.14159265358979*2.00000*(P[VT_junctionLength].value)*
    (P[VT_junctionRadius].value);
  P[VT_SAsl].value          =   3.14159265358979*2.00000*(P[VT_cellRadius].value)*(P[VT_cellLength].value);
  P[VT_FoRT].value          = (P[VT_Frdy].value)/( (P[VT_R].value)*(P[VT_Temp].value));
  P[VT_Qpow].value          = ((P[VT_Temp].value) - 310.000)/10.0000;
  P[VT_Vmyo].value          =  0.650000*(P[VT_Vcell].value);
  P[VT_Vsr].value           =  0.0350000*(P[VT_Vcell].value);
  P[VT_Vsl].value           =  0.0200000*(P[VT_Vcell].value);
  P[VT_Vjunc].value         =  0.0539000*0.0100000*(P[VT_Vcell].value);
  P[VT_Bmax_SLlowsl].value  = (0.0374000*(P[VT_Vmyo].value))/(P[VT_Vsl].value);
  P[VT_Bmax_SLlowj].value   =  ((0.00460000*(P[VT_Vmyo].value))/(P[VT_Vjunc].value))*0.100000;
  P[VT_Bmax_SLhighsl].value = (0.0134000*(P[VT_Vmyo].value))/(P[VT_Vsl].value);
  P[VT_Bmax_SLhighj].value  =  ((0.00165000*(P[VT_Vmyo].value))/(P[VT_Vjunc].value))*0.100000;
  P[VT_Bmax_Csqn].value     = (0.140000*(P[VT_Vmyo].value))/(P[VT_Vsr].value);
  P[VT_ecl].value           =  (1.00000/(P[VT_FoRT].value))*(log(((P[VT_Cli].value)/(P[VT_Clo].value))));
  P[VT_powncx].value        = pow(P[VT_Q10NCX].value, P[VT_Qpow].value);
  P[VT_powcal].value        = pow(P[VT_Q10CaL].value, P[VT_Qpow].value);
  P[VT_powsrcap].value      = pow(P[VT_Q10SRCaP].value, P[VT_Qpow].value);
  P[VT_powslcap].value      = pow(P[VT_Q10SLCaP].value, P[VT_Qpow].value);
  P[VT_powKmPCa].value      = pow(P[VT_KmPCa].value, 1.60000);
  P[VT_Ikiconst].value      = P[VT_GK1].value*(1.0+1.0*P[VT_AF].value)*sqrt(P[VT_K_o].value/5.4);
  P[VT_GNa].value           = P[VT_GNa].value*(1.0-0.1*P[VT_AF].value);
  P[VT_KmNaip].value        = P[VT_KmNaip].value*(1.0-0.25*P[VT_ISO].value);
  P[VT_pNa].value           = (1+0.5*P[VT_ISO].value)*(1-0.5*P[VT_AF].value)*P[VT_pNa].value;
  P[VT_pCa].value           = (1+0.5*P[VT_ISO].value)*(1-0.5*P[VT_AF].value)*P[VT_pCa].value;
  P[VT_pK].value            = (1+0.5*P[VT_ISO].value)*(1-0.5*P[VT_AF].value)*P[VT_pK].value;
  P[VT_IbarNCX].value       = (1+0.4*P[VT_AF].value)*P[VT_IbarNCX].value;
  P[VT_Kmf].value           = (2.5-1.25*P[VT_ISO].value)*P[VT_Kmf].value;
  P[VT_koCa].value          = P[VT_koCa].value*(1.0+2.0*P[VT_AF].value+P[VT_ISO].value*(1.0-P[VT_AF].value));
  P[VT_koff_tncl].value     = (1.+0.5*P[VT_ISO].value)*P[VT_koff_tncl].value;
  P[VT_gks_junc].value      = (1.0+1.0*P[VT_AF].value+2.0*P[VT_ISO].value)*P[VT_gks_junc].value;
  P[VT_gks_sl].value        = (1.0+1.0*P[VT_AF].value+2.0*P[VT_ISO].value)*P[VT_gks_sl].value;
  P[VT_GtoFast].value       = (1.0-0.7*P[VT_AF].value)*P[VT_GtoFast].value;
  P[VT_GNaL].value          = P[VT_AF].value*P[VT_GNaL].value;
  P[VT_Gkur].value          = (1.0-0.5*P[VT_AF].value)*(1.0+2.0*P[VT_ISO].value)* P[VT_Gkur].value *
    (1.0+0.2*P[VT_RA].value);
}  // GrandiEtAlAtriumParameters::Calculate

void GrandiEtAlAtriumParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#GrandiEtAlAtriumParameters - InitTable ..." << endl;
#endif  // if KADEBUG
  ML_CalcType HT = tinc * 1000;
  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);
    mss[Vi]     = 1.00000/(pow((1.00000+(exp((-(56.8600+V)/9.03000)))), 2.00000));
    exptaum[Vi] =
      exp(-HT/
          (0.129200*
           (exp((-(pow(((V+45.7900)/15.5400),
                       2.00000)))))+ 0.0648700*(exp((-(pow(((V - 4.82300)/51.1200), 2.00000)))))));
    xrss[Vi]     = 1.00000/(1.00000+(exp((-(V+10.0000)/5.00000))));
    exptauxr[Vi] =
      exp(-HT/
          ( ( (550.000/(1.00000+(exp(((-22.0000 - V)/9.00000)))))*6.00000)/(1.00000+(exp(((V + 11.0000)/9.00000))))+
            230.000/
            (1.00000+(exp(((V + 40.0000)/20.0000))))));
    xsss[Vi]       = 1.00000/(1.00000+(exp((-(V+40.0*P[VT_ISO].value+3.80000)/14.2500))));
    exptauxs[Vi]   =  exp(-HT/ (990.100/(1.00000+(exp((-(V+40.0*P[VT_ISO].value+2.43600)/14.1200))))));
    xkurss[Vi]     = 1.00000/(1.00000+(exp((-(V + 6.0)/8.6))));
    ykurss[Vi]     = 1.00000/(1.00000+(exp(((V+7.5)/10.0))));
    exptauxkur[Vi] =  exp(-HT/ (9.0/(1.0+exp((V+5.0)/12.0))+0.5));
    exptauykur[Vi] =  exp(-HT/ (590.0/(1.0+exp((V+60.0)/10.0))+3050.0));
    xtoss[Vi]      = 1.00000/(1.00000+(exp((-(V + 1.0)/11.0000))));
    ytoss[Vi]      = 1.00000/(1.00000+(exp(((V+40.5000)/11.50000))));
    exptauxtof[Vi] =  exp(-HT/ (3.50000*(exp((-(pow(((V)/30.0000), 2.00000)))))+1.500000));
    exptauytof[Vi] =  exp(-HT/ (25.635*(exp((-(pow(((V+52.45)/15.8827), 2.00000)))))+24.14));
    dss[Vi]        = 1.00000/(1.00000+(exp((-(V+3.0*P[VT_ISO].value+9.00000)/6.00000))));
    exptaud[Vi]    =
      exp(-HT/
          ((1.00000*dss[Vi]*(1.00000 - (exp((-(V+3.0*P[VT_ISO].value+9.00000)/6.00000)))))/
           (0.0350000*(V+3.0*P[VT_ISO].value+9.00000))));
    fss[Vi] = 1.00000/(1.00000+(exp(((V+3.0*P[VT_ISO].value+30.0000)/7.00000))))+0.200000/
      (1.00000+(exp(((50.0000 - V-3.0*P[VT_ISO].value)/20.0000))));
    exptauf[Vi] =
      exp(-HT/ (1.00000/(0.0197000*(exp((-(pow((0.0337000*(V+3.0*P[VT_ISO].value+25.000)), 2.00000)))))+0.0200000)));
    ML_CalcType ah = (V >= -40.0000 ? 0.00000 :  0.0570000*(exp((-(V+80.0000)/6.80000))));
    ML_CalcType bh =
      (V >=
       -40.0000 ? 0.770000/(0.130000*(1.00000+(exp((-(V+10.6600)/11.1000))))) :  2.70000*(exp((0.0790000*V)))+ 310000.*
       (exp((0.348500*V))));
    exptauh[Vi] = exp(-HT/ (1.00000/(ah+bh)));
    hss[Vi]     = 1.00000/(pow((1.00000+(exp(((V+71.5500)/7.43000)))), 2.00000));
    ML_CalcType aj =
      (V >=
       -40.0000 ? 0.00000 : ( (-25428.0*(exp((0.244400*V))) -  6.94800e-06*(exp((-0.0439100*V))))*(V+37.7800))/
       (1.00000+(exp((0.311000*(V+79.2300))))));
    ML_CalcType bj =
      (V >=
       -40.0000 ? (0.600000*(exp((0.0570000*V))))/
       (1.00000+(exp((-0.100000*(V+32.0000))))) : (0.0242400*(exp((-0.0105200*V))))/
       (1.00000+(exp((-0.137800*(V+40.1400))))));
    exptauj[Vi] = exp(-HT/ (1.00000/(aj+bj)));
    jss[Vi]     = 1.00000/(pow((1.00000+(exp(((V+71.5500)/7.43000)))), 2.00000));
    const ML_CalcType aml = 0.32*(V+47.13)/(1.0-exp(-0.1*(V+47.13)));
    const ML_CalcType bml = 0.08*exp(-V/11);
    mlss[Vi]     = aml/(aml+bml);
    exptauml[Vi] = 1.0/(aml+bml);
    hlss[Vi]     = 1.0/(1.0+exp((V+91.0)/6.1));
    fnak[Vi]     = 1.00000/
      (1.00000+ 0.124500*(exp((-0.100000*V*(P[VT_FoRT].value))))+ 0.0365000*(P[VT_sigma].value)*
       (exp((-V*(P[VT_FoRT].value)))));
    rkr[Vi]    = 1.00000/(1.00000+(exp(((V+74.0000)/24.0000))));
    kp_kp[Vi]  = 1.00000/(1.00000+(exp((7.48800 - V/5.98000))));
    I_Clbk[Vi] =  (P[VT_GClB].value)*(V - (P[VT_ecl].value));
  }
}  // GrandiEtAlAtriumParameters::InitTable
