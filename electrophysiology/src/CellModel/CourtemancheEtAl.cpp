/* -------------------------------------------------------

   CourtemancheEtAl.cpp

   Ver. 1.1.0

   Created:       dw (27.02.2007)
   Last modified: Tobias Gerach (30.06.2022)

   Institute of Biomedical Engineering
   Karlsruhe Institute of Technology (KIT)

   http://www.ibt.kit.edu

   Copyright 2000-2009 - All rights reserved.

   ------------------------------------------------------ */

#include <CourtemancheEtAl.h>

Courtemanche::Courtemanche(CourtemancheParameters *pp) {
  pCmP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(pCmP, NS_CourtemancheParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

Courtemanche::~Courtemanche() {}

#ifdef HETERO

inline bool Courtemanche::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool Courtemanche::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int Courtemanche::GetSize(void) {
  // return ( &fT - &Ca_i + 1 ) * sizeof( ML_CalcType );
  return sizeof(Courtemanche)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(CourtemancheParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char Courtemanche::getSpeed(ML_CalcType adVm) {
  return (unsigned char)(adVm < 1e-6 ? 2 : 1);
}

void Courtemanche::Init() {
#if KADEBUG
  cerr << "Courtemanche::Init" << endl;
#endif  // if KADEBUG
  Na_i = v(VT_Init_Na_i);
  K_i  = v(VT_Init_K_i);
  m    = v(VT_Init_m);
  oa   = v(VT_Init_oa);
  ua   = v(VT_Init_ua);
  d    = v(VT_Init_d);
  f_Ca = v(VT_Init_f_Ca);
  Ca_i = v(VT_Init_Ca_i);
  Ca_up  = v(VT_Init_Ca_up);
  Ca_rel = v(VT_Init_Ca_rel);
  h      = v(VT_Init_h);
  j      = v(VT_Init_j);
  oi     = v(VT_Init_oi);
  ui     = v(VT_Init_ui);
  Xr     = v(VT_Init_Xr);
  Xs     = v(VT_Init_Xs);
  f      = v(VT_Init_f);
  u      = v(VT_Init_u);
  v      = v(VT_Init_v);
  w      = v(VT_Init_w);
  #ifdef TRPN
  Ca_TRPN = v(VT_Init_Ca_TRPN)/v(VT_TRPNmax);
  #endif  // ifdef TRPN
}  // Courtemanche::Init

ML_CalcType Courtemanche::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                               int euler                                            = 2) {
  tinc *= 1000.0;
  const  ML_CalcType V_int = V*1000.0;
  const int Vi             = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);
  const double dNa_i       = 1.0/Na_i;
  const double VmE_Na      = V_int-((v(VT_RTdF))*log((v(VT_Na_o))*dNa_i));
  const double VmE_K       = V_int-((v(VT_RTdF))*log((v(VT_K_o))/K_i));
  const double VmE_Ca      = V_int-((v(VT_RTd2F))*log((v(VT_Ca_o))/(double)Ca_i));
  const double m_m         = pCmP->m_m[Vi];
  m = m_m + (m - m_m) * pCmP->exptau_m[Vi];
  checkGatingVariable(m);
  const double m_h = pCmP->m_h[Vi];
  h = m_h + (h - m_h) * pCmP->exptau_h[Vi];
  const double m_j = pCmP->m_j[Vi];
  j = m_j + (j - m_j) * pCmP->exptau_j[Vi];
  const double m_oa = pCmP->m_oa[Vi];
  oa = m_oa + (oa - m_oa) * pCmP->exptau_oa[Vi];
  const double m_oi = pCmP->m_oi[Vi];
  oi = m_oi + (oi - m_oi) * pCmP->exptau_oi[Vi];
  const double m_ua = pCmP->m_ua[Vi];
  ua = m_ua + (ua - m_ua) * pCmP->exptau_ua[Vi];
  const double m_ui = pCmP->m_ui[Vi];
  ui = m_ui + (ui - m_ui) * pCmP->exptau_ui[Vi];
  const double m_Xr = pCmP->m_Xr[Vi];
  Xr = m_Xr + (Xr - m_Xr) * pCmP->exptau_Xr[Vi];
  const double m_Xs = pCmP->m_Xs[Vi];
  Xs = m_Xs + (Xs - m_Xs) * pCmP->exptau_Xs[Vi];
  const double m_d = pCmP->m_d[Vi];
  d = m_d + (d - m_d) * pCmP->exptau_d[Vi];
  const double m_f = pCmP->m_f[Vi];
  f     = m_f + (f - m_f) * pCmP->exptau_f[Vi];
  f_Ca += (1.0/(1.0+Ca_i*2857.1429)-f_Ca)*0.5*tinc;
  const double m_w = pCmP->m_w[Vi];
  w = m_w + (w - m_w) * pCmP->exptau_w[Vi];

  const double I_Na  = (v(VT_g_Na))*m*m*m*h*j*VmE_Na;
  const double I_to  = (v(VT_g_to))*oa*oa*oa*oi*VmE_K;
  const double I_Kur = pCmP->g_Kur[Vi]*ua*ua*ua*ui*VmE_K;

  const double I_Kr = pCmP->CKr[Vi]*Xr*VmE_K;
  const double I_Ks  = (v(VT_g_Ks))*Xs*Xs*VmE_K;
  const double I_CaL = (v(VT_g_CaL))*d*f*f_Ca*(V_int-(v(VT_E_rL)));

  const double I_K1 = pCmP->CK1[Vi]*VmE_K;
  const double I_bNa    = (v(VT_g_bNa))*VmE_Na;
  const double I_bCa    = (v(VT_g_bCa))*VmE_Ca;
  const double temp     = (v(VT_k_mNai))*dNa_i;
  const double I_NaK    = pCmP->CNaK[Vi]/(1.0+temp*sqrt(temp));
  const double I_pCa    = (v(VT_I_pCamax))*Ca_i/(Ca_i+(v(VT_k_mpCa)));
  const double I_NaCa   = pCmP->CNaCa[Vi]*((Na_i*Na_i*Na_i*(v(VT_Ca_o)))-pCmP->expVm[Vi]*Ca_i);
  const double I_rel    = (v(VT_k_rel))*u*u*v*w*(Ca_rel-Ca_i);
  const double I_up     = (v(VT_I_upmax))*Ca_i/(Ca_i+(v(VT_k_up)));
  const double I_upleak = (v(VT_kupleak))*Ca_up;
  const double I_tr     = (Ca_up-Ca_rel)*(v(VT_dt_tr));

  /// stretch activated channel
#ifdef ISAC
  I_SAC = v(VT_G_sac) * (V_int - v(VT_E_sac)) / (1. + v(VT_K_sac) * exp(-v(VT_alpha) * (stretch - 1.)));
  I_SAC *= v(VT_ISAC_SWITCH);
#endif  // ifdef ISAC

#ifdef TRPN
  ML_CalcType Ca_T50    = v(VT_Ca_T50)+v(VT_beta1)*(std::min(1.2, stretch)-1.);
  ML_CalcType dCaTRPN   = v(VT_k_TRPN)*(pow(1000*Ca_i/Ca_T50, v(VT_n_TRPN))*(1.-Ca_TRPN)-Ca_TRPN);// needs muM with Land parameters
  Ca_TRPN += tinc*dCaTRPN;
#endif // ifdef TRPN

  Ca_up  += tinc*(I_up-I_tr*(v(VT_VreldVup))-I_upleak);
  Ca_rel += tinc*((I_tr-I_rel)/(1.0+(v(VT_csqnkm))/(Ca_rel*Ca_rel+Ca_rel*(v(VT_kmcsqnm2))+(v(VT_kkmcsqn)))));

  Na_i   += tinc*(-3.0*I_NaK-3.0*I_NaCa-I_bNa-I_Na)*(v(VT_CmdFvi));
  K_i    += tinc*(2.0*I_NaK-I_K1-I_to-I_Kur-I_Kr-I_Ks)*(v(VT_CmdFvi));

  ML_CalcType B1 = ((2.0*I_NaCa-I_pCa-I_CaL-I_bCa)*.5*(v(VT_CmdFvi)))+
    ((I_upleak-I_up)*(v(VT_VupdVi))+I_rel*(v(VT_VreldVi)));

#ifdef TRPN
  B1 -= v(VT_TRPNmax)*dCaTRPN;
#endif // ifdef TRPN

#ifdef TRPN
  const double B2 =  1.0+((v(VT_cmdnkm))/((Ca_i*Ca_i)+(Ca_i*(v(VT_kmcmdnm2)))+(v(VT_kkmcmdn))));
#else // ifdef TRPN
  const double B2 =  1.0+((v(VT_trpnkm))/((Ca_i*Ca_i)+(Ca_i*(v(VT_kmtrpnm2)))+(v(VT_kkmtrpn))))+
    ((v(VT_cmdnkm))/((Ca_i*Ca_i)+(Ca_i*(v(VT_kmcmdnm2)))+(v(VT_kkmcmdn))));
#endif // ifdef TRPN
  const double dCa_i = B1 / B2;
  Ca_i += tinc * dCa_i;

  const double Fn  = exp(-((v(VT_Vrel))*I_rel-(0.5*I_CaL-0.2*I_NaCa)*(v(VT_Cmd2F)))*731.5289); // to reproduce Fig. 15 in the original publication, the first part of this equation needs to be divided by 100. Reason: unknown
  const double m_u = 1.0/(1.0+Fn*pCmP->exp250);
  u += (m_u-u)*0.125*tinc;
  const double m_v = 1.0-1.0/(1.0+Fn*pCmP->exp50);
  const double t_v = 1.91+2.09*m_u;
  v += (m_v-v)*tinc/t_v;
#ifdef ISAC
  double I_tot = I_Na+I_CaL+I_pCa+I_K1+I_to+I_Kur+I_Kr+I_Ks+I_bNa+I_bCa+I_NaK+I_NaCa+I_SAC-i_external;
#else // ifdef ISAC
  double I_tot = I_Na+I_CaL+I_pCa+I_K1+I_to+I_Kur+I_Kr+I_Ks+I_bNa+I_bCa+I_NaK+I_NaCa-i_external;
#endif // ifdef ISAC
  return -.001*I_tot*tinc;  // *.001 to revert back from ms to s (conversion was done at the beginning of this function)
}  // Courtemanche::Calc

void Courtemanche::Print(ostream &tempstr, double t,  ML_CalcType V) {
  tempstr << t << ' ' << V << ' ' <<
    m << ' ' << h << ' ' << j << ' ' << oa << ' ' << oi << ' ' << ua << ' ' << ui << ' ' << Xr << ' ' << Xs << ' ' <<
    d << ' ' << f << ' ' << f_Ca << ' ' << u << ' ' << v << ' ' << w << ' ' << Na_i << ' ' << Ca_i << ' ' << K_i <<
    ' ' << Ca_up << ' ' << Ca_rel << ' '
#ifdef TRPN
    << v(VT_TRPNmax)*Ca_TRPN << ' ' // export in mM
#endif // ifdef TRPN
  ;
} // Courtemanche::Print

void Courtemanche::LongPrint(ostream &tempstr, double t,  ML_CalcType V) {
  Print(tempstr, t, V);
  const  ML_CalcType V_int = V*1000.0;
  const int Vi             = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);
  const double dNa_i       = 1.0/Na_i;
  const double VmE_Na      = V_int-((v(VT_RTdF))*log((v(VT_Na_o))*dNa_i));
  const double VmE_K       = V_int-((v(VT_RTdF))*log((v(VT_K_o))/K_i));
  const double VmE_Ca      = V_int-((v(VT_RTd2F))*log((v(VT_Ca_o))/(double)Ca_i));
  const double I_Na        = (v(VT_g_Na))*m*m*m*h*j*VmE_Na;
  const double I_to        = (v(VT_g_to))*oa*oa*oa*oi*VmE_K;
  const double I_Kur       = pCmP->g_Kur[Vi]*ua*ua*ua*ui*VmE_K;
  const double I_Kr        = pCmP->CKr[Vi]*Xr*VmE_K;
  const double I_Ks        = (v(VT_g_Ks))*Xs*Xs*VmE_K;
  const double I_CaL       = (v(VT_g_CaL))*d*f*f_Ca*(V_int-(v(VT_E_rL)));
  const double I_K1        = pCmP->CK1[Vi]*VmE_K;
  const double I_bNa       = (v(VT_g_bNa))*VmE_Na;
  const double I_bCa       = (v(VT_g_bCa))*VmE_Ca;
  const double temp        = (v(VT_k_mNai))*dNa_i;
  const double I_NaK       = pCmP->CNaK[Vi]/(1.0+temp*sqrt(temp));
  const double I_pCa       = (v(VT_I_pCamax))*Ca_i/(Ca_i+(v(VT_k_mpCa)));
  const double I_NaCa      = pCmP->CNaCa[Vi]*((Na_i*Na_i*Na_i*(v(VT_Ca_o)))-pCmP->expVm[Vi]*Ca_i);
  const double I_rel       = (v(VT_k_rel))*u*u*v*w*(Ca_rel-Ca_i);
  const double I_up        = (v(VT_I_upmax))*Ca_i/(Ca_i+(v(VT_k_up)));
  const double I_upleak    = (v(VT_kupleak))*Ca_up;
  const double I_tr        = (Ca_up-Ca_rel)*(v(VT_dt_tr));
  const double Fn          = 1e-12*v(VT_Vrel)*I_rel - 5e-13*v(VT_C_m)*(0.5*I_CaL-0.2*I_NaCa)/v(VT_F);
#ifdef ISAC
  const double I_mem       = I_Na+I_CaL+I_pCa+I_K1+I_to+I_Kur+I_Kr+I_Ks+I_bNa+I_bCa+I_NaK+I_NaCa+I_SAC;
#else // ifdef ISAC
  const double I_mem       = I_Na+I_CaL+I_pCa+I_K1+I_to+I_Kur+I_Kr+I_Ks+I_bNa+I_bCa+I_NaK+I_NaCa;
#endif // ifdef ISAC
  tempstr << I_Na << ' ' << I_K1 << ' ' << I_to << ' ' << I_Kur << ' ' << I_Kr << ' ' << I_Ks << ' ' << I_CaL << ' ' <<
    I_pCa << ' ' <<
    I_NaK << ' ' << I_NaCa << ' ' << I_bNa << ' ' << I_bCa << ' ' << I_mem << ' ' << Fn << ' ' << I_rel << ' ' <<
    I_tr << ' ' << I_up << ' ' << I_upleak << ' '
#ifdef ISAC
    << I_SAC << ' '
#endif // ifdef ISAC
  ;
}  // Courtemanche::LongPrint

void Courtemanche::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] = {
    "m",  "h",  "j", "oa",  "oi",  "ua",  "ui",  "x_r",  "x_s",
    "d",  "f",  "fCa",  "u",  "v",  "w",  "Na_i",  "Ca_i",  "K_i",  "Ca_up",  "Ca_rel"
#ifdef TRPN
    , "Ca_TRPN"
#endif // ifdef TRPN
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
} // Courtemanche::GetParameterNames

void Courtemanche::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] = {
    "I_Na", "I_K1", "I_to", "I_Kur", "I_Kr", "I_Ks", "I_CaL", "I_pCa", "I_NaK", "I_NaCa", "I_bNa", "I_bCa", "I_ion",
    "Fn", "I_rel", "I_tr", "I_up", "I_upLeak"
#ifdef ISAC
    , "I_SAC"
#endif // ifdef ISAC
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
} // Courtemanche::GetLongParameterNames
