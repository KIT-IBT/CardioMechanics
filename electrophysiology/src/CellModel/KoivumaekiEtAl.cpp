/* File: KoivumaekiEtAl.cpp
        Institute of Biomedical Engineering, Karlsruhe Institute of Technology (KIT) */

#include <KoivumaekiEtAl.h>

KoivumaekiEtAl::KoivumaekiEtAl(KoivumaekiEtAlParameters *pp) {
  ptKeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptKeaP, NS_KoivumaekiEtAlParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

KoivumaekiEtAl::~KoivumaekiEtAl() {}

#ifdef HETERO

inline bool KoivumaekiEtAl::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool KoivumaekiEtAl::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int KoivumaekiEtAl::GetSize(void) {
  return (&Ca_SR4 - &m + 1) * sizeof(ML_CalcType);
}

inline unsigned char KoivumaekiEtAl::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void KoivumaekiEtAl::Init() {
#if KADEBUG
  cerr << "#initializing Class: KoivumaekiEtAl ... " << endl;
        #endif  // if KADEBUG
  m          = (v(VT_m_init));
  h1         = (v(VT_h1_init));
  h2         = (v(VT_h2_init));
  d          = (v(VT_d_init));
  f1         = (v(VT_f1_init));
  fca        = (v(VT_fca_init));
  r          = (v(VT_r_init));
  s          = (v(VT_s_init));
  sus_r      = (v(VT_sus_r_init));
  sus_s      = (v(VT_sus_s_init));
  n          = (v(VT_n_init));
  pa         = (v(VT_pa_init));
  y          = (v(VT_y_init));
  RyR_oss    = (v(VT_RyR_oss_init));
  RyR_css    = (v(VT_RyR_css_init));
  RyR_ass    = (v(VT_RyR_ass_init));
  RyR_o1     = (v(VT_RyR_o1_init));
  RyR_c1     = (v(VT_RyR_c1_init));
  RyR_a1     = (v(VT_RyR_a1_init));
  RyR_o2     = (v(VT_RyR_o2_init));
  RyR_c2     = (v(VT_RyR_c2_init));
  RyR_a2     = (v(VT_RyR_a2_init));
  RyR_o3     = (v(VT_RyR_o3_init));
  RyR_c3     = (v(VT_RyR_c3_init));
  RyR_a3     = (v(VT_RyR_a3_init));
  SERCA_Ca1  = (v(VT_SERCA_Ca1_init));
  SERCA_Ca2  = (v(VT_SERCA_Ca2_init));
  SERCA_Ca3  = (v(VT_SERCA_Ca3_init));
  SERCA_Cass = (v(VT_SERCA_Cass_init));
  Na_ss      = (v(VT_Na_ss_init));
  Na_i       = (v(VT_Na_i_init));
  K_i        = (v(VT_K_i_init));
  HTRPNCa1   = (v(VT_HTRPNCa1_init));
  HTRPNCa2   = (v(VT_HTRPNCa2_init));
  HTRPNCa3   = (v(VT_HTRPNCa3_init));
  HTRPNCa4   = (v(VT_HTRPNCa4_init));
  LTRPNCa1   = (v(VT_LTRPNCa1_init));
  LTRPNCa2   = (v(VT_LTRPNCa2_init));
  LTRPNCa3   = (v(VT_LTRPNCa3_init));
  LTRPNCa4   = (v(VT_LTRPNCa4_init));
  Ca_ss      = (v(VT_Ca_ss_init));
  Ca_i1      = (v(VT_Ca_i1_init));
  Ca_i2      = (v(VT_Ca_i2_init));
  Ca_i3      = (v(VT_Ca_i3_init));
  Ca_i4      = (v(VT_Ca_i4_init));
  Ca_SR1     = (v(VT_Ca_SR1_init));
  Ca_SR2     = (v(VT_Ca_SR2_init));
  Ca_SR3     = (v(VT_Ca_SR3_init));
  Ca_SR4     = (v(VT_Ca_SR4_init));
}  // KoivumaekiEtAl::Init

ML_CalcType KoivumaekiEtAl::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                 int euler) {
  ML_CalcType svolt = V*1000;
  ML_CalcType HT    = tinc;

  i_external = i_external*1000;  // convert from nA to pA

  const int Vi = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);

  const ML_CalcType RTdF = v(VT_RTdF);
  const ML_CalcType FdRT = 1.0/RTdF;
  const ML_CalcType Cm   = v(VT_Cm);

  const ML_CalcType K_o    = (v(VT_K_o));
  const ML_CalcType Ca_o   = (v(VT_Ca_o));
  const ML_CalcType E_Na   =  RTdF*log((v(VT_Na_o))/Na_ss);
  const ML_CalcType VmE_Na = (svolt - E_Na);
  const ML_CalcType E_K    =  RTdF*log(K_o/K_i);
  const ML_CalcType VmE_K  = (svolt - E_K);
  const ML_CalcType E_Ca   =  RTdF*0.5*log(Ca_o/Ca_ss);

  const ML_CalcType i_Na = (v(VT_PNa)) * m*m*m * (0.9*h1 + 0.1*h2) * (exp(VmE_Na*FdRT) - 1.0)*
    ptKeaP->CexpINa[Vi];
  const ML_CalcType i_Ca_L       =  (v(VT_gCaL))*d*(fca*f1)*(svolt - (v(VT_ECa_app)));
  const ML_CalcType fca_inf      = 1.0-1.0/(1.0+pow(v(VT_kCa)/Ca_ss, v(VT_kCan)));
  const ML_CalcType i_t          =  (v(VT_gt))*r*s*VmE_K;
  const ML_CalcType i_sus        =  (v(VT_gsus))*sus_r*sus_s*VmE_K;
  const ML_CalcType i_Ks         =  (v(VT_gKs))*n*VmE_K;
  const ML_CalcType i_Kr         =  (v(VT_gKr))*pa*ptKeaP->pip[Vi]*VmE_K;
  const ML_CalcType i_K1         = (v(VT_gK1))*VmE_K*(v(VT_powKo))/(1.00000+exp(1.50000*(VmE_K+3.60000)*FdRT));
  const ML_CalcType i_B_Na       =  (v(VT_gNab))*VmE_Na;
  const ML_CalcType i_B_Ca       =  (v(VT_gCab))*(svolt - E_Ca);
  const ML_CalcType Na_ss3       = Na_ss*Na_ss*Na_ss;
  const ML_CalcType pow_Na_ss_15 = sqrt(Na_ss3);
  const ML_CalcType i_NaK        =
    ( (( (( (v(VT_INaKmax))*K_o)/(K_o+(v(VT_kNaKK))))*pow_Na_ss_15)/(pow_Na_ss_15+(v(VT_pow_kNaKNa))))*(svolt+150.000))/
    (svolt+200.000);
  const ML_CalcType i_NaCa = ( (v(VT_kNaCa))*(Na_ss3*ptKeaP->CI_NaCa1[Vi] - Ca_ss*ptKeaP->CI_NaCa2[Vi]))/
    (1.00000+ (v(VT_dNaCa))*( (v(VT_Na_o3))*Ca_ss+ Na_ss3*Ca_o));
  const ML_CalcType i_CaP = ( (v(VT_ICaPmax))*Ca_ss)/(Ca_ss+(v(VT_kCaP)));
  const ML_CalcType i_fNa = (v(VT_gIf)) * y * 0.2677 * VmE_Na;
  const ML_CalcType i_fK  = (v(VT_gIf)) * y * (1-0.2677) * VmE_K;
  const ML_CalcType i_f   = i_fK + i_fNa;

  const ML_CalcType KdFura2         = (v(VT_KdFura2));
  const ML_CalcType totFura2KdFura2 = (v(VT_totFura2))*KdFura2;
  const ML_CalcType gamma1          = totFura2KdFura2/((Ca_i1 + KdFura2)*(Ca_i1 + KdFura2));
  const ML_CalcType gamma2          = totFura2KdFura2/((Ca_i2 + KdFura2)*(Ca_i2 + KdFura2));
  const ML_CalcType gamma3          = totFura2KdFura2/((Ca_i3 + KdFura2)*(Ca_i3 + KdFura2));
  const ML_CalcType gamma4          = totFura2KdFura2/((Ca_i4 + KdFura2)*(Ca_i4 + KdFura2));

  const ML_CalcType KdCMDN   = (v(VT_KdCMDN));
  const ML_CalcType KdSLlow  = (v(VT_KdSLlow));
  const ML_CalcType KdSLhigh = (v(VT_KdSLhigh));
  const ML_CalcType betass   = 1.0/
    (1.0 + (v(VT_totSLlow))*KdSLlow/((Ca_ss + KdSLlow)*(Ca_ss + KdSLlow)) + (v(VT_totSLhigh))*KdSLhigh/
     ((Ca_ss + KdSLhigh)*(Ca_ss + KdSLhigh)) + (v(VT_totCMDN))*KdCMDN/((Ca_ss + KdCMDN)*(Ca_ss + KdCMDN)) +
     (v(VT_totFura2))*KdFura2/((Ca_ss + KdFura2)*(Ca_ss + KdFura2)));

  const ML_CalcType totCMDNKdCMDN = (v(VT_totCMDN))*KdCMDN;
  const ML_CalcType beta1         = 1.0/ (1.0 + totCMDNKdCMDN/((Ca_i1 + KdCMDN)*(Ca_i1 + KdCMDN)) + gamma1);
  const ML_CalcType beta2         = 1.0/ (1.0 + totCMDNKdCMDN/((Ca_i2 + KdCMDN)*(Ca_i2 + KdCMDN)) + gamma2);
  const ML_CalcType beta3         = 1.0/ (1.0 + totCMDNKdCMDN/((Ca_i3 + KdCMDN)*(Ca_i3 + KdCMDN)) + gamma3);
  const ML_CalcType beta4         = 1.0/ (1.0 + totCMDNKdCMDN/((Ca_i4 + KdCMDN)*(Ca_i4 + KdCMDN)) + gamma4);

  const ML_CalcType KdCSQN        = (v(VT_KdCSQN));
  const ML_CalcType totCSQNKdCSQN = (v(VT_totCSQN))*KdCSQN;
  const ML_CalcType betaSR1       = 1.0/(1.0 + totCSQNKdCSQN/((Ca_SR1 + KdCSQN)*(Ca_SR1 + KdCSQN)));
  const ML_CalcType betaSR2       = 1.0/(1.0 + totCSQNKdCSQN/((Ca_SR2 + KdCSQN)*(Ca_SR2 + KdCSQN)));
  const ML_CalcType betaSR3       = 1.0/(1.0 + totCSQNKdCSQN/((Ca_SR3 + KdCSQN)*(Ca_SR3 + KdCSQN)));
  const ML_CalcType betaSR4       = 1.0/(1.0 + totCSQNKdCSQN/((Ca_SR4 + KdCSQN)*(Ca_SR4 + KdCSQN)));

  const ML_CalcType KdBNa    = (v(VT_KdBNa));
  const ML_CalcType betaNass = 1.0/(1.0 + (v(VT_totBNa))*KdBNa/((Na_ss + KdBNa)*(Na_ss + KdBNa)));

  const ML_CalcType Jj_nj = (v(VT_CJj_nj)) * (Ca_ss-Ca_i4);

  const ML_CalcType k1            = (v(VT_k1));
  const ML_CalcType k2            = (v(VT_k2));
  const ML_CalcType k3            = (v(VT_k3));
  const ML_CalcType k4            = (v(VT_k4));
  const ML_CalcType cpumps        = (v(VT_cpumps));
  const ML_CalcType V_nonjunct1   = (v(VT_V_nonjunct1));
  const ML_CalcType V_nonjunct2   = (v(VT_V_nonjunct2));
  const ML_CalcType V_nonjunct3   = (v(VT_V_nonjunct3));
  const ML_CalcType V_ss          = (v(VT_V_ss))*(v(VT_Dvcell));
  const ML_CalcType J_SERCASR1    = (-k3*Ca_SR1*Ca_SR1*(cpumps-SERCA_Ca1)+k4*SERCA_Ca1)*V_nonjunct1*2.0;
  const ML_CalcType J_bulkSERCA1  = (k1*Ca_i1*Ca_i1*(cpumps-SERCA_Ca1)-k2*SERCA_Ca1)*V_nonjunct1*2.0;
  const ML_CalcType J_SERCASR2    = (-k3*Ca_SR2*Ca_SR2*(cpumps-SERCA_Ca2)+k4*SERCA_Ca2)*V_nonjunct2*2.0;
  const ML_CalcType J_bulkSERCA2  = (k1*Ca_i2*Ca_i2*(cpumps-SERCA_Ca2)-k2*SERCA_Ca2)*V_nonjunct2*2.0;
  const ML_CalcType J_SERCASR3    = (-k3*Ca_SR3*Ca_SR3*(cpumps-SERCA_Ca3)+k4*SERCA_Ca3)*V_nonjunct3*2.0;
  const ML_CalcType J_bulkSERCA3  = (k1*Ca_i3*Ca_i3*(cpumps-SERCA_Ca3)-k2*SERCA_Ca3)*V_nonjunct3*2.0;
  const ML_CalcType J_SERCASRss   = (-k3*Ca_SR4*Ca_SR4*(cpumps-SERCA_Cass)+k4*SERCA_Cass)*V_ss*2.0;
  const ML_CalcType J_bulkSERCAss = (k1*Ca_ss*Ca_ss*(cpumps-SERCA_Cass)-k2*SERCA_Cass)*V_ss*2.0;

  const ML_CalcType RyRSRCass = (1.0 - 1.0/(1.0+exp((Ca_SR4-0.45)/0.15)));
  const ML_CalcType RyRainfss = 0.505-0.427/(1.0 + exp((Ca_ss*1000.0-0.29)/0.082));
  const ML_CalcType RyRoinfss = (1.0 - 1.0/(1.0 + exp((Ca_ss*1000.0-(RyR_ass + 0.22))/0.03)));
  const ML_CalcType RyRcinfss = (1.0/(1.0 + exp((Ca_ss*1000.0-(RyR_ass+0.02))/0.01)));
  const ML_CalcType Jrelss    = 935*V_ss/(v(VT_Dvcell)) * RyR_oss * RyR_css * RyRSRCass * (Ca_SR4 - Ca_ss);

  const ML_CalcType RyRSRCa1 = (1.0 - 1.0/(1.0 + exp((Ca_SR1-0.45)/0.15)));
  const ML_CalcType RyRainf1 = 0.505-0.427/(1.0 + exp((Ca_i1*1000.0-0.29)/0.082));
  const ML_CalcType RyRoinf1 = (1.0 - 1.0/(1.0 + exp((Ca_i1*1000.0-(RyR_a1 + 0.22))/0.03)));
  const ML_CalcType RyRcinf1 = (1.0/(1.0 + exp((Ca_i1*1000.0-(RyR_a1+0.02))/0.01)));
  const ML_CalcType Jrel1    = 2.1*V_nonjunct1/(v(VT_Dvcell)) * RyR_o1 * RyR_c1 * RyRSRCa1 * (Ca_SR1 - Ca_i1);

  const ML_CalcType RyRSRCa2 = (1.0 - 1.0/(1.0 + exp((Ca_SR2-0.45)/0.15)));
  const ML_CalcType RyRainf2 =  0.505-0.427/(1.0 + exp((Ca_i2*1000.0-0.29)/0.082));
  const ML_CalcType RyRoinf2 = (1.0 - 1.0/(1.0 + exp((Ca_i2*1000.0-(RyR_a2 + 0.22))/0.03)));
  const ML_CalcType RyRcinf2 = (1.0/(1.0 + exp((Ca_i2*1000.0-(RyR_a2+0.02))/0.01)));
  const ML_CalcType Jrel2    = 2.1*V_nonjunct2/(v(VT_Dvcell)) * RyR_o2 * RyR_c2 * RyRSRCa2 * (Ca_SR2 - Ca_i2);

  const ML_CalcType RyRSRCa3 = (1.0 - 1.0/(1.0 + exp((Ca_SR3-0.45)/0.15)));
  const ML_CalcType RyRainf3 =  0.505-0.427/(1.0 + exp((Ca_i3*1000.0-0.29)/0.082));
  const ML_CalcType RyRoinf3 = (1.0 - 1.0/(1.0 + exp((Ca_i3*1000.0-(RyR_a3 + 0.22))/0.03)));
  const ML_CalcType RyRcinf3 = (1.0/(1.0 + exp((Ca_i3*1000.0-(RyR_a3+0.02))/0.01)));
  const ML_CalcType Jrel3    = 2.1*V_nonjunct3/(v(VT_Dvcell)) * RyR_o3 * RyR_c3 * RyRSRCa3 * (Ca_SR3 - Ca_i3);

  const ML_CalcType kSRleak     = (v(VT_kSRleak));
  const ML_CalcType JSRCaleak1  = kSRleak * (Ca_SR1 - Ca_i1) * V_nonjunct1/(v(VT_Dvcell));
  const ML_CalcType JSRCaleak2  = kSRleak * (Ca_SR2 - Ca_i2) * V_nonjunct2/(v(VT_Dvcell));
  const ML_CalcType JSRCaleak3  = kSRleak * (Ca_SR3 - Ca_i3) * V_nonjunct3/(v(VT_Dvcell));
  const ML_CalcType JSRCaleakss = kSRleak * (Ca_SR4 - Ca_ss) * V_ss/(v(VT_Dvcell));

  const ML_CalcType JSR_Ca1 = J_SERCASR1 - JSRCaleak1 - Jrel1;
  const ML_CalcType JSR_Ca2 = J_SERCASR2 - JSRCaleak2 - Jrel2;
  const ML_CalcType JSR_Ca3 = J_SERCASR3 - JSRCaleak3 - Jrel3;
  const ML_CalcType JSR_Ca4 = J_SERCASRss - JSRCaleakss - Jrelss;

  const ML_CalcType JNa = (v(VT_CJNa)) * (Na_ss - Na_i);

  // calculating rates part
  const ML_CalcType m_inf     = ptKeaP->m_infinity[Vi];
  const ML_CalcType h_inf     = ptKeaP->h_infinity[Vi];
  const ML_CalcType d_inf     = ptKeaP->d_infinity[Vi];
  const ML_CalcType f_inf     = ptKeaP->f_infinity[Vi];
  const ML_CalcType r_inf     = ptKeaP->r_infinity[Vi];
  const ML_CalcType s_inf     = ptKeaP->s_infinity[Vi];
  const ML_CalcType sus_r_inf = ptKeaP->sus_r_infinity[Vi];
  const ML_CalcType sus_s_inf = ptKeaP->sus_s_infinity[Vi];
  const ML_CalcType n_inf     = ptKeaP->n_infinity[Vi];
  const ML_CalcType pa_inf    = ptKeaP->p_a_infinity[Vi];
  const ML_CalcType y_inf     = ptKeaP->y_infinity[Vi];

  m     = m_inf+(m-m_inf)*ptKeaP->exptau_m[Vi];
  h1    = h_inf+(h1-h_inf)*ptKeaP->exptau_h1[Vi];
  h2    = h_inf+(h2-h_inf)*ptKeaP->exptau_h2[Vi];
  d     = d_inf+(d-d_inf)*ptKeaP->exptau_d[Vi];
  f1    = f_inf+(f1-f_inf)*ptKeaP->exptau_f1[Vi];
  fca   = fca_inf+(fca-fca_inf)*(v(VT_exp_tau_fca));
  r     = r_inf+(r-r_inf)*ptKeaP->exptau_r[Vi];
  s     = s_inf+(s-s_inf)*ptKeaP->exptau_s[Vi];
  sus_r = sus_r_inf+(sus_r-sus_r_inf)*ptKeaP->exptau_sus_r[Vi];
  sus_s = sus_s_inf+(sus_s-sus_s_inf)*ptKeaP->exptau_sus_s[Vi];
  n     = n_inf+(n-n_inf)*ptKeaP->exptau_n[Vi];
  pa    = pa_inf+(pa-pa_inf)*ptKeaP->exptau_pa[Vi];
  y     = y_inf+(y-y_inf)*ptKeaP->exptau_y[Vi];


  const ML_CalcType HT05 = HT*0.5;
  SERCA_Ca1  += HT05*(-J_SERCASR1 + J_bulkSERCA1)/V_nonjunct1;
  SERCA_Ca2  += HT05*(-J_SERCASR2 + J_bulkSERCA2)/V_nonjunct2;
  SERCA_Ca3  += HT05*(-J_SERCASR3 + J_bulkSERCA3)/V_nonjunct3;
  SERCA_Cass += HT05*(-J_SERCASRss + J_bulkSERCAss)/V_ss;

  const ML_CalcType HTRyRtauadapt = HT/(v(VT_RyRtauadapt));
  const ML_CalcType HTRyRtauact   = HT/(v(VT_RyRtauact));
  const ML_CalcType HTRyRtauinact = HT/(v(VT_RyRtauinact));
  RyR_oss += HT*(RyRoinfss-RyR_oss)/(v(VT_RyRtauactss));
  RyR_css += HT*(RyRcinfss-RyR_css)/(v(VT_RyRtauinactss));
  RyR_ass += (RyRainfss-RyR_ass)*HTRyRtauadapt;
  RyR_o1  += (RyRoinf1-RyR_o1)*HTRyRtauact;
  RyR_c1  += (RyRcinf1-RyR_c1)*HTRyRtauinact;
  RyR_a1  += (RyRainf1-RyR_a1)*HTRyRtauadapt;
  RyR_o2  += (RyRoinf2-RyR_o2)*HTRyRtauact;
  RyR_c2  += (RyRcinf2-RyR_c2)*HTRyRtauinact;
  RyR_a2  += (RyRainf2-RyR_a2)*HTRyRtauadapt;
  RyR_o3  += (RyRoinf3-RyR_o3)*HTRyRtauact;
  RyR_c3  += (RyRcinf3-RyR_c3)*HTRyRtauinact;
  RyR_a3  += (RyRainf3-RyR_a3)*HTRyRtauadapt;

  const ML_CalcType Fara = (v(VT_F));
  Na_ss += HT*betaNass*(-JNa/V_ss-(i_Na+i_B_Na+ 3.00000*i_NaCa+ 3.00000*i_NaK+i_fNa)/(V_ss*Fara));
  Na_i  += HT*(JNa/(v(VT_V_nonjunct_Nai)));
  K_i   += HT*(-(i_t+i_sus+i_K1+i_Ks+i_Kr+i_fK - 2.0*i_NaK - i_external)/((v(VT_V_cytosol))*Fara));

  const ML_CalcType kHTRPNon  = (v(VT_kHTRPNon));
  const ML_CalcType kLTRPNon  = (v(VT_kLTRPNon));
  const ML_CalcType kHTRPNoff = (v(VT_kHTRPNoff));
  const ML_CalcType kLTRPNoff = (v(VT_kLTRPNoff));
  const ML_CalcType dHTRPNCa1 = kHTRPNon * Ca_i1 * (1.0 - HTRPNCa1) - kHTRPNoff * HTRPNCa1;
  const ML_CalcType dHTRPNCa2 = kHTRPNon * Ca_i2 * (1.0 - HTRPNCa2) - kHTRPNoff * HTRPNCa2;
  const ML_CalcType dHTRPNCa3 = kHTRPNon * Ca_i3 * (1.0 - HTRPNCa3) - kHTRPNoff * HTRPNCa3;
  const ML_CalcType dHTRPNCa4 = kHTRPNon * Ca_i4 * (1.0 - HTRPNCa4) - kHTRPNoff * HTRPNCa4;
  const ML_CalcType dLTRPNCa1 = kLTRPNon * Ca_i1 * (1.0 - LTRPNCa1) - kLTRPNoff * LTRPNCa1;
  const ML_CalcType dLTRPNCa2 = kLTRPNon * Ca_i2 * (1.0 - LTRPNCa2) - kLTRPNoff * LTRPNCa2;
  const ML_CalcType dLTRPNCa3 = kLTRPNon * Ca_i3 * (1.0 - LTRPNCa3) - kLTRPNoff * LTRPNCa3;
  const ML_CalcType dLTRPNCa4 = kLTRPNon * Ca_i4 * (1.0 - LTRPNCa4) - kLTRPNoff * LTRPNCa4;

  HTRPNCa1 += HT*dHTRPNCa1;
  HTRPNCa2 += HT*dHTRPNCa2;
  HTRPNCa3 += HT*dHTRPNCa3;
  HTRPNCa4 += HT*dHTRPNCa4;
  LTRPNCa1 += HT*dLTRPNCa1;
  LTRPNCa2 += HT*dLTRPNCa2;
  LTRPNCa3 += HT*dLTRPNCa3;
  LTRPNCa4 += HT*dLTRPNCa4;

  const ML_CalcType totHTRPN = (v(VT_totHTRPN));
  const ML_CalcType totLTRPN = (v(VT_totLTRPN));
  const ML_CalcType Jtrpn1   = (totHTRPN*dHTRPNCa1 + totLTRPN*dLTRPNCa1) * V_nonjunct1;
  const ML_CalcType Jtrpn2   = (totHTRPN*dHTRPNCa2 + totLTRPN*dLTRPNCa2) * V_nonjunct2;
  const ML_CalcType Jtrpn3   = (totHTRPN*dHTRPNCa3 + totLTRPN*dLTRPNCa3) * V_nonjunct3;
  const ML_CalcType Jtrpn4   = (totHTRPN*dHTRPNCa4 + totLTRPN*dLTRPNCa4) * (v(VT_V_nonjunct4));

  const ML_CalcType JCa1  = -J_bulkSERCA1 + JSRCaleak1 + Jrel1 - Jtrpn1;
  const ML_CalcType JCa2  = -J_bulkSERCA2 + JSRCaleak2 + Jrel2 - Jtrpn2;
  const ML_CalcType JCa3  = -J_bulkSERCA3 + JSRCaleak3 + Jrel3 - Jtrpn3;
  const ML_CalcType JCa4  = Jj_nj - Jtrpn4;
  const ML_CalcType JCass = -Jj_nj + JSRCaleakss - J_bulkSERCAss + Jrelss;

  Ca_ss += HT*(betass * (JCass/V_ss + (-i_Ca_L - i_B_Ca - i_CaP+ 2.0*i_NaCa) / (2.0*V_ss*Fara) ));

  const ML_CalcType dx     = (v(VT_dx));
  const ML_CalcType d2dx   = 1.0/(2.0*dx);
  const ML_CalcType dxdx   = dx*dx;
  const ML_CalcType d2dxdx = 1.0/(2.0*dxdx);
  const ML_CalcType DCa    = (v(VT_DCa));
  const ML_CalcType DCaBm  = (v(VT_DCaBm));
  Ca_i1 += HT*
    (beta1 * (DCa + gamma1 * DCaBm) * ( (Ca_i2-2.0*Ca_i1+Ca_i1)/dxdx + (Ca_i2-Ca_i1)*d2dxdx) - 2.0*beta1*gamma1*DCaBm/
     (KdCMDN + Ca_i1) * ((Ca_i2-Ca_i1)*d2dx)*((Ca_i2-Ca_i1)*d2dx) + JCa1/V_nonjunct1*beta1);
  Ca_i2 += HT*
    (beta2 * (DCa + gamma2 * DCaBm) * ( (Ca_i3-2.0*Ca_i2+Ca_i1)/dxdx + (Ca_i3-Ca_i1)*d2dxdx/2.0) - 2.0*beta2*gamma2*
     DCaBm/
     (KdCMDN + Ca_i2) * ((Ca_i3-Ca_i1)*d2dx)*((Ca_i3-Ca_i1)*d2dx) + JCa2/V_nonjunct2*beta2);
  Ca_i3 += HT*
    (beta3 * (DCa + gamma3 * DCaBm) * ( (Ca_i4-2.0*Ca_i3+Ca_i2)/dxdx + (Ca_i4-Ca_i2)*d2dxdx/3.0) - 2.0*beta3*gamma3*
     DCaBm/
     (KdCMDN + Ca_i3) * ((Ca_i4-Ca_i2)*d2dx)*((Ca_i4-Ca_i2)*d2dx) + JCa3/V_nonjunct3*beta3);
  Ca_i4 += HT*
    (beta4 * (DCa + gamma4 * DCaBm) * ( (Ca_i4-2.0*Ca_i4+Ca_i3)/dxdx + (Ca_i4-Ca_i3)*d2dxdx/4.0) - 2.0*beta4*gamma4*
     DCaBm/
     (KdCMDN + Ca_i4) * ((Ca_i4-Ca_i3)*d2dx)*((Ca_i4-Ca_i3)*d2dx) + JCa4/(v(VT_V_nonjunct4))*beta4);

  const ML_CalcType DCaSR = (v(VT_DCaSR));
  Ca_SR1 += HT*
    (betaSR1 * DCaSR * ((Ca_SR2-2.0*Ca_SR1+Ca_SR1)/dxdx + (Ca_SR2-Ca_SR1)*d2dxdx) + JSR_Ca1/(v(VT_V_SR1))*betaSR1);
  Ca_SR2 += HT*
    (betaSR2 * DCaSR * ((Ca_SR3-2.0*Ca_SR2+Ca_SR1)/dxdx + (Ca_SR3-Ca_SR1)*d2dxdx/2.0) + JSR_Ca2/(v(VT_V_SR2))*betaSR2);
  Ca_SR3 += HT*
    (betaSR3 * DCaSR * ((Ca_SR4-2.0*Ca_SR3+Ca_SR2)/dxdx + (Ca_SR4-Ca_SR2)*d2dxdx/3.0) + JSR_Ca3/(v(VT_V_SR3))*betaSR3);
  Ca_SR4 += HT*
    (betaSR4 * DCaSR * ((Ca_SR4-2.0*Ca_SR4+Ca_SR3)/dxdx + (Ca_SR4-Ca_SR3)*d2dxdx/4.0) + JSR_Ca4/(v(VT_V_SR4))*betaSR4);

  const ML_CalcType I_tot =
    (-i_external + i_Na + i_Ca_L + i_t + i_sus + i_K1 + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP + i_NaCa + i_f)/
    Cm;
  return tinc*(-I_tot);
}  // KoivumaekiEtAl::Calc

void KoivumaekiEtAl::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<<tArg<< ' ' << V<<' ' <<m<<' ' <<h1<<' ' <<h2<<' ' <<d<<' ' <<f1<<' ' <<fca
         <<' ' <<r<<' ' <<s<<' ' <<sus_r<<' ' <<sus_s<<' ' <<n<<' ' <<pa<<' ' <<y
         <<' ' <<RyR_oss<<' ' <<RyR_css<<' ' <<RyR_ass<<' ' <<RyR_o1<<' ' <<RyR_c1<<' ' <<RyR_a1
         <<' ' <<RyR_o2<<' ' <<RyR_c2<<' ' <<RyR_a2<<' ' <<RyR_o3<<' ' <<RyR_c3<<' ' <<RyR_a3
         <<' ' <<SERCA_Ca1<<' ' <<SERCA_Ca2<<' ' <<SERCA_Ca3<<' ' <<SERCA_Cass
         <<' ' <<Na_ss<<' ' <<Na_i<<' ' <<K_i
         <<' ' <<HTRPNCa1<<' ' <<HTRPNCa2<<' ' <<HTRPNCa3<<' ' <<HTRPNCa4
         <<' ' <<LTRPNCa1<<' ' <<LTRPNCa2<<' ' <<LTRPNCa3<<' ' <<LTRPNCa4
         <<' ' <<Ca_ss<<' ' <<Ca_i1<<' ' <<Ca_i2<<' ' <<Ca_i3<<' ' <<Ca_i4
         <<' ' <<Ca_SR1<<' ' <<Ca_SR2<<' ' <<Ca_SR3<<' ' <<Ca_SR4;
}

void KoivumaekiEtAl::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
  const  ML_CalcType svolt = V*1000.0;
  const int Vi             = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);
  const ML_CalcType RTdF   = v(VT_RTdF);
  const ML_CalcType FdRT   = 1.0/RTdF;
  const ML_CalcType Cm     = v(VT_Cm);

  const ML_CalcType K_o    = (v(VT_K_o));
  const ML_CalcType Ca_o   = (v(VT_Ca_o));
  const ML_CalcType E_Na   =  RTdF*log((v(VT_Na_o))/Na_ss);
  const ML_CalcType VmE_Na = (svolt - E_Na);
  const ML_CalcType E_K    =  RTdF*log(K_o/K_i);
  const ML_CalcType VmE_K  = (svolt - E_K);
  const ML_CalcType E_Ca   =  RTdF*0.5*log(Ca_o/Ca_ss);

  const ML_CalcType i_Na = (v(VT_PNa)) * m*m*m * (0.9*h1 + 0.1*h2) * (exp(VmE_Na*FdRT) - 1.0)*
    ptKeaP->CexpINa[Vi];
  const ML_CalcType i_Ca_L       =  (v(VT_gCaL))*d*(fca*f1)*(svolt - (v(VT_ECa_app)));
  const ML_CalcType fca_inf      = 1.0-1.0/(1.0+pow(v(VT_kCa)/Ca_ss, v(VT_kCan)));
  const ML_CalcType i_t          =  (v(VT_gt))*r*s*VmE_K;
  const ML_CalcType i_sus        =  (v(VT_gsus))*sus_r*sus_s*VmE_K;
  const ML_CalcType i_Ks         =  (v(VT_gKs))*n*VmE_K;
  const ML_CalcType i_Kr         =  (v(VT_gKr))*pa*ptKeaP->pip[Vi]*VmE_K;
  const ML_CalcType i_K1         = (v(VT_gK1))*VmE_K*(v(VT_powKo))/(1.00000+exp(1.50000*(VmE_K+3.60000)*FdRT));
  const ML_CalcType i_B_Na       =  (v(VT_gNab))*VmE_Na;
  const ML_CalcType i_B_Ca       =  (v(VT_gCab))*(svolt - E_Ca);
  const ML_CalcType Na_ss3       = Na_ss*Na_ss*Na_ss;
  const ML_CalcType pow_Na_ss_15 = sqrt(Na_ss3);
  const ML_CalcType i_NaK        =
    ( (( (( (v(VT_INaKmax))*K_o)/(K_o+(v(VT_kNaKK))))*pow_Na_ss_15)/(pow_Na_ss_15+(v(VT_pow_kNaKNa))))*(svolt+150.000))/
    (svolt+200.000);
  const ML_CalcType i_NaCa = ( (v(VT_kNaCa))*(Na_ss3*ptKeaP->CI_NaCa1[Vi] - Ca_ss*ptKeaP->CI_NaCa2[Vi]))/
    (1.00000+ (v(VT_dNaCa))*( (v(VT_Na_o3))*Ca_ss+ Na_ss3*Ca_o));
  const ML_CalcType i_CaP = ( (v(VT_ICaPmax))*Ca_ss)/(Ca_ss+(v(VT_kCaP)));
  const ML_CalcType i_fNa = (v(VT_gIf)) * y * 0.2677 * VmE_Na;
  const ML_CalcType i_fK  = (v(VT_gIf)) * y * (1-0.2677) * VmE_K;
  const ML_CalcType i_f   = i_fK + i_fNa;
  const ML_CalcType I_mem = i_Na + i_Ca_L + i_t + i_sus + i_K1 + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP +
    i_NaCa + i_f;
  tempstr<<' '<<i_Na << ' '<<i_Ca_L << ' '<<i_t << ' '<<i_sus << ' '<<i_Ks << ' '<<i_Kr << ' '<<i_K1 << ' '<<i_B_Na <<
    ' '
         <<i_B_Ca << ' '<<i_NaK << ' '<<i_NaCa << ' '<<i_CaP << ' '<<i_fNa << ' '<<i_fK << ' '<<i_f<< ' '<< I_mem <<
    ' ';
}  // KoivumaekiEtAl::LongPrint

void KoivumaekiEtAl::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 49;
  const string ParaNames[numpara] =
  {"m",         "h1",            "h2",             "d",              "f1",
   "fca",
   "r",
   "s",         "sus_r",         "sus_s",
   "n",         "pa",            "y",
   "RyR_oss",   "RyR_css",       "RyR_ass",        "RyR_o1",         "RyR_c1",
   "RyR_a1",
   "RyR_o2",
   "RyR_c2",    "RyR_a2",        "RyR_o3",         "RyR_c3",         "RyR_a3",
   "SERCA_Ca1", "SERCA_Ca2",     "SERCA_Ca3",      "SERCA_Ca4",      "Na_ss",
   "Na_i",
   "K_i",
   "HTRPNCa1",  "HTRPNCa2",      "HTRPNCa3",       "HTRPNCa4",
   "LTRPNCa1",  "LTRPNCa2",      "LTRPNCa3",       "LTRPNCa4",       "Ca_ss",
   "Ca_i1",
   "Ca_i2",
   "Ca_i3",     "Ca_i4",         "Ca_SR1",         "Ca_SR2",         "Ca_SR3",
   "Ca_SR4"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void KoivumaekiEtAl::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 16;
  const string ParaNames[numpara] =
  {"i_Na", "i_Ca_L", "i_t", "i_sus", "i_Ks", "i_Kr", "i_K1", "i_B_Na", "i_B_Ca", "i_NaK", "i_NaCa", "i_CaP", "i_fNa",
   "i_fK", "i_f",    "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
