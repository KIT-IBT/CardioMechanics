/*
 * File: FabbriEtAl.cpp
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


#include <FabbriEtAl.h>

Fabbri::Fabbri(FabbriParameters *pp) {
  pCmP = pp;
  Init();
}

Fabbri::~Fabbri() {}

inline int Fabbri::GetSize(void) {return sizeof(Fabbri)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(FabbriParameters *);}

void Fabbri::Init() {
#if KADEBUG
  cerr << "Fabbri::Init" << endl;
#endif // if KADEBUG

  // Elektrolyte
  Na_i = v(VT_Init_Nai);
  K_i = v(VT_Init_Ki);
  Ca_i = v(VT_Init_Cai);
  Mg_i = v(VT_Init_Mgi);


  // Gates
  y = v(VT_Init_y);
  m = v(VT_Init_m);
  m_mut = v(VT_Init_m);
  h = v(VT_Init_h);
  h_mut = v(VT_Init_h);
  dL = v(VT_Init_dL);
  fL = v(VT_Init_fL);
  fCa = v(VT_Init_fCa);
  dT = v(VT_Init_dT);
  fT = v(VT_Init_fT);
  r_Kur = v(VT_Init_r_Kur);
  s = v(VT_Init_s);
  r_to = v(VT_Init_r_to);
  q = v(VT_Init_q);
  paS = v(VT_Init_paS);
  paF = v(VT_Init_paF);
  piy = v(VT_Init_piy);
  n = v(VT_Init_n);
  a = v(VT_Init_a);
  x = v(VT_Init_x);


  // SR-Komponenten
  R_Ca_rel = v(VT_Init_R_Ca_rel);
  O_Ca_rel = v(VT_Init_O_Ca_rel);
  I_Ca_rel = v(VT_Init_I_Ca_rel);
  RI_Ca_rel = v(VT_Init_RI_Ca_rel);
  Ca_jsr = v(VT_Init_Ca_jsr);
  Ca_nsr = v(VT_Init_Ca_nsr);
  Ca_sub = v(VT_Init_Ca_sub);
  fTMM = v(VT_Init_fTMM);
  fCMi = v(VT_Init_fCMi);
  fCMs = v(VT_Init_fCMs);
  fTC = v(VT_Init_fTC);
  fTMC = v(VT_Init_fTMC);
  fCQ = v(VT_Init_fCQ);
} // Fabbri::Init

ML_CalcType Fabbri::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                         int euler = 2) {
  // tinc*=1000.0;
  ML_CalcType V_int = V*1000.0;

  // (convert from nA to pA/pF with Cm: 57 pF,
  //  i_external=i_external/(2e-2*Volume()*0.2e6*1e9);
  // hard-coded, because S is currently questionable in TenTusscher2.h)
  // i_external=i_external/57; // convert from nA to pA/pF

  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);

  // Berechnung von RT/F

  const double RTdF = ((v(VT_R)*(v(VT_T)))/(v(VT_F)));

  // V_Nernst
  const double E_Na = RTdF*log((v(VT_Nao))/Na_i);
  const double E_K = RTdF*log((v(VT_Ko))/K_i);
  const double E_Ks = RTdF*log(((v(VT_Ko))+0.12*(v(VT_Nao)))/(K_i+0.12*Na_i));
  const double E_Ca = 0.5*RTdF*log((v(VT_Cao))/Ca_sub);
  const double E_mh = RTdF*log(((v(VT_Nao))+0.12*(v(VT_Ko)))/(Na_i+0.12*K_i));


  // (V-V_Nernst)
  const double VmE_Na = V_int-E_Na;
  const double VmE_K = V_int-E_K;
  const double VmE_Ca = V_int-E_Ca;

  // StrÃ¶me
  //
  // I_f:

  const double m_y = pCmP->m_y[Vi];

  y = m_y + (y - m_y) * pCmP->exptau_y[Vi];
  const double I_fNa = y*(v(VT_g_fNa))*VmE_Na;
  const double I_fK = y*(v(VT_g_fK))*VmE_K;
  const double I_f = I_fK+I_fNa;

  // INaK
  const double I_NaK =
    ((v(VT_i_NaK_max))*
     pow((1+pow(((v(VT_Km_Kp))/(v(VT_Ko))), 1.2)),
         -1)*pow((1+pow(((v(VT_Km_Nap))/Na_i), 1.3)), -1)*pow((1+exp(-(VmE_Na+110)/20)), -1))*v(VT_Iso_INak_increase);

  // INaCa
  const double di = 1+ (Ca_sub/(v(VT_Kci)))*(1+exp(((-1*(v(VT_Qci)))*V_int)/RTdF)+(Na_i/(v(VT_Kcni))))+
    (Na_i/(v(VT_K1ni)))*(1+(Na_i/(v(VT_K2ni)))*(1+(Na_i/(v(VT_K3ni)))));
  const double k43 = Na_i/((v(VT_K3ni))+Na_i);
  const double k12 = ((Ca_sub/(v(VT_Kci)))*exp((-(v(VT_Qci))*V_int)/RTdF))/di;
  const double k14 = ((((Na_i/(v(VT_K1ni)))*Na_i)/(v(VT_K2ni)))*(1+(Na_i/(v(VT_K3ni))))*exp((v(VT_Qn)*V_int)/(2*RTdF)))/
    di;
  const double k41 = pCmP->k41[Vi];
  const double k34 = (v(VT_Nao))/((v(VT_K3no))+(v(VT_Nao)));
  const double k21 = pCmP->k21[Vi];
  const double k23 = pCmP->k23[Vi];
  const double k32 = pCmP->k32[Vi];
  const double x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
  const double x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
  const double x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
  const double x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
  const double I_NaCa = ((v(VT_K_NaCa))*(x2*k21-(x1*k12)))/(x1+x2+x3+x4);
  
  // INa WT
  const double m_m = pCmP->m_m[Vi];
  m = m_m + (m-m_m) *pCmP->exptau_m[Vi];
  const double m_h = pCmP->m_h[Vi];
  h = m_h + (h-m_h) *pCmP->exptau_h[Vi];
  const double I_Na_WT = (v(VT_g_na))*m*m*m*h*(V_int-E_mh);

  // INa Mut
  const double m_m_mut = pCmP->m_m_mut[Vi];
  m_mut = m_m_mut + (m_mut-m_m_mut) *pCmP->exptau_m_mut[Vi];
  const double m_h_mut = pCmP->m_h_mut[Vi];
  h_mut = m_h_mut + (h_mut-m_h_mut) *pCmP->exptau_h_mut[Vi];
  const double I_Na_Mut =
    (v(VT_g_na_mut)*m_mut*m_mut*m_mut*h_mut*(V_int-E_mh) + v(VT_g_na_l)*m_mut*m_mut*m_mut*(V_int-E_mh));

  // INa
  const double I_Na = 0.5 * I_Na_WT + 0.5 * I_Na_Mut; // In case of no mutations: I_Na_Mut = I_Na_WT

  // ICaL
  const double m_dL = pCmP->m_dL[Vi];
  dL = m_dL + (dL-m_dL) *pCmP->exptau_dL[Vi];

  const double m_fL = pCmP->m_fL[Vi];
  fL = m_fL + (fL-m_fL) *pCmP->exptau_fL[Vi];

  const double m_fCa = (v(VT_Km_fCa))/((v(VT_Km_fCa))+Ca_sub);
  const double exptau_fCa = exp(-tinc*((v(VT_alpha_fCa))/(0.001*m_fCa)));
  fCa = m_fCa + (fCa-m_fCa)*exptau_fCa;

  double Ca1 = RTdF*(1-exp((-2*V_int)/RTdF));
  double Ca2 = exp((-2*V_int)/RTdF);
  double CaL3 = RTdF*(1-exp((-V_int)/RTdF));
  double CaL4 = exp((-V_int)/RTdF);

  const double I_siCa = ((2*(v(VT_P_CaL))*V_int)/Ca1)*(Ca_sub-(v(VT_Cao))*Ca2)*dL*fL*fCa;
  const double I_siK = ((0.000365*(v(VT_P_CaL))*V_int)/CaL3)*(K_i-(v(VT_Ko))*CaL4)*dL*fL*fCa;
  const double I_siNa = ((0.0000185*(v(VT_P_CaL))*V_int)/CaL3)*(Na_i-(v(VT_Nao))*CaL4)*dL*fL*fCa;
  const double I_CaL = (I_siCa+I_siK+I_siNa)*v(VT_Iso_ICaL_increase)*(1-v(VT_ACh_iCaL_block));


  // Natrium Dynamik
  if ((v(VT_Vari_Nai)) == 1) {
    double d_Nai = tinc*(((-1)*(I_Na+I_fNa+I_siNa+3*I_NaK+3*I_NaCa))/((((v(VT_V_i))+(v(VT_V_sub))))*(v(VT_F))));
    Na_i += d_Nai;
  } else {
    double d_Nai = 0;
    Na_i += d_Nai;
  }


  // ICaT
  const double m_dT = pCmP->m_dT[Vi];
  dT = m_dT + (dT-m_dT) *pCmP->exptau_dT[Vi];
  const double m_fT = pCmP->m_fT[Vi];
  fT = m_fT + (fT-m_fT) *pCmP->exptau_fT[Vi];

  const double I_CaT = ((2*(v(VT_P_CaT))*V_int)/Ca1)*(Ca_sub-(v(VT_Cao))*Ca2)*dT*fT;

  // IKur
  const double m_r_Kur = pCmP->m_r_Kur[Vi];
  r_Kur = m_r_Kur + (r_Kur-m_r_Kur) *pCmP->exptau_r_Kur[Vi];
  const double m_s = pCmP->m_s[Vi];
  s = m_s + (s-m_s) *pCmP->exptau_s[Vi];
  const double I_Kur = (v(VT_g_Kur))*r_Kur*s*VmE_K;

  // Ito
  const double m_r_to = pCmP->m_r_to[Vi];
  r_to = m_r_to + (r_to-m_r_to) *pCmP->exptau_r_to[Vi];
  const double m_q = pCmP->m_q[Vi];
  q = m_q + (q-m_q) *pCmP->exptau_q[Vi];
  const double I_to = v(VT_g_to)*VmE_K*q*r_to;


  // IKr
  const double m_pa = pCmP->m_pa[Vi];
  paS = m_pa + (paS-m_pa) *pCmP->exptau_paS[Vi];
  paF = m_pa + (paF-m_pa) *pCmP->exptau_paF[Vi];
  const double m_pi = pCmP->m_pi[Vi];
  piy = m_pi + (piy-m_pi) *pCmP->exptau_pi[Vi];

  const double I_Kr = (v(VT_g_Kr))*sqrt(v(VT_Ko)/5.4)*VmE_K*(0.9*paF+0.1*paS)*piy;


  // IKs
  const double m_n = pCmP->m_n[Vi];
  n = m_n + (n-m_n) *pCmP->exptau_n[Vi];

  double R231C = 0;
  if (v(VT_R231C) == 1) {
    double A2      = 0.07;
    double g_ratio =  0.3153;
    double V_WT    = 102.4597;
    double k_WT    = 55.2060;
    double V_R231C = 113.8797;
    double k_R231C = -34.9006;

    R231C = g_ratio*pow((1/(1+exp(-(V_int-V_R231C)/k_R231C))), 2)/pow(((1-A2)/(1+exp(-(V_int-V_WT)/k_WT))+A2), 2);
  } else {
    R231C = 1;
  }
  const double I_Ks = R231C * v(VT_g_Ks)*v(VT_Iso_IKs_increase)*(V_int-E_Ks)*n*n;

  // I_SK
  const double m_x =
    (0.81*pow((1000*Ca_sub), v(VT_n_SK)))/(pow((1000*Ca_sub), v(VT_n_SK))+pow(v(VT_EC50_SK), v(VT_n_SK)));
  const double tau_x = 0.001/(0.047*(1000*Ca_sub)+1.0/76.0);
  const double exptau_x = exp(-tinc/tau_x);
  x = m_x +(x-m_x)*exptau_x;

  const double I_SK = v(VT_g_SK) * VmE_K * x;

  // IKACh
  const double m_a = pCmP->m_a[Vi];
  a = m_a + (a-m_a) *pCmP->exptau_a[Vi];
  double I_KACh = 0;
  if ((v(VT_ACh)) > 0) {
    I_KACh = (v(VT_g_KACh))*VmE_K*(1+exp((V_int+20)/20))*a;
  } else {
    I_KACh = 0;
  }

  // Kalium Dynamik
  if ((v(VT_Vari_Ki)) == 1) {
    double d_Ki = tinc * (-1) *((I_Kur+I_to+I_Kr+I_Ks+I_fK+I_siK+I_SK-2*I_NaK)/((v(VT_V_i)+v(VT_V_sub))*v(VT_F)));
    K_i += d_Ki;
  } else {
    double d_Ki = 0;
    K_i += d_Ki;
  }

  // SR_Release
  const double j_SRCarel = (v(VT_ks))*O_Ca_rel*(Ca_jsr-Ca_sub);
  const double diff = Ca_jsr-Ca_sub;
  double MaxSR = (v(VT_MaxSR));
  double EC50 = (v(VT_EC50_SR));
  double HSR = (v(VT_HSR));
  double Bruch = (v(VT_EC50_SR))/Ca_jsr;
  double downk_CaSR = (1+pow(Bruch, (v(VT_HSR))));
  const double k_CaSR = (v(VT_MaxSR))-(((v(VT_MaxSR))-(v(VT_MinSR)))/(1+pow((v(VT_EC50_SR))/Ca_jsr, (v(VT_HSR)))));
  const double koSRCa = (v(VT_koCa))/(k_CaSR);
  const double kiSRCa = (v(VT_kiCa))*k_CaSR;
  R_Ca_rel += tinc*((v(VT_kim))*RI_Ca_rel-kiSRCa*Ca_sub*R_Ca_rel-(koSRCa*Ca_sub*Ca_sub*R_Ca_rel-(v(VT_kom))*O_Ca_rel));
  O_Ca_rel += tinc*(koSRCa*Ca_sub*Ca_sub*R_Ca_rel-(v(VT_kom))*O_Ca_rel-(kiSRCa*Ca_sub*O_Ca_rel-(v(VT_kim))*I_Ca_rel));
  I_Ca_rel += tinc*(kiSRCa*Ca_sub*O_Ca_rel-(v(VT_kim))*I_Ca_rel-((v(VT_kom))*I_Ca_rel-koSRCa*Ca_sub*Ca_sub*RI_Ca_rel));
  RI_Ca_rel += tinc*
    ((v(VT_kom))*I_Ca_rel-koSRCa*Ca_sub*Ca_sub*RI_Ca_rel-((v(VT_kim))*RI_Ca_rel-kiSRCa*Ca_sub*R_Ca_rel));

  // IntrazellulÃ¤re Ca_StrÃ¶me:
  double b_up = (0.7*v(VT_ACh)/(0.00009+v(VT_ACh)))-v(VT_Iso_b_up);

  const double P_up = 5; // in mM/s
  const double j_Ca_dif = (Ca_sub-Ca_i)/(v(VT_tau_dif_Ca));
  const double j_up = P_up*(1-b_up)/(1+exp(((-1*Ca_i)+(v(VT_K_up)))/(v(VT_slope_up))));
  const double j_tr = (Ca_nsr-Ca_jsr)/(v(VT_tau_tr));

  // Ca-Buffering:
  const double delta_fTC = (v(VT_kf_TC))*Ca_i*(1-fTC)-((v(VT_kb_TC))*fTC);
  const double delta_fTMC = (v(VT_kf_TMC))*Ca_i*(1-(fTMC+fTMM))-((v(VT_kb_TMC))*fTMC);
  const double delta_fTMM = (v(VT_kf_TMM))*Mg_i*(1-(fTMC+fTMM))-((v(VT_kb_TMM))*fTMM);
  const double delta_fCMi = (v(VT_kf_CM))*Ca_i*(1-fCMi)-((v(VT_kb_CM))*fCMi);
  const double delta_fCMs = (v(VT_kf_CM))*Ca_sub*(1-fCMs)-((v(VT_kb_CM))*fCMs);
  const double delta_fCQ = (v(VT_kf_CQ))*Ca_jsr*(1-fCQ)-((v(VT_kb_CQ))*fCQ);
  fTC += tinc*delta_fTC;
  fTMC += tinc*delta_fTMC;
  fTMM += tinc*delta_fTMM;
  fCMi += tinc*delta_fCMi;
  fCMs += tinc*delta_fCMs;
  fCQ += tinc*delta_fCQ;

  // Ca-Dynamik:
  double V_nsr = (v(VT_V_nsr));
  double V_i = (v(VT_V_i));
  double CM_tot = (v(VT_CM_tot));
  double TC_tot = (v(VT_TC_tot));
  double TMC_tot = (v(VT_TMC_tot));
  double V_jsr = (v(VT_V_jsr));
  Ca_i += tinc *
    (((1*(j_Ca_dif*(v(VT_V_sub))-(j_up*V_nsr)))/V_i)-(CM_tot*delta_fCMi+TC_tot*delta_fTC+TMC_tot*delta_fTMC));
  Ca_sub += tinc*
    (((j_SRCarel*V_jsr)/(v(VT_V_sub)))-
     (((I_siCa+I_CaT-(2*I_NaCa))/(2*(v(VT_F))*(v(VT_V_sub))))+j_Ca_dif+CM_tot*delta_fCMs));
  Ca_nsr += tinc*(j_up-((j_tr*V_jsr)/V_nsr));
  Ca_jsr += tinc*(j_tr-(j_SRCarel+(v(VT_CQ_tot))*delta_fCQ));


  double I_tot = (I_f+I_Kr+I_Ks+I_to+I_NaK+I_NaCa+I_Na+I_CaL+I_CaT+I_KACh+I_Kur+I_SK);

  // convert I_tot from nA to pA/pF (as expected from the outsided)
  I_tot *= (1000.0/(v(VT_C_m)));

  // add external current (defined in pA/pF)
  I_tot -= i_external;
  return -I_tot*tinc;
} // Fabbri::Calc

void Fabbri::Print(ostream &tempstr, double t,  ML_CalcType V) {
  tempstr<<t<<' '<<V<<' '
         <<Na_i<<' '<<K_i<<' '<<Ca_i<<' '
         <<Mg_i<<' '<<y<<' '<<m<<' '
         <<h<<' '<<dL<<' '<<fL<<' '
         <<fCa<<' '<<dT<<' '<<fT<<' '
         <<r_Kur<<' '<<s<<' '<<r_to<<' '
         <<q<<' '<<paS<<' '<<paF<<' '
         <<piy<<' '<<n<<' '<<a<<' '<<x<<' '
         <<R_Ca_rel<<' '<<O_Ca_rel<<' '<<I_Ca_rel<<' '<<RI_Ca_rel<<' '
         <<Ca_jsr<<' '<<Ca_nsr<<' '<<Ca_sub<<' '
         <<fTMM<<' '<<fCMi<<' '<<fCMs<<' '
         <<fTC<<' '<<fTMC<<' '<<fCQ<<' ';
}

void Fabbri::LongPrint(ostream &tempstr, double t,  ML_CalcType V) {
  Print(tempstr, t, V);
  const ML_CalcType V_int = V*1000.0;
  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);

  const double RTdF = ((v(VT_R)*(v(VT_T)))/(v(VT_F)));

  // V_Nernst
  double E_Na = RTdF*log((v(VT_Nao))/Na_i);
  double E_K = RTdF*log((v(VT_Ko))/K_i);
  double E_Ca = 0.5*RTdF*log((v(VT_Cao))/Ca_sub);


  // (V-V_Nernst)
  const double VmE_Na = V_int-E_Na;
  const double VmE_K = V_int-E_K;
  const double VmE_Ca = V_int-E_Ca;

  // I_f
  const double I_fNa = y*(v(VT_g_fNa))*VmE_Na;
  const double I_fK = y*(v(VT_g_fK))*VmE_K;
  const double I_f = I_fK+I_fNa;

  // I_NaK
  double i_Nak_max = (v(VT_i_NaK_max)*v(VT_Iso_INak_increase));
  double I_nak1 = i_Nak_max*pow((1+pow(((v(VT_Km_Kp))/(v(VT_Ko))), 1.2)), -1);
  double I_nak2 = pow((1+pow(((v(VT_Km_Nap))/Na_i), 1.3)), -1);
  double I_nak3 = pow((1+exp(-(VmE_Na+110)/20)), -1);

  const double I_NaK = I_nak1*I_nak2*I_nak3;

  // I_NaCa
  const double di = 1+ (Ca_sub/(v(VT_Kci)))*(1+exp(((-1*(v(VT_Qci)))*V_int)/RTdF)+(Na_i/(v(VT_Kcni))))+
    (Na_i/(v(VT_K1ni)))*(1+(Na_i/(v(VT_K2ni)))*(1+(Na_i/(v(VT_K3ni)))));
  const double k43 = Na_i/((v(VT_K3ni))+Na_i);
  const double k12 = ((Ca_sub/(v(VT_Kci)))*exp((-(v(VT_Qci))*V_int)/RTdF))/di;
  const double k14 = ((((Na_i/(v(VT_K1ni)))*Na_i)/(v(VT_K2ni)))*(1+(Na_i/(v(VT_K3ni))))*exp((v(VT_Qn)*V_int)/(2*RTdF)))/
    di;
  const double k41 = pCmP->k41[Vi];
  const double k34 = v(VT_k34); // (v(VT_Nao))/((v(VT_K3no))+(v(VT_Nao)));
  const double k21 = pCmP->k21[Vi];
  const double k23 = pCmP->k23[Vi];
  const double k32 = pCmP->k32[Vi];
  const double x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
  const double x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
  const double x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
  const double x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
  const double I_NaCa = ((v(VT_K_NaCa))*(x2*k21-(x1*k12)))/(x1+x2+x3+x4);

  // I_Na
  const double E_mh = RTdF*log(((v(VT_Nao))+0.12*(v(VT_Ko)))/(Na_i+0.12*K_i));
  const double I_Na_WT = (v(VT_g_na))*m*m*m*h*(V_int-E_mh);
  const double I_Na_Mut =
    (v(VT_g_na_mut)*m_mut*m_mut*m_mut*h_mut*(V_int-E_mh) + v(VT_g_na_l)*m_mut*m_mut*m_mut*(V_int-E_mh));
  const double I_Na = 0.5 * I_Na_WT + 0.5 * I_Na_Mut;

  // const double I_Na = (v(VT_g_na))*m*m*m*h*(V_int-E_mh);

  // I_CaL
  double Ca1 = RTdF*(1-exp((-2*V_int)/RTdF));
  double Ca2 = exp((-2*V_int)/RTdF);
  double CaL3 = RTdF*(1-exp((-V_int)/RTdF));
  double CaL4 = exp((-V_int)/RTdF);


  const double I_siCa = ((2*(v(VT_P_CaL))*V_int)/Ca1)*(Ca_sub-(v(VT_Cao))*Ca2)*dL*fL*fCa;
  const double I_siK = ((0.000365*(v(VT_P_CaL))*V_int)/CaL3)*(K_i-(v(VT_Ko))*CaL4)*dL*fL*fCa;
  const double I_siNa = ((0.0000185*(v(VT_P_CaL))*V_int)/CaL3)*(Na_i-(v(VT_Nao))*CaL4)*dL*fL*fCa;
  const double I_CaL = (I_siCa+I_siK+I_siNa)*v(VT_Iso_ICaL_increase)*(1-v(VT_ACh_iCaL_block));

  const double I_CaT = ((2*(v(VT_P_CaT))*V_int)/Ca1)*(Ca_sub-(v(VT_Cao))*Ca2)*dT*fT;
  const double I_Kur = ((v(VT_g_Kur))*r_Kur*s*VmE_K);
  const double I_to = ((v(VT_g_to))*VmE_K*q*r_to);
  const double I_Kr = ((v(VT_g_Kr))*sqrt(v(VT_Ko)/5.4)*VmE_K*(0.9*paF+0.1*paS)*piy);

  // I_Ks
  const double E_Ks = RTdF*log(((v(VT_Ko))+0.12*(v(VT_Nao)))/(K_i+0.12*Na_i));
  double R231C = 1;
  if (v(VT_R231C) == 1) {
    double A2      = 0.07;
    double g_ratio =  0.3153;
    double V_WT    = 102.4597;
    double k_WT    = 55.2060;
    double V_R231C = 113.8797;
    double k_R231C = -34.9006;

    R231C = g_ratio*pow((1/(1+exp(-(V_int-V_R231C)/k_R231C))), 2)/pow(((1-A2)/(1+exp(-(V_int-V_WT)/k_WT))+A2), 2);
  }
  const double I_Ks = R231C*v(VT_g_Ks)*v(VT_Iso_IKs_increase)*(V_int-E_Ks)*n*n;

  // I_SK
  const double I_SK = v(VT_g_SK) * VmE_K * x;

  // I_KACh
  double I_KACh = 0;
  if ((v(VT_ACh)) > 0) {
    I_KACh = (v(VT_g_KACh))*VmE_K*1+exp((V_int+20)/20)*a;
  } else {
    I_KACh = 0;
  }

  // IntrazellulÃ¤re Ca_StrÃ¶me:
  const double j_SRCarel = (v(VT_ks))*O_Ca_rel*(Ca_jsr-Ca_sub);

  double b_up = 0.7*v(VT_ACh)/(0.00009+v(VT_ACh))-v(VT_Iso_b_up);

  const double P_up = 5; // in mM/s
  const double j_Ca_dif = (Ca_sub-Ca_i)/(v(VT_tau_dif_Ca));
  const double j_up = P_up*(1-b_up)/(1+exp(((-1*Ca_i)+(v(VT_K_up)))/(v(VT_slope_up))));
  const double j_tr = (Ca_nsr-Ca_jsr)/(v(VT_tau_tr));


  double MaxSR = (v(VT_MaxSR));
  double EC50 = (v(VT_EC50_SR));
  double HSR = (v(VT_HSR));
  double Bruch = (v(VT_EC50_SR))/Ca_jsr;
  double downk_CaSR = (1+pow(Bruch, (v(VT_HSR))));
  const double k_CaSR = (v(VT_MaxSR))-(((v(VT_MaxSR))-(v(VT_MinSR)))/(1+pow((v(VT_EC50_SR))/Ca_jsr, (v(VT_HSR)))));

  const double I_mem = I_f+I_Kr+I_Ks+I_to+I_NaK+I_NaCa+I_Na+I_CaL+I_CaT+I_KACh+I_Kur+I_SK; // in nA

  tempstr<<I_f<<' '<<I_NaK<<' '
         <<I_NaCa<<' '<<I_Na<<' '
         <<I_CaL<<' '<<I_CaT<<' '
         <<I_Kur<<' '<<I_to<<' '
         <<I_Kr<<' '<<I_Ks<<' '
         <<I_SK<<' '<<I_KACh<<' '
         <<j_up<<' '<<j_SRCarel<<' '
         <<I_mem<<' ';
} // Fabbri::LongPrint

void Fabbri::GetParameterNames(vector<string> &getpara) {
  const int numpara = 35;
  const string ParaNames[numpara] =
  {"Na_i", "K_i", "Ca_i", "Mg_i", "y", "m", "h", "dL", "fL", "fCa", "dT", "fT", "r_Kur",
   "s", "r_to", "q", "paS", "paF", "piy", "n", "a", "x", "R_Ca_rel", "O_Ca_rel",
   "I_Ca_rel", "RI_Ca_rel", "Ca_jsr", "Ca_nsr",
   "Ca_sub", "fTMM", "fCMi", "fCMs", "fTC", "fTMC", "fCQ"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void Fabbri::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara = 15;
  const string ParaNames[numpara] =
  {"I_f", "I_NaK", "I_NaCa", "I_Na", "I_CaL", "I_CaT", "I_Kur", "I_to", "I_Kr", "I_Ks",
   "I_SK", "I_KACh", "J_UP", "j_SRCarel", "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
