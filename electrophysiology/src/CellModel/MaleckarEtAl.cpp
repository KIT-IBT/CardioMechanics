/*
 * File: MaleckarEtAl.cpp
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


#include <MaleckarEtAl.h>

MaleckarEtAl::MaleckarEtAl(MaleckarEtAlParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_MaleckarEtAlParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

MaleckarEtAl::~MaleckarEtAl() {}

#ifdef HETERO

inline bool MaleckarEtAl::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool MaleckarEtAl::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int MaleckarEtAl::GetSize(void) {
  return (&yKur - &Ca_i + 1) * sizeof(ML_CalcType);
}

inline unsigned char MaleckarEtAl::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void MaleckarEtAl::Init() {
#if KADEBUG
  cerr << "#initializing Class: MaleckarEtAl ... " << endl;
        #endif  // if KADEBUG
  Na_c    = (v(VT_Na_c_init));
  Na_i    = (v(VT_Na_i_init));
  m       = (v(VT_m_init));
  h1      = (v(VT_h1_init));
  h2      = (v(VT_h2_init));
  Ca_d    = (v(VT_Ca_d_init));
  d_L     = (v(VT_d_L_init));
  f_L1    = (v(VT_f_L1_init));
  f_L2    = (v(VT_f_L2_init));
  K_c     = (v(VT_K_c_init));
  K_i     = (v(VT_K_i_init));
  r       = (v(VT_r_init));
  s       = (v(VT_s_init));
  a_ur    = (v(VT_a_ur_init));
  i_ur    = (v(VT_i_ur_init));
  n       = (v(VT_n_init));
  pa      = (v(VT_pa_init));
  Ca_c    = (v(VT_Ca_c_init));
  Ca_i    = (v(VT_Ca_i_init));
  O_C     = (v(VT_O_C_init));
  O_TC    = (v(VT_O_TC_init));
  O_TMgC  = (v(VT_O_TMgC_init));
  O_TMgMg = (v(VT_O_TMgMg_init));
  O       = (v(VT_O_init));
  Ca_rel  = (v(VT_Ca_rel_init));
  Ca_up   = (v(VT_Ca_up_init));
  O_Calse = (v(VT_O_Calse_init));
  F1      = (v(VT_F1_init));
  F2      = (v(VT_F2_init));
  yKur    = v(VT_Init_yKur);
}  // MaleckarEtAl::Init

ML_CalcType MaleckarEtAl::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
  ML_CalcType svolt = V*1000;
  ML_CalcType HT    = tinc;

  i_external = i_external*1000;  // convert from nA to pA
  const int Vi = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);

  // calculating algebraic part
  const ML_CalcType RTdF          = v(VT_RTdF);
  const ML_CalcType FdRT          = 1.0/RTdF;
  const ML_CalcType dVol_iF       = v(VT_dVol_iF);
  const ML_CalcType dVol_cF       = v(VT_dVol_cF);
  const ML_CalcType phi_Na_en     = v(VT_phi_Na_en);
  const ML_CalcType J_O_TMgMg     =  2000.00*(v(VT_Mg_i))*((1.00000 - O_TMgC) - O_TMgMg) -  666.000*O_TMgMg;
  const ML_CalcType r_Ca_d_term   = Ca_d/(Ca_d+(v(VT_k_rel_d)));
  const ML_CalcType r_Ca_d_factor =  r_Ca_d_term*r_Ca_d_term*r_Ca_d_term*r_Ca_d_term;
  const ML_CalcType r_Ca_i_term   = Ca_i/(Ca_i+(v(VT_k_rel_i)));
  const ML_CalcType r_Ca_i_factor =  r_Ca_i_term*r_Ca_i_term*r_Ca_i_term*r_Ca_i_term;
  const ML_CalcType r_act         =  203.800*(r_Ca_i_factor+r_Ca_d_factor);
  const ML_CalcType r_inact       = 33.9600+ 339.600*r_Ca_i_factor;
  const ML_CalcType E_K           =  RTdF*log(K_c/K_i);
  const ML_CalcType VmE_K         = (svolt - E_K);
  const ML_CalcType i_t           =  (v(VT_g_t))*r*s*VmE_K;
  const ML_CalcType i_Kur         =  (v(VT_g_kur))*yKur*a_ur*i_ur*VmE_K;
  const ML_CalcType i_K1          =
    ( (v(VT_g_K1))*(pow(K_c, 0.445700))*VmE_K)/(1.00000+(exp(((1.50000*(VmE_K+3.60000)*FdRT)))));
  const ML_CalcType i_Kr        =  (v(VT_g_Kr))*pa*ptTeaP->pip[Vi]*VmE_K;
  const ML_CalcType i_Ks        =  (v(VT_g_Ks))*n*VmE_K;
  const ML_CalcType pow_Na_i_15 = pow(Na_i, 1.50000);
  const ML_CalcType i_NaK       =
    ( (( (( (v(VT_i_NaK_max))*K_c)/(K_c+(v(VT_K_NaK_K))))*pow_Na_i_15)/(pow_Na_i_15+(v(VT_pow_K_NaK_Na_15))))*
      (svolt+150.000))/(svolt+200.000);
  const ML_CalcType E_Na   =  RTdF*log(Na_c/Na_i);
  const ML_CalcType VmE_Na = (svolt - E_Na);
  const ML_CalcType i_Na   = (v(VT_P_Na))*
    ( ((m*m*m*(0.900000*h1+ 0.100000*h2)*Na_c))*((exp(((VmE_Na*FdRT)))) - 1.00000))*ptTeaP->CexpINa[Vi];
  const ML_CalcType i_B_Na =  (v(VT_g_B_Na))*VmE_Na;
  const ML_CalcType i_NaCa =
    ( (v(VT_K_NaCa))*(Na_i*Na_i*Na_i*Ca_c*ptTeaP->CI_NaCa1[Vi] - Na_c*Na_c*Na_c*Ca_i*ptTeaP->CI_NaCa2[Vi]))/
    (1.00000+ (v(VT_d_NaCa))*(Na_c*Na_c*Na_c*Ca_i+ Na_i*Na_i*Na_i*Ca_c));
  const ML_CalcType f_Ca     = Ca_d/(Ca_d+(v(VT_k_Ca)));
  const ML_CalcType i_Ca_L   =  (v(VT_g_Ca_L))*d_L*(f_Ca*f_L1+ (1.00000 - f_Ca)*f_L2)*(svolt - (v(VT_E_Ca_app)));
  const ML_CalcType E_Ca     =  RTdF*0.5*log(Ca_c/Ca_i);
  const ML_CalcType i_B_Ca   =  (v(VT_g_B_Ca))*(svolt - E_Ca);
  const ML_CalcType i_CaP    = ( (v(VT_i_CaP_max))*Ca_i)/(Ca_i+(v(VT_k_CaP)));
  const ML_CalcType i_di     = (Ca_d - Ca_i)*(v(VT_Ci_di));
  const ML_CalcType i_KACh   =  ptTeaP->KIKACH[Vi]*VmE_K;
  const ML_CalcType J_O_C    =  200000.*Ca_i*(1.00000 - O_C) -  476.000*O_C;
  const ML_CalcType J_O_TC   =  78400.0*Ca_i*(1.00000 - O_TC) -  392.000*O_TC;
  const ML_CalcType J_O_TMgC =  200000.*Ca_i*((1.00000 - O_TMgC) - O_TMgMg) -  6.60000*O_TMgC;
  const ML_CalcType J_O      =  0.0800000*J_O_TC+ 0.160000*J_O_TMgC+ 0.0450000*J_O_C;
  const ML_CalcType i_up     = ( (v(VT_I_up_max))*(Ca_i/(v(VT_k_cyca)) - (v(VT_k_xcs2ds))*Ca_up))/
    ((Ca_i/(v(VT_k_cyca))+1.0)+( (v(VT_k_xcs))*(Ca_up/(v(VT_k_srca))+1.0)));
  const ML_CalcType i_rel_f2     = F2/(F2+0.250000);
  const ML_CalcType i_rel_factor =  i_rel_f2*i_rel_f2;
  const ML_CalcType i_rel        =  (v(VT_alpha_rel))*i_rel_factor*(Ca_rel - Ca_i);
  const ML_CalcType i_tr         = (Ca_up - Ca_rel)*(v(VT_Ci_tr));
  const ML_CalcType J_O_Calse    =  480.000*Ca_rel*(1.00000 - O_Calse) -  400.000*O_Calse;

  yKur = ptTeaP->yKur_infinity[Vi]+(yKur-ptTeaP->yKur_infinity[Vi])*ptTeaP->exptau_yKur[Vi];

  // calculating rates part
  O_TMgMg += HT*(J_O_TMgMg);

  // r += HT*( (ptTeaP->r_infinity[Vi] - r)/ptTeaP->tau_r[Vi]);
  r = ptTeaP->r_infinity[Vi]+(r-ptTeaP->r_infinity[Vi])*ptTeaP->exptau_r[Vi];

  // a_ur += HT*( (ptTeaP->a_ur_infinity[Vi] - a_ur)/ptTeaP->tau_a_ur[Vi]);
  a_ur = ptTeaP->a_ur_infinity[Vi]+(a_ur-ptTeaP->a_ur_infinity[Vi])*ptTeaP->exptau_a_ur[Vi];

  // i_ur += HT*( (ptTeaP->i_ur_infinity[Vi] - i_ur)/ptTeaP->tau_i_ur[Vi]);
  i_ur = ptTeaP->i_ur_infinity[Vi]+(i_ur-ptTeaP->i_ur_infinity[Vi])*ptTeaP->exptau_i_ur[Vi];

  // m += HT*( (ptTeaP->m_infinity[Vi] - m)/tau_m);
  m = ptTeaP->m_infinity[Vi]+(m-ptTeaP->m_infinity[Vi])*ptTeaP->exptau_m[Vi];

  // h1 += HT*( (ptTeaP->h_infinity[Vi] - h1)/tau_h1);
  h1 = ptTeaP->h_infinity[Vi]+(h1-ptTeaP->h_infinity[Vi])*ptTeaP->exptau_h1[Vi];

  // h2 += HT*( (ptTeaP->h_infinity[Vi] - h2)/tau_h2);
  h2 = ptTeaP->h_infinity[Vi]+(h2-ptTeaP->h_infinity[Vi])*ptTeaP->exptau_h2[Vi];

  // d_L += HT*( (ptTeaP->d_L_infinity[Vi] - d_L)/tau_d_L);
  d_L = ptTeaP->d_L_infinity[Vi]+(d_L-ptTeaP->d_L_infinity[Vi])*ptTeaP->exptau_d_L[Vi];

  // f_L1 += HT*( (ptTeaP->f_L_infinity[Vi] - f_L1)/tau_f_L1);
  f_L1 = ptTeaP->f_L_infinity[Vi]+(f_L1-ptTeaP->f_L_infinity[Vi])*ptTeaP->exptau_f_L1[Vi];

  // f_L2 += HT*( (ptTeaP->f_L_infinity[Vi] - f_L2)/tau_f_L2);
  f_L2 = ptTeaP->f_L_infinity[Vi]+(f_L2-ptTeaP->f_L_infinity[Vi])*ptTeaP->exptau_f_L2[Vi];

  // s += HT*( (ptTeaP->s_infinity[Vi] - s)/tau_s);
  s = ptTeaP->s_infinity[Vi]+(s-ptTeaP->s_infinity[Vi])*ptTeaP->exptau_s[Vi];

  // n += HT*( (ptTeaP->n_infinity[Vi] - n)/tau_n);
  n = ptTeaP->n_infinity[Vi]+(n-ptTeaP->n_infinity[Vi])*ptTeaP->exptau_n[Vi];

  // pa += HT*( (ptTeaP->p_a_infinity[Vi] - pa)/tau_pa);
  pa       = ptTeaP->p_a_infinity[Vi]+(pa-ptTeaP->p_a_infinity[Vi])*ptTeaP->exptau_pa[Vi];
  F1      += HT*(  (v(VT_r_recov))*((1.00000 - F1) - F2) -  r_act*F1);
  F2      += HT*(r_act*F1 -  r_inact*F2);
  K_i     += HT*(-(((i_t+i_Kur+i_K1+i_Ks+i_Kr) -  2.00000*i_NaK) - i_external)*dVol_iF);
  K_c     += HT*( ((v(VT_K_b)) - K_c)/(v(VT_tau_K))+((i_t+i_Kur+i_K1+i_Ks+i_Kr) -  2.00000*i_NaK)*dVol_cF);
  Na_i    += HT*(-(i_Na+i_B_Na+ 3.00000*i_NaCa+ 3.00000*i_NaK+phi_Na_en)*dVol_iF);
  Ca_c    += HT*( ((v(VT_Ca_b)) - Ca_c)/(v(VT_tau_Ca))+((i_Ca_L+i_B_Ca+i_CaP) -  2.00000*i_NaCa)*dVol_cF*0.5);
  Na_c    += HT*( ((v(VT_Na_b)) - Na_c)/(v(VT_tau_Na))+(i_Na+i_B_Na+ 3.00000*i_NaCa+ 3.00000*i_NaK+phi_Na_en)*dVol_cF);
  Ca_d    += HT*(-(i_Ca_L+i_di)*(v(VT_dVol_dF)));
  O_C     += HT*(J_O_C);
  O_TC    += HT*(J_O_TC);
  O_TMgC  += HT*(J_O_TMgC);
  O       += HT*(J_O);
  Ca_i    += HT*(-((i_B_Ca+i_CaP+i_up) - (i_di+i_rel+ 2.00000*i_NaCa))*dVol_iF*0.5 - J_O);
  Ca_up   += HT*( (i_up - i_tr)*(v(VT_dVol_upF)));
  O_Calse += HT*(J_O_Calse);
  Ca_rel  += HT*( (i_tr - i_rel)*(v(VT_dVol_relF)) -  31.0000*J_O_Calse);
  const ML_CalcType I_tot =
    (-i_external + i_Na + i_Ca_L + i_t + i_Kur + i_K1 + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP + i_NaCa +
     i_KACh)/
    v(VT_Cm);
  return tinc*(-I_tot);
}  // MaleckarEtAl::Calc

void MaleckarEtAl::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<<tArg<< ' ' << V<<' ' <<Na_c<<' ' <<Na_i<<' ' <<m<<' ' <<h1<<' ' <<h2<<' ' <<Ca_d<<' ' <<d_L<<' ' <<f_L1<<
    ' ' <<f_L2
         <<' ' <<K_c<<' ' <<K_i<<' ' <<r<<' ' <<s<<' ' <<a_ur<<' ' <<i_ur<<' ' <<n<<' ' <<pa<<' ' <<Ca_c
         <<' ' <<Ca_i<<' ' <<O_C<<' ' <<O_TC<<' ' <<O_TMgC<<' ' <<O_TMgMg<<' ' <<O<<' ' <<Ca_rel<<' ' <<Ca_up<<' ' <<
    O_Calse
         <<' ' <<F1<<' ' <<F2<<' '<<yKur;
}

void MaleckarEtAl::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
  const  ML_CalcType svolt        = V*1000.0;
  const int Vi                    = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);
  const ML_CalcType J_O_TMgMg     =  2000.00*(v(VT_Mg_i))*((1.00000 - O_TMgC) - O_TMgMg) -  666.000*O_TMgMg;
  const ML_CalcType r_Ca_d_term   = Ca_d/(Ca_d+(v(VT_k_rel_d)));
  const ML_CalcType r_Ca_d_factor =  r_Ca_d_term*r_Ca_d_term*r_Ca_d_term*r_Ca_d_term;
  const ML_CalcType r_Ca_i_term   = Ca_i/(Ca_i+(v(VT_k_rel_i)));
  const ML_CalcType r_Ca_i_factor =  r_Ca_i_term*r_Ca_i_term*r_Ca_i_term*r_Ca_i_term;
  const ML_CalcType r_act         =  203.800*(r_Ca_i_factor+r_Ca_d_factor);
  const ML_CalcType r_inact       = 33.9600+ 339.600*r_Ca_i_factor;
  const ML_CalcType E_K           =  (v(VT_RTdF))*(log((K_c/K_i)));
  const ML_CalcType i_t           =  (v(VT_g_t))*r*s*(svolt - E_K);
  const ML_CalcType i_Kur         =  (v(VT_g_kur))*yKur*a_ur*i_ur*(svolt - E_K);
  const ML_CalcType i_K1          =
    ( (v(VT_g_K1))*
      (pow((K_c/1.00000),
           0.445700))*(svolt - E_K))/
    (1.00000+(exp(((1.50000*((svolt - E_K)+3.60000)*(v(VT_F)))/( (v(VT_R))*(v(VT_T)))))));
  const ML_CalcType i_Kr        =  (v(VT_g_Kr))*pa*ptTeaP->pip[Vi]*(svolt - E_K);
  const ML_CalcType i_Ks        =  (v(VT_g_Ks))*n*(svolt - E_K);
  const ML_CalcType pow_Na_i_15 = pow(Na_i, 1.50000);
  const ML_CalcType i_NaK       =
    ( (( (( (v(VT_i_NaK_max))*K_c)/(K_c+(v(VT_K_NaK_K))))*pow_Na_i_15)/(pow_Na_i_15+(v(VT_pow_K_NaK_Na_15))))*
      (svolt+150.000))/(svolt+200.000);

  // const ML_CalcType i_Stim = (VOI - ptTeaP->past[Vi]>=(v(VT_stim_offset))&&VOI -
  // ptTeaP->past[Vi]<=(v(VT_stim_offset))+(v(VT_stim_duration)) ? (v(VT_stim_amplitude))*1.00000 : 0.00000)*1.00000;
  const ML_CalcType E_Na =  (v(VT_RTdF))*(log((Na_c/Na_i)));
  const ML_CalcType i_Na =
    ( (( (v(VT_P_Na))*m*m*m*(0.900000*h1+ 0.100000*h2)*Na_c*svolt*(v(VT_F))*(v(VT_F)))/( (v(VT_R))*(v(VT_T))))*
      ((exp((( (svolt - E_Na)*(v(VT_F)))/( (v(VT_R))*(v(VT_T)))))) - 1.00000))/
    ((exp(((svolt*(v(VT_F)))/( (v(VT_R))*(v(VT_T)))))) - 1.00000);
  const ML_CalcType i_B_Na =  (v(VT_g_B_Na))*(svolt - E_Na);
  const ML_CalcType i_NaCa =
    ( (v(VT_K_NaCa))*
      (Na_i*Na_i*Na_i*Ca_c*(exp((( (v(VT_F))*svolt*(v(VT_gamma_Na)))/( (v(VT_R))*(v(VT_T)))))) -  Na_c*Na_c*Na_c*Ca_i*
       (exp((( ((v(VT_gamma_Na)) - 1.00000)*svolt*(v(VT_F)))/( (v(VT_R))*(v(VT_T))))))))/
    (1.00000+ ((v(VT_d_NaCa))*1.00000)*(Na_c*Na_c*Na_c*Ca_i+ Na_i*Na_i*Na_i*Ca_c));
  const ML_CalcType f_Ca     = Ca_d/(Ca_d+(v(VT_k_Ca)));
  const ML_CalcType i_Ca_L   =  (v(VT_g_Ca_L))*d_L*(f_Ca*f_L1+ (1.00000 - f_Ca)*f_L2)*(svolt - (v(VT_E_Ca_app)));
  const ML_CalcType E_Ca     =  (v(VT_RTd2F))*(log((Ca_c/Ca_i)));
  const ML_CalcType i_B_Ca   =  (v(VT_g_B_Ca))*(svolt - E_Ca);
  const ML_CalcType i_CaP    = ( (v(VT_i_CaP_max))*Ca_i)/(Ca_i+(v(VT_k_CaP)));
  const ML_CalcType i_di     = ( (Ca_d - Ca_i)*2.00000*(v(VT_Vol_d))*(v(VT_F)))/(v(VT_tau_di));
  const ML_CalcType i_KACh   =  ptTeaP->KIKACH[Vi]*(svolt - E_K);
  const ML_CalcType J_O_C    =  200000.*Ca_i*(1.00000 - O_C) -  476.000*O_C;
  const ML_CalcType J_O_TC   =  78400.0*Ca_i*(1.00000 - O_TC) -  392.000*O_TC;
  const ML_CalcType J_O_TMgC =  200000.*Ca_i*((1.00000 - O_TMgC) - O_TMgMg) -  6.60000*O_TMgC;
  const ML_CalcType J_O      =  0.0800000*J_O_TC+ 0.160000*J_O_TMgC+ 0.0450000*J_O_C;
  const ML_CalcType i_up     =
    ( (v(VT_I_up_max))*(Ca_i/(v(VT_k_cyca)) - ( (v(VT_k_xcs))*(v(VT_k_xcs))*Ca_up)/(v(VT_k_srca))))/
    ((Ca_i+(v(VT_k_cyca)))/(v(VT_k_cyca))+( (v(VT_k_xcs))*(Ca_up+(v(VT_k_srca))))/(v(VT_k_srca)));
  const ML_CalcType i_rel_f2     = F2/(F2+0.250000);
  const ML_CalcType i_rel_factor =  i_rel_f2*i_rel_f2;
  const ML_CalcType i_rel        =  (v(VT_alpha_rel))*i_rel_factor*(Ca_rel - Ca_i);
  const ML_CalcType i_tr         = ( (Ca_up - Ca_rel)*2.00000*(v(VT_Vol_rel))*(v(VT_F)))/(v(VT_tau_tr));
  const ML_CalcType J_O_Calse    =  480.000*Ca_rel*(1.00000 - O_Calse) -  400.000*O_Calse;
  const ML_CalcType I_mem        = i_Na + i_Ca_L + i_t + i_Kur + i_K1 + i_Kr + i_Ks + i_B_Na + i_B_Ca + i_NaK + i_CaP +
    i_NaCa + i_KACh + i_di + i_up + i_rel + i_tr;
  tempstr<<' '<<i_Na << ' '<<i_Ca_L << ' '<<i_t << ' '<<i_Kur << ' '<<i_K1 << ' '<<i_Kr << ' '<<i_Ks << ' '<<i_B_Na <<
    ' '
         <<i_B_Ca << ' '<<i_NaK << ' '<<i_CaP << ' '<<i_NaCa << ' '<<i_KACh << ' '<<i_di << ' '<<i_up << ' '
         <<i_rel << ' '<<i_tr << ' '<< I_mem << ' ';
}  // MaleckarEtAl::LongPrint

void MaleckarEtAl::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 30;
  const string ParaNames[numpara] =
  {"Na_c",    "Na_i",       "m",           "h1",           "h2",             "Ca_d",          "d_L",          "f_L1",
   "f_L2",
   "K_c",
   "K_i",
   "r",
   "s",       "a_ur",
   "i_ur",    "n",          "pa",
   "Ca_c",    "Ca_i",       "O_C",         "O_TC",         "O_TMgC",         "O_TMgMg",       "O",            "Ca_rel",
   "Ca_up",
   "O_Calse",
   "F1",
   "F2",
   "yKur"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void MaleckarEtAl::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 18;
  const string ParaNames[numpara] =
  {"i_Na",   "i_Ca_L",   "i_t",       "i_Kur",       "i_K1",       "i_Kr",       "i_Ks",       "i_B_Na",       "i_B_Ca",
   "i_NaK",
   "i_CaP",
   "i_NaCa",
   "i_KACh",
   "i_di",   "i_up",     "i_rel",     "i_tr",        "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
