/*
 * File: TenTusscher.cpp
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


#include <TenTusscher.h>

TenTusscherEtAl::TenTusscherEtAl(TenTusscherEtAlParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_TenTusscherEtAlParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

TenTusscherEtAl::~TenTusscherEtAl() {}

#ifdef HETERO

inline bool TenTusscherEtAl::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  //    cerr<<"add value in "<<PS<<endl;
  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool TenTusscherEtAl::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int TenTusscherEtAl::GetSize(void) {
  return sizeof(TenTusscherEtAl)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(TenTusscherEtAlParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char TenTusscherEtAl::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void TenTusscherEtAl::Init() {
#if KADEBUG
  cerr << "#initializing Class: TenTusscherEtAl ... " << endl;
#endif  // if KADEBUG
  m    = (v(VT_m_init));
  h    = (v(VT_h_init));
  j    = (v(VT_j_init));
  xr1  = (v(VT_xr1_init));
  xr2  = (v(VT_xr2_init));
  xs   = (v(VT_xs_init));
  r    = (v(VT_r_init));
  s    = (v(VT_s_init));
  d    = (v(VT_d_init));
  f    = (v(VT_f_init));
  fCa  = (v(VT_fCa_init));
  g    = (v(VT_g_init));
  Ca_i = (v(VT_Cai_init));
  CaSR = (v(VT_CaSR_init));
  Na_i = (v(VT_Nai_init));
  K_i  = fabs((v(VT_Ki_init)));
}

ML_CalcType TenTusscherEtAl::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch,
                                  int euler) {
  //    printf("TenTusscher: Calc
  // tinc=0.000000\tV=0.000000\ti_external=0.000000\tstretch=0.000000\teuler=0.000000\n",tinc,V,i_external,stretch,euler);
  ML_CalcType svolt = V*1000;  // membrane voltage in mV
  ML_CalcType HT    = tinc*1000; // timestep in ms
  const int   Vi    = (int)(DivisionTab*(RangeTabhalf+svolt)+.5); // array position
  // Needed to compute currents
  const  ML_CalcType EK      = (v(VT_RTONF))*(log(((v(VT_K_o))/K_i)));
  const  ML_CalcType ENa     = (v(VT_RTONF))*(log(((v(VT_Na_o))/Na_i)));
  const  ML_CalcType EKs     = (v(VT_RTONF))*(log((v(VT_KopKNaNao))/(K_i+(v(VT_pKNa))*Na_i)));
  const  ML_CalcType ECa     = 0.5*(v(VT_RTONF))*(log(((v(VT_Ca_o))/Ca_i)));
  const  ML_CalcType AK1     = 0.1/(1.+exp(0.06*(svolt-EK-200)));
  const  ML_CalcType BK1     = (3.*exp(0.0002*(svolt-EK+100))+exp(0.1*(svolt-EK-10)))/(1.+exp(-0.5*(svolt-EK)));
  const  ML_CalcType rec_iK1 = AK1/(AK1+BK1);
  const  ML_CalcType I_to    = (v(VT_g_to))*r*s*(svolt-EK);
  const  ML_CalcType I_Kr    = (v(VT_g_Kr))*sqrt((v(VT_K_o))/5.4)*xr1*xr2*(svolt-EK);
  const  ML_CalcType I_Ks    = (v(VT_g_Ks))*xs*xs*(svolt-EKs);
  const  ML_CalcType I_K1    = (v(VT_g_K1))*rec_iK1*(svolt-EK);
  const  ML_CalcType I_pCa   = (v(VT_g_pCa))*Ca_i/((v(VT_KpCa))+Ca_i);
  const  ML_CalcType I_pK    = (v(VT_g_pK))*ptTeaP->rec_ipK[Vi]*(svolt-EK);
  const  ML_CalcType I_bNa   = (v(VT_g_bNa))*(svolt-ENa);
  const  ML_CalcType I_bCa   = (v(VT_g_bCa))*(svolt-ECa);

  // Compute scalar currents
  const  ML_CalcType I_Na   = (v(VT_g_Na))*m*m*m*h*j*(svolt-ENa);
  const  ML_CalcType I_CaL  = d*f*fCa*(ptTeaP->CaL_P1[Vi]*Ca_i+ptTeaP->CaL_P2[Vi]);
  const  ML_CalcType I_NaCa = (ptTeaP->NaCa_P1[Vi])*Na_i*Na_i*Na_i-Ca_i*(ptTeaP->NaCa_P2[Vi]);
  const  ML_CalcType I_NaK  = (ptTeaP->NaK_P1[Vi])*(Na_i/(Na_i+(v(VT_KmNa))));
  const  ML_CalcType I_tot  = I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaK+I_NaCa+I_pCa+I_pK-i_external;

  // update concentrations
  const  ML_CalcType Caisquare     = Ca_i*Ca_i;
  const  ML_CalcType CaSRsquare    = CaSR*CaSR;
  const  ML_CalcType CaCurrent     = -(I_CaL+I_bCa+I_pCa-2*I_NaCa)*(v(VT_inverseVcF2C));
  const  ML_CalcType A             = 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
  const  ML_CalcType I_rel         = A*d*g;
  const  ML_CalcType I_leak        = 0.00008*(CaSR-Ca_i);
  const  ML_CalcType SERCA         = (v(VT_Vmaxup))/(1.+((v(VT_Kupsquare))/Caisquare));
  const  ML_CalcType HTCaSRCurrent = HT*(SERCA-I_rel-I_leak);
  const  ML_CalcType CaCSQN        = (v(VT_Bufsr))*CaSR/(CaSR+(v(VT_Kbufsr)));
  const  ML_CalcType dCaSR         = (v(VT_VcdVsr))*HTCaSRCurrent;
  const  ML_CalcType CaSum_1       = CaCSQN+dCaSR+CaSR;
  const  ML_CalcType bjsr          = (v(VT_BufsrPKbufsr))-CaSum_1;
  const  ML_CalcType cjsr          = (v(VT_Kbufsr))*CaSum_1;

  CaSR = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
  const  ML_CalcType CaBuf   = (v(VT_Bufc))*Ca_i/(Ca_i+(v(VT_Kbufc)));
  const  ML_CalcType dCai    = HT*CaCurrent-HTCaSRCurrent;
  const  ML_CalcType CaSum_2 = CaBuf+dCai+Ca_i;
  const  ML_CalcType bc      = (v(VT_BufcPKbufc))-CaSum_2;
  const  ML_CalcType cc      = (v(VT_Kbufc))*CaSum_2;
  Ca_i = (sqrt(bc*bc+4*cc)-bc)/2;
  const  ML_CalcType dNai = -(I_Na+I_bNa+3*I_NaK+3*I_NaCa)*(v(VT_inverseVcFC));
  Na_i += HT*dNai;
  const  ML_CalcType dKi = -(I_to+I_Kr+I_Ks+I_K1-2*I_NaK+I_pK-i_external)*(v(VT_inverseVcFC));
  K_i += HT*dKi;
  const  ML_CalcType FCa_INF =
    (1./
     (1.+
      pow((Ca_i*3076.923076923077), 8))+0.1/(1.+exp(10000*Ca_i-5))+0.20/(1.+exp(1250*Ca_i-0.9375))+0.23)*.684931506849;
  const  ML_CalcType G_INF =
    (Ca_i < .00035 ? 1./(1.+pow((Ca_i*2857.142857142857), 6)) : 1./(1.+pow((Ca_i*2857.142857142857), 16)));

  // set update gates prerequisites
  const  ML_CalcType M_INF   = ptTeaP->m_inf[Vi];
  const  ML_CalcType H_INF   = ptTeaP->h_inf[Vi];
  const  ML_CalcType J_INF   = ptTeaP->j_inf[Vi];
  const  ML_CalcType Xr1_INF = ptTeaP->Xr1_inf[Vi];
  const  ML_CalcType Xr2_INF = ptTeaP->Xr2_inf[Vi];
  const  ML_CalcType Xs_INF  = ptTeaP->Xs_inf[Vi];
  const  ML_CalcType S_INF   = ptTeaP->s_inf[Vi];
  const  ML_CalcType R_INF   = ptTeaP->r_inf[Vi];
  const  ML_CalcType D_INF   = ptTeaP->d_inf[Vi];
  const  ML_CalcType F_INF   = ptTeaP->f_inf[Vi];

  // update gates
  m   = M_INF-(M_INF-m)*ptTeaP->exptau_m[Vi];
  h   = H_INF-(H_INF-h)*ptTeaP->exptau_h[Vi];
  j   = J_INF-(J_INF-j)*ptTeaP->exptau_j[Vi];
  xr1 = Xr1_INF-(Xr1_INF-xr1)*ptTeaP->exptau_Xr1[Vi];
  xr2 = Xr2_INF-(Xr2_INF-xr2)*ptTeaP->exptau_Xr2[Vi];
  xs  = Xs_INF-(Xs_INF-xs)*ptTeaP->exptau_Xs[Vi];
  s   = S_INF-(S_INF-s)*ptTeaP->exptau_s[Vi];
  r   = R_INF-(R_INF-r)*ptTeaP->exptau_r[Vi];

  const  ML_CalcType fcaold    = fCa;
  const  ML_CalcType gold      = g;
  const  ML_CalcType exptaufca = exp(-HT/(v(VT_taufca)));
  const  ML_CalcType exptaug   = exp(-HT/(v(VT_taug)));
  d   = D_INF-(D_INF-d)*ptTeaP->exptau_d[Vi];
  f   = F_INF-(F_INF-f)*ptTeaP->exptau_f[Vi];
  fCa = FCa_INF-(FCa_INF-fCa)*exptaufca;
  g   = G_INF-(G_INF-g)*exptaug;
  if ((fCa > fcaold) && ( (svolt) > -37) )
    fCa = fcaold;
  if ((g > gold) && ( (svolt) > -37) )
    g = gold;
  /*    printf("Vectors:\n");printVector(&v_Ito_IKr_IKs_IK1);
      printVector(&v_IpCa_IpK_IbNa_IbCa);
      printf("return: 0.000000\n",tinc*(-I_tot));*/
  return tinc*(-I_tot);
}  // TenTusscherEtAl::Calc

void TenTusscherEtAl::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<<tArg<<' '<<V<<' '
         <<m<<' '<<h<<' '<<j<<' '<<d<<' '
         <<f<<' '<<fCa<<' '<<g<<' '<<xr1<<' '<<xr2<<' '
         <<xs<<' '<<r<<' '<<s<<' '<<Ca_i<<' '<<CaSR<<' '<<Na_i<<' '<<K_i;
}

void TenTusscherEtAl::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
  const  ML_CalcType svolt   = V*1000.0;
  const int Vi               = (int)(DivisionTab*(RangeTabhalf+svolt)+.5);
  const  ML_CalcType EK      = (v(VT_RTONF))*(log(((v(VT_K_o))/K_i)));
  const  ML_CalcType ENa     = (v(VT_RTONF))*(log(((v(VT_Na_o))/Na_i)));
  const  ML_CalcType EKs     = (v(VT_RTONF))*(log((v(VT_KopKNaNao))/(K_i+(v(VT_pKNa))*Na_i)));
  const  ML_CalcType ECa     = 0.5*(v(VT_RTONF))*(log(((v(VT_Ca_o))/Ca_i)));
  const  ML_CalcType AK1     = 0.1/(1.+exp(0.06*(svolt-EK-200)));
  const  ML_CalcType BK1     = (3.*exp(0.0002*(svolt-EK+100))+exp(0.1*(svolt-EK-10)))/(1.+exp(-0.5*(svolt-EK)));
  const  ML_CalcType rec_iK1 = AK1/(AK1+BK1);
  const  ML_CalcType I_pCa   = (v(VT_g_pCa))*Ca_i/((v(VT_KpCa))+Ca_i);
  const  ML_CalcType I_pK    = (v(VT_g_pK))*ptTeaP->rec_ipK[Vi]*(svolt-EK);
  const  ML_CalcType I_Na    = (v(VT_g_Na))*m*m*m*h*j*(svolt-ENa);
  const  ML_CalcType I_CaL   = d*f*fCa*(ptTeaP->CaL_P1[Vi]*Ca_i+ptTeaP->CaL_P2[Vi]);
  const  ML_CalcType I_NaCa  = (ptTeaP->NaCa_P1[Vi])*Na_i*Na_i*Na_i-Ca_i*(ptTeaP->NaCa_P2[Vi]);
  const  ML_CalcType I_NaK   = (ptTeaP->NaK_P1[Vi])*(Na_i/(Na_i+(v(VT_KmNa))));
  const  ML_CalcType I_to    = (v(VT_g_to))*r*s*(svolt-EK);
  const  ML_CalcType I_Kr    = (v(VT_g_Kr))*sqrt((v(VT_K_o))/5.4)*xr1*xr2*(svolt-EK);
  const  ML_CalcType I_Ks    = (v(VT_g_Ks))*xs*xs*(svolt-EKs);
  const  ML_CalcType I_K1    = (v(VT_g_K1))*rec_iK1*(svolt-EK);
  const  ML_CalcType I_bNa   = (v(VT_g_bNa))*(svolt-ENa);
  const  ML_CalcType I_bCa   = (v(VT_g_bCa))*(svolt-ECa);

  // const  ML_CalcType Caisquare=Ca_i*Ca_i;
  const  ML_CalcType CaSRsquare = CaSR*CaSR;

  // const  ML_CalcType CaCurrent=-(I_CaL+I_bCa+I_pCa-2*I_NaCa)*(v(VT_inverseVcF2C));
  const  ML_CalcType A      = 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
  const  ML_CalcType I_rel  = A*d*g;
  const  ML_CalcType I_leak = 0.00008*(CaSR-Ca_i);
  const  ML_CalcType I_mem  = I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaK+I_NaCa+I_pCa+I_pK;
  tempstr<<' '<<I_Na/0.185<<' '<<I_CaL/0.185<<' '<<I_bCa/0.185<<' '<<I_pCa/0.185<<' '<<I_to/0.185<<' '<<I_Ks/0.185<<
    ' '<<I_Kr/0.185<<' '<<I_K1/0.185
         <<' '<<I_pK/0.185<<' '<<I_bNa/0.185<<' '<<I_NaK/0.185<<' '<<I_NaCa/0.185<<' '<<I_rel
         <<' '<<I_leak <<' '<<I_mem << ' ';
}  // TenTusscherEtAl::LongPrint

void TenTusscherEtAl::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 16;
  const string ParaNames[numpara] =
  {"m", "h", "j", "d", "f", "fCa", "g", "Xr1", "Xr2", "Xs", "r", "s", "Ca_i", "Ca_SR", "Na_i", "K_i"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void TenTusscherEtAl::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 15;
  const string ParaNames[numpara] =
  {"I_Na",   "I_CaL",   "I_bCa", "I_pCa", "I_to", "I_Ks", "I_Kr", "I_K1", "I_pK", "I_bNa", "INaK", "I_NaCa", "I_rel",
   "I_leak", "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
