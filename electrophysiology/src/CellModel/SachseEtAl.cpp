/*
 * File: SachseEtAl.cpp
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



#include <SachseEtAl.h>

void SachseEtAl::Init() {
  Ki = v(VT_Ki);

  C0Shaker = v(VT_C0Shaker);
  C1Shaker = v(VT_C1Shaker);
  C2Shaker = v(VT_C2Shaker);
  C3Shaker = v(VT_C3Shaker);
  C4Shaker = v(VT_C4Shaker);
  OShaker  = v(VT_OShaker);
}

void SachseEtAl::Print(ostream &tempstr, double t, ML_CalcType V) {
  tempstr << t<< ' '  // 1
          << V<< ' '
          << Ki << ' '
          << C0Shaker << ' '
          << C1Shaker << ' '
          << C2Shaker << ' '
          << C3Shaker << ' '
          << C4Shaker << ' '
          << OShaker << ' '; // 9
}

void SachseEtAl::LongPrint(ostream &tempstr, double t, ML_CalcType V) {
  Print(tempstr, t, V);

  if (fabs(V) < 1e-7)
    V = 1e-7;

  const double Tx     = v(VT_Tx);
  const double Cm     = v(VT_Cm);
  const double Vfibro = v(VT_Vfibro);

  const double F         = ElphyModelConstants::F*1000;
  const double R         = ElphyModelConstants::R;
  const double RTdF      = R*Tx/F;
  const double FdRT      = 1./RTdF;
  const double VFdRT     = V*F/(R*Tx);
  const double dexpVFdRT = exp(-VFdRT);
  const double eVdkT     = ElphyModelConstants::qe *V/ElphyModelConstants::k/Tx;

  const double Ko = v(VT_Ko);

  const double EK      = RTdF*log(Ko/Ki);
  const double I_GHK_K = VFdRT*F*(Ki-Ko*dexpVFdRT)/(1.-dexpVFdRT);

  const double I_Shaker = v(VT_PShaker)*I_GHK_K*OShaker;

  tempstr << ' ' << I_Shaker;
  tempstr << ' ' << I_Shaker/v(VT_Cm);

  const double Kir   = 1./(v(VT_aKir)+exp(v(VT_bKir)*FdRT*(V-EK)));
  const double I_Kir = v(VT_GKir)*Kir*sqrt(Ko)*(V-EK);
  tempstr << ' ' << I_Kir;
  tempstr << ' ' << I_Kir/v(VT_Cm);

  ML_CalcType I_b = v(VT_Gb)*(V-v(VT_Eb));
  tempstr << ' ' << I_b;
  tempstr << ' ' << I_b/v(VT_Cm);

  ML_CalcType I_stretch = v(VT_Gstretch)*(V-v(VT_Estretch));
  tempstr << ' ' << I_stretch;
  tempstr << ' ' << I_stretch/v(VT_Cm);

  // Current through BK channel (static model from Horrigan&Aldrich, JGP 2002)
  const double L = v(VT_L0BK)*exp(v(VT_zLBK)*eVdkT);
  const double J = v(VT_J0BK)*exp(-v(VT_zJBK)*eVdkT);
  const double K = v(VT_Cai)/v(VT_KDBK);
  const double C = v(VT_CBK);
  const double D = v(VT_DBK);
  const double E = v(VT_EBK);
  double Po      = L*pow(1+K*C+J*D+J*K*C*D*E, 4);
  Po = Po/(Po+pow(1+J+K+J*K*E, 4));

  const double I_BK = v(VT_GBK)*Po*(V-EK);
  tempstr << ' ' << I_BK;
  tempstr << ' ' << I_BK/v(VT_Cm);

  ML_CalcType I_m = I_Shaker+I_Kir+I_b+I_stretch+I_BK;
  tempstr << ' ' << I_m;
  tempstr << ' ' << I_m/v(VT_Cm);
}  // SachseEtAl::LongPrint

void SachseEtAl::GetParameterNames(vector<string> &getpara) {
  const char *ParaNames[] = {"Ki [M]", "C0Shaker", "C1Shaker", "C2Shaker", "C3Shaker", "C4Shaker", "OShaker"};

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}

void SachseEtAl::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const char *ParaNames[] =
  {"I_Shaker [A]",    "I_Shaker [A/F]",    "I_Kir [A]",        "I_Kir [A/F]",  "I_b [A]", "I_b [A/F]", "I_Stretch [A]",
   "I_Stretch [A/F]", "I_BK [A]",          "I_BK [A/F]",       "I_m [A]",      "I_m [A/F]"};
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}

ML_CalcType SachseEtAl::Calc(double tinc, ML_CalcType V, ML_CalcType I_Stim, ML_CalcType stretch, int euler) {
  if (fabs(V) < 1e-7)
    V = 1e-7;

  const double Tx     = v(VT_Tx);
  const double Cm     = v(VT_Cm);
  const double Vfibro = v(VT_Vfibro);

  const double F         = ElphyModelConstants::F*1000;
  const double R         = ElphyModelConstants::R;
  const double RTdF      = R*Tx/F;
  const double FdRT      = 1./RTdF;
  const double VRTdF     = V*RTdF;
  const double VFdRT     = V*F/(R*Tx);
  const double dexpVFdRT = exp(-VFdRT);
  const double eVdkT     = ElphyModelConstants::qe *V/ElphyModelConstants::k/Tx;

  const double Ko = v(VT_Ko);

  const double EK      = RTdF*log(Ko/Ki);
  const double I_GHK_K = VFdRT*F*(Ki-Ko*dexpVFdRT)/(1.-dexpVFdRT);

  //   cerr << EK << endl;

  const double Qb = pow(2., (Tx-295.)/10.);

  const double rv  = v(VT_Shakerkv)*exp(VFdRT*v(VT_Shakerzv))*Qb;
  const double rvm = v(VT_Shakerkvm)*exp(VFdRT*v(VT_Shakerzvm))*Qb;
  const double ro  = v(VT_Shakerko)*Qb;
  const double rom = v(VT_Shakerkom)*exp(VFdRT*v(VT_Shakerzom))*Qb;

  const double dC0Shaker = -4.*rv*C0Shaker+        rvm*C1Shaker;
  const double dC1Shaker = 4.*rv*C0Shaker-(3.*rv+rvm)*C1Shaker+        2.*rvm*C2Shaker;
  const double dC2Shaker =                      3.*rv*C1Shaker-(2.*rv+2.*rvm)*C2Shaker+     3.*rvm*C3Shaker;
  const double dC3Shaker =                                              2.*rv*C2Shaker-(3.*rvm+rv)*C3Shaker+4.*rvm*
    C4Shaker;
  const double dC4Shaker =                                                                      rv*C3Shaker-(4.*rvm+ro)*
    C4Shaker+rom*OShaker;
  const double dOShaker =                                                                                           ro*
    C4Shaker-rom*OShaker;

  C0Shaker += tinc*dC0Shaker;
  if (C0Shaker < 0.)
    C0Shaker = 0.; else if (C0Shaker > 1.)
    C0Shaker = 1.;
  C1Shaker += tinc*dC1Shaker;
  if (C1Shaker < 0.)
    C1Shaker = 0.; else if (C1Shaker > 1.)
    C1Shaker = 1.;
  C2Shaker += tinc*dC2Shaker;
  if (C2Shaker < 0.)
    C2Shaker = 0.; else if (C2Shaker > 1.)
    C2Shaker = 1.;
  C3Shaker += tinc*dC3Shaker;
  if (C3Shaker < 0.)
    C3Shaker = 0.; else if (C3Shaker > 1.)
    C3Shaker = 1.;
  C4Shaker += tinc*dC4Shaker;
  if (C4Shaker < 0.)
    C4Shaker = 0.; else if (C4Shaker > 1.)
    C4Shaker = 1.;
  OShaker += tinc*dOShaker;
  if (OShaker < 0.)
    OShaker = 0.; else if (OShaker > 1.)
    OShaker = 1.;

  const double I_Shaker = v(VT_PShaker)*I_GHK_K*OShaker;

  const double Kir   = 1./(v(VT_aKir)+exp(v(VT_bKir)*FdRT*(V-EK)));
  const double I_Kir = v(VT_GKir)*Kir*sqrt(Ko*1e3)*(V-EK);  // bug Ko in mM -fs

  //    Ki-=-tinc*(I_Shaker+I_Kir-I_Stim)/(Volume()*F);

  // Unspecific background current
  const double I_b = v(VT_Gb)*(V-v(VT_Eb));

  // Current through stretch activated ion channels
  const double I_stretch = v(VT_Gstretch)*(V-v(VT_Estretch));

  // Current through BK channel (static model from Horrigan&Aldrich, JGP 2002)
  const double L = v(VT_L0BK)*exp(v(VT_zLBK)*eVdkT);
  const double J = v(VT_J0BK)*exp(-v(VT_zJBK)*eVdkT);
  const double K = v(VT_Cai)/v(VT_KDBK);
  const double C = v(VT_CBK);
  const double D = v(VT_DBK);
  const double E = v(VT_EBK);
  double Po      = L*pow(1+K*C+J*D+J*K*C*D*E, 4);
  Po = Po/(Po+pow(1+J+K+J*K*E, 4));

  const double I_BK = v(VT_GBK)*Po*(V-EK);

  // correction of stimulus current
  I_Stim *= 1e-9;

  return -tinc*(I_Shaker+I_Kir+I_b+I_stretch+I_BK-I_Stim)/Cm;
}  // SachseEtAl::Calc

void SachseEtAl::GetStatus(double *p) const {
  p[0] = Ki;
  p[1] = C0Shaker;
  p[2] = C1Shaker;
  p[3] = C2Shaker;
  p[4] = C3Shaker;
  p[5] = C4Shaker;
  p[6] = OShaker;
}

void SachseEtAl::SetStatus(const double *p) {
  Ki       = p[0];
  C0Shaker = p[1];
  C1Shaker = p[2];
  C2Shaker = p[3];
  C3Shaker = p[4];
  C4Shaker = p[5];
  OShaker  = p[6];
}
