/*
 * File: BeelerReuter.cpp
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


#include <BeelerReuter.h>

BeelerReuter::BeelerReuter(BeelerReuterParameters *pp) {
  pBRP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(pBRP, NS_BeelerReuterParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

#ifdef HETERO

inline bool BeelerReuter::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool BeelerReuter::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

void BeelerReuter::Init() {
#if KADEBUG
  cerr << "initializing Class: BeelerReuter ... " << endl;
#endif  // if KADEBUG
  Ca_i = v(VT_Init_Ca_i);
  m    = v(VT_Init_m);
  h    = v(VT_Init_h);
  j    = v(VT_Init_j);
  d    = v(VT_Init_d);
  f    = v(VT_Init_f);
  x1   = v(VT_Init_x1);
}

void BeelerReuter::Print(ostream &tempstr, double t, ML_CalcType V) {
  tempstr<<t<<' '<<V<<' '
         <<m<<' '<<h<<' '<<j<<' '
         <<d<<' '<<f<<' '
         <<x1<<' '<<Ca_i<<' ';
}

void BeelerReuter::LongPrint(ostream &tempstr, double t, ML_CalcType V) {
  Print(tempstr, t, V);
  const ML_CalcType V_int = V*1000.0;
  const int Vi            = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);
  const double I_s        = v(VT_g_s)*d*f*(V_int+(82.3+13.0287*log(Ca_i)));
  const double I_K1       = pBRP->i_K_1[Vi];
  const double I_x1       = pBRP->i_x_1[Vi]*x1;
  const double I_Na       = v(VT_g_Na)*m*m*m*h*j+v(VT_g_NaC)*(V_int-v(VT_E_Na));
  tempstr<<I_s<<' '<<I_K1<<' '
         <<I_x1<<' '<<I_Na<<' ';
}

void BeelerReuter::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 7;
  const string ParaNames[numpara] = {"m", "h", "j", "d", "f", "x1", "Ca_i"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void BeelerReuter::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 4;
  const string ParaNames[numpara] = {"I_s", "I_K1", "I_x1", "I_Na"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

ML_CalcType BeelerReuter::Calc(double tinc, ML_CalcType V, ML_CalcType i_external = .0, ML_CalcType stretch = 1.,
                               int euler                                          = 1) {
  tinc *= 1000.0;
  const ML_CalcType V_int = V*1000.0;
  const int Vi            = (int)(DivisionTab*(RangeTabhalf+V_int)+.5);

  // x1+=tinc*(pBRP->a_x1[Vi]-pBRP->b_x1[Vi]*x1);
  const double x1_inf = pBRP->x1_inf[Vi];
  x1 = x1_inf + (x1 - x1_inf) * pBRP->exptau_x1[Vi];

  const double m_inf = pBRP->m_inf[Vi];
  m = m_inf + (m - m_inf) * pBRP->exptau_m[Vi];

  const double h_inf = pBRP->h_inf[Vi];
  h = h_inf + (h - h_inf) * pBRP->exptau_h[Vi];

  const double j_inf = pBRP->j_inf[Vi];
  j = j_inf + (j - j_inf) * pBRP->exptau_j[Vi];

  const double d_inf = pBRP->d_inf[Vi];
  d = d_inf + (d - d_inf) * pBRP->exptau_d[Vi];

  const double f_inf = pBRP->f_inf[Vi];
  f = f_inf + (f - f_inf) * pBRP->exptau_f[Vi];

  const double i_s = v(VT_g_s)*d*f*(V_int+(82.3+13.0287*log(Ca_i)));
  Ca_i += tinc*(-1e-7*i_s+.07*(1e-7-Ca_i));
  return -tinc*
         (pBRP->i_K_1[Vi]+pBRP->i_x_1[Vi]*x1+(v(VT_g_Na)*m*m*m*h*j+v(VT_g_NaC))*(V_int-v(VT_E_Na))+i_s-i_external)*
         (v(VT_dC_m));
}  // BeelerReuter::Calc

void BeelerReuter::GetStatus(double *p) const {
  p[0] = Ca_i;
  p[1] = m;
  p[2] = h;
  p[3] = j;
  p[4] = d;
  p[5] = f;
  p[6] = x1;
}

void BeelerReuter::SetStatus(const double *p) {
  Ca_i = p[0];
  m    = p[1];
  h    = p[2];
  j    = p[3];
  d    = p[4];
  f    = p[5];
  x1   = p[6];
}
