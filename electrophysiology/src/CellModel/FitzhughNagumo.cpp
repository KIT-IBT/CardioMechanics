/*
 * File: FitzhughNagumo.cpp
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


#include <FitzhughNagumo.h>

FitzhughNagumo::FitzhughNagumo(FitzhughNagumoParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_FitzhughNagumoParameters::vtLast);
#endif // ifdef HETERO
  Init();
}

FitzhughNagumo::~FitzhughNagumo() {}

#ifdef HETERO

inline bool FitzhughNagumo::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else // ifdef HETERO

inline bool FitzhughNagumo::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif // ifdef HETERO

inline int FitzhughNagumo::GetSize(void) {
  return sizeof(FitzhughNagumo)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(FitzhughNagumoParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif // ifdef HETERO
  ;
}

inline unsigned char FitzhughNagumo::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void FitzhughNagumo::Init() {
#if KADEBUG
  cerr << "#initializing Class: FitzhughNagumo ... " << endl;
        #endif // if KADEBUG

  v = (v(VT_v_init));
  w = (v(VT_w_init));
}

ML_CalcType FitzhughNagumo::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                                 int euler = 2) {
  tinc *= 1000.0; // second to millisecond conversion

  const ML_CalcType J_stim = i_external;
  const ML_CalcType J_tot = J_stim + v*(v-v(VT_alpha))*(1-v)-w;

  v += J_tot*tinc;
  w += (v(VT_epsilon)*(v-v(VT_gamma)*w))*tinc;

  return .001 * J_tot; // needs calibration term
}

void FitzhughNagumo::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<< tArg << ' ' << V << ' ' << v << ' ' << w << ' ';
}

void FitzhughNagumo::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
}

void FitzhughNagumo::GetParameterNames(vector<string> &getpara) {
  const int numpara = 2;
  const string ParaNames[numpara] = {"v", "w"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void FitzhughNagumo::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
}
