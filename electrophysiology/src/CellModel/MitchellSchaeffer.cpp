/* File: MitchellSchaeffer.cpp
        automatically created by CellML2Elphymodel.pl
        Institute of Biomedical Engineering, Universität Karlsruhe (TH) */

#include <MitchellSchaeffer.h>

MitchellSchaeffer::MitchellSchaeffer(MitchellSchaefferParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_MitchellSchaefferParameters::vtLast);
#endif // ifdef HETERO
  Init();
}

MitchellSchaeffer::~MitchellSchaeffer() {}

#ifdef HETERO

inline bool MitchellSchaeffer::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else // ifdef HETERO

inline bool MitchellSchaeffer::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif // ifdef HETERO

inline int MitchellSchaeffer::GetSize(void) {
  return sizeof(MitchellSchaeffer)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(MitchellSchaefferParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif // ifdef HETERO
  ;
}

inline unsigned char MitchellSchaeffer::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void MitchellSchaeffer::Init() {
#if KADEBUG
  cerr << "#initializing Class: MitchellSchaeffer ... " << endl;
        #endif // if KADEBUG

  u = (v(VT_u_init));
  h = (v(VT_h_init));
}

ML_CalcType MitchellSchaeffer::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,
                                    ML_CalcType stretch = 1.0, int euler = 2) {
  tinc *= 1000.0; // second to millisecond conversion
  const ML_CalcType J_stim = i_external;


  const ML_CalcType J_in = (h * u * u * (1 - u))/v(VT_tau_in);
  const ML_CalcType J_out = -(u / v(VT_tau_out));

  const ML_CalcType J_tot = J_in + J_out + J_stim;

  u += J_tot * tinc;
  h += (u < v(VT_V_gate) ? ((1 - h)/v(VT_tau_open))*tinc : (-h/v(VT_tau_close))*tinc);

  return .001 * J_tot;
}

void MitchellSchaeffer::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr<<tArg<< ' ' << V << ' ' << u << ' ' << h << ' ';
}

void MitchellSchaeffer::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);
}

void MitchellSchaeffer::GetParameterNames(vector<string> &getpara) {
  const int numpara = 2;
  const string ParaNames[numpara] = {"u", "h"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void MitchellSchaeffer::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
}
