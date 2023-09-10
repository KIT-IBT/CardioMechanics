/**@file HybridModel.cpp
 * @brief <Brief (one-line) description here.>
 *
 * Please see the wiki for details on how to fill out this header:
 * https://intern.ibt.uni-karlsruhe.de/wiki/Document_a_IBT_C%2B%2B_tool
 *
 * @version 1.0.0
 *
 * @date Created <Your Name> (yyyy-mm-dd)
 *
 * @author Your Name\n
 *         Institute of Biomedical Engineering\n
 *         Karlsruhe Institute of Technology (KIT)\n
 *         http://www.ibt.kit.edu\n
 *         Copyright yyyy - All rights reserved.
 *
 * @see ...
 */

#include <HybridModel.h>

HybridModel::HybridModel(HybridModelParameters *pp) {
  hmp = pp;
  Init();
}

void HybridModel::Init() {
  AMATP   = v(VT_AMATPa);
  MATP    = v(VT_MATPa);
  MADPP   = v(VT_MADPPa);
  AwMADPP = v(VT_AwMADPPa);
  AsMADPP = v(VT_AsMADPPa);
  AsMADP  = v(VT_AsMADPa);
  AMADP   = v(VT_AMADPa);
  TCa     = v(VT_TCaa);
  TMon    = v(VT_TMona);

  MADP = v(VT_MADPa);
  M    = v(VT_Ma);
}

ML_CalcType HybridModel::Calc(double tinc, ML_CalcType stretch, ML_CalcType velocity, ML_CalcType &Ca, int euler = 1) {
  ML_CalcType dTCa = tinc*
    (-v(VT_k_on)*(2-MADP-M-AMATP-MATP-MADPP-AwMADPP)*
     pow(stretch, ML_CalcType(v(VT_TCa_stretch)))*Ca*(1-TCa)+v(VT_k_off)*TCa);

  TCa  -= dTCa;
  Ca   += dTCa*(v(VT_TCaMax)-v(VT_TCaMin));
  TMon += tinc*(v(VT_tm_on)*pow(1+(v(VT_TMon_coop)+stretch)*TMon, v(VT_TMon_pow))*TCa*(1-TMon)-v(VT_tm_off)*TMon);
  double k_5            = v(VT_k_5)*TMon;
  ML_CalcType velFactor = pow(fabs(velocity), v(VT_N_v))/(pow(fabs(velocity), v(VT_N_v))+pow(v(VT_v50), v(VT_N_v)));
  ML_CalcType k_7       = v(VT_k_7)*(v(VT_k7_base)-v(VT_k7_stretch)*stretch+fabs(velocity))/
    (1+v(VT_k7_force)*
     (Overlap(stretch, hmp->getOverlapID(),
              hmp->getOverlapParameters())*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP)*v(VT_dFmax)));
  return ForceEulerNEuler(euler, tinc/euler, k_5, stretch, velFactor, k_7);
}

ML_CalcType HybridModel::CalcTrop(double tinc, ML_CalcType stretch, ML_CalcType velocity, ML_CalcType TnCa,
                                  int euler = 1) {
  TCa = (TnCa-v(VT_TCaMin))/(v(VT_TCaMax)-v(VT_TCaMin));
  if (TCa < 0.)
    TCa = 0.; else if (TCa > 1.)
    TCa = 1.;
  TMon += tinc*(v(VT_tm_on)*pow(1+(v(VT_TMon_coop)+stretch)*TMon, v(VT_TMon_pow))*TCa*(1-TMon)-v(VT_tm_off)*TMon);
  double k_5            = v(VT_k_5)*TMon;
  ML_CalcType velFactor = pow(fabs(velocity), v(VT_N_v))/(pow(fabs(velocity), v(VT_N_v))+pow(v(VT_v50), v(VT_N_v)));
  ML_CalcType k_7       = v(VT_k_7)*(v(VT_k7_base)-v(VT_k7_stretch)*stretch+fabs(velocity))/
    (1+v(VT_k7_force)*
     (Overlap(stretch, hmp->getOverlapID(),
              hmp->getOverlapParameters())*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP)*v(VT_dFmax)));
  return ForceEulerNEuler(euler, tinc/euler, k_5, stretch, velFactor, k_7);
}

inline ML_CalcType HybridModel::ForceEulerNEuler(int r, ML_CalcType tinc, ML_CalcType k_5, ML_CalcType stretch,
                                                 ML_CalcType velFactor, ML_CalcType k_7) {
  while (r--) {
    double transition1 = v(VT_k_1)*v(VT_ATP)*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP-AsMADP-AMADP)-v(VT_k_m1)*AMATP;
    double transition2 = v(VT_k_2)*(1+v(VT_detach_vel)*velFactor)*AMATP;
    double transition3 = v(VT_k_3)*MATP-v(VT_k_m3)*MADPP;
    double transition4 = v(VT_k_4)*MADPP-v(VT_k_m4)*(1+v(VT_detach_vel)*velFactor)*AwMADPP;
    double transition5 = k_5*(v(VT_k5_stretch)*stretch+0.4)*
      pow((1+v(VT_k5_xb)*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP)), 2)*AwMADPP-v(VT_k_m5)*AsMADPP;
    double transition6 = v(VT_k_6)*AsMADPP-v(VT_k_m6)*AsMADP;
    double transition7 = k_7*AsMADP;
    double transition8 = v(VT_k_8)*AMADP-v(VT_k_m8)*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP-AsMADP-AMADP);


    AMATP   += tinc*(transition1-transition2);
    MATP    += tinc*(transition2-transition3+v(VT_ATP)*v(VT_k_14)*M);
    MADPP   += tinc*(transition3-transition4+velFactor*v(VT_k_13)*AsMADPP);
    AwMADPP += tinc*(transition4-transition5);
    AsMADPP += tinc*(transition5-transition6-velFactor*v(VT_k_13)*AsMADPP);
    AsMADP  += tinc*(transition6-transition7-velFactor*v(VT_k_11)*AsMADP);
    AMADP   += tinc*(transition7-transition8-velFactor*v(VT_k_10)*AMADP);
    MADP    += tinc*(velFactor*v(VT_k_11)*AsMADP+velFactor*v(VT_k_10)*AMADP-v(VT_k_12)*MADP);
    M       += tinc*
      (v(VT_k_12)*MADP+velFactor*v(VT_k_9)*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP-AsMADP-AMADP)-v(VT_ATP)*
       v(VT_k_14)*
       M);
  }

  return v(VT_F_physiol)*
         Overlap(stretch, hmp->getOverlapID(),
                 hmp->getOverlapParameters())*(1-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP)*v(VT_dFmax);
} // HybridModel::ForceEulerNEuler

void HybridModel::Print(ostream &tempstr) {
  tempstr<<1.0-MADP-M-AMATP-MATP-MADPP-AwMADPP-AsMADPP-AsMADP-AMADP<<' '<<AMATP<<' '
         <<MATP<<' '<<MADPP<<' '
         <<AwMADPP<<' '<<AsMADPP<<' '
         <<AsMADP<<' '<<AMADP<<' '
         <<MADP<<' '<<M<<' '
         <<1.0-TMon<<' '<<TMon<<' '
         <<1.0-TCa<<' '<<TCa<<' ';
}

void HybridModel::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 14;
  const string ParaNames[numpara] =
  {"AM",     "AMATP",      "MATP",                                      "MADPP",       "AwMADPP",       "AsMADPP",
   "AsMADP",
   "AMADP",  "MADP",       "M",                                         "TMoff",       "TMon",          "T",
   "TCa"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
