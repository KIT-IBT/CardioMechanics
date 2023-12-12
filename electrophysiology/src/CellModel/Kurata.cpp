/*
 * File: Kurata.cpp
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


#include <Kurata.h>

Kurata::Kurata(KurataParameters *pp) {
  pCmP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(pCmP, NS_KurataParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

Kurata::~Kurata() {}

#ifdef HETERO

inline bool Kurata::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool Kurata::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int Kurata::GetSize(void) {
  return sizeof(Kurata)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(KurataParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char Kurata::getSpeed(ML_CalcType adVm) {
  return (unsigned char)(adVm < 1e-6 ? 2 : 1);
}

inline ML_CalcType Kurata::GetVm() {
  return v(VT_Init_Vm);
}

void Kurata::Init() {
#if KADEBUG
  cerr << "Kurata::Init" << endl;
#endif  // if KADEBUG

  // Init values
  gdL   = v(VT_Init_gdL);
  gfL   = v(VT_Init_gfL);
  gpa   = v(VT_Init_gpa);
  gn    = v(VT_Init_gn);
  gq    = v(VT_Init_gq);
  gh    = v(VT_Init_gh);
  Cai   = v(VT_Init_Cai);
  Carel = v(VT_Init_Carel);
  Caup  = v(VT_Init_Caup);
  gdR   = v(VT_Init_gdR);
  gfR   = v(VT_Init_gfR);
  Rtc   = v(VT_Init_Rtc);
  Nai   = v(VT_Init_Nai);
  Ki    = v(VT_Init_Ki);
}

ML_CalcType Kurata::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external = .0,  ML_CalcType stretch = 1.,
                         int euler                                            = 2) {
  tinc *= 1000.0;                        // time is in ms (milliseconds, 10^3)
  const  ML_CalcType V_int = V * 1000.0;  // voltage is in mV (millivolt, 10^3)
  // for accessing the table
  const int Vi       = (int)(DivisionTab * (RangeTabhalf + V_int) + .5);
  ML_CalcType I_stim = i_external;
  /**********************************************************************
  * Equilibrium Potentials
  **********************************************************************/
  const double ENa = v(VT_RTdF) * log(v(VT_Nao) / Nai);
  const double EK  = v(VT_RTdF) * log(v(VT_Ko) / Ki);
  const double ECa = v(VT_RTdF) * log(v(VT_Cao) / Cai) / 2;
  /**********************************************************************
   * Gating Variables
   ***********************************************************************/

  // only one here, for the rest, see table
  const double gfCas = v(VT_KmfCa) / (v(VT_KmfCa) + Cai);  // (8)
  /**********************************************************************
   * Ionic Currents
   ***********************************************************************/

  // iCa,L - L-type Ca2+ channel current
  const double iCaL = v(VT_gCaLmax) * (V_int - v(VT_ECaL) ) * gdL * gfL * gfCas;  // (3)
  // iKr - Rapidly activating delayed-rectifier K+ current
  const double EKr = EK;
  const double iKr = v(VT_gKrmax) * (V_int - EKr) * gpa * pCmP->gpis[Vi];  // (9) gpa and gpis are gating variables,
                                                                           // dependent from time
  // iKs - Slowly activating delayed-rectifier K+ current
  const double EKs = v(VT_RTdF) * log( (v(VT_PNaKs) * v(VT_Nao) + v(VT_Ko) ) / (v(VT_PNaKs) * Nai + Ki) );  // (16)
  const double iKs = v(VT_gKsmax) * (V_int - EKs) * (gn * gn);  // (15)
  // ito - Transient outward current
  const double Eto = v(VT_RTdF) * log( (v(VT_PNato) * v(VT_Nao) + v(VT_Ko) ) / (v(VT_PNato) * Nai + Ki) );  // (20)
  const double ito = v(VT_gtomax) * (V_int - Eto) * pCmP->grs[Vi] * gq;  // (19)
  // iNa - Na+ channel current
  const double Emh = v(VT_RTdF) * log( (v(VT_Nao) + v(VT_PKNa) * v(VT_Ko) ) / (Nai + v(VT_PKNa) * Ki) );  // (29)
  const double iNa = v(VT_gNamax) * (V_int - Emh) * pCmP->gms[Vi]*pCmP->gms[Vi]*pCmP->gms[Vi] * gh*gh;  // (28)
  // iK1 - Inward-rectifier K+ channels current
  const double gk1a = .1 / (1 + exp(.06 * (V_int - EK - 200) ) );  // (36)
  const double gk1b = (3 * exp(.0002 * (V_int - EK + 100) ) + exp(.1 * (V_int - EK - 10) ) ) /
    (1 + exp(-(.5) * (V_int - EK) ) );  // (37)
  const double iK1 = v(VT_gK1) * (V_int - EK) * gk1a / (gk1a + gk1b);  // (35)
  // background Na+/Ca2+ currents
  const double ibNa = v(VT_gbNa) * (V_int - ENa);  // (38)
  const double ibCa = v(VT_gbCa) * (V_int - ECa);  // (39)
  const double ibK  = v(VT_gbK) * (V_int - EK);

  // iNaK - Na+ K+ pump current
  const double iNaK = v(VT_iNaKmax) * (v(VT_Ko) / (v(VT_Ko) + v(VT_KmKp) ) ) /
    (1 +
     pow( (v(VT_KmNap) / Nai),
          v(VT_nNa) ) ) /
    (1 + .1245 * exp(-.1 * v(VT_FdRT) * V_int) + .0365 * v(VT_RhoNaK) * exp(-(v(VT_FdRT) ) * V_int) );  // (40)
  // iNaCa - Na+/Ca2+ exchanger current
  const double iNaCa = v(VT_kNaCa) /
    ( (v(VT_KmNaex))*(v(VT_KmNaex))*(v(VT_KmNaex)) + (v(VT_Nao))*(v(VT_Nao))*(v(VT_Nao)) ) /
    (v(VT_KmCaex) + v(VT_Cao) ) *
    (v(VT_Cao) * Nai*Nai*Nai * exp(v(VT_rNaCa) * v(VT_FdRT) * V_int) - (v(VT_Nao))*(v(VT_Nao))*(v(VT_Nao)) * Cai *
     exp( (v(VT_rNaCa) - 1) * v(VT_FdRT) * V_int) ) / (1 + v(VT_ksat) * exp( (v(VT_rNaCa) - 1) * v(VT_FdRT) * V_int) );  //
                                                                                                                         //
                                                                                                                         // (41)
  // ipCa - Sarcolemmal
  // Ca2+ pump current
  const double ipCa = v(VT_ipCamax) / (1 + (v(VT_KmCap) / Cai) );  // (42)
  // total membrane
  // current
  double itotal = iCaL + iKr + iKs + ito + iNa + iK1 + ibNa + ibCa + ibK + iNaK + iNaCa + ipCa;  // (43) ibK=0 in
                                                                                                 // control
  // net ion fluxes
  const double jNanet = (iNa + ibNa + 3 * iNaK + 3 * iNaCa) * v(VT_CmdF);  // (44)
  const double jKnet  = (iKr + iKs + ito + iK1 + ibK - 2 * iNaK) * v(VT_CmdF); // (45) ibK=0 in control
  const double jCanet = (iCaL + ibCa - 2 * iNaCa + ipCa) * v(VT_CmdTwoF);  // (46)
  // Ca2+-induced Ca2+ release via Ryanodine receptor in SR (mM/msec)
  // const double gdRt = 4 * 4; // (51)
  // const double gfRt = 4 * 4; // Kurata.doc says 4 * 4; // (51)
  const double gdRt = 4 * 4;
  const double gfRt = 4 * 4;
  const double jrel = (v(VT_Prel) / (1 + exp( (iCaL + ibCa - 2 * iNaCa + ipCa + 5) / .9) ) ) * (gdR * gfR)*(gdR * gfR)*
    (gdR * gfR) * (Carel - Cai);  // (47) + (48)
  // Ca2+ uptake via Ca2+-pump in SR (mM/msec)
  const double jup = v(VT_Pup) * (Cai * Cai) / ( (Cai * Cai) + (v(VT_Kup) * v(VT_Kup) ) );  // (52)
  // Ca2+ transfer/leak in SR (mM/msec)
  const double jtr   = (Caup - Carel) / v(VT_Ttr); // (53)
  const double jleak = v(VT_Pleak) * (Caup - Cai);  // (54)
  // Ca buffering flux
  const double Ftc = v(VT_kfTC) * Cai * (1 - Rtc) - v(VT_kbTC) * Rtc;  // (57)
  const double Bcm = 1 / (1 + v(VT_ConcCM) * v(VT_KdCM) / (v(VT_KdCM) + Cai)*(v(VT_KdCM) + Cai) );  // (58)
  const double Bcq = 1 / (1 + v(VT_ConcCQ) * v(VT_KdCQ) / (v(VT_KdCQ) + Carel)*(v(VT_KdCQ) + Carel) );  // (59)
  /**********************************************************************
   * Final Calculation
   ***********************************************************************/

  // (1) dVm/dt => see return
  // Gating Variables
  // (2) ddL/dt
  gdL += ( (pCmP->gdLs[Vi] - gdL) / pCmP->gdLt[Vi]) * tinc;  // (56)
  // (3) dfL/dt
  gfL += ( (pCmP->gfLs[Vi] - gfL) / pCmP->gfLt[Vi]) * tinc;  // (56)
  // (4) dpa/dt
  gpa += ( (pCmP->gpas[Vi] - gpa) / pCmP->gpat[Vi]) * tinc;  // (56)
  // (5) dn/dt
  gn += ( (pCmP->gns[Vi] - gn) / pCmP->gnt[Vi]) * tinc;  // (56)
  // (6) dq/dt
  gq += ( (pCmP->gqs[Vi] - gq) / pCmP->gqt[Vi]) * tinc;  // (56)
  // (7) dh/dt
  gh += ( (pCmP->ghs[Vi] - gh) / pCmP->ght[Vi]) * tinc;  // (56)
  // Ion Concentrations
  // (8) d[Ca]i/dt
  Cai +=
    (Bcm * ( (-jCanet + jrel * v(VT_Vrel) - jup * v(VT_Vup) + jleak * v(VT_Vup) ) / v(VT_Vi) - v(VT_ConcTC) * Ftc) ) *
    tinc;  // (60)
           // (9) d[Ca]rel/dt
  Carel += (Bcq * (jtr - jrel) ) * tinc;  // (61)
  // (10) d[Ca]up/dt
  Caup += (jup - jtr * v(VT_Vrel) / v(VT_Vup) - jleak) * tinc;  // (62)
  // (11) ddR/dt
  gdR += ( (pCmP->gdRs[Vi] - gdR) / gdRt) * tinc;  // (56)
  // (12) dfR/dt
  gfR += ( (pCmP->gfRs[Vi] - gfR) / gfRt) * tinc;  // (56)
  // (13) dfTnCa/dt
  Rtc += Ftc * tinc;

  // (14) d[Na]i/dt
  Nai += (-jNanet / v(VT_Vi) ) * tinc;  // (63)
  // (15) d[K]i/dt
  Ki += ( (I_stim * v(VT_C_m) / v(VT_F) - jKnet) / v(VT_Vi) ) * tinc;  // (64)
  return .001 * (I_stim - itotal) * tinc;  // (55)
}  // Kurata::Calc

void Kurata::Print(ostream &tempstr, double t,  ML_CalcType V) {
  tempstr << t << ' ' << V << ' ' << Cai << ' '  << Carel << ' ' << Caup << ' ' << Nai << ' ' << Ki << ' ';
}

void Kurata::LongPrint(ostream &tempstr, double t,  ML_CalcType V) {
  Print(tempstr, t, V);
  const  ML_CalcType V_int = V * 1000.0;  // voltage is in mV (millivolt, 10^3)
  // for accessing the table
  const int Vi = (int)(DivisionTab * (RangeTabhalf + V_int) + .5);
  /**********************************************************************
  * Equilibrium Potentials
  **********************************************************************/
  const double ENa = v(VT_RTdF) * log(v(VT_Nao) / Nai);
  const double EK  = v(VT_RTdF) * log(v(VT_Ko) / Ki);
  const double ECa = v(VT_RTdF) * log(v(VT_Cao) / Cai) / 2;
  /**********************************************************************
   * Gating Variables
   ***********************************************************************/

  // only one here, for the rest, see table
  const double gfCas = v(VT_KmfCa) / (v(VT_KmfCa) + Cai);  // (8)
  /**********************************************************************
   * Ionic Currents
   ***********************************************************************/

  // iCa,L - L-type Ca2+ channel current
  const double iCaL = v(VT_gCaLmax) * (V_int - v(VT_ECaL) ) * gdL * gfL * gfCas;  // (3)
  // iKr - Rapidly activating delayed-rectifier K+ current
  const double EKr = EK;
  const double iKr = v(VT_gKrmax) * (V_int - EKr) * gpa * pCmP->gpis[Vi];  // (9) gpa and gpis are gating variables,
                                                                           // dependent from time
  // iKs - Slowly activating delayed-rectifier K+ current
  const double EKs = v(VT_RTdF) * log( (v(VT_PNaKs) * v(VT_Nao) + v(VT_Ko) ) / (v(VT_PNaKs) * Nai + Ki) );  // (16)
  const double iKs = v(VT_gKsmax) * (V_int - EKs) * (gn * gn);  // (15)
  // ito - Transient outward current
  const double Eto = v(VT_RTdF) * log( (v(VT_PNato) * v(VT_Nao) + v(VT_Ko) ) / (v(VT_PNato) * Nai + Ki) );  // (20)
  const double ito = v(VT_gtomax) * (V_int - Eto) * pCmP->grs[Vi] * gq;  // (19)
  // iNa - Na+ channel current
  const double Emh = v(VT_RTdF) * log( (v(VT_Nao) + v(VT_PKNa) * v(VT_Ko) ) / (Nai + v(VT_PKNa) * Ki) );  // (29)
  const double iNa = v(VT_gNamax) * (V_int - Emh) * pCmP->gms[Vi]*pCmP->gms[Vi]*pCmP->gms[Vi] * gh*gh;  // (28)
  // iK1 - Inward-rectifier K+ channels current
  const double gk1a = .1 / (1 + exp(.06 * (V_int - EK - 200) ) );  // (36)
  const double gk1b = (3 * exp(.0002 * (V_int - EK + 100) ) + exp(.1 * (V_int - EK - 10) ) ) /
    (1 + exp(-(.5) * (V_int - EK) ) );  // (37)
  const double iK1 = v(VT_gK1) * (V_int - EK) * gk1a / (gk1a + gk1b);  // (35)
  // background Na+/Ca2+ currents
  const double ibNa = v(VT_gbNa) * (V_int - ENa);  // (38)
  const double ibCa = v(VT_gbCa) * (V_int - ECa);  // (39)
  const double ibK  = v(VT_gbK) * (V_int - EK);

  // iNaK - Na+ K+ pump current
  const double iNaK = v(VT_iNaKmax) * (v(VT_Ko) / (v(VT_Ko) + v(VT_KmKp) ) ) /
    (1 +
     pow( (v(VT_KmNap) / Nai),
          v(VT_nNa) ) ) /
    (1 + .1245 * exp(-.1 * v(VT_FdRT) * V_int) + .0365 * v(VT_RhoNaK) * exp(-(v(VT_FdRT) ) * V_int) );  // (40)
  // iNaCa - Na+/Ca2+ exchanger current
  const double iNaCa = v(VT_kNaCa) /
    ( (v(VT_KmNaex))*(v(VT_KmNaex))*(v(VT_KmNaex)) + (v(VT_Nao))*(v(VT_Nao))*(v(VT_Nao)) ) /
    (v(VT_KmCaex) + v(VT_Cao) ) *
    (v(VT_Cao) * Nai*Nai*Nai * exp(v(VT_rNaCa) * v(VT_FdRT) * V_int) - (v(VT_Nao))*(v(VT_Nao))*(v(VT_Nao)) * Cai *
     exp( (v(VT_rNaCa) - 1) * v(VT_FdRT) * V_int) ) / (1 + v(VT_ksat) * exp( (v(VT_rNaCa) - 1) * v(VT_FdRT) * V_int) );  //
                                                                                                                         //
                                                                                                                         // (41)
  // ipCa - Sarcolemmal Ca2+ pump current
  const double ipCa  = v(VT_ipCamax) / (1 + (v(VT_KmCap) / Cai) ); // (42)
  const double I_mem = iCaL + iKr + iKs + ito + iNa + iK1 + ibNa + ibCa + ibK + iNaK + iNaCa + ipCa;

  tempstr << iCaL << ' ' << iKr << ' ' << iKs << ' ' << ito << ' ' << iNa << ' ' << iK1 << ' ' << ibNa << ' ' << ibCa <<
    ' ' << ibK << ' '
          << iNaK << ' ' << iNaCa << ' ' << ipCa << ' ';
}  // Kurata::LongPrint

void Kurata::GetParameterNames(vector<string> &getpara) {
  const int numpara               = 5;
  const string ParaNames[numpara] = {"Cai", "Carel", "Caup", "Nai", "Ki"};

  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}

void Kurata::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const int numpara               = 13;
  const string ParaNames[numpara] =
  {"iCaL", "iKr", "iKs", "ito", "iNa", "iK1", "ibNa", "ibCa", "ibK", "iNaK", "iNaCa", "ipCa", "I_mem"};
  for (int i = 0; i < numpara; i++)
    getpara.push_back(ParaNames[i]);
}
