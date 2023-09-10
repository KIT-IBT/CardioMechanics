/*      File: KurataParameters.cpp
    automatically created by ExtractParameterClass.pl - done by dw (22.03.2007)
    Institute of Biomedical Engineering, Universitt Karlsruhe (TH)
        send comments to dw@ibt.uka.de  */


#include <KurataParameters.h>

KurataParameters::KurataParameters(const char *initFile) {
  // Konstruktor
  P = new Parameter[vtLast];
  Init(initFile);
}

KurataParameters::~KurataParameters() {
  // Destruktor
}

void KurataParameters::Init(const char *initFile) {
#if KADEBUG
  cerr << "Loading the Kurata parameters from " << initFile << " ...\n";
#endif  // if KADEBUG

  // init values
  // constant values used in Calculate()


  // Initialization of the Parameters ...
  P[VT_F].name          =   "F";
  P[VT_R].name          =   "R";
  P[VT_Tx].name         =  "Tx";
  P[VT_nvar].name       =        "nvar";
  P[VT_Nao].name        = "Nao";
  P[VT_Ko].name         =  "Ko";
  P[VT_Cao].name        = "Cao";
  P[VT_C_m].name        = "C_m";
  P[VT_Vcell].name      =       "Vcell";
  P[VT_Vi].name         =  "Vi";
  P[VT_Vrel].name       =        "Vrel";
  P[VT_Vup].name        = "Vup";
  P[VT_gCaLmax].name    =     "gCaLmax";
  P[VT_KmfCa].name      =       "KmfCa";
  P[VT_ECaL].name       =        "ECaL";
  P[VT_gKrmax].name     =      "gKrmax";
  P[VT_gKsmax].name     =      "gKsmax";
  P[VT_PNaKs].name      =       "PNaKs";
  P[VT_gtomax].name     =      "gtomax";
  P[VT_PNato].name      =       "PNato";
  P[VT_gNamax].name     =      "gNamax";
  P[VT_PKNa].name       =        "PKNa";
  P[VT_gK1].name        = "gK1";
  P[VT_gbNa].name       =        "gbNa";
  P[VT_gbCa].name       =        "gbCa";
  P[VT_gbK].name        = "gbK";
  P[VT_iNaKmax].name    =     "iNaKmax";
  P[VT_KmNap].name      =       "KmNap";
  P[VT_nNa].name        = "nNa";
  P[VT_KmKp].name       =        "KmKp";
  P[VT_kNaCa].name      =       "kNaCa";
  P[VT_KmNaex].name     =      "KmNaex";
  P[VT_KmCaex].name     =      "KmCaex";
  P[VT_rNaCa].name      =       "rNaCa";
  P[VT_ksat].name       =        "ksat";
  P[VT_ipCamax].name    =     "ipCamax";
  P[VT_KmCap].name      =       "KmCap";
  P[VT_Prel].name       =        "Prel";
  P[VT_nrel].name       =        "nrel";
  P[VT_Pup].name        = "Pup";
  P[VT_Kup].name        = "Kup";
  P[VT_Pleak].name      =       "Pleak";
  P[VT_Ttr].name        = "Ttr";
  P[VT_ConcTC].name     =      "ConcTC";
  P[VT_kfTC].name       =        "kfTC";
  P[VT_kbTC].name       =        "kbTC";
  P[VT_ConcCM].name     =      "ConcCM";
  P[VT_KdCM].name       =        "KdCM";
  P[VT_ConcCQ].name     =      "ConcCQ";
  P[VT_KdCQ].name       =        "KdCQ";
  P[VT_Init_Vm].name    =     "Init_Vm";
  P[VT_Init_gdL].name   =    "Init_gdL";
  P[VT_Init_gfL].name   =    "Init_gfL";
  P[VT_Init_gpa].name   =    "Init_gpa";
  P[VT_Init_gn].name    =     "Init_gn";
  P[VT_Init_gq].name    =     "Init_gq";
  P[VT_Init_gh].name    =     "Init_gh";
  P[VT_Init_Cai].name   =    "Init_Cai";
  P[VT_Init_Carel].name =  "Init_Carel";
  P[VT_Init_Caup].name  =   "Init_Caup";
  P[VT_Init_gdR].name   =    "Init_gdR";
  P[VT_Init_gfR].name   =    "Init_gfR";
  P[VT_Init_Rtc].name   =    "Init_Rtc";
  P[VT_Init_Nai].name   =    "Init_Nai";
  P[VT_Init_Ki].name    =     "Init_Ki";
  P[VT_RTdF].name       =        "RTdF";
  P[VT_FdRT].name       =        "FdRT";
  P[VT_RhoNaK].name     =      "RhoNaK";
  P[VT_CmdF].name       =        "CmdF";
  P[VT_CmdTwoF].name    =     "CmdTwoF";
  P[VT_Amp].name        =         "Amp";

  P[VT_RTdF].readFromFile    = false;
  P[VT_FdRT].readFromFile    = false;
  P[VT_RhoNaK].readFromFile  = false;
  P[VT_CmdF].readFromFile    = false;
  P[VT_CmdTwoF].readFromFile = false;

  ParameterLoader EPL(initFile, EMT_Kurata);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile) {
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
    }

    // End Initialization of the Parameters ...
    else {
#if KADEBUG
#endif
    }
  Calculate();
  InitTable();
}  // KurataParameters::Init

void KurataParameters::Calculate() {
#if KADEBUG
  cerr << "KurataParameters - Calculate ..." << endl;
#endif  // if KADEBUG
  P[VT_RTdF].value = P[VT_R].value * P[VT_Tx].value / P[VT_F].value; // multiply RTdF with 1000 since it used
                                                                     // together with mV later
  P[VT_FdRT].value    = 1 / P[VT_RTdF].value;
  P[VT_RhoNaK].value  = (exp(P[VT_Nao].value / 67.3) - 1) / 7;
  P[VT_CmdF].value    = P[VT_C_m].value / P[VT_F].value;
  P[VT_CmdTwoF].value = P[VT_C_m].value / (2 * P[VT_F].value);
  if (PrintParameterMode == PrintParameterModeOn)
    PrintParameters();
}

void KurataParameters::InitTable() {
  for (double V = -RangeTabhalf + .0001; V < RangeTabhalf; V += dDivisionTab) {
    int Vi = (int)(DivisionTab * (RangeTabhalf + V) + .5);
    /**********************************************************************
     * Gating Variables
     ***********************************************************************/

    // dL/fL/fCa (iCa,L) - L-type Ca2+ channel current
    gdLs[Vi] = 1 / (1 + exp(-(V + 7.64) / 6.32) );  // (4)
    gfLs[Vi] = 1 / (1 + exp( (V + 24.6) / 6.9) );  // (5)
    gdLt[Vi] = (1.4 / (1 + exp(-(V + 35) / 13) ) + .25) * 1.4 / (1 + exp( (V + 5) / 5) ) + 1 /
      (1 + exp(-(V - 50) / 20) );  // (6)
    gfLt[Vi] = 17.925 / (.1389 * exp(-( (.0358 * (V - 10.9) ) * (.0358 * (V - 10.9) ) ) ) + .0519);  // (7)
    // pa/pi (iKr) - Rapidly activating delayed-rectifier K+ current
    gpas[Vi] = 1 / (1 + exp(-(V + 14) / 7.7) );  // (10)
    double gpaa = .0003 * (V + 14) / (1 - exp(-(V + 14) / 5) );  // (13)
    double gpab = .000073898 * (V - 3.4328) / (exp( (V - 3.4328) / 5.1237) - 1);  // (14)
    gpat[Vi] = 1 / (gpaa + gpab);  // (12)
    gpis[Vi] = 1 / (1 + exp( (V + 15) / 22.4) );  // (11)
    // n (iKs) - Slowly activating delayed-rectifier K+ current
    gns[Vi] = 1 / sqrt(1 + exp(-(V - 9.4) / 11.8) );  // (17)
    gnt[Vi] = 555 / (1 + exp(-(V + 22) / 11.3) ) + 129;  // (18)
    // r/q (ito) - Transient outward current
    double gra = .5266 * exp(-.0166 * (V - 42.2912) ) / (1 + exp(-.0943 * (V - 42.2912) ) );  // (22) alpha
    double grb = (.5149 * exp(-.1344 * (V - 5.0027) ) + .00005186 * V) / (1 + exp(-.1348 * (V - .00005186) ) );  // (23)
                                                                                                                 // beta
    grs[Vi] = gra / (gra + grb);  // n_unendlich
    double gqa = (.0721 * exp(-(.173) * (V + 34.2531) ) + .00005612 * V) / (1 + exp(-.1732 * (V + 34.2531) ) );  // (26)
    double gqb = (.0767 * exp(-(.00000000166) * (V + 34.0235) ) + .0001215 * V) / (1 + exp(-(.1604) * (V + 34.0235) ) );  //
                                                                                                                          //
                                                                                                                          // (27)
    gqs[Vi] = gqa / (gqa + gqb);  // (24)
    gqt[Vi] = 1 / (gqa + gqb);  // (25) tau
    // m/h (iNa)
    double gma = .32 * (V + 47.13) / (1 - exp(-(V + 47.13) / 10) );  // (31)
    double gmb = .08 * exp(-V / 11);  // (32)
    gms[Vi]  = gma / (gma + gmb); // (30)
    ghs[Vi]  = .5 * (1 - tanh(7.74 + .12 * V) ); // (33)
    ght[Vi]  = .25 + 2.24 * (1 - tanh(7.74 + .12 * V) ) / (1 - tanh(.07 * (V + 92.4) ) ); // (34)
    gdRs[Vi] = 1 / (1 + exp(-(V + 7.64) / 6.32) );  // (49) gdRs=gdLs;
    gfRs[Vi] = 1 / (1 + exp( (V + 24.6) / 6.9) );  // (50) gfRs=gfLs;
  }
}  // KurataParameters::InitTable

void KurataParameters::PrintParameters() {
  // print the parameter to the stdout
  cout<<"KurataParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= " << P[i].value << endl;
  }
}
