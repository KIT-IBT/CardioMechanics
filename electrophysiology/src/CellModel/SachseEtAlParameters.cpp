/*! \file SachseEtAlParameters.cpp
   \brief Implementation of  fibroblast model
   This model uses only SI units.
   Rat ventricular fibroblast Ann Biomed Eng Jan 2008
   Human synovial fibrobast 2013
   \author fs, CVRTI - University of Utah, USA
 */


#include <SachseEtAlParameters.h>

SachseEtAlParameters::SachseEtAlParameters(const char *initFile) {
  P = new Parameter[vtLast];
  Init(initFile);
}

void SachseEtAlParameters::PrintParameters() {
  cout << "SachseEtAlParameter:" << endl;

  for (int i = vtFirst; i < vtLast; i++)
    cout << P[i].name << "\t = " << P[i].value << endl;
}

void SachseEtAlParameters::Init(const char *initFile) {
#if KADEBUG
  cerr << "SachseEtAlParameters:init " << initFile << endl;
#endif  // if KADEBUG

  P[VT_Tx].name        = "Tx";
  P[VT_Vm].name        = "Vm";
  P[VT_Cm].name        = "Cm";
  P[VT_Vfibro].name    = "Vfibro";
  P[VT_Ko].name        = "Ko";
  P[VT_Ki].name        = "Ki";
  P[VT_Cao].name       = "Cao";
  P[VT_Cai].name       = "Cai";
  P[VT_C0Shaker].name  = "C0Shaker";
  P[VT_C1Shaker].name  = "C1Shaker";
  P[VT_C2Shaker].name  = "C2Shaker";
  P[VT_C3Shaker].name  = "C3Shaker";
  P[VT_C4Shaker].name  = "C4Shaker";
  P[VT_OShaker].name   = "OShaker";
  P[VT_PShaker].name   = "PShaker";
  P[VT_Shakerkvm].name = "kvm";
  P[VT_Shakerkv].name  = "kv";
  P[VT_Shakerko].name  = "ko";
  P[VT_Shakerkom].name = "kom";
  P[VT_Shakerzv].name  = "zv";
  P[VT_Shakerzvm].name = "zvm";
  P[VT_Shakerzom].name = "zom";
  P[VT_GKir].name      = "GKir";
  P[VT_aKir].name      = "aKir";
  P[VT_bKir].name      = "bKir";
  P[VT_Gb].name        = "Gb";
  P[VT_Eb].name        = "Eb";
  P[VT_Gstretch].name  = "Gstretch";
  P[VT_Estretch].name  = "Estretch";
  P[VT_GBK].name       = "GBK";
  P[VT_L0BK].name      = "L0BK";
  P[VT_zLBK].name      = "zLBK";
  P[VT_J0BK].name      = "J0BK";
  P[VT_zJBK].name      = "zJBK";
  P[VT_KDBK].name      = "KDBK";
  P[VT_CBK].name       = "CBK";
  P[VT_DBK].name       = "DBK";
  P[VT_EBK].name       = "EBK";
  P[VT_Amp].name       = "Amp";

  ParameterLoader EPL(initFile, EMT_SachseEtAl);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);

  Calculate();
}  // SachseEtAlParameters::Init
