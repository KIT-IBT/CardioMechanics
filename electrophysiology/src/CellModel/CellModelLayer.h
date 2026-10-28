/*
 * File: CellModelLayer.h
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


#ifndef CELLMODELLAYER_H
#define CELLMODELLAYER_H

#include <ForceDummy.h>
#include <HybridModel.h>
#include <Land17.h>

#include <ElphyDummy.h>
#include <BeelerReuter.h>
#include <CourtemancheEtAl.h>
#include <TenTusscher.h>
#include <TenTusscher2.h>
#include <SachseEtAl.h>
#include <MaleckarEtAl.h>
#include <KoivumaekiEtAl.h>
#include <GrandiEtAlVentricle.h>
#include <GrandiEtAlAtrium.h>
#include <FabbriEtAl.h>
#include <Himeno.h>
#include <OHaraRudy.h>
#include <OHaraRudyIso.h>
#include <Tomek.h>
#include <TWorld.h>
#include <MitchellSchaeffer.h>
#include <FitzhughNagumo.h>

template<class CalcType>

inline void initForceParameters(vbForceParameters<CalcType> **pfp, const char *initFVFile, ForceModelType fmt) {
  if (*pfp) {
    delete *pfp;
    *pfp = NULL;
  }
  switch (fmt) {
    case FMT_Hybrid:
      *pfp = new HybridModelParameters(initFVFile); return;

    case FMT_Land17:
      *pfp = new Land17Parameters(initFVFile); return;

    case FMT_Dummy:
      return;

    case FMT_Last:
    default:
      throw kaBaseException("Init force parameters failed"); return;
  }  // switch
}  // initForceParameters

template<class CalcType>

inline void initForceModel(vbForceModel<CalcType> **pfm, vbForceParameters<CalcType> *pfp, ForceModelType fmt) {
  if (*pfm) {
    delete *pfm;
    *pfm = NULL;
  }
  switch (fmt) {
    case FMT_Dummy:
      *pfm = new ForceDummy<CalcType>();
      break;
    case FMT_Hybrid:
      *pfm = new HybridModel((HybridModelParameters *)pfp);
      break;
    case FMT_Land17:
      *pfm = new Land17((Land17Parameters *)pfp);
      break;
    default:
      throw kaBaseException("Undefined force model");
      break;
  }  // switch
}  // initForceModel

template<class CalcType>

inline void initElphyParameters(vbElphyParameters<CalcType> **pep, const char *initEVFile, ElphyModelType emt,
                                ML_CalcType tinc) {
  if (*pep) {
    delete *pep;
    *pep = NULL;
  }
  switch (emt) {
    case EMT_BeelerReuter:
      *pep = new BeelerReuterParameters(initEVFile, tinc); return;

    case EMT_CourtemancheEtAl:
      *pep = new CourtemancheParameters(initEVFile, tinc); return;

    case EMT_TenTusscher:
      *pep = new TenTusscherEtAlParameters(initEVFile, tinc); return;

    case EMT_TenTusscher2:
      *pep = new TenTusscher2Parameters(initEVFile, tinc); return;

    case EMT_SachseEtAl:
      *pep = new SachseEtAlParameters(initEVFile); return;

    case EMT_MaleckarEtAl:
      *pep = new MaleckarEtAlParameters(initEVFile, tinc); return;

    case EMT_KoivumaekiEtAl:
      *pep = new KoivumaekiEtAlParameters(initEVFile, tinc); return;

    case EMT_GrandiEtAlVentricle:
      *pep = new GrandiEtAlVentricleParameters(initEVFile, tinc); return;

    case EMT_GrandiEtAlAtrium:
      *pep = new GrandiEtAlAtriumParameters(initEVFile, tinc); return;

    case EMT_FabbriEtAl:
      *pep = new FabbriParameters(initEVFile, tinc); return;

    case EMT_HimenoEtAl:
      *pep = new HimenoParameters(initEVFile, tinc); return;

    case EMT_OHaraRudy:
      *pep = new OHaraRudyParameters(initEVFile, tinc); return;

    case EMT_OHaraRudyIso:
      *pep = new OHaraRudyIsoParameters(initEVFile, tinc); return;

    case EMT_Tomek:
      *pep = new TomekParameters(initEVFile, tinc); return;

    case EMT_TWorld:
      *pep = new TWorldParameters(initEVFile, tinc); return;

    case EMT_MitchellSchaeffer:
      *pep = new MitchellSchaefferParameters(initEVFile); return;

    case EMT_FitzhughNagumo:
      *pep = new FitzhughNagumoParameters(initEVFile); return;

    case EMT_Dummy:
      return;

    case EMT_Last:
    default:
      throw kaBaseException("Init elphy parameters failed"); return;
  }  // switch
}  // initElphyParameters

template<class CalcType>

inline void initElphyModel(vbElphyModel<CalcType> **pem, vbElphyParameters<CalcType> *pep, ElphyModelType emt,
                           bool opt = false) {
  if (*pem) {
    delete *pem;
    *pem = NULL;
  }
  switch (emt) {
    case EMT_Dummy:
      *pem = new ElphyDummy<CalcType>();
      break;

    case EMT_BeelerReuter:
      *pem = new BeelerReuter((BeelerReuterParameters *)pep);
      break;

    case EMT_CourtemancheEtAl:
      *pem = new Courtemanche((CourtemancheParameters *)pep);
      break;

    case EMT_TenTusscher:
      *pem = new TenTusscherEtAl((TenTusscherEtAlParameters *)pep);
      break;

    case EMT_TenTusscher2:
      *pem = new TenTusscherEtAl2((TenTusscher2Parameters *)pep);
      break;

    case EMT_SachseEtAl:
      *pem = new SachseEtAl((SachseEtAlParameters *)pep);
      break;

    case EMT_MaleckarEtAl:
      *pem = new MaleckarEtAl((MaleckarEtAlParameters *)pep);
      break;

    case EMT_KoivumaekiEtAl:
      *pem = new KoivumaekiEtAl((KoivumaekiEtAlParameters *)pep);
      break;

    case EMT_GrandiEtAlVentricle:
      *pem = new GrandiEtAlVentricle((GrandiEtAlVentricleParameters *)pep);
      break;

    case EMT_GrandiEtAlAtrium:
      *pem = new GrandiEtAlAtrium((GrandiEtAlAtriumParameters *)pep);
      break;

    case EMT_FabbriEtAl:
      *pem = new Fabbri((FabbriParameters *)pep);
      break;

    case EMT_HimenoEtAl:
      *pem = new Himeno((HimenoParameters *)pep);
      break;

    case EMT_OHaraRudy:
      *pem = new OHaraRudy((OHaraRudyParameters *)pep);
      break;

    case EMT_OHaraRudyIso:
      *pem = new OHaraRudyIso((OHaraRudyIsoParameters *)pep);
      break;

    case EMT_Tomek:
      *pem = new Tomek((TomekParameters *)pep);
      break;

    case EMT_TWorld:
      *pem = new TWorld((TWorldParameters *)pep);
      break;

    case EMT_MitchellSchaeffer:
      *pem = new MitchellSchaeffer((MitchellSchaefferParameters *)pep);
      break;

    case EMT_FitzhughNagumo:
      *pem = new FitzhughNagumo((FitzhughNagumoParameters *)pep);
      break;

    default:
      throw kaBaseException("Undefined elphy model");
  }  // switch
}  // initElphyModel

template<class CalcType>

inline void PreCalcModels(vbElphyModel<CalcType> **PCEM, vbForceModel<CalcType> **PCFM, CalcType bclmin,
                          CalcType *Vm_loc, CalcType *Force_loc, double tinc, CalcType basis = 0.96) {
#if KADEBUG
  cerr<<"PreCalcModels\n";
#endif  // if KADEBUG
  int nmrimp = 4;
  vector<CalcType> imptime;
  imptime.resize(nmrimp);
  double calcend = 4.0*bclmin;
  if (bclmin < 1.0) {
    nmrimp += (int)(log(bclmin)/log(basis)+1.0);
    imptime.resize(nmrimp);
    double addtime = 0.0;
    for (int i = 0; i < nmrimp-4; i++) {
      imptime[i] = addtime;
      addtime   += pow(basis, (CalcType)i);
    }
    calcend = addtime+4.0*bclmin;
  }

  imptime[nmrimp-4] = calcend-4.0*bclmin;
  imptime[nmrimp-3] = calcend-3.0*bclmin;
  imptime[nmrimp-2] = calcend-2.0*bclmin;
  imptime[nmrimp-1] = calcend-bclmin;

  CalcType istim = (*PCEM)->GetAmplitude(), dVm = 0.0;
  const CalcType defStimTime = (*PCEM)->GetStimTime();
  *Vm_loc = (*PCEM)->GetVm();
  for (double t = 0; t < calcend; t += tinc) {
    CalcType i_external = 0.0;
    for (int j = 0; j < nmrimp; j++) {
      i_external += ((t >= imptime[j] && t <= imptime[j]+defStimTime) ? istim : 0.0);
    }
    dVm = (*PCEM)->Calc(tinc, *Vm_loc, i_external);
    CalcType calcium_loc = (*PCEM)->GetCai()*1000.0;
    *Force_loc = (*PCFM)->Calc(tinc, 1.0, .0, calcium_loc);
    *Vm_loc   += dVm;
  }
#if KADEBUG
  cerr<<"PreCalcModels\n";
#endif  // if KADEBUG
}  // PreCalcModels

#endif  // ifndef CELLMODELLAYER_H
