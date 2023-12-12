/*
 * File: acltCellModel.h
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


#ifndef ACLTCELLMODEL_H
#define ACLTCELLMODEL_H

#include <CellModelLayer.h>
#include <PETScLSE.h>

class CMS {
 public:
  CMS();
  ~CMS() {clear();}

  void clear();

  int material;                    //!< material class
  string emd;                      //!< The cell model corresponding to material
  string fmd;                      //!< The force model corresponding to material
  vbElphyParameters<double> *pep;  //! the used elphy model parameter, will be initialized before usage
  vbForceParameters<double> *pfp;  //! the used force model parameter, will be initialized before usage
  vbElphyModel<double> *prepem;    //! precalculation elphy model
  vbForceModel<double> *prepfm;    //! precalculation force model
  ElphyModelType emt;              //! elphymodeltype (from Cellmodelbasis.h: enum CellModelType)
  ForceModelType fmt;              //! forcemodeltype (from Cellmodelbasis.h: enum CellModelType)
  double VmInit;                   //! Vm after precalculation of the used elphy models
  double ForceInit;                //! Force after precalculation of the used force models
};


class CellModelStruct {
 public:
  CellModelStruct() {}

  ~CellModelStruct() {}

  vector<CMS> cellmodels;  //!< One cell model struct for each tissue class
  void LoadCellModels(string);
  void InitParameters(double);
  void setElphyModel(vbElphyModel<double> **, int);
  void setForceModel(vbForceModel<double> **, int);
};


class EMCoupling {
 public:
  //! function pointer for the following two functions
  typedef double (*CalculateCoupling)(vbElphyModel<double> *, vbForceModel<double> *, double, double, double, double &,
                                      double, double);
  typedef double (*GetParameters)(vbElphyModel<double> *, string &, int &);
  typedef double (*GetValues)(vbElphyModel<double> *, int, double &, double, double);

  CalculateCoupling *CouplingMethod;
  GetParameters *GetResultParameters;
  GetValues *GetResultValues;

  EMCoupling() {
    CouplingMethod      = NULL;
    GetResultParameters = NULL;
    GetResultValues     = NULL;
  }

  ~EMCoupling() {Delete();}

  void Delete();
  void New(unsigned long);
  void Set(int, vbElphyModel<double> *, vbForceModel<double> *, bool);

  /*!
   * In the function the calculation of electrophysiology is done without force calculation. This function belongs to a
   * function pointer
   */
  static double NoCoupling(vbElphyModel<double> *pem, vbForceModel<double> *pfm, double tinc, double stretch,
                           double velocity, double &dVm, double Vm, double i_external) {
    dVm = pem->Calc(tinc, Vm, i_external, stretch, 1);
    double Calcium = pem->GetCai()*1000.0;
    return Calcium;
  }

  static double GetResPara(vbElphyModel<double> *pem, string &resName, int &resPos) {
    // Funktion liest den String der Parameternamen des Zellmodells aus und sucht nach der einem angegebenen Parameter
    // und dessen Position im String
    bool foundParameter = false;

    vector<string> NameString;
    NameString.push_back("t");
    NameString.push_back("V_m");
    pem->GetLongParameterNames(NameString);
    int length = NameString.size();
    for (int i = 0; i < length; i++) {
      int compareNames = NameString[i].compare(resName);
      if (compareNames == 0) {
        foundParameter = true;
        resName        = NameString[i];
        resPos         = i;
        break;
      }
    }

    // gibt 1 zurÃ¼ck wenn der Parameter gefunden wurde, ansonsten 0
    if (foundParameter == true)
      return 1.;
    else
      return 0.;
  }  // GetResPara

  static double GetResVal(vbElphyModel<double> *pem, int resPos, double &resValue, double calctime, double Vm) {
    // Funktion liest den String der Parameterwerte des Zellmodells aus und sucht den Wert an der angegebenen Position
    // im String
    double tempValue = 0;
    ostringstream ValueString(ostringstream::out);

    pem->LongPrint(ValueString, calctime, Vm);
    ValueString<<endl;
    ValueString << '\0';
    istringstream ParameterValue(ValueString.str().c_str());
    for (int k = 0; k < resPos; k++) {
      ParameterValue >> tempValue;
    }
    ParameterValue >> resValue;
    return .0;
  }

  /*!
   * In the function the calculation of electrophysiology and force with Calcium bound to Troponin coupling is done.
   * This function belongs to a function pointer
   */
  static double TroponinCoupling(vbElphyModel<double> *pem, vbForceModel<double> *pfm, double tinc, double stretch,
                                 double velocity, double &dVm, double Vm, double i_external) {
    dVm = pem->Calc(tinc, Vm, i_external, stretch, 1);
    double TroponinCa = pem->GetTCa()*1000.0;
    return pfm->CalcTrop(tinc, stretch, velocity, TroponinCa);
  }

  /*!
   * In the function the calculation of electrophysiology and force with Calcium bound to Troponin coupling is done.
   * This function belongs to a function pointer
   */
  static double CalciumCoupling(vbElphyModel<double> *pem, vbForceModel<double> *pfm, double tinc, double stretch,
                                double velocity, double &dVm, double Vm, double i_external) {
    dVm = pem->Calc(tinc, Vm, i_external, stretch, 1);
    double Calcium = pem->GetCai()*1000.0;
    double Force   = pfm->Calc(tinc, stretch, velocity, Calcium);

    // pem->SetCai(Calcium*.001);
    return Force;
  }
};  // class EMCoupling


class HCM {
 public:
  HCM() {filename = ""; variable = "";}

  ~HCM() {}

  string filename;  //!< The name of the heterogeneous file
  string variable;  //!< Variable to which the heterogeneous file should be assigned
};


class HeterogeneCellModel {
 public:
  HeterogeneCellModel() {}

  ~HeterogeneCellModel() {}

  vector<HCM> hetCM;
  void Load(string);
  void Set(vbElphyModel<double> **, PetscInt);
};


#endif  // ifndef ACLTCELLMODEL_H
