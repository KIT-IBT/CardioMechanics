/*
 * File: acCELLerate.h
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


#ifndef ACCELERATE_H
#define ACCELERATE_H

#include <kaMachineOS.h>
#include <timer.h>

typedef int32_t IndexType;

#include "PETScLSE.h"
#include "acltCellModel.h"
#include "acltConditions.h"
#include "acltProjectFile.h"
#include "acltSensors.h"
#include "acltTime.h"
#include <map>

enum struct CurrentScheme {
  ICI,  /* Cell models at nodes, integrated using mass matrix */
  SVI,  /* Cell models at Gauss points, state variables interpolated from nodes */
  Godunow,  /* Cell models at nodes, Godunow splitting, default */

  // Strang /* Cell models at nodes, symmetric Strang splitting (evaluated twice) */
};


namespace ResultVars {
enum {
  Ii = 0, Ie, If, Ve, Vem, Vf, Vm, AT, numResultVars
};
}

// TODO godunow/strang splitting could be implemented as theta-scheme, however,
// theta != .5 / theta != 0 does not usually make sense...
// But well, that's also true for the other theta scheme...

/*!
   Classs to handle all cardiac simulations.
 */

class acCELLerate : public ACLTConditions, public ACLTSensors, public CellModelStruct {
  PetscInt StartCells;                //!< Start point of all cell-based vectors for the processor working on a part of
  //!< the problem
  PetscInt EndCells;                  //!< End point of all cell-based vectors for the processor working on a part of
  //!< the problem
  PetscInt StartNodes;                //!< Start point of all node-based vectors for the processor working on a part of
  //!< the problem
  PetscInt EndNodes;                  //!< End point of all node-based vectors for the processor working on a part of
  //!< the problem
  PetscInt StartExtraNodes;           //!< Start point of local portion of extracellular problem
  PetscInt EndExtraNodes;             //!< End point of local portion of extracellular problem
  vbElphyModel<double> **pElphyIntra;  //!< Pointer on the intracellular electrophysiological equations
  vbElphyModel<double> **pElphyFibro;  //!< Pointer on the fibroblast electrophysiological equations
  vbForceModel<double> **pForceIntra;  //!< Pointer on the intracellular force equations
  vbForceModel<double> **pForceFibro;  //!< Pointer on the fibroblast force equations

  EMCoupling EMIntra;                 //!< Electromechanical coupling for intracellular space
  EMCoupling EMFibro;                 //!< Electromechanical coupling for fibroblast space

  HeterogeneCellModel HeteroCMIntra;  //!< Heterogeneous cell model for intracellular space
  HeterogeneCellModel HeteroCMFibro;  //!< Heterogeneous cell model for fibroblast space

  string resultprefix;                //!< prefix of result
  string protocolfile;                //!< name of protocol file
  string backuploadname;              //!< name of the backup file that is loaded
  string backupsavename;              //!< name of the backup file that is saved
  string Material;                    //!< material vector name
  string MaterialFibro;               //!< material fibro vector name
  string matrixextraname;             //!< filename of the extracellular matrix
  string acCELLerateBUVersion;        //!< backup version string and number
  acltTime currenttime;  //!< current time of simulation [s]
  acltTime starttime;                 //!< current time of simulation [s]
  acltTime calclen;                   //!< length of the calculation [s]
  acltTime dtcell;                    //!< time step for cell model calculation [s]
  acltTime dtintra;                   //!< time step for intracellular current calculation [s]
  acltTime dtextra;                   //!< time step for extracellular current calculation [s]
  acltTime beginsave;                 //!< start time for saving the results [s]
  acltTime dtsave;                    //!< time step for saving the results [s]
  acltTime dtbackup;                  //!< time step for backup [s]
  double vm_def_range;                //!< Upper limit where transmembrane voltage is defined [V]
  PetscBool verbose;                  //!< verbose mode
  PetscBool combinedExtra;            //!< Flag if Extra Matrix was combined alreay

  string resultName[64];              //!< Array with the result names
  int resultPos[64][256];             //!< [maxParameters][maxMaterial]
  double resultValue[64];             //!< Result Values for each calculation step
  bool saveresultLattice[64];         //!< Also save lattices for the result Parameters?
  bool setResults;                    //!< Save result Parameters?
  bool exportResults;
  int numberParameters;               //!< Number Parameters
  int maxParameters;                  //!< Maximum Parameters
  string results;                     //!< Result string in aclt file
  bool parametersinSensorfile;        //!< CellmodelParameters in Sensorfile?
  const char *saveLattice_Name[8];    //!< Names of Parameters in the Sensorfiles
  bool saveLattice[8];                //!< Save Lattices of the Parameters in the Sensorfiles

  FILE *fprot;                        //!< Pointer to the protocol file
  char prottext[256];                 //!< Name of the protocol file

  PETScLSE Intra;                     //!< An object of PETScLSE for the intracellular space to hold and calculate on
  //!< vectors and matrices
  PETScLSE Extra;                     //!< An object of PETScLSE for the extracellular space to hold and calculate on
  //!< vectors and matrices
  PETScLSE Fibro;                     //!< An object of PETScLSE for the fibroblast space to hold and calculate on
  //!< vectors and matrices
  PETScLSE CurrentDensity;            //!< LSE solver for intercellular currents when using mass matrix

  PetscScalar CmMyo;                  //!< Myocyte membrane capacitynce per unit area [F/m^2]
  PetscScalar BetaMyo;                //!< Myocyte surface to volume ratio [1/m]
  PetscScalar BetaMyoFib;             //!< Number of Myo-Fibro gap junctions [1/m^3]
  PetscScalar RMyoFib;                //!< Resistor of single myo-fibro gap junction [Ohm]
  PetscScalar VolumeMyo;              //!< Relative volume of intracellular space of myocytes
  PetscScalar VolumeFibro;            //!< Relative volume of intracellular space of fibroblasts
  PetscScalar VolumeExtra;            //!< Relative volume of extracellular space
  Vec VecCmMyo;                       //!< Vector for myocyte membrane capacitynce per unit area [F/m^2]
  Vec VecBetaMyo;                     //!< Vector for myocyte surface to volume ratio [1/m]
  Vec VecVolMyo;                      //!< A vector to hold the distribution of intracellular volume
  Vec VecVolFibro;                    //!< A vector to hold the distribution of extracellular volume
  Vec VecVolExtra;                    //!< A vector to hold the distribution of fibroblasts volume
  Vec VecSaveTmp;                     //!< Temporary vector to hold results to be written in case of IS
  CurrentScheme currentScheme;
  IS IntraIndexSet;
  std::string IntraIndexSetFile;
  std::map<PetscInt, PetscInt> GlobalToIntraMap;

  bool implicit;              //!< true if implicit calculation is selected
  int jacobi;                 //!< number of Jacobi iterations during implicit calculation
  double thetaIntra;          //!< Theta for intra solving scheme [0.0..1.0]
  PetscScalar tinc;           //!< time step for implicit calculation [s]
  Mat MatImp;                 //!< A Matrix (I+dt*A) for implicit calculation
  Mat IntraRHS;               //!< Right-hand side matrix for intracellular domain
  Mat IntraMass;              //!< Mass matrix for intracellular domain
  Mat ExtraRHS;               //!< Right-hand side matrix for extracellular domain
  Vec VecDiag;                //!< A Vector that contains the diagonal elements of (I+dt*A)
  Vec VecVCell;               //!< A Vector to hold the change of voltage dVm of the cell calculations
  Vec VecICell;
  Vec VecRes;                 //!< A Vector to hold the residuum of the Jacobi-method
  Vec VecVRes;                //!< A Vector to hold the iterative solution of dVm of the Jacobi-method
  Mat GaussFwd, GaussInt;     //!< Intracellular current interpolation matrices

  Vec ActivationTime;         //!< Vector to hold activation times
  double ActivationThreshold;  //!< Threshold for cell activation

  Vec VecResultValues;        //!< A vector to hold the result values
  Vec MaterialV;              //!< A vector to hold the material class of each voxel

  Vec Stim;                   //!< Vector holding extracellular potential boundary condition information
  Vec Force;                  //!< A vector to hold the force distribution of intracellular volume

  int domains;                //!< 1 = Monodomain; 2 = Bidomain; 3 = Tridomain
  int mpirank;                //!< The rank of the executing process
  int mpisize;                //!< The amount of processes that are involved in the calculation
  timer t;                    //!< timer to calculate the time for different steps inside the code
  timer t_mono_matrix, t_mono_cond, t_mono_cell;
  time_t t_start, t_now;

  acltProjectFile aPF;        //!< Object of acCELLerate project file to handle loading project information
  const char *VecFnTemplate;  //!< Template for output vector filenames

  bool forceset;  //!< true if force is selected
  void initTimeSteps(acltTime &, acltTime &, acltTime &, acltTime &, acltTime * = NULL);

 public:
  acCELLerate(void);       //!< Constructor to initialize
  virtual ~acCELLerate();  //!< Destructor to free

  //! Function to print debug information synchronized for all involved processes
  inline void PrintDebug(const char *format, ...) {
#if PSDEBUG > 0
    char errstr[256];
    va_list args;
    va_start(args, format);
    vsnprintf(errstr, sizeof(errstr), format, args);
    va_end(args);

    PetscSynchronizedFPrintf(PETSC_COMM_WORLD, stderr, "acCELLerate::%s rank %d\n", errstr, mpirank);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#else  // if PSDEBUG > 0
#endif  // if PSDEBUG > 0
  }

  //! Function to print debug information synchronized for all involved processes with a higher level of debug
  //! information
  inline void PrintDebug2(const char *format, ...) {
#if PSDEBUG == 2
    PrintDebug(format);
#else  // if PSDEBUG == 2
#endif  // if PSDEBUG == 2
  }

  void WriteProtocol(string, int type = 0);  //!< Writing information to terminal and/or protocol file
  void LoadProject(const char *);

  void InitMono(Vec = NULL);
  void InitBi();
  void InitTri();
  void BuildIntraOperators(void);         //!< Build intra/mono operators from stiffness&mass matrices
  void SetIntraMatrices(Mat, Mat = NULL);  //!< Set intra/mono stiffness and mass matrix
  void SetConditionStatus(ConditionType ct1, ConditionType ct2, bool redo);
  void ApplyConditionsToVec(PetscScalar *px, PetscScalar *pb, ConditionType ct1, ConditionType ct2);
  void ApplyConditions(PetscScalar *px, PetscScalar *pb, ConditionType ct1, ConditionType ct2, bool redo);
  void ApplyConditionsExtra();
  void InitSensors();
  bool IsFloat(const char *s);
  void SaveSensors(PetscScalar *pix, PetscScalar *pib, PetscScalar *pex = NULL, PetscScalar *peb = NULL,
                   PetscScalar *pfx                                     = NULL, PetscScalar *pfb = NULL);
  void SaveParameterSensors(PetscScalar *pres, string currentPara);
  void SaveResultIntra();
  void SaveResultExtra();
  void SaveResultFibro();
  void SaveResultParameters(string);

  void SetLoadBackup(string);
  void SetSaveBackup(string);
  bool ReadBackup(vbElphyModel<double> **, vbForceModel<double> **, PETScLSE *, bool firstdomain = true);  //!< Read
  //!< Backup if it matches with input data
  void   WriteBackup(vbElphyModel<double> **, vbForceModel<double> **, PETScLSE *, bool firstdomain = true);
  void   SendInteger(int);
  int    ReceiveInteger();
  void   SendDouble(double);
  double ReceiveDouble();
  void SendTime(acltTime);  // TODO const&?
  acltTime ReceiveTime();
  void     Run();
  void     CalcCellModels(PetscScalar const *, PetscScalar const *, PetscScalar *, PetscScalar const * = NULL,
                          PetscScalar const *                                                          = NULL,
                          PetscScalar *                                                                = NULL);
  void SaveCellModelVariables(PetscScalar const *, acltTime);
  void IntraStep(acltTime &, PetscScalar * = NULL, PetscScalar * = NULL);
  void MonoDomain(void);
  void MonoDomain(acltTime, Vec = NULL, Vec = NULL);
  void BiDomain();
  void TriDomain();

  Vec GetVmVec() {return Intra.X;}

  Vec GetCurrentVec() {return Intra.B;}

  Vec GetPhiEVec() {return Extra.X;}

  Vec GetForceVec() const {return Force;}

  Mat GetSystemMatrix() {return Intra.A;}

  Mat GetMassMatrix() {return IntraRHS;}

  void PrintTime(void);
};  // class acCELLerate

#endif  // ifndef ACCELERATE_H

/*! \page acCELLerate acCELLerate
   Calculates 1-3 domain models of electrical conduction in cardiac tissue.

   \section SYNOPSIS_acCELLerate SYNOPSIS
   acCELLerate \<project file\>\n

   \section OPTIONS_acCELLerate OPTIONS
   \param "<project file>" The procejt file. Just run acCELLerate without options to get the possible entries with
      explanation

   \section DESCRIPTION_acCELLerate DESCRIPTION
   To be done...

   \section SOURCE_acCELLerate SOURCE
   acCELLerate.cpp acCELLerate.h

   \section SEEALSO_acCELLerate SEE ALSO
   \ref xCELLent

   \section CHANGELOG_acCELLerate CHANGELOG
   V1.0.0 - 31.07.2008 (Gunnar Seemann): Starting with changelog\n
   V1.0.1 - 20.04.2009 (Meike Karl): Bugfix in acCELLerate::ReadBackup, now backup file is checked if compatible with
      input data. Attention: acCELLerate executed with more than one process, exception handling in ReadBackup does not
      work properly (no output written). \n
   V1.1.0 - 04.12.2009 (Christian Rombach): Model Parameter are now possible to exctract each. Semi-implicite solver was
      added. Conditions and sensors have new convention of names. Sensors is now in array format (not in Lattice format
      anymore)\n
 */
