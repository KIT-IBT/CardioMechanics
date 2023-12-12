/*
 * File: acCELLerate.cpp
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


#include "acCELLerate.h"
#include <algorithm>


static const KSPType KSPTYPE  = KSPGMRES;
static const PCType  PCTYPE   = PCBJACOBI;
static const PetscReal RTOL   = 1e-6;
static const PetscReal ATOL   = 1e-6;
static const PetscReal DTOL   = 10000;
static const PetscInt  MAXITS = 1000;

/*!
 Constructor to initialize MPI parameters and to set all values on standard
 */
acCELLerate::acCELLerate() {
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank);
    MPI_Comm_size(PETSC_COMM_WORLD, &mpisize);
    
    PrintDebug("acCELLerate()");
    
    pElphyIntra = NULL;
    pElphyFibro = NULL;
    pForceIntra = NULL;
    pForceFibro = NULL;
    
    MatImp            = NULL;
    VecDiag           = NULL;
    VecVCell          = NULL;
    VecICell          = NULL;
    VecRes            = NULL;
    VecVRes           = NULL;
    IntraRHS          = NULL;
    IntraMass         = NULL;
    ExtraRHS          = NULL;
    GaussFwd          = NULL;
    GaussInt          = NULL;
    IntraIndexSet     = NULL;
    IntraIndexSetFile = "";
    VecSaveTmp        = NULL;
    MaterialV         = NULL;
    
    ActivationTime      = NULL;
    ActivationThreshold = .0;
    
    thetaIntra = -42.;
    
    CmMyo      = .0; // Default: cell model specific
    BetaMyo    = .0; // Default: cell model specific
    VecCmMyo   = NULL;
    VecBetaMyo = NULL;
    
    RMyoFib     = 1e7;
    BetaMyoFib  = 1e7;
    VolumeMyo   = 1.;
    VolumeExtra = 1.;
    VolumeFibro = 0.0;
    
    resultprefix         = "Res";
    backuploadname       = "";
    acCELLerateBUVersion = "acCELLerateBackupVersion1.0";
    backupsavename       = "";
    protocolfile         = "";
    fprot                = NULL;
    exportResults        = true;
    
    calclen     = 1_s;
    dtcell      = 1e-5_s;
    dtintra     = 1e-5_s;
    dtextra     = 1e-5_s;
    beginsave   = 0_s;
    dtsave      = 1e-3_s;
    dtbackup    = 1e-1_s;
    currenttime = 0_s;
    starttime   = 0_s;
    
    setResults             = false;
    parametersinSensorfile = false;
    numberParameters       = 0;
    maxParameters          = 64;
    vm_def_range           = RangeTabhalf/1000.0;
    
    saveLattice_Name[ResultVars::Ii]  = "Ii";
    saveLattice_Name[ResultVars::Ie]  = "Ie";
    saveLattice_Name[ResultVars::If]  = "If";
    saveLattice_Name[ResultVars::Ve]  = "Ve";
    saveLattice_Name[ResultVars::Vem] = "Vem";
    saveLattice_Name[ResultVars::Vf]  = "Vf";
    saveLattice_Name[ResultVars::Vm]  = "Vm";
    saveLattice_Name[ResultVars::AT]  = "AT";
    
    domains       = 1;
    verbose       = PETSC_FALSE;
    combinedExtra = PETSC_FALSE;
    forceset      = false;
    implicit      = false;
    currentScheme = CurrentScheme::Godunow;
    
    // VecFnTemplate = "%s_%011.6lf_%s.vec";
    VecFnTemplate = "%s_%s_%s.vec";
    
    PrintDebug("acCELLerate finished");
}

/*!
 Destructor to free memory and close protocol file
 */
acCELLerate::~acCELLerate() {
    PrintDebug("~acCELLerate()");
    
    delete[] pElphyIntra;
    delete[] pElphyFibro;
    delete[] pForceIntra;
    delete[] pForceFibro;
    
    if (fprot)
        fclose(fprot);
    
    PetscErrorCode ierr;
    if (MatImp) {
        ierr = MatDestroy(&MatImp); CHKERRQ(ierr);
    }
    if (VecDiag) {
        ierr = VecDestroy(&VecDiag); CHKERRQ(ierr);
    }
    if (VecVCell) {
        ierr = VecDestroy(&VecVCell); CHKERRQ(ierr);
    }
    if (VecRes) {
        ierr = VecDestroy(&VecRes); CHKERRQ(ierr);
    }
    if (VecVRes) {
        ierr = VecDestroy(&VecVRes); CHKERRQ(ierr);
    }
    if (VecCmMyo) {
        ierr = VecDestroy(&VecCmMyo); CHKERRQ(ierr);
    }
    if (VecBetaMyo) {
        ierr = VecDestroy(&VecBetaMyo); CHKERRQ(ierr);
    }
    
    if (IntraMass) {
        ierr = MatDestroy(&IntraMass); CHKERRQ(ierr);
    }
    if (IntraRHS) {
        ierr = MatDestroy(&IntraRHS); CHKERRQ(ierr);
    }
    if (ExtraRHS) {
        ierr = MatDestroy(&ExtraRHS); CHKERRQ(ierr);
    }
    if (IntraIndexSet) {
        ierr = ISDestroy(&IntraIndexSet); CHKERRQ(ierr);
        if (VecSaveTmp) {
            ierr = VecDestroy(&VecSaveTmp); CHKERRQ(ierr);
        }
    }
    
    if (GaussFwd) {
        ierr = MatDestroy(&GaussFwd); CHKERRQ(ierr);
    }
    if (GaussInt) {
        ierr = MatDestroy(&GaussInt); CHKERRQ(ierr);
    }
    
    if (ActivationTime) {
        ierr = VecDestroy(&ActivationTime); CHKERRQ(ierr);
    }
}

/*!
 Writing the protocol information to terminal and/or protocol file (only performed for main process)
 \param ptext protocol text to be printed
 \param type only of type > 0 this information is printed to terminal
 */
void acCELLerate::WriteProtocol(string ptext, int type) {
    if (mpirank == 0) {
        if (verbose || type) {
            int len = 0;
            printf("\r%s%n", ptext.c_str(), &len);
            for (int i = len; i < 60; ++i) {
                putchar(' ');
            }
            putchar('\n');
        }
        if (fprot)
            fprintf(fprot, "%s\n", ptext.c_str());
    }
}

void acCELLerate::LoadProject(const char *projectFile) {
    PrintDebug("LoadProject()");
    PetscErrorCode ierr;
    
    ifstream proj(projectFile);
    if (!proj.is_open())
        throw kaBaseException("Project file %s does not exist", projectFile);
    
    string token, args;
    while (!proj.eof()) {
        proj>>token;
        getline(proj, args);
#if KADEBUG
        if (!mpirank) {
            std::cerr << "Token: " << token;
        }
#endif  // if KADEBUG
        if (token[0] == '#')
            continue;
        while (std::isspace(args[0]))
            args = args.substr(1);
        
        // trim comments
        args = args.substr(0, args.find_first_of('#'));
        
        // trim trailing whitespace
        args.erase(std::find_if(args.rbegin(), args.rend(), std::not1(std::ptr_fun<int, int>(
                                                                                             std::isspace))).base(), args.end());
        
#if KADEBUG
        if (!mpirank) {
            std::cerr <<  " - Args: " << args << "\n";
        }
#endif  // if KADEBUG
        
        switch (aPF.getIndex(token)) {
            case 0:
                resultprefix = args; break;
            case 1:
                ACLTConditions::Load(args.c_str()); break;
            case 2:
                ACLTSensors::Load(args.c_str()); break;
            case 3:
                calclen = args; break;
            case 4:
                dtcell = args; break;
            case 5:
                dtextra = args; break;
            case 6:
                dtintra = args; break;
            case 7:
                beginsave = args; break;
            case 8:
                dtsave = args; break;
            case 9:
                CellModelStruct::LoadCellModels(args); break;
            case 10:
                this->Material = args; break;
            case 11:
                this->MaterialFibro = args; break;
            case 12:
                Intra.LoadA(args.c_str()); break;
            case 13:
                if (!combinedExtra) {matrixextraname = args; Extra.LoadA(matrixextraname.c_str());}
                break;
            case 14:
                matrixextraname = args; Extra.LoadA(matrixextraname.c_str()); combinedExtra = PETSC_TRUE; break;
            case 15:
                Fibro.LoadA(args.c_str()); break;
            case 16:
                BetaMyoFib = atof(args.c_str()); break;
            case 17:
                RMyoFib = atof(args.c_str()); break;
            case 18:
                if (IsFloat(args.c_str())) {VolumeMyo = atof(args.c_str());} else {
                    LoadVec(args.c_str(), VecVolMyo, ot_bin); VolumeMyo = 0.0;
                }
                break;
            case 19:
                if (IsFloat(args.c_str())) {VolumeExtra = atof(args.c_str());} else {
                    LoadVec(args.c_str(), VecVolExtra, ot_bin); VolumeExtra = 0.0;
                }
                break;
            case 20:
                if (IsFloat(args.c_str())) {VolumeFibro = atof(args.c_str());} else {
                    LoadVec(args.c_str(), VecVolFibro, ot_bin); VolumeFibro = 0.0;
                }
                break;
            case 21:
                if (!args.compare("Bi"))
                    domains = 2;
                else if (!args.compare("Tri"))
                    domains = 3;
                break;
            case 22:
                verbose = PETSC_TRUE; break;
            case 23:
                protocolfile = args; break;
            case 24:
                HeteroCMIntra.Load(args); break;
            case 25:
                HeteroCMFibro.Load(args); break;
            case 26:
                SetLoadBackup(args); break;
            case 27:
                SetSaveBackup(args); break;
            case 28:
                dtbackup = args; break;
            case 29:
                results = args.c_str(); setResults = true; break;
            case 30:
                jacobi = atoi(args.c_str()); implicit = true; break;
            case 31:
                forceset = true; break;
            case 32: {  // MassMatrixIntra
                        // Set IntraRHS to mass matrix, will be used in BuildMonoOperators()
                        // to create theta-scheme LHS and RHS matrices.
                PetscViewer fd = NULL;
                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, args.c_str(),
                                             FILE_MODE_READ, &fd); CHKERRQ(ierr);
                ierr = MatCreate(PETSC_COMM_WORLD, &IntraRHS); CHKERRQ(ierr);
                ierr = MatSetType(IntraRHS, MATMPIAIJ); CHKERRQ(ierr);
                ierr = MatLoad(IntraRHS, fd); CHKERRQ(ierr);
                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
            }
                break;
            case 33:  // ThetaIntra: Theta for intra PDE solver scheme
                thetaIntra = atof(args.c_str());
                if ((thetaIntra < .0) || (thetaIntra > 1.0) ) {
                    throw kaBaseException("ThetaIntra must be between 0.0 and 1.0!");
                }
                break;
            case 34:  // C_m
                if (IsFloat(args.c_str())) {CmMyo = atof(args.c_str());} else {
                    LoadVec(args.c_str(), VecCmMyo, ot_bin); CmMyo = 0.0;
                }
                break;
            case 35:  // Beta
                if (IsFloat(args.c_str())) {BetaMyo = atof(args.c_str());} else {
                    LoadVec(args.c_str(), VecBetaMyo, ot_bin); BetaMyo = 0.0;
                }
                break;
            case 36: {  // Gauss interpolation and integration matrices
                PetscViewer fd  = NULL;
                const auto  i   = args.find(" ");
                std::string n2c = args.substr(0, i);
                std::string c2n = args.substr(i+1);
                
                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, n2c.c_str(),
                                             FILE_MODE_READ, &fd); CHKERRQ(ierr);
                ierr = MatCreate(PETSC_COMM_WORLD, &GaussFwd); CHKERRQ(ierr);
                ierr = MatSetType(GaussFwd, MATMPIAIJ); CHKERRQ(ierr);
                ierr = MatLoad(GaussFwd, fd); CHKERRQ(ierr);
                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
                
                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, c2n.c_str(),
                                             FILE_MODE_READ, &fd); CHKERRQ(ierr);
                ierr = MatCreate(PETSC_COMM_WORLD, &GaussInt); CHKERRQ(ierr);
                ierr = MatSetType(GaussInt, MATMPIAIJ); CHKERRQ(ierr);
                ierr = MatLoad(GaussInt, fd); CHKERRQ(ierr);
                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
            }
                break;
            case 37:
                ActivationThreshold = std::atof(args.c_str());
                break;
            case 38:  /* Splitting / current interpolation method */
                if (args == "ICI") {
                    currentScheme = CurrentScheme::ICI;
                } else if (args == "SVI") {
                    currentScheme = CurrentScheme::SVI;
                } else if (args == "Godunow") {
                    currentScheme = CurrentScheme::Godunow;
                } else {
                    throw kaBaseException("Unknown current interpolation method: %s", args.c_str());
                }
                break;
            case 39:  // IntraIndexSet
                IntraIndexSetFile = args;
                break;
            case 40:  /* Compress */
                
                // VecFnTemplate = "%s_%011.6lf_%s.vec.gz";
                VecFnTemplate = "%s_%s_%s.vec.gz";
                break;
            case 41: // disable Res export
                exportResults = false;
                break;
            case -1:
            default:
                WriteProtocol("acCELLerate::LoadProject Unknown option: "+token+" in project file "+projectFile, 1); break;
        }  // switch
    }
    
    /* Assign default value to thetaIntra depending on FD or FE case */
    if (thetaIntra == -42.) {
        if (IntraRHS) {
            thetaIntra = .5;
        } else {
            thetaIntra = .0;
        }
    }
    
    /* Sanity checks */
    
    if ((dtcell > dtintra) || ((domains == 2) && (dtcell > dtextra))) {
        throw kaBaseException("DTIntra must be less than or equal to DTIntra and DTExtra.");
    }
    
    if ((dtintra % dtcell != 0_s) || ((domains == 2) && (dtextra % dtcell != 0_s))) {
        throw kaBaseException("DTIntra and DTExtra must be multiples of DTCell.");
    }
    
    switch (currentScheme) {
        case CurrentScheme::SVI:
            if (!GaussFwd || !GaussInt) {
                throw kaBaseException("Gauss matrices needed for selected current scheme.");
            }
            break;
        case CurrentScheme::ICI:
            if (!IntraRHS) {
                throw kaBaseException("ICI requires Mass matrix.");
            }
        default:
            if (GaussFwd || GaussInt) {
                WriteProtocol("acCELLerate::LoadProject: WARNING: Gauss matrices will be ignored with current scheme.", 1);
                if (GaussFwd) {
                    ierr = MatDestroy(&GaussFwd); CHKERRQ(ierr);
                }
                if (GaussInt) {
                    ierr = MatDestroy(&GaussInt); CHKERRQ(ierr);
                }
                GaussFwd = GaussInt = NULL;
            }
    }
    
    if (domains == 3) {
        if ((currentScheme != CurrentScheme::Godunow) || IntraRHS || (thetaIntra != .0)) {
            throw kaBaseException("TriDomain can not be used with MassMatrix, Theta != 0. or SVI/ICI.");
        }
    }
    
    CellModelStruct::InitParameters((double)dtcell);
    
    PrintDebug("LoadProject() finished");
}  // acCELLerate::LoadProject

void acCELLerate::BuildIntraOperators() {
    PrintDebug("BuildIntraOperators()");
    PetscErrorCode ierr;
    
    if (VolumeMyo != 0.0) {
        ierr = MatScale(Intra.A, VolumeMyo); CHKERRQ(ierr);
    } else {
        Intra.MatVecScale(VecVolMyo);
    }
    
    if (implicit) {
        PetscScalar *pvolc, *pvol;
        PetscInt Istart, Iend;
        PetscInt ncols;
        const PetscInt *cols;
        const PetscScalar *vals;
        tinc = PetscScalar(dtcell);  // dt
        if (MatImp) {
            ierr = MatDestroy(&MatImp); CHKERRQ(ierr);
        }
        ierr = MatConvert(Intra.A, MATSAME, MAT_INITIAL_MATRIX, &MatImp); CHKERRQ(ierr);  // copy A
        if (!VecDiag) {
            ierr = VecDuplicate(Intra.X, &VecDiag); CHKERRQ(ierr);
        }
        ierr = VecSet(VecDiag, 0.0); CHKERRQ(ierr);
        
        // create vector for volume conversion
        ierr = VecGetArray(VecDiag, &pvolc); CHKERRQ(ierr);
        if (VolumeMyo == 0.0) {
            ierr = VecGetArray(VecVolMyo, &pvol); CHKERRQ(ierr);
        }
        for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
            PetscScalar volume_local = 0.0;
            if (VolumeMyo == 0.0)
                volume_local = (pvol[Ii-StartCells] < 1e-9 ? 1e-9 : pvol[Ii-StartCells]);
            else
                volume_local = VolumeMyo;
            
            pvolc[Ii-StartCells] = pElphyIntra[Ii-StartCells]->Volume()/volume_local*1e9;  // convert A/m^3 to nA/cell
        }
        if (VolumeMyo == 0.0) {
            ierr = VecRestoreArray(VecVolMyo, &pvol); CHKERRQ(ierr);
        }
        
        // multiply rows of A with volume factor
        ierr = MatGetOwnershipRange(Intra.A, &Istart, &Iend); CHKERRQ(ierr);
        int i = Istart;
        for (; i < Iend; i++) {
            ierr = MatGetRow(Intra.A, i, &ncols, &cols, &vals); CHKERRQ(ierr);
            vector<PetscInt> indizes(ncols);
            vector<PetscScalar> values(ncols);
            for (int j = 0; j < ncols; j++) {
                indizes[j] = cols[j];
                values[j]  = vals[j]*(*(pvolc+i-Istart));
                if (cols[j] == i) {
                    if (values[j] < 1e-9)
                        values[j] = 1e-9;
                }
            }
            ierr = MatRestoreRow(Intra.A, i, &ncols, &cols, &vals); CHKERRQ(ierr);
            
            // in rows: volume/volume_local*1e9*A
            ierr = MatSetValues(MatImp, 1, &i, ncols, indizes.data(), values.data(), INSERT_VALUES); CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(MatImp, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(MatImp, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        
        ierr = VecRestoreArray(VecDiag, &pvolc); CHKERRQ(ierr);
        
        ierr = MatScale(MatImp, tinc); CHKERRQ(ierr);                      // dt*A
        
        ierr = VecSet(VecDiag, 1.0); CHKERRQ(ierr);                        // I
        ierr = MatDiagonalSet(MatImp, VecDiag, ADD_VALUES); CHKERRQ(ierr);  // I+dt*A
        ierr = MatGetDiagonal(MatImp, VecDiag); CHKERRQ(ierr);             // Diagonal of I+dt*A
        ierr = VecReciprocal(VecDiag); CHKERRQ(ierr);                      // Diagonal^-1
        ierr = MatScale(MatImp, -1.0); CHKERRQ(ierr);                      // -(I+dt*A)
        if (!VecRes) {
            ierr = VecDuplicate(Intra.X, &VecRes); CHKERRQ(ierr);
        }
        ierr = VecSet(VecRes, 0.0); CHKERRQ(ierr);
    }
    
    // Initialization for intra PDE theta solver scheme
    
    // NOTATION
    // K: Stiffness matrix
    // M: Mass matrix
    // C_m: Cell membrane capacitance per unit area
    // beta: Cell surface to volume ratio
    // dt: Intra time inecrement
    
    if (this->IntraRHS || (thetaIntra != .0)) {
        if ((domains > 1) || saveLattice[ResultVars::Ii]) {
            // Save intra stiffness matrix for bidomain or Ii calculation
            if (ExtraRHS) {
                ierr = MatDestroy(&ExtraRHS); CHKERRQ(ierr);
            }
            ierr = MatDuplicate(Intra.A, MAT_COPY_VALUES, &ExtraRHS); CHKERRQ(ierr);
        }
        
        // Scale M with (C_m*beta)/dt
        // C_m and beta can be cell model specific
        Vec betaCm_dt;
        ierr = VecDuplicate(Intra.X, &betaCm_dt); CHKERRQ(ierr);
        ierr = VecSet(betaCm_dt, 1.0 / (double)dtintra);
        
        /* If capacitance has been set (vec or scalar), use it */
        if (VecCmMyo) {
            ierr = VecPointwiseMult(betaCm_dt, betaCm_dt, VecCmMyo); CHKERRQ(ierr);
        } else if (CmMyo != .0) {
            ierr = VecScale(betaCm_dt, CmMyo); CHKERRQ(ierr);
        } else {  /* Otherwise use values from the cell model */
            PetscInt istart, iend;
            PetscScalar *v;
            ierr = VecGetOwnershipRange(betaCm_dt, &istart, &iend); CHKERRQ(ierr);
            ierr = VecGetArray(betaCm_dt, &v); CHKERRQ(ierr);
            
            for (PetscInt i = 0; i < (iend - istart); ++i)
                v[i] *= pElphyIntra[i]->GetCm();
            
            ierr = VecRestoreArray(betaCm_dt, &v); CHKERRQ(ierr);
        }
        
        /* The same applies to teh surface-area-to-volume-ratio */
        if (VecBetaMyo) {
            ierr = VecPointwiseMult(betaCm_dt, betaCm_dt, VecBetaMyo); CHKERRQ(ierr);
        } else if (BetaMyo != .0) {
            ierr = VecScale(betaCm_dt, BetaMyo); CHKERRQ(ierr);
        } else {
            PetscInt istart, iend;
            PetscScalar *v;
            ierr = VecGetOwnershipRange(betaCm_dt, &istart, &iend); CHKERRQ(ierr);
            ierr = VecGetArray(betaCm_dt, &v); CHKERRQ(ierr);
            
            for (PetscInt i = 0; i < (iend - istart); ++i)
                v[i] *= pElphyIntra[i]->GetBeta();
            
            ierr = VecRestoreArray(betaCm_dt, &v); CHKERRQ(ierr);
        }
        
        ierr = VecAssemblyBegin(betaCm_dt); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(betaCm_dt); CHKERRQ(ierr);
        
        if (this->IntraRHS) {  /* a mass matrix has been set */
            if (saveLattice[ResultVars::Ii] || saveLattice[ResultVars::Ie]) {
                /* need to save Ii, so we need an unscaled mass
                 matrix */
                if (CurrentDensity.A) {
                    ierr = MatDestroy(&CurrentDensity.A); CHKERRQ(ierr);
                    ierr = VecDestroy(&CurrentDensity.X); CHKERRQ(ierr);
                    ierr = VecDestroy(&CurrentDensity.B); CHKERRQ(ierr);
                }
                ierr = MatDuplicate(IntraRHS, MAT_COPY_VALUES, &CurrentDensity.A); CHKERRQ(ierr);
                CurrentDensity.Set(Intra.xLattice, Intra.yLattice, Intra.zLattice);
                CurrentDensity.CreateB();
                CurrentDensity.CreateX();
                CurrentDensity.InitSolve(KSPTYPE, PCTYPE, RTOL, ATOL, DTOL, MAXITS);
            }
            
            // IntraRHS := (C_m*beta)/dt M - (1-theta) * K
            ierr = MatDiagonalScale(IntraRHS, betaCm_dt, NULL); CHKERRQ(ierr);
            if (currentScheme == CurrentScheme::ICI) {
                /* We need the so-scaled mass matrix (beta*C_m/dt*M) for ICI. */
                if (thetaIntra == .0) {
                    ierr      = PetscObjectReference((PetscObject)Intra.A); CHKERRQ(ierr);
                    IntraMass = Intra.A;
                } else if (thetaIntra == 1.0) {
                    ierr      = PetscObjectReference((PetscObject)IntraRHS); CHKERRQ(ierr);
                    IntraMass = IntraRHS;
                } else {
                    if (IntraMass) {
                        ierr = MatDestroy(&IntraMass); CHKERRQ(ierr);
                    }
                    ierr = MatDuplicate(IntraRHS, MAT_COPY_VALUES, &IntraMass); CHKERRQ(ierr);
                }
            }
            ierr = MatAXPY(this->IntraRHS, -(1.0-thetaIntra), Intra.A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        } else {  // assume A = Inv(M)*K and use identity matrix as mass matrix
                  // IntraRHS := A = K
            if (IntraRHS) {
                ierr = MatDestroy(&IntraRHS); CHKERRQ(ierr);
            }
            ierr = MatDuplicate(Intra.A, MAT_COPY_VALUES, &this->IntraRHS); CHKERRQ(ierr);
            
            // IntraRHS := -(1-theta) * K
            ierr = MatScale(IntraRHS, -(1.0-thetaIntra)); CHKERRQ(ierr);
            
            // IntraRHS = (Cm*beta)/dt * I - (1-theta) * K
            // ierr = VecScale(betaCm_dt, -1.0); CHKERRQ(ierr);
            ierr = MatDiagonalSet(IntraRHS, betaCm_dt, ADD_VALUES); CHKERRQ(ierr);
            
            // IntraRHS := -IntraRHS = -(-(beta*Cm/dt)+A) = (beta*Cm/dt)-A = (beta*Cm/dt) - (1-theta) * K
            // ierr = MatScale(IntraRHS, -1.0); CHKERRQ(ierr);
        }
        
        ierr = VecDestroy(&betaCm_dt); CHKERRQ(ierr);
        
        // Intra.A = dt/(C_m*beta) * K + IntraRHS == K + a*M - (1-theta)*K == theta*K + a*M
        ierr = MatAXPY(Intra.A, 1.0, IntraRHS, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        
        Intra.InitSolve(KSPTYPE, PCTYPE, RTOL, ATOL, DTOL, MAXITS);
        
        // Now, Intra.A  == C_m*beta/dt*M +     theta*K
        // and  IntraRHS == C_m*beta/dt*M - (1-theta)*K
        // and IntraMass == C_m*beta/dt*M
    } else if (domains > 1) {  /* and NOT (this->IntraRHS || thetaIntra != .0)  */
        ierr     = PetscObjectReference((PetscObject)Intra.A); CHKERRQ(ierr);
        ExtraRHS = Intra.A;
    }
    
    // TODO the steps from initbi are missing...
}  // acCELLerate::BuildIntraOperators

void acCELLerate::SetIntraMatrices(Mat stiffness, Mat mass) {
    PrintDebug("SetIntraMatrices()");
    PetscErrorCode ierr;
    if (Intra.A && (Intra.A != stiffness)) {
        ierr = MatDestroy(&Intra.A); CHKERRQ(ierr);
    }
    Intra.A = stiffness;
    
    if (IntraRHS && (IntraRHS != mass)) {
        ierr = MatDestroy(&IntraRHS); CHKERRQ(ierr);
    }
    IntraRHS = mass;
    
    BuildIntraOperators();
}

void acCELLerate::InitMono(Vec materials) {
    PrintDebug("InitMono()");
    
    PetscErrorCode ierr;
    
    Intra.CreateX();
    Intra.CreateB();
    
    ierr = VecGetOwnershipRange(Intra.X, &StartNodes, &EndNodes); CHKERRQ(ierr);
    
    PetscInt problemSize = domains == 1 ? Intra.size : Extra.size;
    ACLTConditions::Parse(problemSize);
    
    if (domains == 1) {
        for (auto pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end();) {
            if ((pSC->cflag == CT_II) || (pSC->cflag == CT_UI)) {
                if ((pSC->position < StartNodes) || (pSC->position >= EndNodes)) {
                    pSC = ACLTConditions::SC.erase(pSC);
                    continue;
                }
            } else {
                throw kaBaseException("Condition type at pos " PETSCINT_FORMAT " not implemented.", pSC->position);
            }
            ++pSC;
        }
    }
    
    if (materials) {
        // if InitMono is called with a material vector, ignore material vec in .aclt if given
        MaterialV = materials;
    } else if (Material.size()) {
        LoadVec(Material.c_str(), MaterialV, ot_bin);
    } else {
        throw kaBaseException("No Material vector was specified.");
    }
    
    ierr = VecDuplicate(MaterialV, &Force); CHKERRQ(ierr);
    
    PetscInt lenMaterialV;
    ierr = VecGetSize(MaterialV, &lenMaterialV); CHKERRQ(ierr);
    
    if (((lenMaterialV != Intra.size) && !GaussFwd)) // TODO check size against gauss matrix
        throw kaBaseException("Sizes of LSE and material vector do not fit");
    
    ierr = VecGetOwnershipRange(MaterialV, &StartCells, &EndCells); CHKERRQ(ierr);
    
    pElphyIntra = new vbElphyModel<double> *[EndCells-StartCells];
    assert(pElphyIntra);
    pForceIntra = new vbForceModel<double> *[EndCells-StartCells];
    assert(pForceIntra);
    
    EMIntra.New(EndCells-StartCells);
    
    PetscScalar *pm;
    ierr = VecGetArray(MaterialV, &pm); CHKERRQ(ierr);
    
    ierr = VecDuplicate(MaterialV, &VecVCell); CHKERRQ(ierr);
    ierr = VecSet(VecVCell, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(MaterialV, &VecICell); CHKERRQ(ierr);
    ierr = VecSet(VecICell, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(MaterialV, &VecVRes); CHKERRQ(ierr);
    ierr = VecSet(VecVRes, 0.0); CHKERRQ(ierr);
    
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        PetscInt lindex = Ii-StartCells;
        int modelIndex  = (int)pm[lindex]; // mwk: modelIndex == material class?
        pElphyIntra[lindex] = NULL;
        CellModelStruct::setElphyModel(&pElphyIntra[lindex], modelIndex);
        pForceIntra[lindex] = NULL;
        CellModelStruct::setForceModel(&pForceIntra[lindex], modelIndex);  // mwk: force model read-in
        PetscScalar Vm = pElphyIntra[lindex]->GetVm();
        ierr = VecSetValues(VecVCell, 1, &Ii, &Vm, INSERT_VALUES); CHKERRQ(ierr);
        EMIntra.Set((int)lindex, pElphyIntra[lindex], pForceIntra[lindex], forceset);
    }
    
    ierr = VecAssemblyBegin(VecVCell); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(VecVCell); CHKERRQ(ierr);
    if (GaussInt) {
        // ierr = VecSet(Intra.X, pElphyIntra[0]->GetVm()); CHKERRQ(ierr);
        
        // Solve M*V_node = P*V_gauss
        PETScLSE lse;
        lse.A = IntraRHS;  /* Still the mass matrix */
        
        ierr = VecDuplicate(Intra.X, &lse.B); CHKERRQ(ierr);
        ierr = MatMult(GaussInt, VecVCell, lse.B);
        
        lse.InitSolve(KSPTYPE, PCTYPE, RTOL, ATOL, DTOL, MAXITS);
        lse.Solve(lse.B, Intra.X);
        lse.A = NULL;  /* Otherwise matrix is destroyed */
    } else {
        ierr = VecCopy(VecVCell, Intra.X); CHKERRQ(ierr);
    }
    
    HeteroCMIntra.Set(pElphyIntra, EndCells-StartCells);
    
    if (backuploadname.size())
        ReadBackup(pElphyIntra, pForceIntra, &Intra);
    
    // Lattices von Ii, Ve... werden standardmÃ¤Ãig nicht gespeichert
    for (int i = 0; i < ResultVars::numResultVars; i++) {
        saveLattice[i] = false;
    }
    
    // Vm wird gespeichert
    int ri = 0;
    bool compareType;
    bool comparePara;
    bool nextPara;
    
    if (setResults == true) {
        int EndacltResult = 0;
        string acltResult;
        string acltResult_old;
        istringstream readresult(results);
        
        // Schleife geht die angegebenen Parameter im aclt file durch und sucht die Position ihres Auftretens im Zellmodell
        do {
            compareType = false;
            if (ri >= maxParameters)
                throw kaBaseException("Too much Resultparameters in acltfile. Exiting...");
            readresult >> acltResult;
            EndacltResult  = acltResult.compare(acltResult_old);
            acltResult_old = acltResult;
            if (EndacltResult == 0)
                continue;
            resultName[ri] = acltResult;
            nextPara       = false;
            
            // FÃ¼r jedes Voxel wird die Materialklasse/verwendetes Zellmodell festgestellt und die Position des Parameters im
            // Outputstring fÃ¼r dieses Modell gespeichert
            for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
                PetscInt lindex = Ii-StartCells;
                
                // rm entspricht dem Vektor der Zellmodelltypen
                int rm                         = (int)pm[lindex];
                double CallGetResultParameters =
                (EMIntra.GetResultParameters[lindex])(pElphyIntra[lindex], resultName[ri], resultPos[ri][rm]);
                
                // falls Parameter nicht im Zellmodell gefunden wird
                if (CallGetResultParameters == 0.) {
                    // check ob Ii, Ie ...
                    for (int i = 0; i < ResultVars::numResultVars; i++) {
                        if (!resultName[ri].compare(saveLattice_Name[i])) {
                            // Abfrage ob die Standardlattices im aclt File gespeichert werden sollen
                            saveLattice[i] = true;
                            compareType    = true;
                            break;
                        }
                    }
                    
                    // unbekannter Parameter, Resultwerte werden spÃ¤ter auf "0" gesetzt
                    if (compareType == false) {
                        resultPos[ri][rm] = -1;
                        nextPara          = true;
                    }
                }
                
                // Parameter im Zellmodell gefunden....
                else if (CallGetResultParameters == 1.) {
                    saveresultLattice[ri] = true;
                    nextPara              = true;
                }
            }
            if (nextPara == true)
                ri++;
        } while (EndacltResult != 0);
        numberParameters = ri;
    }
    
    if (numberParameters > 0) {
        // Vektor fÃ¼r Resultwerte erstellen
        ierr = VecDuplicate(Intra.X, &VecResultValues); CHKERRQ(ierr);
    } else {
        setResults = false;
    }
    
    if (saveLattice[ResultVars::AT]) {
        ierr = VecDuplicate(Intra.X, &ActivationTime); CHKERRQ(ierr);
        ierr = VecSet(ActivationTime, -1.0); CHKERRQ(ierr);
    }
    
    MPI_Barrier(PETSC_COMM_WORLD);
    numberParameters = ri;
    
    ierr = VecRestoreArray(MaterialV, &pm); CHKERRQ(ierr);
    
    /* Initialize Sensors */
    ACLTSensors::Parse(problemSize);
    if (domains == 1) {
        StartExtraNodes = StartNodes;
        EndExtraNodes   = EndNodes;
    }  /* Otherwise, these have already been set by InitBi(). */
    InitSensors();
    
    if (Intra.A) // At least a stiffness matrix has already been set from project file
        BuildIntraOperators();
    
    PrintDebug("InitMono finished");
}  // acCELLerate::InitMono

void acCELLerate::InitBi() {
    PrintDebug("InitBi()");
    
    PetscErrorCode ierr;
    
    if (VolumeExtra != 0.0) {
        ierr = MatScale(Extra.A, VolumeExtra); CHKERRQ(ierr);
    } else {
        Extra.MatVecScale(VecVolExtra);
    }
    
    Extra.CreateX();
    Extra.CreateB();
    
    ierr = VecGetOwnershipRange(Extra.X, &StartExtraNodes, &EndExtraNodes); CHKERRQ(ierr);
    
    ierr = VecDuplicate(Extra.B, &Stim); CHKERRQ(ierr);
    
    InitMono();
    
    /* Read IndexSet file (map of extracellular nodes to intracellular nodes in
     * case of unequal domain sizes). Unfortunately, a reader is not yet
     * implemented in PETSc, so we ahve to do the reading manually. */
    if (!IntraIndexSetFile.empty()) {
        if (saveLattice[ResultVars::Ie]) {
            /* If we want to save Ie and have an IS, we need a storage vecor */
            ierr = VecDuplicate(Extra.B, &VecSaveTmp); CHKERRQ(ierr);
            ierr = VecSet(VecSaveTmp, .0); CHKERRQ(ierr);
        }
        PetscViewer fd;
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, IntraIndexSetFile.c_str(),
                                     FILE_MODE_READ, &fd); CHKERRQ(ierr);
        PetscInt header[2];
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        ierr = PetscViewerBinaryRead(fd, header, 2, NULL, PETSC_INT); CHKERRQ(ierr);
#else  // if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        ierr = PetscViewerBinaryRead(fd, header, 2, PETSC_INT); CHKERRQ(ierr);
#endif  // if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        if (header[0] != IS_FILE_CLASSID) {
            throw kaBaseException("File %s is not an index set.", IntraIndexSetFile.c_str());
        } else if (header[1] > Extra.size) {
            throw kaBaseException("IndexSet %s is too large (%lld > %lld).",
                                  IntraIndexSetFile.c_str(), header[1], Extra.size);
        }
        PetscInt *idcs;
        ierr = PetscMalloc(header[1] * sizeof(PetscInt), &idcs); CHKERRQ(ierr);
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        ierr = PetscViewerBinaryRead(fd, idcs, header[1],  /*&n*/ NULL, PETSC_INT); CHKERRQ(ierr);
        
        // if (n < header[1]) {
        //  throw kaBaseException("Could not read index set %s: %lld of %lld entries read.",
        //      IntraIndexSetFile.c_str(), n, header[1]);
        // }
#else  // if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        ierr = PetscViewerBinaryRead(fd, idcs, header[1], PETSC_INT); CHKERRQ(ierr);
#endif  // if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR > 5
        ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
        
        GlobalToIntraMap.clear();
        for (PetscInt i = 0; i < header[1]; ++i) {
            GlobalToIntraMap[idcs[i]] = i;
        }
        
        /* Map conditions. We assume that all conditions use the node numbering of
         * the extracellular mesh. Erase all intra conditions that are not located
         * in the thread-local portion of the intra nodes. */
        ConditionContainerType::iterator pSC;
        for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); /* no increment because erase() */) {
            if ((pSC->cflag == CT_II) || (pSC->cflag == CT_UI)) {
                const auto &newpos = GlobalToIntraMap.find(pSC->position);
                if (newpos == GlobalToIntraMap.end()) {
                    throw kaBaseException("Condition position " PETSCINT_FORMAT " is not in intracellular domain.",
                                          pSC->position);
                } else if ((newpos->second < StartNodes) || (newpos->second >= EndNodes)) {
                    pSC = ACLTConditions::SC.erase(pSC);
                    continue;
                }
                pSC->position = newpos->second;
            } else if ((pSC->cflag == CT_IE) || (pSC->cflag == CT_UE)) {
                if ((pSC->position < StartExtraNodes) || (pSC->position >= EndExtraNodes)) {
                    pSC = ACLTConditions::SC.erase(pSC);
                    continue;
                }
            } else {
                throw kaBaseException("Condition type at pos " PETSCINT_FORMAT " not implemented.", pSC->position);
            }
            ++pSC;
        }
        
        if (!combinedExtra) {
            const PetscInt *rhs_cols;
            PetscInt rhs_ncols;
            const PetscScalar *rhs_vals;
            
            for (PetscInt i = StartNodes; i < EndNodes; ++i) {
                ierr = MatGetRow(ExtraRHS, i, &rhs_ncols, &rhs_cols, &rhs_vals); CHKERRQ(ierr);
                
                /* Add local RHS values to corresponding global entry */
                for (int j = 0; j < rhs_ncols; ++j) {
                    PetscScalar v;
                    ierr = MatSetValue(Extra.A, idcs[i], idcs[rhs_cols[j]], rhs_vals[j], ADD_VALUES); CHKERRQ(ierr);
                }
                
                ierr = MatRestoreRow(ExtraRHS, i, &rhs_ncols, &rhs_cols, &rhs_vals); CHKERRQ(ierr);
            }
            
            ierr          = MatAssemblyBegin(Extra.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr          = MatAssemblyEnd(Extra.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            combinedExtra = PETSC_TRUE;
        }
        
        IS is2;  /* Create local index set first */
        ierr = ISCreateGeneral(PETSC_COMM_SELF, EndCells-StartCells, idcs+StartCells,
                               PETSC_COPY_VALUES, &is2); CHKERRQ(ierr);
        ierr = ISOnComm(is2, PETSC_COMM_WORLD, PETSC_USE_POINTER, &IntraIndexSet); CHKERRQ(ierr);
    }
    
    if (!combinedExtra) {
        ierr = MatAXPY(Extra.A, 1.0, ExtraRHS, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    }
    
    /* TODO move somewhere else, delete ns when finished */
    MatNullSpace ns;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &ns); CHKERRQ(ierr);
    ierr = MatSetNullSpace(Extra.A, ns); CHKERRQ(ierr);
    if (CurrentDensity.A) {
        ierr = MatSetNullSpace(CurrentDensity.A, ns); CHKERRQ(ierr);
    }
    
    // ierr = MatSetTransposeNullSpace(Extra.A, ns); CHKERRQ(ierr);
    // ??? the documentation has this method but its not there...?
    
    /* For some reason, FECC matrices are not symmetric. BJACOBI works well for them. Eisenstat works much better for the
     * matrices of the UnstructuredGrid Bidomain Matrix Generator. */
    Extra.InitSolve(KSPTYPE, IntraRHS ? PCTYPE : PCBJACOBI, RTOL, ATOL, DTOL, MAXITS);
    
    PrintDebug("InitBi finished");
}  // acCELLerate::InitBi

void acCELLerate::InitTri() {
    PrintDebug("InitTri()");
    
    PetscErrorCode ierr;
    
    if (VolumeFibro != 0.0) {
        ierr = MatScale(Fibro.A, VolumeFibro); CHKERRQ(ierr);
    } else {
        Fibro.MatVecScale(VecVolExtra);
    }
    
    Fibro.CreateX();
    Fibro.CreateB();
    
    Vec MaterialFV;
    LoadVec(MaterialFibro.c_str(), MaterialFV, ot_bin);
    
    // TODO assert len(MaterialFV) == IntraSize
    ierr = VecGetOwnershipRange(MaterialFV, &StartCells, &EndCells); CHKERRQ(ierr);
    
    pElphyFibro = new vbElphyModel<double> *[EndCells-StartCells];
    assert(pElphyFibro);
    pForceFibro = new vbForceModel<double> *[EndCells-StartCells];
    assert(pForceFibro);
    
    EMFibro.New(EndCells-StartCells);
    
    PetscScalar *pm;
    ierr = VecGetArray(MaterialFV, &pm); CHKERRQ(ierr);
    
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        PetscInt lindex = Ii-StartCells;
        int modelIndex  = pm[lindex];
        pElphyFibro[lindex] = NULL;
        CellModelStruct::setElphyModel(&pElphyFibro[lindex], modelIndex);
        pForceFibro[lindex] = NULL;
        CellModelStruct::setForceModel(&pForceFibro[lindex], modelIndex);
        
        PetscScalar Vm = pElphyFibro[lindex]->GetVm();
        ierr = VecSetValues(Fibro.X, 1, &Ii, &Vm, INSERT_VALUES); CHKERRQ(ierr);
        
        EMFibro.Set((int)lindex, pElphyFibro[lindex], pForceFibro[lindex], forceset);
    }
    HeteroCMFibro.Set(pElphyFibro, EndCells-StartCells);
    
    ierr = VecRestoreArray(MaterialFV, &pm); CHKERRQ(ierr);
    ierr = VecDestroy(&MaterialFV); CHKERRQ(ierr);
    
    InitBi();
    if (!combinedExtra)
        MatAXPY(Extra.A, 1.0, Fibro.A, DIFFERENT_NONZERO_PATTERN);
    
    PrintDebug("InitTri finished");
}  // acCELLerate::InitTri

inline bool acCELLerate::IsFloat(const char *s) {
    char *endptr;
    double d = strtod(s, &endptr);
    
    if ((endptr == s) && (d == 0.0))
        return false;
    
    return true;
}

inline void acCELLerate::SetConditionStatus(ConditionType ct1, ConditionType ct2, bool redo) {
    PrintDebug2("SetConditionStatus()");
    
    ConditionContainerType::iterator pSC;
    
    for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); pSC++)
        if ((pSC->cflag == ct1) || (pSC->cflag == ct2)) {
            double cyclelength = (double)pSC->cl;
            
            // TODO this must be possible more clean!!!
            int cycles  = (int)(((double)currenttime-(double)pSC->toff)/cyclelength);
            double time = (double)currenttime-cycles*cyclelength-(double)pSC->toff;
            
            if ((0.0 <= time) && (time <= (double)pSC->sl) && (currenttime >= pSC->toff)) {
                switch (pSC->cs) {
                    case CS_applied:
                    case CS_apply:
                        pSC->cs = (redo ? CS_apply : CS_applied);
                        break;
                    default:
                        pSC->cs = CS_apply;
                        break;
                }
            } else {
                switch (pSC->cs) {
                    case CS_applied:
                    case CS_apply:
                        pSC->cs = CS_remove;
                        break;
                    default:
                        pSC->cs = CS_nothing_to_do;
                        break;
                }
            }
        }
}  // acCELLerate::SetConditionStatus

inline void acCELLerate::ApplyConditionsToVec(PetscScalar *px, PetscScalar *pb, ConditionType ct1, ConditionType ct2) {
    PrintDebug2("ApplyConditions()");
    
    ConditionContainerType::iterator pSC;
    
    for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); pSC++)
        if ((pSC->cs == CS_apply) && ((pSC->cflag == ct1) || (pSC->cflag == ct2) )) {
            if (pSC->cflag == ct1)
                pb[pSC->position-StartNodes] += pSC->val;
            if (pSC->cflag == ct2)
                px[pSC->position-StartNodes] = pSC->val;
        }
}

inline void acCELLerate::ApplyConditions(PetscScalar *px, PetscScalar *pb, ConditionType ct1, ConditionType ct2,
                                         bool redo) {
    PrintDebug2("ApplyConditions()");
    
    SetConditionStatus(ct1, ct2, redo);
    ApplyConditionsToVec(px, pb, ct1, ct2);
    
    PrintDebug2("ApplyConditions finished");
}

inline void acCELLerate::ApplyConditionsExtra() {
    PrintDebug2("ApplyConditionsExtra()");
    
    SetConditionStatus(CT_UE, CT_IE, false);
    
    PetscErrorCode ierr;
    bool potentialBoundaryChanged = false;
    
    vector<BoundaryStorage> BS;
    
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    
    ConditionContainerType::iterator pSC;
    for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); pSC++) {
        if (pSC->cflag == CT_IE) {
            if ((pSC->cs == CS_apply) || (pSC->cs == CS_applied))
                Extra.SetCurrent(pSC->position, pSC->val);
        }
        if (pSC->cflag == CT_UE) {
            if ((pSC->cs == CS_apply) || (pSC->cs == CS_applied)) {
                BoundaryStorage lBS;
                lBS.index = pSC->position;
                lBS.value = pSC->val;
                BS.push_back(lBS);
            }
            if ((pSC->cs == CS_apply) || (pSC->cs == CS_remove))
                potentialBoundaryChanged = true;
        }
    }
    
    if (potentialBoundaryChanged) {
        for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); pSC++) {
            if (pSC->cflag == CT_UE) {
                if (pSC->GetColumnEntries(&ncols, &cols, &vals)) {
                    ierr = MatSetValues(Extra.A, ncols, cols, 1, &pSC->position, vals, INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValues(Extra.A, 1, &pSC->position, ncols, cols, vals, INSERT_VALUES); CHKERRQ(ierr);
                }
            }
        }
        ierr = MatAssemblyBegin(Extra.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Extra.A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        
        for (pSC = ACLTConditions::SC.begin(); pSC != ACLTConditions::SC.end(); pSC++) {
            if (pSC->cflag == CT_UE) {
                if (pSC->cs == CS_apply) {
                    Extra.AddColumnVec(&Stim, pSC->position, pSC->val);
                    ierr = MatGetRow(Extra.A, pSC->position, &ncols, &cols, &vals); CHKERRQ(ierr);
                    pSC->SetColumnEntries(ncols, const_cast<PetscInt *>(cols), (double *)vals);
                    ierr    = MatRestoreRow(Extra.A, pSC->position, &ncols, &cols, &vals); CHKERRQ(ierr);
                    pSC->cs = CS_applied;
                }
                if (pSC->cs == CS_remove) {
                    pSC->GetColumnEntries(&ncols, &cols, &vals);
                    vector<PetscScalar> values(ncols);
                    for (int j = 0; j < ncols; j++)
                        values[j] = -1.0*vals[j]*pSC->val;
                    ierr = VecSetValues(Stim, ncols, cols, values.data(), ADD_VALUES); CHKERRQ(ierr);
                    pSC->FreeColumnEntries();
                    pSC->cs = CS_nothing_to_do;
                }
            }
        }
        ierr = VecAssemblyBegin(Stim); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(Stim); CHKERRQ(ierr);
    }
    
    PetscInt *indizes;
    PetscScalar *values;
    
    if (BS.size()) {
        ierr = VecAXPY(Extra.B, -1.0, Stim); CHKERRQ(ierr);
        
        MPI_Barrier(PETSC_COMM_WORLD);
        indizes = new PetscInt[BS.size()];
        values  = new PetscScalar[BS.size()];
        
        vector<BoundaryStorage>::iterator iBS, bBS = BS.begin(), eBS = BS.end();
        int i = 0;
        for (iBS = bBS; iBS != eBS; iBS++, i++) {
            indizes[i] = (*iBS).index;
            values[i]  = (*iBS).value;
        }
        if (mpirank == 0) {
            ierr = VecSetValues(Extra.B, BS.size(), indizes, values, INSERT_VALUES); CHKERRQ(ierr);
        }
        
        ierr = VecAssemblyBegin(Extra.B); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(Extra.B); CHKERRQ(ierr);
        
        if (potentialBoundaryChanged) {
            Extra.MatZeroColumns(BS.size(), indizes);
            ierr = MatZeroRows(Extra.A, BS.size(), indizes, 1.0, 0, 0); CHKERRQ(ierr);
        }
        delete indizes;
        delete values;
    }
    PrintDebug2("ApplyConditionsExtra finished");
}  // acCELLerate::ApplyConditionsExtra

inline void acCELLerate::InitSensors() {
    PrintDebug2("InitSensors()");
    
    PetscInt ri = numberParameters;
    bool comparePara, compareType;
    
    PetscScalar *pm;
    PetscErrorCode ierr = VecGetArray(MaterialV, &pm); CHKERRQ(ierr);
    
    for (auto pSC = ACLTSensors::SC.begin(); pSC != ACLTSensors::SC.end(); /*no increment because of erase*/) {
        PetscInt row = pSC->SensorIndex;
        if ((pSC->cflag == ST_II) || (pSC->cflag == ST_UI)) {
            if (IntraIndexSet) {
                try {
                    pSC->SensorIndex = GlobalToIntraMap.at(pSC->SensorIndex);
                } catch (...) {
                    throw kaBaseException("Node %d for intracellular sensor is not in intracellular domain.", row);
                }
            }
            
            if ((pSC->SensorIndex >= StartNodes) && (pSC->SensorIndex < EndNodes)) {
                pSC->Init();
            } else {
                pSC = ACLTSensors::SC.erase(pSC);
                continue;
            }
        } else if ((pSC->cflag == ST_IE) || (pSC->cflag == ST_UE)) {
            if ((row >= StartExtraNodes) && (row < EndExtraNodes)) {
                pSC->Init();
            } else {
                pSC = ACLTSensors::SC.erase(pSC);
                continue;
            }
        }
        
        if (pSC->cflag == ST_PARA) {
            comparePara    = false;
            compareType    = false;
            resultName[ri] = pSC->SensorParameterType;
            
            if ((row >= StartCells) && (row < EndCells)) {
                pSC->Init();
                int rm                         = (int)pm[row-StartCells];
                double CallGetResultParameters =
                (EMIntra.GetResultParameters[row-StartCells])(pElphyIntra[row-StartCells], resultName[ri], resultPos[ri][rm]);
                if (CallGetResultParameters == 1.) {
                    parametersinSensorfile = true;
                    if (ri > 0) {
                        for (int i = 0; i < ri; i++) {
                            if (!resultName[ri].compare(resultName[i])) {
                                comparePara = true;
                            }
                        }
                    }
                    
                    // Wenn Parameter aus Zellmodell im Sensorfile vorhanden wird er an das Parameterarray angehÃ¤ngt
                    if (comparePara == false) {
                        saveresultLattice[ri] = false;
                        ri++;
                    }
                } else {
                    // Standardparameter?
                    for (int i = 0; i < ResultVars::numResultVars; i++) {
                        if (!resultName[ri].compare(saveLattice_Name[i]))
                            compareType = true;
                    }
                    if (compareType == false) {
                        throw kaBaseException("Parameter %s not found in the Cellmodel (%i) at position %i. Exiting...",
                                              resultName[ri].c_str(), rm, row-StartCells);
                    }
                }
            } else {  /* Sensor not in local cells */
                pSC = ACLTSensors::SC.erase(pSC);
                continue;
            }
        }
        
        ++pSC;
    }
    
    MPI_Barrier(PETSC_COMM_WORLD);
    numberParameters = ri;
    
    ierr = VecRestoreArray(MaterialV, &pm); CHKERRQ(ierr);
    
    PrintDebug2("InitSensors finished");
}  // acCELLerate::InitSensors

inline void acCELLerate::SaveSensors(PetscScalar *pix, PetscScalar *pib, PetscScalar *pex, PetscScalar *peb,
                                     PetscScalar *pfx, PetscScalar *pfb) {
    PrintDebug2("SaveSensors(%p %p %p %p %p %p)", pix, pib, pex, peb, pfx, pfb);
    
    for (auto pSC = ACLTSensors::SC.begin(); pSC != ACLTSensors::SC.end(); pSC++) {
        if (currenttime < pSC->beginsave)
            continue;
        pSC->beginsave += pSC->dtsave;
        PetscInt row = pSC->SensorIndex;
        
        switch (pSC->cflag) {
            case ST_II:
                if (pib)
                    pSC->Save(currenttime, pib[row-StartNodes]); break;
            case ST_UI:
                if (pix)
                    pSC->Save(currenttime, pix[row-StartNodes]); break;
            case ST_IF:
                if (pfb)
                    pSC->Save(currenttime, pfb[row-StartNodes]); break;
            case ST_UF:
                if (pfx)
                    pSC->Save(currenttime, pfx[row-StartNodes]); break;
            case ST_IE:
                if (peb)
                    pSC->Save(currenttime, peb[row-StartExtraNodes]); break;
            case ST_UE:
                if (pex)
                    pSC->Save(currenttime, pex[row-StartExtraNodes]); break;
            case ST_PARA:
                break;
            case ST_UEM:
                if (pex) {
                    PetscScalar r;
                    PetscErrorCode ierr = VecSum(Extra.X, &r); CHKERRQ(ierr);
                    r /= Extra.size;
                    pSC->Save(currenttime, pex[pSC->SensorIndex-StartNodes]-r);
                }
                break;
            default:
                throw kaBaseException("Sensor type invalid");
        }  // switch
    }
    
    PrintDebug2("SaveSensors finished");
}  // acCELLerate::SaveSensors

void acCELLerate::SaveResultIntra() {
    PrintDebug2("SaveResultIntra()");
    
    char tmp[256];
    if (saveLattice[ResultVars::Vm] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::Vm]);
        Intra.SaveX(tmp);
    }
    if (saveLattice[ResultVars::Ii] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::Ii]);
        if (IntraRHS) {
            /* we have to solve [K Phi_e +] K V_m = M I_i for I_i*/
            PetscErrorCode ierr = MatMult(ExtraRHS, Intra.X, CurrentDensity.B); CHKERRQ(ierr);
            
            if (domains > 1) {
                Vec tmpvec = Extra.X;
                if (IntraIndexSet) {
                    ierr = VecGetSubVector(Extra.X, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                }
                ierr = MatMultAdd(ExtraRHS, tmpvec, CurrentDensity.B, CurrentDensity.B); CHKERRQ(ierr);
                if (IntraIndexSet) {
                    ierr = VecRestoreSubVector(Extra.X, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                }
            }
            
            CurrentDensity.Solve(CurrentDensity.B, CurrentDensity.X);
            
            KSPConvergedReason cr;
            ierr = KSPGetConvergedReason(CurrentDensity.ksp, &cr); CHKERRQ(ierr);
            if (cr < 0) {
                ierr = PetscFPrintf(PETSC_COMM_WORLD, stderr, "%s: %s LSE did not converge: %s\n",
                                    currenttime.toStr(4, 6).c_str(), saveLattice_Name[ResultVars::Ii], KSPConvergedReasons[cr]);
                CHKERRQ(ierr);
            }
            
            SaveVecMatInfo(tmp, CurrentDensity.X, Intra.xot, Intra.m, Intra.xLattice, Intra.yLattice, Intra.zLattice);
        } else {
            Intra.SaveB(tmp);
        }
    }
    if (forceset) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(), "Force");
        SaveVecMatInfo(tmp, Force, Intra.xot, Intra.m, Intra.xLattice, Intra.yLattice, Intra.zLattice);
    }
}  // acCELLerate::SaveResultIntra

void acCELLerate::SaveResultParameters(string currentParameter) {
    PrintDebug2("SaveResultParameters()");
    char tmp[256];
    sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(), currentParameter.c_str());
    SaveVecMatInfo(tmp, VecResultValues, ot_bin, Intra.m, Intra.xLattice, Intra.yLattice, Intra.zLattice);
}

void acCELLerate::SaveResultExtra() {
    PrintDebug2("SaveResultExtra()");
    
    char tmp[256];
    if (saveLattice[ResultVars::Ve] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::Ve]);
        Extra.SaveX(tmp);
    }
    if (saveLattice[ResultVars::Ie] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::Ie]);
        if (IntraRHS) {
            /* We have to solve K V_m = M I_e */
            KSPConvergedReason cr;
            PetscErrorCode ierr;
            if (IntraIndexSet) {
                Vec tmpvec;
                ierr = VecGetSubVector(Extra.B, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                CurrentDensity.Solve(tmpvec, CurrentDensity.X);
                ierr = VecRestoreSubVector(Extra.B, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                ierr = VecGetSubVector(VecSaveTmp, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                ierr = VecCopy(CurrentDensity.X, tmpvec); CHKERRQ(ierr);
                ierr = VecRestoreSubVector(VecSaveTmp, IntraIndexSet, &tmpvec); CHKERRQ(ierr);
                SaveVecMatInfo(tmp, VecSaveTmp, Extra.xot, Extra.m, Extra.xLattice, Extra.yLattice, Extra.zLattice);
            } else {
                CurrentDensity.Solve(Extra.B, CurrentDensity.X);
                SaveVecMatInfo(tmp, CurrentDensity.X, Extra.xot, Extra.m, Extra.xLattice, Extra.yLattice, Extra.zLattice);
            }
            ierr = KSPGetConvergedReason(CurrentDensity.ksp, &cr); CHKERRQ(ierr);
            if (cr < 0) {
                ierr = PetscFPrintf(PETSC_COMM_WORLD, stderr, "%s: %s LSE did not converge: %s\n",
                                    currenttime.toStr(4, 6).c_str(), saveLattice_Name[ResultVars::Ie], KSPConvergedReasons[cr]);
                CHKERRQ(ierr);
            }
        } else {
            Extra.SaveB(tmp);
        }
    }
}  // acCELLerate::SaveResultExtra

void acCELLerate::SaveResultFibro() {
    PrintDebug2("SaveResultFibro()");
    
    char tmp[256];
    if (saveLattice[ResultVars::Vf] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::Vf]);
        Fibro.SaveX(tmp);
    }
    if (saveLattice[ResultVars::If] == true) {
        sprintf(tmp, VecFnTemplate, resultprefix.c_str(), currenttime.toStr(4, 6).c_str(),
                saveLattice_Name[ResultVars::If]);
        Fibro.SaveB(tmp);
    }
}

void acCELLerate::SetLoadBackup(string bckname) {
    if (bckname.size()) {
        backuploadname = bckname;
    } else {
        backuploadname = resultprefix;
        backuploadname.append(".abp");
    }
}

void acCELLerate::SetSaveBackup(string bckname) {
    if (bckname.size()) {
        backupsavename = bckname;
    } else {
        backupsavename = resultprefix;
        backupsavename.append(".abp");
    }
}

bool acCELLerate::ReadBackup(vbElphyModel<double> **ElphyModel, vbForceModel<double> **ForceModel, PETScLSE *Vvector,
                             bool firstdomain) {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] ReadBackup(...)\n", mpirank);  // PrintDebug("ReadBackup(...)");
#endif  // if PSDEBUG > 0
    
    // Calculate theoretical size of backup file without header
    long int objectSizeTemp = 0;
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        int my_rank;
        PetscInt lindex = Ii-StartCells;
        objectSizeTemp += (ElphyModel[lindex]->GetSize());
        objectSizeTemp += sizeof(ML_CalcType);
    }
    long int sumOfObjectSizes = 0;
    MPI_Reduce(&objectSizeTemp, &sumOfObjectSizes, 1, MPI_LONG, MPI_SUM, 0, PETSC_COMM_WORLD);  // Collect and sum up
                                                                                                // objectSizeTemp, store in sumOfObjectSizes
    
    int buSize, buSecondLineLength;
    int seekpos = ReceiveInteger();  // between seekpos=ReceiveInteger() here and SendInteger(seekpos) beneath, all
                                     // functions using MPI_Recv are
                                     // leading to a deadlock with more than one CPU! Avoid coding these functions (e.g. VecCreate()...) in this section!
    acltTime timetemp = ReceiveTime();
    fpos_t   buPosition1, buPosition2;
    
    FILE *bu = fopen(backuploadname.c_str(), "r");
    if (!bu)
        throw kaBaseException("Opening backupfile %s failed. Exiting...", backuploadname.c_str());
    
    fseek(bu, seekpos, SEEK_SET);
    
    if ((mpirank == 0) && firstdomain) {
        char tempstr[256];
        fscanf(bu, "%s\n", tempstr);
        const char *aBUv = acCELLerateBUVersion.c_str();
        if (strcasecmp(aBUv, tempstr) != 0) {
            // This exception leads to a deadlock when using MORE THAN ONE PROCESS caused by missing SendInteger()
            // and SendDouble() needed for Receive...() from above! No Error is displayed!
            throw kaBaseException("Backup file %s is incompatible to acCELLerate backup version 1.0. Exiting...",
                                  backuploadname.c_str());
        }
        
        int indextemp;
        fgetpos(bu, &buPosition1);
        fscanf(bu, "%[^\n]", tempstr);
        buSecondLineLength = strlen(tempstr);  // Stringlength needed for calculation of theoretical size of backup file
                                               // header
        fsetpos(bu, &buPosition1);
        char ts[64];
        fscanf(bu, "%63s %d\n", ts, &indextemp);
        timetemp = ts;
        if (indextemp != Intra.size) {
            throw kaBaseException("Dimension of Intra vector is not compatible with Backup file %s. Exiting...",
                                  backuploadname.c_str());
        }
        
        // get size of backupfile
        fgetpos(bu, &buPosition2);
        fseek(bu, 0, SEEK_END);
        buSize = ftell(bu);
    }
    
    if ((mpirank == 0) && firstdomain) {
        // Calculate theoretical size of backup file: sum of object sizes + header + endl sizes (1 respectively)
        long int calcSize = sumOfObjectSizes + acCELLerateBUVersion.size()+1 + buSecondLineLength+1;
        
        if (buSize != calcSize)
            throw kaBaseException(" Backup file %s does not match to project file. Exiting...", backuploadname.c_str());
        
        fsetpos(bu, &buPosition2);
    }
    
    struct stat file;
    if (stat(backuploadname.c_str(), &file))
        printf("\n bu file size: %i \n", (int)file.st_size);
    
    starttime   = timetemp;
    currenttime = timetemp;
    
    PetscErrorCode ierr;
    PetscScalar   *pV, *pVe;
    ierr = VecGetArray(Vvector->X, &pV); CHKERRQ(ierr);
    
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        PetscInt lindex = Ii-StartCells;
        if (firstdomain)
            kaRead(&pV[lindex], sizeof(PetscScalar), 1, bu);
        
        ElphyModel[lindex]->ReadStatus(bu);
        
        // fread(NULL, 0, 1, bu);
        // ForceModel[lindex]->ReadStatus(bu); // mwk: force model cannot be read from initializeXCLT backup yet in yet.
    }
    ierr = VecRestoreArray(Vvector->X, &pV); CHKERRQ(ierr);
    
    seekpos = ftell(bu);
    fclose(bu);
    
    SendInteger(seekpos);  //! !!!!!!!! seekpos muss wieder zurueckgegeben werden bei mehreren Domaenen!!!!!
    SendTime(timetemp);
    MPI_Barrier(PETSC_COMM_WORLD);
    
    return true;
}  // acCELLerate::ReadBackup

void acCELLerate::WriteBackup(vbElphyModel<double> **ElphyModel, vbForceModel<double> **ForceModel, PETScLSE *Vvector,
                              bool firstdomain) {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] WriteBackup(...)\n", mpirank);  // PrintDebug("WriteBackup(...)");
#endif  // if PSDEBUG > 0
    
    int seekpos = ReceiveInteger();
    FILE *bu    = fopen(backupsavename.c_str(), (mpirank == 0 && firstdomain) ? "w" : "r+");
    if (!bu)
        throw kaBaseException("Opening backupfile %s failed. Exiting...", backupsavename.c_str());
    
    fseek(bu, seekpos, SEEK_SET);
    
    if ((mpirank == 0) && firstdomain) {
        fprintf(bu, "%s\n", acCELLerateBUVersion.c_str());
        fprintf(bu, "%s %d\n", currenttime.toStr().c_str(), Intra.size);
    }
    
    PetscErrorCode ierr;
    PetscScalar   *pV, *pVe;
    ierr = VecGetArray(Vvector->X, &pV); CHKERRQ(ierr);
    
    // ierr=VecGetArray(Vevector->X, &pVe);CHKERRQ(ierr); // mwk 20.07.2012: to write Ve in Backup as well. If desired,
    // also change calcSize above.
    
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        PetscInt lindex = Ii-StartCells;
        if (firstdomain) {
            kaWrite(&pV[lindex], sizeof(PetscScalar), 1, bu);
            
            // kaWrite(&pVe[lindex], sizeof(PetscScalar), 1, bu); // mwk 20.07.2012: to write Ve in Backup as well If desired,
            // also change calcSize above.
        }
        
        ElphyModel[lindex]->WriteStatus(bu);
        
        // ForceModel[lindex]->WriteStatus(bu); // mwk 20.07.2012: ForceModel should not be written, otherwise acCELLerate
        // cannot read in Matrix.
    }
    ierr = VecRestoreArray(Vvector->X, &pV); CHKERRQ(ierr);
    
    seekpos = ftell(bu);
    fclose(bu);
    SendInteger(seekpos);
    MPI_Barrier(PETSC_COMM_WORLD);
}  // acCELLerate::WriteBackup

void acCELLerate::SendInteger(int valuetosend) {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] SendInteger(%d)\n", mpirank, valuetosend);
#endif  // if PSDEBUG > 0
    if (mpirank != mpisize-1)
        MPI_Send(&valuetosend, 1, MPI_INT, mpirank+1, 0, PETSC_COMM_WORLD);
}

int acCELLerate::ReceiveInteger() {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] ReceiveInteger()\n", mpirank);  // PrintDebug("ReceiveInteger()");
#endif  // if PSDEBUG > 0
    int valuetoreceive = 0;
    if (mpirank != 0) {
        MPI_Status stat;
        MPI_Recv(&valuetoreceive, 1, MPI_INT, mpirank-1, 0, PETSC_COMM_WORLD, &stat);
    }
    return valuetoreceive;
}

void acCELLerate::SendDouble(double valuetosend) {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] SendDouble(%lf)\n", mpirank, valuetosend);
#endif  // if PSDEBUG > 0
    
    if (mpirank != mpisize-1)
        MPI_Send(&valuetosend, 1, MPI_DOUBLE, mpirank+1, 0, PETSC_COMM_WORLD);
}

double acCELLerate::ReceiveDouble() {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] ReceiveDouble()\n", mpirank);  // PrintDebug("ReceiveDouble()");
#endif  // if PSDEBUG > 0
    double valuetoreceive = 0.0;
    if (mpirank != 0) {
        MPI_Status stat;
        MPI_Recv(&valuetoreceive, 1, MPI_DOUBLE, mpirank-1, 0, PETSC_COMM_WORLD, &stat);
    }
    return valuetoreceive;
}

void acCELLerate::SendTime(acltTime valuetosend) {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] SendTime(%s)\n", mpirank, valuetosend.toStr().c_str());  //  PrintDebug("SendDouble(%lf)",
                                                                                                // valuetosend);
#endif  // if PSDEBUG > 0
    
    if (mpirank != mpisize-1)
        MPI_Send(&valuetosend, sizeof valuetosend, MPI_BYTE, mpirank+1, 0, PETSC_COMM_WORLD);
}

acltTime acCELLerate::ReceiveTime() {
#if PSDEBUG > 0
    PetscPrintf(PETSC_COMM_SELF, "[%d] ReceiveTime()\n", mpirank);  //  PrintDebug("ReceiveDouble()");
#endif  // if PSDEBUG > 0
    acltTime valuetoreceive;
    if (mpirank != 0) {
        MPI_Status stat;
        MPI_Recv(&valuetoreceive, sizeof valuetoreceive, MPI_BYTE, mpirank-1, 0, PETSC_COMM_WORLD, &stat);
    }
    return valuetoreceive;
}

void acCELLerate::Run() {
    PrintDebug2("Run()");
    
    if (protocolfile.size() && (mpirank == 0)) {
        fprot = fopen(protocolfile.c_str(), "w");
        if (!fprot)
            throw kaBaseException("Opening file %s for output", protocolfile.c_str());
    }
    
    switch (domains) {
        case 1:
            InitMono();
            MonoDomain(calclen);
            break;
        case 2:
            InitBi();
            BiDomain();
            break;
        case 3:
            InitTri();
            TriDomain();
            break;
    }
    
    if (ActivationTime) {
        std::string tmp(resultprefix + "_" + saveLattice_Name[ResultVars::AT] + ".vec");
        const char *fn = tmp.c_str();
        SaveVec(fn, ActivationTime, ot_bin);
    }
}  // acCELLerate::Run

void acCELLerate::CalcCellModels(PetscScalar const *pix, PetscScalar const *pib, PetscScalar *pidv,
                                 PetscScalar const *ps, PetscScalar const *pv, PetscScalar *pif) {
    // falls weder im acltfile noch im sensorfile Zellmodellparameter auftauchen
    for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
        PetscInt lindex = Ii-StartCells;
        
        double force =
        (EMIntra.CouplingMethod[lindex])(pElphyIntra[lindex], pForceIntra[lindex], (double)dtcell,
                                         (ps ? ps[lindex] : 1.0), (pv ? pv[lindex] : 0.0),
                                         pidv[lindex], pix[lindex], pib[lindex]);
        if (pif)
            pif[lindex] = force;
    }
}

void acCELLerate::SaveCellModelVariables(PetscScalar const *pix, acltTime saveat) {
    PetscErrorCode ierr;
    
    // pres -> Pointer auf Resultvektor
    PetscScalar *pres;
    
    if (setResults || parametersinSensorfile) {
        PetscScalar *pm;
        
        // pm -> pointer auf Zellmodelltypvektor
        ierr = VecGetArray(MaterialV, &pm); CHKERRQ(ierr);
        
        // fÃ¼r jeden Parameter...
        for (int j = 0; j < numberParameters; j++) {
            if (setResults == true) {
                ierr = VecGetArray(VecResultValues, &pres);
                CHKERRQ(ierr);
            }
            
            for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
                PetscInt lindex = Ii-StartCells;
                int m           = (int)pm[lindex];
                
                // Suche Parameterwert und speichere das Lattice falls entsprechendes Zeitintervall erreicht und der Parameter
                // im acltfile vor kam (Lattices der Parameter im
                // Sensorfile werden nciht gespeichert)
                if ((currenttime >= saveat) && (saveresultLattice[j] == true) ) {
                    if (resultPos[j][m] == -1) {
                        resultValue[j] = 0;
                    } else {
                        double CallGetResultValues =
                        (EMIntra.GetResultValues[lindex])(pElphyIntra[lindex], resultPos[j][m], resultValue[j],
                                                          (double)currenttime,
                                                          pix[lindex]);
                    }
                    pres[lindex] = resultValue[j];
                }
                
                // Suche Parameter der Sensorfiles und speichere die Werte in den Sensorfiles
                if (parametersinSensorfile) {
                    for (auto pSCPS = ACLTSensors::SC.begin(); pSCPS != ACLTSensors::SC.end(); ++pSCPS) {
                        if ((Ii == pSCPS->SensorIndex)
                            && (pSCPS->cflag == ST_PARA)
                            && !pSCPS->SensorParameterType.compare(resultName[j])
                            && (currenttime >=  pSCPS->beginsave)) {
                            double CallGetResultValues_sensor =
                            (EMIntra.GetResultValues[lindex])(pElphyIntra[lindex], resultPos[j][m], resultValue[j],
                                                              (double)currenttime, pix[lindex]);
                            pSCPS->beginsave += pSCPS->dtsave;
                            pSCPS->Save(currenttime, resultValue[j]);
                        }
                    }
                }
            }
            
            if (setResults) {
                ierr = VecRestoreArray(VecResultValues, &pres);
                CHKERRQ(ierr);
            }
            
            if ((currenttime >= saveat) && (saveresultLattice[j] == true) ) {
                SaveResultParameters(resultName[j]);
            }
        }
    }
}  // acCELLerate::SaveCellModelVariables

void acCELLerate::IntraStep(acltTime &intraat, PetscScalar *stretch, PetscScalar *velocity) {
    PetscErrorCode ierr;
    PetscScalar   *pix, *pib, *pvvi, *piv, *pidv, *pif;
    
    t_mono_matrix.start();
    bool createIntraBVector = (currenttime >= intraat);
    if (createIntraBVector) {
        intraat += dtintra;
        
        if (!this->IntraRHS) {  /* e.g. FD or FECC or FE w/ lumped M, i.e., where A = Inv(M)*K */
            ierr = MatMult(Intra.A, Intra.X, Intra.B); CHKERRQ(ierr);
            if (Extra.X) {       // bidomain
                if (IntraIndexSet) {
                    Vec PhiEIntra;
                    ierr = VecGetSubVector(Extra.X, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
                    ierr = MatMultAdd(Intra.A, PhiEIntra, Intra.B, Intra.B); CHKERRQ(ierr);
                    ierr = VecRestoreSubVector(Extra.X, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
                } else {
                    ierr = MatMultAdd(Intra.A, Extra.X, Intra.B, Intra.B); CHKERRQ(ierr);
                }
            }
        } else {
            ierr = VecSet(Intra.B, .0); CHKERRQ(ierr);
        }
    }
    t_mono_matrix.stop();
    
    if (GaussFwd) {
        ierr = MatMult(GaussFwd, Intra.X, VecVRes); CHKERRQ(ierr);
        ierr = VecGetArray(VecVRes, &pix); CHKERRQ(ierr);
    } else {
        ierr = VecGetArray(Intra.X, &pix); CHKERRQ(ierr);
    }
    
    ierr = VecGetArray(Intra.B, &pib); CHKERRQ(ierr);
    ierr = VecGetArray(VecVCell, &piv); CHKERRQ(ierr);
    
    if (Force) {
        ierr = VecGetArray(Force, &pif); CHKERRQ(ierr);
    }
    
    if (createIntraBVector && !this->IntraRHS) {
        if (VolumeMyo == 0.0) {ierr = VecGetArray(VecVolMyo, &pvvi); CHKERRQ(ierr);}
        for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
            PetscScalar volume_local = 0.0;
            if (VolumeMyo == 0.0)
                volume_local = (pvvi[Ii-StartCells] < 1e-9 ? 1e-9 : pvvi[Ii-StartCells]);
            else
                volume_local = VolumeMyo;
            
            pib[Ii-StartCells] *= -pElphyIntra[Ii-StartCells]->Volume()/volume_local*1e9;  // convert A/m^3 to nA/cell
        }
        if (VolumeMyo == 0.0) {ierr = VecRestoreArray(VecVolMyo, &pvvi); CHKERRQ(ierr);}
    }
    
    t_mono_cond.start();
    ApplyConditions(pix, pib, CT_II, CT_UI, createIntraBVector);
    if (GaussFwd) {
        ierr = VecRestoreArray(Intra.B, &pib); CHKERRQ(ierr);
        ierr = MatMult(GaussFwd, Intra.B, VecICell); CHKERRQ(ierr);
        ierr = VecGetArray(VecICell, &pib); CHKERRQ(ierr);
    }
    t_mono_cond.stop();
    
    t_mono_cell.start();
    CalcCellModels(pix, pib, piv, stretch, velocity, pif);
    t_mono_cell.stop();
    
    ierr = VecRestoreArray(VecVCell, &piv); CHKERRQ(ierr);
    
    if (GaussFwd) {
        ierr = VecRestoreArray(VecICell, &pib); CHKERRQ(ierr);
        ierr = VecRestoreArray(VecVRes, &pix); CHKERRQ(ierr);
    } else {
        ierr = VecRestoreArray(Intra.B, &pib); CHKERRQ(ierr);
    }
    
    if (pif) {
        ierr = VecRestoreArray(Force, &pif); CHKERRQ(ierr);
    }
    
    // TODO: only if not gauss?
    ierr = VecCopy(VecVCell, VecVRes); CHKERRQ(ierr);
    if (implicit) {
        for (int jacit = 1; jacit < jacobi; jacit++) {  // damped Jacobi method
            ierr = MatMultAdd(MatImp, VecVRes, VecVCell, VecRes); CHKERRQ(ierr);
            ierr = VecPointwiseMult(VecRes, VecRes, VecDiag); CHKERRQ(ierr);
            ierr = VecAXPY(VecVRes, 0.5, VecRes); CHKERRQ(ierr);
        }
    }
    
    if (currentScheme == CurrentScheme::Godunow) {
        ierr = VecGetArray(VecVRes, &pidv); CHKERRQ(ierr);
        for (PetscInt Ii = StartNodes; Ii < EndNodes; Ii++) {
            PetscInt lindex = Ii-StartNodes;
            pix[lindex] += pidv[lindex];
            
            if (::isnan(pix[lindex])) {
                if (mpisize != 1)
                    SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)!", Ii);
                else
                    throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)", Ii);
            } else if ((pix[lindex] > vm_def_range) || (pix[lindex] < -vm_def_range)) {
                if (mpisize != 1)
                    SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is out of range!", Ii);
                else
                    throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is out of range", Ii);
            }
        }
        ierr = VecRestoreArray(VecVRes, &pidv); CHKERRQ(ierr);
    }
    
    t_mono_cell.stop();
    
    if (this->IntraRHS) {
        t_mono_matrix.start();
        
        // b = (M - (1-theta)*dt/(C_m*beta)*K) * x
        ierr = MatMult(this->IntraRHS, Intra.X, Intra.B); CHKERRQ(ierr);  // b = Mx
        
        if (GaussInt) {
            /* b = b + P*dv */
            ierr = MatMultAdd(GaussInt, VecVRes, Intra.B, Intra.B); CHKERRQ(ierr);
        } else {
            ierr = VecRestoreArray(Intra.X, &pix); CHKERRQ(ierr);
        }
        
        if (domains > 1) {  // bidomain
            /* b = b - (dt/C_m*beta)*K_i * Phi_e */
            if (IntraIndexSet) {
                Vec PhiEIntra;
                ierr = VecGetSubVector(Extra.X, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
                ierr = MatMultAdd(ExtraRHS, PhiEIntra, Intra.B, Intra.B); CHKERRQ(ierr);
                ierr = VecRestoreSubVector(Extra.X, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
            } else {
                ierr = MatMultAdd(ExtraRHS, Extra.X, Intra.B, Intra.B); CHKERRQ(ierr);
            }
        }
        
        if (currentScheme == CurrentScheme::ICI) {
            // b = b + beta M * C_m * dV/dt = b - beta M * I_ion
            ierr = MatMultAdd(IntraMass, VecVRes, Intra.B, Intra.B); CHKERRQ(ierr);
        }
        
        // Solve (M + theta*dt/(C_m*beta)*K) * x = b
        Intra.Solve(Intra.B, Intra.X);
        
        KSPConvergedReason cr;
        ierr = KSPGetConvergedReason(Intra.ksp, &cr); CHKERRQ(ierr);
        if (cr < 0) {
            ierr = PetscFPrintf(PETSC_COMM_WORLD, stderr, "%s: Intra LSE did not converge: %s\n",
                                currenttime.toStr(4, 6).c_str(), KSPConvergedReasons[cr]); CHKERRQ(ierr);
        }
        
        t_mono_matrix.stop();
        ierr = VecGetArray(Intra.X, &pix); CHKERRQ(ierr);
    }
    
    if (ActivationTime) {
        PetscScalar *at;
        ierr = VecGetArray(ActivationTime, &at); CHKERRQ(ierr);
        for (PetscInt Ii = 0; Ii < EndNodes - StartNodes; ++Ii) {
            if ((at[Ii] == -1.0) && (pix[Ii] >= ActivationThreshold))
                at[Ii] = (double)currenttime;
        }
        ierr = VecRestoreArray(ActivationTime, &at); CHKERRQ(ierr);
    }
    
    ierr = VecRestoreArray(Intra.X, &pix); CHKERRQ(ierr);
}  // acCELLerate::IntraStep

void acCELLerate::PrintTime() {
    if (!mpirank) {
        time_t now;
        time(&now);
        if (difftime(now, t_now) > .1) {
            t_now = now;
            double dtime = difftime(t_now, t_start);
            double prog  = ((double)(currenttime - starttime)+1e-6) / (double)(calclen - starttime);
            fprintf(stdout, "\r%s (%5.1lf%%, %.1fm remaining) ", currenttime.toStr(4, 6).c_str(), prog*100,
                    (dtime/prog-dtime)/60.0);
            fflush(stdout);
        }
    }
}

void acCELLerate::initTimeSteps(acltTime &currenttime, acltTime &intraat, acltTime &saveat, acltTime &backupat,
                                acltTime *extraat) {
    if (starttime != 0_s) {  /* A backup has been loaded */
        if (beginsave > starttime) {
            saveat = beginsave;
        } else {
            saveat = starttime + dtsave - ((starttime - beginsave) % dtsave);
        }
        backupat     = starttime + dtbackup - (starttime % dtbackup);
        intraat      = starttime + dtintra - (starttime % dtintra);
        currenttime += dtcell;
        if (extraat)
            *extraat = starttime + dtextra - (starttime % dtextra);
    } else {
        saveat   = beginsave + currenttime;
        backupat = starttime + dtbackup;
        intraat  = currenttime;
        if (extraat)
            *extraat = currenttime;
    }
    
#if KADEBUG
    std::cerr << "Start time:           " << starttime << std::endl;
    std::cerr << "Next intra step:      " << intraat << std::endl;
    std::cerr << "Next cell model step: " << currenttime << std::endl;
    std::cerr << "Next save time:       " << saveat << std::endl;
    std::cerr << "Next backup time:     " << backupat << std::endl;
#endif  // if KADEBUG
} // acCELLerate::initTimeSteps

void acCELLerate::MonoDomain(acltTime tend, Vec strVec, Vec velVec) {
    PrintDebug("MonoDomain()");
    
    time_t start, now;
    
    acltTime saveat, intraat, backupat;
    initTimeSteps(currenttime, intraat, saveat, backupat);
    
    PetscErrorCode ierr;
    PetscScalar   *pix, *pib;
    PetscScalar   *stretch = NULL, *velocity = NULL;
    
    Vec strInterp = NULL, velInterp = NULL;
    
    if (strVec) {
        if (GaussFwd) {  /* interpolate stretch vec if needed */
            PetscInt siz, siz2;
            ierr = VecGetSize(strVec, &siz); CHKERRQ(ierr);
            ierr = VecGetSize(VecVRes, &siz2); CHKERRQ(ierr);
            if (siz == siz2) {
                ierr = VecGetArray(strVec, &stretch); CHKERRQ(ierr);
            } else if (siz == Intra.size) {
                ierr = VecDuplicate(VecVRes, &strInterp); CHKERRQ(ierr);
                ierr = MatMult(GaussFwd, strVec, strInterp); CHKERRQ(ierr);
                ierr = VecGetArray(strInterp, &stretch); CHKERRQ(ierr);
            } else {
                throw kaBaseException("Stretch vector size does not match.");
            }
        } else {
            ierr = VecGetArray(strVec, &stretch); CHKERRQ(ierr);
        }
    }
    if (velVec) {
        if (GaussFwd) {  /* interpolate velocity vec if needed */
            PetscInt siz, siz2;
            ierr = VecGetSize(velVec, &siz); CHKERRQ(ierr);
            ierr = VecGetSize(VecVRes, &siz2); CHKERRQ(ierr);
            if (siz == siz2) {
                ierr = VecGetArray(velVec, &velocity); CHKERRQ(ierr);
            } else if (siz == Intra.size) {
                ierr = VecDuplicate(VecVRes, &velInterp); CHKERRQ(ierr);
                ierr = MatMult(GaussFwd, velVec, velInterp); CHKERRQ(ierr);
                ierr = VecGetArray(velInterp, &velocity); CHKERRQ(ierr);
            } else {
                throw kaBaseException("Velocity vector size does not match.");
            }
        } else {
            ierr = VecGetArray(velVec, &velocity); CHKERRQ(ierr);
        }
    }
    
    if (!mpirank) {
        time(&t_start);
        t_now = t_start;
    }
    
    // TODO timestpe = min(dtintra, dtcell)
    // also: what to do when dtintra, dtcell not integer multiples of one another
    for (; currenttime <= tend; currenttime += dtcell) {
        PrintTime();
        
        if (backupsavename.size()) {
            if (currenttime >= backupat) {
                backupat += dtbackup;
                WriteBackup(pElphyIntra, pForceIntra, &Intra);
                int n = 0;
                sprintf(prottext, "Wrote backup at %s%n", currenttime.toStr(4, 6).c_str(), &n);
                if (backupat <= tend) {
                    sprintf(prottext+n, ", next backup at %s.", backupat.toStr(4, 6).c_str());
                } else {
                    sprintf(prottext+n, ".");
                }
                WriteProtocol(prottext, 1);
            }
        }
        
        IntraStep(intraat, stretch, velocity);
        
        ierr = VecGetArray(Intra.X, &pix); CHKERRQ(ierr);
        ierr = VecGetArray(Intra.B, &pib); CHKERRQ(ierr);
        
        if (setResults || parametersinSensorfile)
            SaveCellModelVariables(pix, saveat);
        
        SaveSensors(pix, pib);
        ierr = VecRestoreArray(Intra.B, &pib); CHKERRQ(ierr);
        ierr = VecRestoreArray(Intra.X, &pix); CHKERRQ(ierr);
        
        if ((currenttime >= saveat) && exportResults) {
            saveat += dtsave;
            SaveResultIntra();
            sprintf(prottext, "Saved time step %s.", currenttime.toStr(4, 6).c_str());
            WriteProtocol(prottext, 1);
        }
    }
    
    if (strInterp) {
        ierr = VecRestoreArray(strInterp, &stretch); CHKERRQ(ierr);
        ierr = VecDestroy(&strInterp); CHKERRQ(ierr);
    } else if (stretch) {
        ierr = VecRestoreArray(strVec, &stretch); CHKERRQ(ierr);
    }
    
    if (velInterp) {
        ierr = VecRestoreArray(velInterp, &velocity); CHKERRQ(ierr);
        ierr = VecDestroy(&velInterp); CHKERRQ(ierr);
    } else if (velocity) {
        ierr = VecRestoreArray(velVec, &velocity); CHKERRQ(ierr);
    }
    
#ifndef ACLT_NO_MAIN
    ierr = VecDestroy(&Force);
    CHKERRQ(ierr);
#endif // ifndef ACLT_NO_MAIN
    
    if (ActivationTime && (tend == calclen)) {
        std::string tmp(resultprefix + "_" + saveLattice_Name[ResultVars::AT] + ".vec");
        const char *fn = tmp.c_str();
        SaveVec(fn, ActivationTime, ot_bin);
    }
    
    if (!mpirank) {
        cout << "\rTotal time consumption for matrix multiplication  " << t_mono_matrix << endl;
        cout << "Total time consumption for applying conditions    " << t_mono_cond << endl;
        cout << "Total time consumption for cell model calculation " << t_mono_cell << endl;
    }
}  // acCELLerate::MonoDomain

void acCELLerate::BiDomain() {
    PrintDebug("BiDomain()");
    
    PetscErrorCode ierr;
    PetscScalar   *pex, *peb, *pix, *pib, *pvvi, *pres;
    acltTime saveat, intraat, extraat, backupat;
    initTimeSteps(currenttime, intraat, saveat, backupat, &extraat);
    
    if (!mpirank) {
        time(&t_start);
        t_now = t_start;
    }
    
    for (; currenttime <= calclen; currenttime += dtcell) {
        PrintTime();
        
        if (backupsavename.size()) {
            if (currenttime >= backupat) {
                backupat += dtbackup;
                WriteBackup(pElphyIntra, pForceIntra, &Intra);
            }
        }
        
        if (currenttime >= extraat) {
            extraat += dtextra;
            
            if (IntraIndexSet) {
                Vec PhiEIntra;
                ierr = VecGetSubVector(Extra.B, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
                ierr = MatMult(ExtraRHS, Intra.X, PhiEIntra); CHKERRQ(ierr);
                ierr = VecRestoreSubVector(Extra.B, IntraIndexSet, &PhiEIntra); CHKERRQ(ierr);
            } else {
                ierr = MatMult(ExtraRHS, Intra.X, Extra.B); CHKERRQ(ierr);
            }
            
            VecScale(Extra.B, -1.0);
            
            ApplyConditionsExtra();
            Extra.Solve(Extra.B, Extra.X);
            
            KSPConvergedReason cr;
            ierr = KSPGetConvergedReason(Extra.ksp, &cr); CHKERRQ(ierr);
            if (cr < 0) {
                ierr = PetscFPrintf(PETSC_COMM_WORLD, stderr, "%s: Extra LSE did not converge: %s\n",
                                    currenttime.toStr(4, 6).c_str(), KSPConvergedReasons[cr]); CHKERRQ(ierr);
            }
        }
        
        IntraStep(intraat);
        
        ierr = VecGetArray(Intra.X, &pix); CHKERRQ(ierr);
        ierr = VecGetArray(Intra.B, &pib); CHKERRQ(ierr);
        ierr = VecGetArray(Extra.X, &pex); CHKERRQ(ierr);
        ierr = VecGetArray(Extra.B, &peb); CHKERRQ(ierr);
        
        if (setResults || parametersinSensorfile)
            SaveCellModelVariables(pix, saveat);
        
        SaveSensors(pix, pib, pex, peb);
        
        ierr = VecRestoreArray(Extra.X, &pex); CHKERRQ(ierr);
        ierr = VecRestoreArray(Extra.B, &peb); CHKERRQ(ierr);
        ierr = VecRestoreArray(Intra.X, &pix); CHKERRQ(ierr);
        ierr = VecRestoreArray(Intra.B, &pib); CHKERRQ(ierr);
        
        if (currenttime >= saveat) {
            saveat += dtsave;
            SaveResultIntra();
            SaveResultExtra();
        }
    }
    ierr = VecDestroy(&Stim); CHKERRQ(ierr);
}  // acCELLerate::BiDomain

void acCELLerate::TriDomain() {
    PrintDebug("TriDomain()");
    
    acltTime saveat, intraat, extraat, backupat;
    initTimeSteps(currenttime, intraat, saveat, backupat, &extraat);
    
    PetscErrorCode ierr;
    PetscScalar   *pex, *peb, *pix, *pib, *pfx, *pfb, *pvvi, *pvvf, *pres;
    
    for (/*currenttime*/; currenttime <= calclen; currenttime += dtcell) {
        sprintf(prottext, "Time calculated: %s", currenttime.toStr(4, 6).c_str());
        WriteProtocol(prottext);
        
        if (currenttime >= extraat) {
            extraat += dtextra;
            ierr     = MatMult(Intra.A, Intra.X, Intra.B); CHKERRQ(ierr);
            ierr     = MatMult(Fibro.A, Fibro.X, Fibro.B); CHKERRQ(ierr);
            ierr     = VecWAXPY(Extra.B, 1.0, Intra.B, Fibro.B); CHKERRQ(ierr);
            VecScale(Extra.B, -1.0);
            ApplyConditionsExtra();
            Extra.Solve(Extra.B, Extra.X);
        }
        
        bool createIntraBVector = (currenttime >= intraat);
        if (createIntraBVector) {
            intraat += dtintra;
            ierr     = MatMultAdd(Intra.A, Extra.X, Intra.B, Intra.B); CHKERRQ(ierr);
            ierr     = MatMultAdd(Fibro.A, Extra.X, Fibro.B, Fibro.B); CHKERRQ(ierr);
        }
        
        ierr = VecGetArray(Extra.X, &pex); CHKERRQ(ierr);
        ierr = VecGetArray(Extra.B, &peb); CHKERRQ(ierr);
        ierr = VecGetArray(Intra.X, &pix); CHKERRQ(ierr);
        ierr = VecGetArray(Intra.B, &pib); CHKERRQ(ierr);
        ierr = VecGetArray(Fibro.X, &pfx); CHKERRQ(ierr);
        ierr = VecGetArray(Fibro.B, &pfb); CHKERRQ(ierr);
        
        if (VolumeMyo == 0.0) {
            ierr = VecGetArray(VecVolMyo, &pvvi); CHKERRQ(ierr);
            ierr = VecGetArray(VecVolFibro, &pvvf); CHKERRQ(ierr);
        }
        
        if (createIntraBVector) {
            for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
                PetscInt lindex = Ii-StartCells;
                
                PetscScalar volumeintra_local = 0.0;
                if (VolumeMyo == 0.0)
                    volumeintra_local = (pvvi[lindex] < 1e-9 ? 1e-9 : pvvi[lindex]);
                else
                    volumeintra_local = VolumeMyo;
                
                PetscScalar volumefibro_local = 0.0;
                if (VolumeFibro == 0.0)
                    volumefibro_local = (pvvf[lindex] < 1e-9 ? 1e-9 : pvvf[lindex]);
                else
                    volumefibro_local = VolumeFibro;
                
                PetscScalar IMyoFib = BetaMyoFib*(pfx[lindex]-pix[lindex])/RMyoFib;  // phi_i-phi_f==V_i-V_f!
                pib[lindex] -= IMyoFib;
                pfb[lindex] += IMyoFib;
                pib[lindex] *= -pElphyIntra[lindex]->Volume()/volumeintra_local*1e9;  // convert A/m^3 to nA/cell
                pfb[lindex] *= -pElphyFibro[lindex]->Volume()/volumefibro_local*1e9;  // convert A/m^3 to nA/cell
            }
        }
        
        ApplyConditions(pix, pib, CT_II, CT_UI, createIntraBVector);
        ApplyConditions(pfx, pfb, CT_IF, CT_UF, createIntraBVector);
        
        
        if ((setResults == false) && (parametersinSensorfile == false) ) {
            for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
                PetscInt lindex = Ii-StartCells;
                double   dVm    = .0;
                double   force  =
                (EMIntra.CouplingMethod[lindex])(pElphyIntra[lindex], pForceIntra[lindex], (double)dtcell, 0.0, 0.0, dVm,
                                                 pix[lindex], pib[lindex]);
                pix[lindex] += dVm;
                force        =
                (EMFibro.CouplingMethod[lindex])(pElphyFibro[lindex], pForceFibro[lindex], (double)dtcell, 0.0, 0.0, dVm,
                                                 pfx[lindex], pfb[lindex]);
                pfx[lindex] += dVm;
                pix[lindex] += dVm;
                
                if (::isnan(pix[lindex])) {
                    if (mpisize != 1)
                        SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)!", Ii); else
                            throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)", Ii);
                } else if ((pix[lindex] > vm_def_range) || (pix[lindex] < -vm_def_range) ) {
                    if (mpisize != 1)
                        SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is out of range!", Ii); else
                            throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is out of range", Ii);
                }
            }
        } else {
            if (setResults == true) {
                ierr = VecGetArray(VecResultValues, &pres); CHKERRQ(ierr);
            }
            PetscScalar *pm;
            ierr = VecGetArray(MaterialV, &pm); CHKERRQ(ierr);
            for (int j = 0; j < numberParameters; j++) {
                for (PetscInt Ii = StartCells; Ii < EndCells; Ii++) {
                    PetscInt lindex = Ii-StartCells;
                    int m           = (int)pm[lindex];
                    if (j == 0) {
                        double dVm   = .0;
                        double force =
                        (EMIntra.CouplingMethod[lindex])(pElphyIntra[lindex], pForceIntra[lindex], (double)dtcell, 0.0, 0.0, dVm,
                                                         pix[lindex], pib[lindex]);
                        pix[lindex] += dVm;
                        force        =
                        (EMFibro.CouplingMethod[lindex])(pElphyFibro[lindex], pForceFibro[lindex], (double)dtcell, 0.0, 0.0, dVm,
                                                         pfx[lindex], pfb[lindex]);
                        pfx[lindex] += dVm;
                        pix[lindex] += dVm;
                        
                        if (::isnan(pix[lindex])) {
                            if (mpisize != 1)
                                SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)!", Ii); else
                                    throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is not a number (nan)", Ii);
                        } else if ((pix[lindex] > vm_def_range) || (pix[lindex] < -vm_def_range) ) {
                            if (mpisize != 1)
                                SETERRQ1(1, "Transmembrane voltage at " PETSCINT_FORMAT " is out of range!", Ii); else
                                    throw kaBaseException("Transmembrane voltage at " PETSCINT_FORMAT " is out of range", Ii);
                        }
                    } else {}
                    
                    if ((currenttime >= saveat) && (saveresultLattice[j] == true) ) {
                        if (resultPos[j][m] == -1) {
                            resultValue[j] = 0;
                        } else {
                            double CallGetResultValues =
                            (EMIntra.GetResultValues[lindex])(pElphyIntra[lindex], resultPos[j][m], resultValue[j],
                                                              (double)currenttime, pix[lindex]);
                        }
                        pres[lindex] = resultValue[j];
                    }
                    
                    if (parametersinSensorfile) {
                        for (auto pSCPS = ACLTSensors::SC.begin(); pSCPS != ACLTSensors::SC.end(); ++pSCPS) {
                            if (!pSCPS->SensorParameterType.compare(resultName[j]) && (lindex == (pSCPS->SensorIndex-StartCells)) &&
                                (pSCPS->SensorIndex >= StartCells) && (pSCPS->SensorIndex < EndCells) &&
                                ((double)currenttime >= (double)pSCPS->beginsave)) {
                                double CallGetResultValues_sensor =
                                (EMIntra.GetResultValues[lindex])(pElphyIntra[lindex], resultPos[j][m], resultValue[j],
                                                                  (double)currenttime, pix[lindex]);
                                pSCPS->beginsave += pSCPS->dtsave;
                                pSCPS->Save(currenttime, resultValue[j]);
                            }
                        }
                    }
                }
                if (setResults == true)
                    ierr = VecRestoreArray(VecResultValues, &pres); CHKERRQ(ierr);
                
                if ((currenttime >= saveat) && (saveresultLattice[j] == true) ) {
                    SaveResultParameters(resultName[j]);
                }
            }
        }
        
        
        SaveSensors(pix, pib, pex, peb, pfx, pfb);
        
        ierr = VecRestoreArray(Extra.X, &pex); CHKERRQ(ierr);
        ierr = VecRestoreArray(Extra.B, &peb); CHKERRQ(ierr);
        ierr = VecRestoreArray(Intra.X, &pix); CHKERRQ(ierr);
        ierr = VecRestoreArray(Intra.B, &pib); CHKERRQ(ierr);
        ierr = VecRestoreArray(Fibro.X, &pfx); CHKERRQ(ierr);
        ierr = VecRestoreArray(Fibro.B, &pfb); CHKERRQ(ierr);
        if (VolumeMyo == 0.0) {
            ierr = VecRestoreArray(VecVolMyo, &pvvi); CHKERRQ(ierr);
            ierr = VecRestoreArray(VecVolFibro, &pvvf); CHKERRQ(ierr);
        }
        
        if (currenttime >= saveat) {
            saveat += dtsave;
            SaveResultIntra();
            SaveResultExtra();
            SaveResultFibro();
        }
    }
    ierr = VecDestroy(&Stim); CHKERRQ(ierr);
}  // acCELLerate::TriDomain

PrintParameterModus PrintParameterMode = PrintParameterModeOff;

#ifndef ACLT_NO_MAIN

int main(int argc, char *argv[]) {
    PetscErrorCode ierr;
    
    PetscInitialize(&argc, &argv, (char *)0, "");
    
    if (argc < 2) {
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            cerr << argv[0] << " <project file> " << endl;
            cerr << "\tPossible entries in the project file are:" << endl;
            acltProjectFile PF;
            PF.printHelp();
        }
        MPI_Barrier(PETSC_COMM_WORLD);
        ierr = PetscFinalize(); CHKERRQ(ierr);
        exit(-1);
    }
    try {
        acCELLerate ACT;
        ACT.LoadProject(argv[1]);
        ACT.Run();
    } catch (kaBaseException &e) {
        cerr << argv[0] << " Error: " << e << endl;
        ierr = PetscFinalize(); CHKERRQ(ierr);
        exit(-1);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}  // main

#endif  // ifndef ACLT_NO_MAIN
