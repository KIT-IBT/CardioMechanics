//
//  DCCtrlPETSc.h
//  CardioMechanics
//
//  Created by Thomas Fritz on 29.10.11.
//  Copyright (c) 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
//

#ifndef DC_CONTROLLER_PETSC
#define DC_CONTROLLER_PETSC

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>

#include <petscksp.h>
#include "DCCtrl.h"
class DCCtrl;
class DCVectorBase;
class DCMatrixBase;

class DCCtrlPETSc : public DCCtrl
{
public:
    DCCtrlPETSc(){}
    virtual ~DCCtrlPETSc()
    {
        PetscFinalize();
    }
    
    
    // Petsc specific functions, access with DCCtrl::Petsc::...
    static void CreateVector(PetscInt n,PetscInt N, Vec* v);
    static void CreateSeqVector(PetscInt n,Vec* v);
    static void CreateMatrix(PetscInt n, PetscInt m,PetscInt N, PetscInt M, PetscInt numNonZeros, Mat* v);
    static void CreateSeqMatrix(PetscInt N, PetscInt M, PetscInt numNonZeros, Mat* v);
    static MPI_Comm Comm(){return PETSC_COMM_WORLD;}
    
    static std::string SNESReasonToString(int i);
    static std::string KSPReasonToString(int i);
    
protected:
private:
    
    int InitImpl(int argc, char** argv) override;
    
    unsigned int GetProcessIDImpl() override;
    unsigned int GetNumberOfProcessesImpl() override;
    
    bool IsProcessZeroImpl() override;
    bool IsParallelImpl() override;
    
    DCVectorBase* NewVectorImpl() override;
    DCMatrixBase* NewMatrixImpl() override;
    
    void GatherToZeroImpl(std::vector<int>& a) override;
    void GatherToZeroImpl(std::vector<long>& a)override;
    void GatherToZeroImpl(std::vector<float>& a)override;
    void GatherToZeroImpl(std::vector<double>& a)override;
    void GatherToZeroImpl(std::vector<std::string>& a)override;
    
    void CopyFromZeroToAllImpl(std::vector<int>& a)override;
    void CopyFromZeroToAllImpl(std::vector<long>& a)override;
    void CopyFromZeroToAllImpl(std::vector<float>& a)override;
    void CopyFromZeroToAllImpl(std::vector<double>& a)override;
    
    void WeightedAverageImpl(double& localAverage, double& localWeight, double& globalAverage)override;
    
    typedef DCCtrl  Base;
    
    PetscMPIInt mpiRank_  = 0;
    PetscMPIInt mpiSize_  = 0;
    
    static std::map<int,std::string> snesReasonStr_;
    static std::map<int,std::string> kspReasonStr_;
};

typedef DCCtrlPETSc Petsc;

#endif
