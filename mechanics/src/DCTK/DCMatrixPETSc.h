/*
 * File: DCMatrixPETSc.h
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


#ifndef DC_MATRIX_PETSC
#define DC_MATRIX_PETSC

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <petscksp.h>
#include "DCMatrix.h"

class DCCtrl;
class DCMatrixBase;

class DCMatrixPETSc : public DCMatrixBase
{
public:
    DCMatrixPETSc(){}
    virtual ~DCMatrixPETSc();
    Mat GetMat(){return mat_;}
protected:
private:
    
    void SetType(DCMatrixType type);
    DCMatrixType GetType();
    
    void SetNumberOfLocalRows(PetscInt numLocalRows) override;
    void SetNumberOfGlobalRows(PetscInt numGlobalRows) override;
    void SetNumberRows(PetscInt numLocalRows, PetscInt numGlobalRows) override;
    
    void SetNumberOfLocalCols(PetscInt numLocalCols) override;
    void SetNumberOfGlobalCols(PetscInt numGlobalCols) override;
    void SetNumberOfCols(PetscInt numLocalCols, PetscInt numGlobalCols)override;
    
    void SetNumberOfNonzeros(PetscInt numNzDiag, PetscInt numNzOffDiag)override;
    void SetNonzeros(PetscInt* nonZerosDiag, PetscInt* nonZerosDiagOff)override;
    
    void DuplicateTo(DCMatrixBase* m);
    void DuplicateFrom(DCMatrixBase* m);
    
    void CopyValuesTo(DCMatrixBase* m);
    void CopyValuesFrom(DCMatrixBase* m);
    
    void Assemble();
    void Build();
    
    PetscInt numLocalRows_ = 0;
    PetscInt numGlobalRows_  = 0;
    
    PetscInt numLocalCols_ = 0;
    PetscInt numGlobalCols_ = 0;
    
    PetscInt numNzDiag_ = 0;
    PetscInt numNzOffDiag_ = 0;
    
    PetscInt* nonZerosDiag_ = 0;
    PetscInt* nonZerosOffDiag_ = 0;
    
    DCMatrixType type_ = MAT_PARALLEL;
    
    bool hasBeenBuilt_ = false;
    
    Mat mat_ = 0;
    
    typedef DCCtrl Ctrl;
};



#endif
