/*
 * File: DCVectorPETSc.h
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


#ifndef DC_VECTOR_PETSC
#define DC_VECTOR_PETSC

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <petscksp.h>

#include "DCVector.h"

class DCVectorPETSc : public DCVectorBase
{
public:
    DCVectorPETSc(){}
    virtual ~DCVectorPETSc();
    Vec GetVec(){return vec_;}
protected:
private:
    
    void SetType(DCVectorType type);
    DCVectorType GetType();
    
    void SetLocalSize(PetscInt size);
    void SetGlobalSize(PetscInt size);
    PetscInt GetLocalSize();
    PetscInt GetGlobalSize();
    void SetGhostIndices(PetscInt numGhostIndices, const PetscInt* ghostIndices);
    
    void DuplicateTo(DCVectorBase* v);
    void CopyValuesTo(DCVectorBase* v);
    
    void Assemble();
    void Build();
    
    void InsertValues(PetscInt numValues, const PetscInt* indices, const PetscScalar* values);
    void InsertValuesLocal(PetscInt numValues, const PetscInt* indices, const PetscScalar* values);
    
    void AddValues(PetscInt numValues, const PetscInt* indices, const PetscScalar* values);
    void AddValuesLocal(PetscInt numValues, const PetscInt* indices, const PetscScalar* values);
    
    void GetValues(PetscInt numValues, const PetscInt* indices, PetscScalar* values);
    void GetValuesLocal(PetscInt numValues, const PetscInt* indices, PetscScalar* values);
    
    void SetLocalToGlobalMapping(PetscInt* mapping);
    void ConvertLocalToGlobalIndices(PetscInt numIndices, PetscInt* indices);
    
    PetscInt localSize_ = 0;
    PetscInt globalSize_ = 0;
    DCVectorType type_ = VEC_PARALLEL;
    
    PetscInt numGhostIndices_ = 0;
    const PetscInt* ghostIndices_ = 0;
    
    bool isBuild_ = false;
    bool isGhost_ = false;
    
    Vec vec_ = 0;
    Vec localForm_ = 0;
    
    typedef DCCtrl Ctrl;
};



#endif
