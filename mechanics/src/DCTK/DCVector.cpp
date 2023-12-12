/*
 * File: DCVector.cpp
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


#include "DCCtrl.h"
#include "DCVector.h"
#include "DCVectorPETSc.h"

DCVector::DCVector()
{
    v_ = DCCtrl::NewVector();
}

DCVector::DCVector(const TInt globalSize)
{
    v_ = DCCtrl::NewVector();
    v_->SetGlobalSize(globalSize);
}

DCVector::DCVector(const TInt localSize, const TInt globalSize)
{
    v_ = DCCtrl::NewVector();
    v_->SetGlobalSize(globalSize);
    v_->SetLocalSize(globalSize);
}

DCVector::~DCVector()
{
    delete v_;
}


DCVector::DCVector(const DCVector& v) : v_(0)
{
    v.v_->DuplicateTo(v_);
    v.v_->CopyValuesTo(v_);
}

void DCVector::operator=(const DCVector& v)
{
    v.v_->DuplicateTo(v_);
    v.v_->CopyValuesTo(v_);
}


Vec DCVector::Petsc()
{
    DCVectorPETSc* p = dynamic_cast<DCVectorPETSc*>(v_);
    if(p)
        return p->GetVec();
    else
        throw std::runtime_error("Vec DCVector::Petsc(): DCCtrl has not been initialized with Petsc");
}
