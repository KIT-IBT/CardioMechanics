/*! \file acltCellModel.cpp
 \brief
 
 \author gs, IBT - Universitaet Karlsruhe Dec 2007. Copyright 2007 IBT Universitaet Karlsruhe (TH). All rights
 reserved.
 */

#ifndef CELLMODELSTRUCT_H
#define CELLMODELSTRUCT_H

#include <acltCellModel.h>


// ***********CMS**************

CMS::CMS() {
#if KADEBUG
    cerr<<"CMS::CMS()"<<endl;
#endif  // if KADEBUG
    pep    = NULL;
    prepem = NULL;
    pfp    = NULL;
    prepfm = NULL;
    clear();
}

void CMS::clear() {
#if KADEBUG
    cerr<<"CMS::clear()"<<endl;
#endif  // if KADEBUG
    emd = "";
    fmd = "";
    
    emt = EMT_Dummy;
    if (pep)
        delete pep;
    pep = NULL;
    if (prepem)
        delete prepem;
    prepem = NULL;
    
    fmt = FMT_Dummy;
    if (pfp)
        delete pfp;
    pfp = NULL;
    if (prepfm)
        delete prepfm;
    prepfm = NULL;
    
    VmInit    = 0.0;
    ForceInit = 0.0;
    material  = 0;
}  // CMS::clear

// ***********CellModelStruct**************

void CellModelStruct::LoadCellModels(string cmfilename) {
#if KADEBUG
    cerr<<"CellModelStruct::LoadCellModels("<< cmfilename <<")"<<endl;
#endif  // if KADEBUG
    
    ifstream cm(cmfilename.c_str());
    if (!cm.is_open())
        throw kaBaseException("Cell model file %s does not exist", cmfilename.c_str());
    
    CMS cms;
    
    while (!cm.eof()) {
        char c = cm.peek();
        if (c == '#') {
            do
                cm.get(c);
            while (!cm.eof() && c != '\n');
            continue;
        }
        
        cm >> cms.material >> cms.emd >> cms.fmd;
        cellmodels.push_back(cms);
#if KADEBUG
        cerr<<"Read line: " << cms.material << " " << cms.emd << " " << cms.fmd << endl;
#endif  // if KADEBUG
    }
    cm.close();
    
    char temptext[256];
    vector<CMS>::iterator iCM, bCM = cellmodels.begin(), eCM = cellmodels.end();
    for (iCM = bCM; iCM != eCM; iCM++) {
        int modtyp = 0;
        if ((*iCM).emd.size() > 0) {
            for (modtyp = (int)EMT_Dummy; modtyp <= (int)EMT_Last; modtyp++) {
                sprintf(temptext, "%d", modtyp);
                if (((*iCM).emd == FindElphyModelDescriptor((ElphyModelType)modtyp)) || ((*iCM).emd == temptext)) {
                    (*iCM).emt = (ElphyModelType)modtyp;
                    (*iCM).emd = FindElphyModelFileName((ElphyModelType)modtyp);
                    break;
                }
                if (modtyp == (int)EMT_Last)
                    (*iCM).emt = FindElphyModelFileType((*iCM).emd.c_str());
            }
        }
        if ((*iCM).fmd.size() > 0) {
            for (modtyp = (int)FMT_Dummy; modtyp <= (int)FMT_Last; modtyp++) {
                sprintf(temptext, "%d", (int)modtyp);
                if (( (*iCM).fmd == FindForceModelDescriptor((ForceModelType)modtyp)) || ( (*iCM).fmd == temptext) ) {
                    (*iCM).fmt = (ForceModelType)modtyp;
                    (*iCM).fmd = FindForceModelFileName((ForceModelType)modtyp);
                    break;
                }
                if (modtyp == (int)FMT_Last)
                    (*iCM).fmt = FindForceModelFileType((*iCM).fmd.c_str());
            }
        }
    }
    
    
#if KADEBUG
    cerr<<"CellModelStruct::LoadCellModels finished"<<endl;
#endif  // if KADEBUG
}  // CellModelStruct::LoadCellModels

void CellModelStruct::InitParameters(double tinc) {
#if KADEBUG
    cerr<<"CellModelStruct::InitParameters("<< tinc <<")"<<endl;
#endif  // if KADEBUG
    
    vector<CMS>::iterator iCM, bCM = cellmodels.begin(), eCM = cellmodels.end();
    for (iCM = bCM; iCM != eCM; iCM++) {
        initElphyParameters<double>(&(*iCM).pep, (*iCM).emd.c_str(), (*iCM).emt, tinc);
        initForceParameters<double>(&(*iCM).pfp, (*iCM).fmd.c_str(), (*iCM).fmt);
    }
    
#if KADEBUG
    cerr<<"CellModelStruct::InitParameters finished"<<endl;
#endif  // if KADEBUG
}

void CellModelStruct::setElphyModel(vbElphyModel<double> **pem, int mat) {
#if KADEBUG == 2
    cerr<<"CellModelStruct::initElphyModel()"<<endl;
#endif  // if KADEBUG == 2
    vector<CMS>::iterator iCM, bCM = cellmodels.begin(), eCM = cellmodels.end();
    for (iCM = bCM; iCM != eCM; iCM++) {
        if ((*iCM).material == mat) {
            initElphyModel<double>(pem, (*iCM).pep, (*iCM).emt);
            return;
        }
    }
    throw kaBaseException("CellModelStruct::initElphyModel: No material defined in cell model file for material %d", mat);
}

void CellModelStruct::setForceModel(vbForceModel<double> **pfm, int mat) {
#if KADEBUG == 2
    cerr<<"CellModelStruct::initForceModel()"<<endl;
#endif  // if KADEBUG == 2
    vector<CMS>::iterator iCM, bCM = cellmodels.begin(), eCM = cellmodels.end();
    for (iCM = bCM; iCM != eCM; iCM++) {
        if ((*iCM).material == mat) {
            initForceModel<double>(pfm, (*iCM).pfp, (*iCM).fmt);
            return;
        }
    }
    throw kaBaseException("CellModelStruct::initForceModel: No material defined in cell model file for material %d", mat);
}

// ***********EMCoupling**************

void EMCoupling::Delete() {
    if (CouplingMethod) {
        delete CouplingMethod;
        CouplingMethod = NULL;
    }
    if (GetResultParameters) {
        delete GetResultParameters;
        GetResultParameters = NULL;
    }
    if (GetResultValues) {
        delete GetResultValues;
        GetResultValues = NULL;
    }
}

void EMCoupling::New(unsigned long size) {
    Delete();
    CouplingMethod      = new CalculateCoupling[size];
    GetResultParameters = new GetParameters[size];
    GetResultValues     = new GetValues[size];
}

void EMCoupling::Set(int index, vbElphyModel<double> *pem, vbForceModel<double> *pfm, bool forceset) {
    GetResultParameters[index] = &GetResPara;
    GetResultValues[index]     = &GetResVal;
    if (forceset) {
        if (pem->OutIsTCa() && pfm->InIsTCa())
            CouplingMethod[index] = &TroponinCoupling;
        else
            CouplingMethod[index] = &CalciumCoupling;
    } else {
        CouplingMethod[index] = &NoCoupling;
    }
}

// ***********HeterogeneCellModel**************

void HeterogeneCellModel::Load(string hcmfilename) {
#if KADEBUG
    cerr<<"HeterogeneCellModel::Load "<<hcmfilename<<endl;
#endif  // if KADEBUG
    ifstream hcmf(hcmfilename.c_str());
    if (!hcmf.is_open())
        throw kaBaseException("HeterogeneCellModel::Load: cannot open File: %s \n", hcmfilename.c_str());
    
    string token, args;
    while (!hcmf.eof()) {
        hcmf>>token;
        getline(hcmf, args);
        while (args[0] == ' ' || args[0] == '\n')
            args = args.substr(1);
        
        ifstream hfile(token.c_str());
        if (hfile.is_open()) {
            HCM tempHCM;
            tempHCM.filename = token;
            tempHCM.variable = args;
            hetCM.push_back(tempHCM);
        } else {
            throw kaBaseException("HeterogeneCellModel::Load: Unknown filename: %s \n", token.c_str());
        }
    }
}  // HeterogeneCellModel::Load

void HeterogeneCellModel::Set(vbElphyModel<double> **pem, PetscInt numCells) {
#if KADEBUG
    cerr<<"HeterogeneCellModel::Set "<<endl;
#endif  // if KADEBUG
    vector<HCM>::iterator iHCM, bHCM = hetCM.begin(), eHCM = hetCM.end();
    for (iHCM = bHCM; iHCM != eHCM; iHCM++) {
        PetscErrorCode ierr;
        PetscInt startCell, endCell;
        
        Vec hetCMV;
        LoadVec( (*iHCM).filename.c_str(), hetCMV, ot_bin);
        
        ierr = VecGetOwnershipRange(hetCMV, &startCell, &endCell); CHKERRQ(ierr);
        
        if (numCells != endCell-startCell) {
            SETERRQ1(-42, "Heterogeneous vector has wrong size, local size %d, expected %d.", endCell-startCell, numCells);
        }
        
        PetscScalar *pm;
        ierr = VecGetArray(hetCMV, &pm); CHKERRQ(ierr);
        
        string hcmvariable = (*iHCM).variable;
        PetscInt endfor    = endCell-startCell;
        for (PetscInt lindex = 0; lindex < endfor; lindex++)
            pem[lindex]->AddHeteroValue(hcmvariable, pm[lindex]);
        
        ierr = VecRestoreArray(hetCMV, &pm); CHKERRQ(ierr);
        ierr = VecDestroy(&hetCMV); CHKERRQ(ierr);
    }
}  // HeterogeneCellModel::Set

#endif  // ifndef CELLMODELSTRUCT_H
