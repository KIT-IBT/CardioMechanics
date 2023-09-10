/*! \file acltConditions.h
 \brief Managing conditions for acCELLerate simulations
 
 \author gs, IBT - Universitaet Karlsruhe (TH) on Dec 2007, ported from MultiDomainConditions.h
 */


#ifndef ACLTCONDITIONS_H
#define ACLTCONDITIONS_H

enum ConditionType { CT_II = 1, CT_IE = 2, CT_IF = 4, CT_UI = 8, CT_UE = 16, CT_UF = 32 };
enum ConditionStatus { CS_nothing_to_do, CS_apply, CS_applied, CS_remove };

#include "PETScConditions.h"
#include "acltTime.h"

//! Base class for conditions in the bidomain case

class BaseConditionACLT {
public:
    BaseConditionACLT() {
        cs            = CS_nothing_to_do;
        ColumnIndizes = NULL;
        ColumnValues  = NULL;
    }
    
    virtual ~BaseConditionACLT() {}
    
    acltTime cl;   //!< cyclelength of stimulus
    acltTime toff;  //!< temporal offset of stimulus
    acltTime sl;   //!< stimulus length
    double val;    //!< potential or current
    
    ConditionStatus cs;
    
    PetscInt numEntries;      //!< Amount of entries in a matrix row for this boundary condition
    PetscInt *ColumnIndizes;  //!< Column index of each entry
    double *ColumnValues;     //!< Value of each entry for column index
    inline void Init();
};

void BaseConditionACLT::Init() {
    numEntries = -1;
    if (ColumnIndizes)
        delete ColumnIndizes;
    ColumnIndizes = NULL;
    
    if (ColumnValues)
        delete ColumnValues;
    ColumnValues = NULL;
}

template<class ValTyp, class BaseCondition>
class SingleCondition : public BaseCondition {
public:
    SingleCondition() {BaseCondition::Init();}
    
    ~SingleCondition() {BaseCondition::Init();}
    
    PetscInt position;
    ConditionType cflag;
    void ScanLine(const char *buf);
    void SetColumnEntries(PetscInt, PetscInt *, double *);
    bool GetColumnEntries(PetscInt *, const PetscInt **, const double **);
    
    void FreeColumnEntries() {BaseCondition::Init();}
};

typedef PETScConditions<double, BaseConditionACLT, ConditionType> ACLTConditions;


/*!
 * the function is here overwritten, because BaseCondition is larger than in
 * the normal case. the function sets the coordinates, the amplitude, length,
 * cyclelength, temporal offset and the kind of impulse
 */

template<class ValTyp, class BaseCondition>
void SingleCondition<ValTyp, BaseCondition>::ScanLine(const char *buf) {
    // TODO switch to input stream
    double C;
    char   C1[64], C2[64], C3[64];
    PetscInt pos;
    char type[256];
    
    int rc = sscanf(buf, PETSCINT_FORMAT " %lf %63s %63s %63s %255s", &pos, &C, C1, C2, C3, type);
    
#if KADEBUG
    fprintf(stderr, "ScanLine '" PETSCINT_FORMAT "' '%lf' '%s' '%s' '%s' '%s'\n", pos, C, C1, C2, C3, type);
#endif  // if KADEBUG
    
    if (rc != 6)
        throw ConditionError("Condition line includes not enough words");
    
    position = pos;
    
    if (!strcasecmp(type, "Ie"))
        cflag = CT_IE;
    else if (!strcasecmp(type, "Ii"))
        cflag = CT_II;
    else if (!strcasecmp(type, "If"))
        cflag = CT_IF;
    else if (!strcasecmp(type, "Ve"))
        cflag = CT_UE;
    else if (!strcasecmp(type, "Vm"))
        cflag = CT_UI;
    else if (!strcasecmp(type, "Vf"))
        cflag = CT_UF;
    else
        throw ConditionError("Condition type is unknown");
    
    this->val  = C;
    this->cl   = C1;
    this->sl   = C2;
    this->toff = C3;
}  // SingleCondition<ValTyp, BaseCondition>::ScanLine

template<class ValTyp, class BaseCondition>
void SingleCondition<ValTyp, BaseCondition>::SetColumnEntries(PetscInt nEntries, PetscInt *cIndizes, double *cValues) {
    this->numEntries    = nEntries;
    this->ColumnIndizes = new PetscInt[nEntries];
    this->ColumnValues  = new double[nEntries];
    
    for (PetscInt i = 0; i < nEntries; i++) {
        this->ColumnIndizes[i] = *cIndizes++;
        this->ColumnValues[i]  = *cValues++;
    }
}

template<class ValTyp, class BaseCondition>
bool SingleCondition<ValTyp, BaseCondition>::GetColumnEntries(PetscInt *nEntry, const PetscInt **cIndizes,
                                                              const double **cValues) {
    if (this->numEntries == -1)
        return false;
    
    *nEntry   = this->numEntries;
    *cIndizes = (const PetscInt *)this->ColumnIndizes;
    *cValues  = (const double *)this->ColumnValues;
    
    return true;
}

#endif  // ifndef ACLTCONDITIONS_H
