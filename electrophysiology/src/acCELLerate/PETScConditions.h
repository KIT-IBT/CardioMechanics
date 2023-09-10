/*
 **      Name
 **              PETScConditions.h
 **
 **
 **      History
 **              18.12.07 -gs ported from Conditions.h
 **
 **  created at IBT - Universitaet Karlsruhe (TH)
 */

#ifndef PETSCCONDITIONS_H
#define PETSCCONDITIONS_H

#include <list>
#include <string>

const int MaxConditions = 400000;

typedef kaBaseException ConditionError;

template<class ValTyp, class BaseCondition> class SingleCondition;

template<class ValTyp, class BaseCondition, class ConditionType>
class PETScConditions {
protected:
    std::string file;
    typedef std::list<SingleCondition<ValTyp, BaseCondition>> ConditionContainerType;
    ConditionContainerType SC;
    
public:
    PETScConditions();
    virtual ~PETScConditions();
    void Load(const char *name);
    void Parse(PetscInt);
    virtual void Parse(PetscInt xLattice, PetscInt yLattice, PetscInt zLattice);
    
    inline int GetnumOfConditions() const {return SC.size();}
};


template<class ValTyp, class BaseCondition, class ConditionType>
PETScConditions<ValTyp, BaseCondition, ConditionType>::PETScConditions() {
    // SC.reserve(MaxConditions); /* Initialize vector with estimated size to usually avoid reallocation */
}

template<class ValTyp, class BaseCondition, class ConditionType>
PETScConditions<ValTyp, BaseCondition, ConditionType>::~PETScConditions()
{}

template<class ValTyp, class BaseCondition, class ConditionType>
void PETScConditions<ValTyp, BaseCondition, ConditionType>::Load(const char *name) {
#if KADEBUG
    fprintf(stderr, "PETScConditions::Load %s\n", name);
#endif  // if KADEBUG
    
    file = name;
}

template<class ValTyp, class BaseCondition, class ConditionType>
void PETScConditions<ValTyp, BaseCondition, ConditionType>::Parse(PetscInt maxpoints) {
#if KADEBUG
    fprintf(stderr, "PETScConditions::Parse(%ld)\n", (long)maxpoints);
#endif  // if KADEBUG
    
    if (file.empty())
        return;
    
    ifstream conditions(file.c_str());
    if (!conditions.good()) {
        throw kaBaseException("Unable to read conditions from file %s.", file.c_str());
    }
    string line = "";
    while (std::getline(conditions, line) ) {
        SingleCondition<ValTyp, BaseCondition> s;
        s.ScanLine(line.c_str());
        
        if ((s.position < 0) || (s.position >= maxpoints) )
            throw ConditionError("PETScCondition position %d is outside of range of data %d", s.position, maxpoints);
        
#if KADEBUG > 1
        fprintf(stderr, "PETScConditions::Parse %d.) %d\n", this->GetnumOfConditions(), (int)s.position);
#endif  // if KADEBUG
        
        SC.push_back(s);
    }
    conditions.close();
    
#if KADEBUG
    fprintf(stderr, "PETScConditions::Init finished Conditions# %d\n", this->GetnumOfConditions());
#endif  // if KADEBUG
}  // PETScConditions<ValTyp, BaseCondition, ConditionType>::Parse

#if !defined(ACLTCONDITIONS_H) && !defined(PETSCLSECONDITIONS_H)

template<class ValTyp, class BaseCondition, class ConditionType>
void PETScConditions<ValTyp, BaseCondition, ConditionType>::Parse(PetscInt xLattice, PetscInt yLattice,
                                                                  PetscInt zLattice) {
# if KADEBUG
    fprintf(stderr, "Conditions::Parse(%ld,%ld,%ld)\n", (long)xLattice, (long)yLattice, (long)zLattice);
# endif  // if KADEBUG
    
    ifstream conditions(file.c_str());
    string   line;
    while (std::getline(conditions, line) ) {
        SingleCondition<ValTyp, BaseCondition> s;
        s.ScanLine(line.c_str());
        
        if ((s.cx < 0) || (s.cx >= xLattice) || (s.cy < 0) || (s.cy >= yLattice) || (s.cz < 0) || (s.cz >= zLattice) )
            throw ConditionError("Condition coordinates %d %d %d are outside of lattice with dimension %d x %d x %d ", s.cx,
                                 s.cy, s.cz, xLattice, yLattice, zLattice);
        
        if (s.dx > xLattice)
            s.dx = xLattice;
        if (s.dy > yLattice)
            s.dy = yLattice;
        if (s.dz > zLattice)
            s.dz = zLattice;
        
# if KADEBUG
        fprintf(stderr, "Conditions::Init %d.) %d - %d %d - %d %d - %d\n",
                this->GetnumOfConditions(), s.cx, s.dx, s.cy, s.dy, s.cz, s.dz);
# endif  // if KADEBUG
        
        SC.push_back(s);
    }
    conditions.close();
    
# if KADEBUG
    fprintf(stderr, "Conditions::Init finished Conditions# %d\n", this->GetnumOfConditions());
# endif  // if KADEBUG
}  // PETScConditions<ValTyp, BaseCondition, ConditionType>::Parse

#else  // if !defined(ACLTCONDITIONS_H) && !defined(PETSCLSECONDITIONS_H)

template<class ValTyp, class BaseCondition, class ConditionType>
void PETScConditions<ValTyp, BaseCondition, ConditionType>::Parse(PetscInt xLattice, PetscInt yLattice,
                                                                  PetscInt zLattice)
{}

#endif  // if !defined(ACLTCONDITIONS_H) && !defined(PETSCLSECONDITIONS_H)

#endif  // ifndef PETSCCONDITIONS_H
