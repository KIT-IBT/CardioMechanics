/*
 * File: PETScLSEConditions.h
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



#ifndef PETSCLSECONDITIONS_H
#define PETSCLSECONDITIONS_H

enum ConditionType { CT_I, CT_U };

#include "PETScConditions.h"

//! Base class for conditions

class BaseConditionPETScLSE {
 public:
  BaseConditionPETScLSE() {}

  virtual ~BaseConditionPETScLSE() {}

  double val;  //!< potential or current
  void Init();
};

void BaseConditionPETScLSE::Init() {
  val = 0;
}

template<class ValTyp, class BaseCondition>
class SingleCondition : public BaseCondition {
 public:
  SingleCondition() {BaseCondition::Init();}

  ~SingleCondition() {BaseCondition::Init();}

  long position;
  ConditionType cflag;
  void ScanLine(const char *buf);

  void FreeColumnEntries() {BaseConditionPETScLSE::Init();}
};

typedef PETScConditions<double, BaseConditionPETScLSE, ConditionType> PETScLSEConditions;


/*!
 * Different implementation from acltConditions.h because possible ConditionTypes are different
 */

template<class ValTyp, class BaseCondition>
void SingleCondition<ValTyp, BaseCondition>::ScanLine(const char *buf) {
  double C, C1, C2, C3;
  long   pos;
  char   type[256];

  int rc = sscanf(buf, "%ld %lf %s", &pos, &C, type);

#if KADEBUG
  fprintf(stderr, "ScanLine '%ld' '%lf' '%s'\n", pos, C, type);
#endif  // if KADEBUG

  if (rc != 3)
    throw ConditionError("Wrong word count (!=3) in condition line");

  position = pos;

  if (type[0] == 'I')
    cflag = CT_I;
  else if (type[0] == 'V')
    cflag = CT_U;
  else
    throw ConditionError("Condition type is unknown");

  this->val = C;
}

#endif  // ifndef PETSCLSECONDITIONS_H
