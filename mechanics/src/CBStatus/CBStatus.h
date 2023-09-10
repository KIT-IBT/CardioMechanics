/*
 *  CBStatus.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 17.08.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_STATUS_H
#define CB_STATUS_H

#include <map>

enum class CBStatus
{
    DACCORD = 6,
    WAITING = 5,
    NOT_INITIALIZED      = 4,
    PREPARING_SIMULATION = 3,
    SUCCESS              = 2,
    RUNNING              = 1,
    NOTHING_DONE         = 0,
    INTERRUPTED          = -1,
    REPEAT               = -2,
    FAILED               = -3,
    INFINITIVE           = -4,
    NOT_A_NUMBER         = -5,
    CORRUPT_ELEMENT      = -6,
    CRITICAL_STRESS      = -7,
    SYSTEM_WRONG_STATE   = -8
};
std::string CBStatusToStr(CBStatus s);
#endif
