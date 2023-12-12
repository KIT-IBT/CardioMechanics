/*
 * File: CBStatus.h
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
