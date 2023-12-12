/*
 * File: CBStatus.cpp
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


#include <stdio.h>
#include <map>
#include <string>
#include "CBStatus.h"

std::string CBStatusToStr(CBStatus s)
{
    std::map<CBStatus,std::string> str =
    {
        {CBStatus::DACCORD,"DACCORD"},
        {CBStatus::WAITING,"WAITING"},
        {CBStatus::NOT_INITIALIZED,"NOT_INITIALIZED"},
        {CBStatus::PREPARING_SIMULATION,"INITIALIZE"},
        {CBStatus::SUCCESS,"SUCCESS"},
        {CBStatus::RUNNING,"RUNNING"},
        {CBStatus::NOTHING_DONE,"NOTHING_DONE"},
        {CBStatus::INTERRUPTED,"INTERRUPTED"},
        {CBStatus::REPEAT,"REPEAT"},
        {CBStatus::FAILED,"FAILED"},
        {CBStatus::INFINITIVE,"INFINITIVE"},
        {CBStatus::NOT_A_NUMBER,"NOT_A_NUMBER"},
        {CBStatus::CORRUPT_ELEMENT,"CORRUPT_ELEMENT"},
        {CBStatus::CRITICAL_STRESS,"CRITICAL_STRESS"},
        {CBStatus::SYSTEM_WRONG_STATE,"SYSTEM_WRONG_STATE"}
    };
    return str[s];
}
