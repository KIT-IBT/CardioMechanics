/*
 * File: CBDataCtrl.cpp
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

#include "CBDataCtrl.h"

void CBDataCtrl::Init(TInt size)
{
    val_.resize(size,0);
    lastVal_.resize(size,0);
}
