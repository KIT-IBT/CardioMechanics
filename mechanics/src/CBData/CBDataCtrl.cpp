//
//  CBDataCtrl.cpp
//  CardioMechanics_unstable
//
//  Created by Thomas Fritz on 12.03.13.
//
//

#include "CBDataCtrl.h"

void CBDataCtrl::Init(TInt size)
{
    val_.resize(size,0);
    lastVal_.resize(size,0);
}
