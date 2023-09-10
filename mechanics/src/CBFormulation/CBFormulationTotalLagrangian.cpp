/*
 *  CBFormulationTotalLagrangian.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 14.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBFormulationTotalLagrangian.h"
#include "CBSolver.h"

CBStatus CBFormulationTotalLagrangian::CalcNodalForces()
{
    CBStatus rc=CBStatus::NOTHING_DONE;
    for(auto& it : Base::solver_->GetSolidElementVector())
    {
        rc = it->CalcNodalForces();
        if(rc != CBStatus::SUCCESS)
            return(rc);
    }
    return(rc);
}

CBStatus CBFormulationTotalLagrangian::CalcNodalForcesJacobian()
{
    CBStatus rc=CBStatus::NOTHING_DONE;
    for(auto& it : Base::solver_->GetSolidElementVector())
    {
        rc = it->CalcNodalForcesJacobian();
        if(rc != CBStatus::SUCCESS)
            return(rc);
    }
    return(rc);
}
