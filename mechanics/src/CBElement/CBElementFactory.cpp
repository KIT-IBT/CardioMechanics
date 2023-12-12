/*
 * File: CBElementFactory.cpp
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


#include "CBElementFactory.h"

#include "CBElement.h"
#include "CBElementSolid.h"
#include "CBElementSolidT4.h"
#include "CBElementSolidT10.h"
#include "CBElementSurface.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"
#include "CBElementCavity.h"
#include "CBElementContactMaster.h"
#include "CBElementContactSlave.h"
#include "CBElementContactRobin.h"


CBElementFactory::CBElementFactory() {
    producers_["T4"]      = CBElementSolidT4::New;
    producers_["T10"]     = CBElementSolidT10::New;
    producers_["CAVITY"]  = CBElementSurfaceT3::New;
    producers_["T3"]      = CBElementSurfaceT3::New;
    producers_["T6"]      = CBElementSurfaceT6::New;
    producers_["CONTACT_ROBIN"] = CBElementContactRobin::New;
    producers_["CONTACT_MASTER"]  = CBElementContactMaster::New;
    producers_["CONTACT_SLAVE"]   = CBElementContactSlave::New;
}

CBElement *CBElementFactory::New(std::string elementType) {
    CBElement *result = nullptr;
    
    auto it = producers_.find(elementType);
    
    if (it == producers_.end()) {
        throw std::runtime_error(std::string("Unkown element type: [" + elementType + "] You might have to extend CBElementFactory \n Available elements are: \n")
                                 + "\t T4\t:4-node Iso-P1 tetrahedral element,\n"
                                 + "\t T10\t:10-node Iso-P2 tetrahedral element,\n"
                                 + "\t CONTACT_MASTER\t:3-node triangle surface element for contact problems\n"
                                 + "\t CONTACT_SLAVE\t:3-node triangle surface element for contact problems\n"
                                 + "\t CONTACT_ROBIN\t:3-node triangle surface element for Robin boundary condition\n"
                                 + "\t CAVITY\t:3-node triangle surface element for circulatory system plugin\n"
                                 + "\t T3\t: 3-node triangle surface element\n"
                                 + "\t T6\t: 6-node triangle surface element\n");
    } else {
        result = it->second();
    }
    return result;
}
