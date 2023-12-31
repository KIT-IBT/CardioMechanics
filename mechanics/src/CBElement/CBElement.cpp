/*
 * File: CBElement.cpp
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




#include "CBElement.h"
#include "CBElementAdapter.h"
#include "CBElementSurface.h"

CBElement::CBElement(CBElement& other) {
    index_ = other.index_;
    localIndex_ = other.localIndex_;
    materialIndex_ = other.materialIndex_;
    neighbors_ = other.neighbors_;
    adapter_ = other.adapter_;
    parameters_ = other.parameters_;
    material_ = other.material_;
    properties = other.properties;
}


void CBElement::PrintNodesIndices()
{
    std::cout << GetNodeIndex(0);
    for(unsigned int i = 0; i < GetNumberOfNodesIndices(); i++)
        std::cout << " " << GetNodeIndex(i);
}

bool CBElement::CheckFiberLength() {
    return false;
}


void CBElement::SetProperty(ElementProperty property) {
    properties |= property;
}

bool CBElement::IsProperty(ElementProperty property) {
    return (properties && property);
}

void CBElement::DeleteProperty(ElementProperty property) {
    properties &= ~property;
}

