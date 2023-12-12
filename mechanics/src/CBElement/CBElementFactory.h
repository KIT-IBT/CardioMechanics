/*
 * File: CBElementFactory.h
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

#ifndef CB_ELEMENT_FACTORY_H
#define CB_ELEMENT_FACTORY_H

#include <vector>
#include <string>
#include <list>
#include <set>
#include <map>
#include <memory>
#include <functional>

#include "CBElement.h"

class CBElementFactory {
public:
    using FactoryFunction = std::function<CBElement*()>;
    
    CBElementFactory();
    ~CBElementFactory() {}
    
    CBElement* New(std::string elementType);
    
protected:
private:
    std::map<std::string, FactoryFunction> producers_;
};

#endif
