/*
 *  CBElementFactory.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 04.06.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
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
