/*
 * File: CBModelLoaderTetgen.h
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

#ifndef CB_MODEL_LOADER_TETGEN
#define CB_MODEL_LOADER_TETGEN

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>

#include "Matrix3.h"

#include "CBModelLoader.h"
#include "DCCtrl.h"

class CBModelLoaderTetgen : public CBModelLoader
{
public:
    CBModelLoaderTetgen(ParameterMap* parameters) : Base(parameters){};
    virtual ~CBModelLoaderTetgen(){}
    CBModel* Load(std::string Tag = "Mesh");
    CBModel* Load(std::string nodesFilename, std::string elementsFilename, std::string basesFilename, std::string elementType, std::string surfacesFilename, TFloat unit, bool determineNeighbors=true, bool enforceOrthonormalBases=false);
protected:
    virtual void LoadNodes(std::string, TFloat);
    virtual void LoadElements(std::string, std::string, std::string, std::vector<CBElement*>&, bool);
    virtual void LoadElementsNeighbors();
    virtual void LoadSurfaces(std::string, std::vector<CBElement*>&);
    
private:
    typedef CBModelLoader   Base;
    TInt numNodes_;
    TInt numElements_;
    TInt nodesPerElement_;
};

#endif
