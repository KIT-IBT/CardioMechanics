/*
 * File: CBModelLoaderTetgen.cpp
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


#include "CBModelLoaderTetgen.h"
#include "sstream"
#include "cstring"

#include <unordered_set>
#include <algorithm>     // intersection

CBModel *CBModelLoaderTetgen::Load(std::string Tag) {
    std::string nodesFilename    = Base::parameters_->Get<std::string>(Tag+".Tetgen.Nodes");
    std::string elementsFilename = Base::parameters_->Get<std::string>(Tag+".Tetgen.Elements");
    std::string basesFilename    = Base::parameters_->Get<std::string>(Tag+".Tetgen.Bases", "");
    std::string elementType      = Base::parameters_->Get<std::string>(Tag+".Type");
    std::string surfacesFilename = Base::parameters_->Get<std::string>(Tag+".Tetgen.Surfaces", "");
    TFloat   unit                = Base::parameters_->Get<TFloat>(Tag+".Tetgen.Unit", 1.0);
    bool determineNeighbors      = Base::parameters_->Get<bool>(Tag+".Tetgen.DetermineNeighbors", true);
    bool enforceOrthonormalBases = Base::parameters_->Get<bool>(Tag+".Tetgen.EnforceOrthonormalBases", false);
    
    this->Load(nodesFilename, elementsFilename, basesFilename, elementType, surfacesFilename, unit, determineNeighbors,
               enforceOrthonormalBases);
    
    return Base::model_;
}

CBModel *CBModelLoaderTetgen::Load(std::string nodesFilename, std::string elementsFilename, std::string basesFilename,
                                   std::string elementType, std::string surfacesFilename, TFloat unit,
                                   bool determineNeighbors, bool enforceOrthonormalBases) {
    std::vector<CBElement *> newElements;
    
    LoadNodes(nodesFilename, unit);
    LoadElements(elementsFilename, basesFilename, elementType, newElements, enforceOrthonormalBases);
    LoadSurfaces(surfacesFilename, newElements);
    
    model_->SetElements(newElements);
    model_->CheckElementsIndexes();
    model_->ExtractSolidElements();
    Base::model_->ExtractSurfaceElements();
    
    Base::model_->DetermineNeighbors();
    

    
    return Base::model_;
}

void CBModelLoaderTetgen::LoadNodes(std::string Filename, TFloat unit) {
    // Format Mesh.node:
    // <amount of nodes> <Dimension (must be 3)> <Boundary Conditions> (0: no boundary condition, 1: Fixation)> <BoundaryMarker (must be 0) >
    // <Node Index> <x> <y> <z> <Boundary Condition in decimal number: Binary: 001 (dec: 1) -> x-direction fixes, 010 (dec: 2) -> y-direction fixes, 100 (dec: 4) -> z-direction fixed, 111 (dec: 1+2+4 = 7) all direction fixed >
    
    std::ifstream file(Filename.c_str());
    
    if (!file.good()) {
        file.close();
        throw std::runtime_error("CBModelLoaderTetgen::LoadNodes(): File " + Filename + " is not accessable !");
    } else {
        TInt         numNodes;
        TInt         currentNode;
        unsigned int dim;
        unsigned int attr;
        unsigned int boundaryMarker;
        
        file >> numNodes;
        file >> dim;
        file >> attr;
        file >> boundaryMarker;
        
        numNodes_ = numNodes;
        std::vector<Vector3<TFloat>> nodes;
        std::vector<TInt> boundaryConditions;
        
        nodes.reserve(numNodes);
        boundaryConditions.reserve(numNodes);
        
        if (boundaryMarker != 0) {
            throw std::runtime_error(
                                     "CBModelLoaderTetgen::LoadNodes(): Boundary Marker (col 4 in node header) not yet implemented ! Maybe you don't have 4 numbers in the header, or the last one is not zero?");
        }
        if (attr > 1) {
            throw std::runtime_error(
                                     "CBModelLoaderTetgen::LoadNodes(): Just one attribute is implemented yet (dirichlet boundary conditions) !!!");
        }
        if (dim != 3)
            throw std::runtime_error("CBModelLoaderTetgen::LoadNodes(): Only for 3D model !");
        
        std::string line;
        std::string currentNodeAsString;
        TFloat x, y, z;
        int i = 0;
        while (getline(file, line)) {
            std::istringstream linestream(line);
            linestream >> currentNodeAsString;
            
            // check for comment lines and empty lines
            if ( (std::string(&currentNodeAsString[0]) != std::string("#")) &&
                (std::string(&currentNodeAsString[0]) != std::string("")) ) {
                currentNode = std::atoi(currentNodeAsString.c_str());
                if (i+1 != currentNode) {
                    std::cout << line << std::endl;
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadNodes(): Node indices are not consistent! Do they increase one by one, beginning from 1?");
                }
                
                // read x, y, z
                if (!(linestream >> x >> y >> z) ) {
                    std::cout << line << std::endl;
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadNodes(): Something went wrong with reading a line from the nodes file " + Filename +
                                             "! Did you check if all lines have a sufficient number of coordinates?");
                }
                Vector3<TFloat> node(unit * x, unit *y, unit *z);
                nodes.push_back(node);
                
                // set fixation, read bc if available
                TInt bc;
                if (attr > 0) {
                    if (!(linestream >> bc) ) {
                        std::cout << line << std::endl;
                        throw std::runtime_error(
                                                 "CBModelLoaderTetgen::LoadNodes(): Something went wrong with reading a line from the nodes file " + Filename +
                                                 "! Do all lines contain fixation information?");
                    }
                    boundaryConditions.push_back(bc);
                } else {
                    boundaryConditions.push_back(0);
                }
                
                i++;
            }
        }
        if (i != numNodes) {
            throw std::runtime_error(
                                     "CBModelLoaderTetgen::LoadNodes(): Error, the expected number of nodes does not correspond to the number of lines in " + Filename +
                                     " !");
        }
        
        model_->SetNodes(std::move(nodes));
        model_->SetBoundaryConditions(std::move(boundaryConditions));
    }
    
    file.close();
} // CBModelLoaderTetgen::LoadNodes

void CBModelLoaderTetgen::LoadElements(std::string elementsFilename, std::string basesFilename, std::string elementType,
                                       std::vector<CBElement *> &newElements, bool enforceOrthonormalBases) { // CBElement<T> Class Factory
                                                                                                              // Format Mesh.ele
                                                                                                              // <amount of elements> <Nodes per Element> (must be 4 or 10)> <Attributes> (must be 1: Material class)>
                                                                                                              // <Element Index> <node 1> <node 2> ... <node n> <material index>
    
    std::ifstream elementsFile;
    std::ifstream basesFile;
    bool     anisotropy = false;
    
    if (basesFilename != "") {
        // Format Mesh.bases
        // <amount of elements>
        // <element index> <fiber_x> <fiber_y> <fiber_z> <sheet_x> <sheet_y> <sheet_z> <sheetnormal_x> <sheetnormal_y> <sheetnormal_z>
        basesFile.open(basesFilename.c_str());
        if (!basesFile.good()) {
            basesFile.close();
            throw std::runtime_error("CBModelLoaderTetgen::LoadElements(): File " + basesFilename + " is not accessable !");
        } else {
            anisotropy = true; // bases filename exists and can be opened [lb451]
        }
    }
    
    elementsFile.open(elementsFilename.c_str());
    if (!elementsFile.good()) {
        elementsFile.close();
        throw std::runtime_error("CBModelLoaderTetgen::LoadElements(): File " + elementsFilename + " is not accessable !");
    } else {
        TInt         numElements         = 0;
        TInt         numBases            = 0;
        TInt         numQuadraturePoints = 0;
        TInt         currentElement      = 0;
        unsigned int nodesPerElement     = 0;
        unsigned int attr                = 0;
        
        
        elementsFile >> numElements;
        elementsFile >> nodesPerElement;
        elementsFile >> attr;
        
        numElements_ = numElements;
        nodesPerElement_ = nodesPerElement;
        
        if (anisotropy) {
            std::string str;
            getline(basesFile, str);
            std::stringstream ss(str);
            ss >> numBases;
            ss >> numQuadraturePoints;
            
            if (numBases != numElements) {
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadElements(): Elements file " + elementsFilename + " and bases file " + basesFilename +
                                         " don't fit to each other !\n");
            }
        }
        
        if ((nodesPerElement != 4) && (nodesPerElement != 10) )
            throw std::runtime_error("CBModelLoaderTetgen::LoadElements(): Unkown Volume Element Type !");
        
        if (attr != 1)
            throw std::runtime_error("CBModelLoaderTetgen::LoadElements(): Exactly one attribute (Material index) is needed !");
        
        for (TInt i = 0; i < numElements; i++) {
            if (elementsFile.eof()) {
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadElements(): File " + elementsFilename +
                                         " is corrupt. Check number of elements !");
            }
            if (basesFile.eof()) {
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadElements(): File " + basesFilename +
                                         " is corrupt. Check number of elements !");
            }
            elementsFile >> currentElement;
            CBElement *element;
            
            element = Base::model_->GetElementFactory()->New(elementType);
            
            if (element->GetNumberOfNodesIndices() != nodesPerElement) {
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadElements(): Mesh Type: " + elementType + " does not fit to mesh file: " + elementsFilename +
                                         "!");
            }
            
            element->SetIndex(currentElement - 1);
            
            for (unsigned int j = 0; j < nodesPerElement; j++) {
                TInt node;
                elementsFile >> node;
                element->SetNodeIndex(j, node - 1);
            }
            if (attr == 1) {
                unsigned int a;
                elementsFile >> a;
                element->SetMaterialIndex(a);
            }
            
            Matrix3<TFloat> basis;
            
            if (anisotropy) {
                TInt currentElementBases = 0;
                std::string line;
                getline(basesFile, line);
                std::stringstream ss(line);
                ss >> currentElementBases;
                TFloat a = 0;
                std::vector<TFloat> vals;
                
                while (ss >> a)
                    vals.push_back(a);
                
                if (vals.size() == 0)
                    continue;
                
                bool spherical = false;
                
                if (vals.size() == 2*numQuadraturePoints)
                    spherical = true;
                else if (vals.size() == 9 * numQuadraturePoints)
                    spherical = false;
                else
                    std::runtime_error("CBModelLoaderTetgen::LoadElements(): Bases file: " + basesFilename+ " is corrupt !");
                
                if ( (numQuadraturePoints >= element->GetNumberOfQuadraturePoints()) || (numQuadraturePoints == 1)) {
                    // for(int n=0; n <  numQuadraturePoints; n++) {
                    for (int n = 0; n < std::min(numQuadraturePoints, element->GetNumberOfQuadraturePoints()); n++) {
                        Vector3<TFloat> b1;
                        Vector3<TFloat> b2;
                        Vector3<TFloat> b3;
                        
                        if (!spherical) {
                            b1 = {vals.at(n*9+0), vals.at(n*9+1), vals.at(n*9+2)};
                            b2 = {vals.at(n*9+3), vals.at(n*9+4), vals.at(n*9+5)};
                            b3 = {vals.at(n*9+6), vals.at(n*9+7), vals.at(n*9+8)};
                        } else {
                            TFloat theta = vals.at(2*numQuadraturePoints*n);
                            TFloat phi = vals.at(2*numQuadraturePoints*n+1);
                            b1 = {sin(theta)*cos(phi),  sin(theta)*sin(phi), cos(theta)};
                            b2 = {-sin(phi),  cos(phi), 0};
                            b3 = {-cos(theta)*cos(phi), -cos(theta)*sin(phi), sin(theta)};
                        }
                        if (enforceOrthonormalBases) {
                            /// do Gram Schmidt orthonormalization
                            b1.Normalize();
                            b2 = b2 - (b1*b2) * b1;
                            b2.Normalize();
                            b3 = b3 - (b1*b3) * b1 - (b2*b3) * b2;
                            b3.Normalize();
                        }
                        
                        Matrix3<TFloat> basis(b1, b2, b3);
                        
                        if (basis.Det() < 0) {
                            basis(0, 2) *= -1;
                            basis(1, 2) *= -1;
                            basis(2, 2) *= -1;
                        }
                        TFloat  precision = parameters_->Get<TFloat>("Solver.Precision", 1e-8);
                        if ((basis.Det() < 1 - precision) &&  (basis.Det() > 1 + precision)) {
                            std::cout << "CBModelLoaderTetgen::LoadElements(): Bases file: Element: " << element->GetIndex() <<
                            " QP: " << n <<
                            " Does not have a high enough precision, try enforcing orthonormal basis matrix <EnforceOrthonormalBase>true</EnforceOrthonormalBase>";
                        }
                        element->SetBasisAtQuadraturePoint(n, basis);
                    }
                    
                    // both, 1 and 5 fit to t10 elements and work with them [lb451]
                    if ((numQuadraturePoints == 1) && (element->GetNumberOfQuadraturePoints() > numQuadraturePoints) ) {
                        for (int n = 1; n < element->GetNumberOfQuadraturePoints(); n++)
                            element->SetBasisAtQuadraturePoint(n, *(element->GetBasisAtQuadraturePoint(0)));
                    }
                } else {
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadElements(): Number of quadrature points doesn't fit to element type and isn't one, which would be okay ...");
                }
            } else {
                basis.SetToIdentityMatrix();
                element->SetBasis(basis);
                for (int n = 0; n < element->GetNumberOfQuadraturePoints(); n++)
                    element->SetBasisAtQuadraturePoint(n, basis);
            }
            element->SetParameters(parameters_);
            newElements.push_back(element);
        }
    }
    
    for (auto n = 0; n < newElements.size(); n++) {
        if (newElements.at(n)->CheckFiberLength()) {
            std::cout << "Warning from CBModelLoaderTetgen::LoadElement on process " << DCCtrl::GetProcessID() <<
            ": Found inappropriate fiber length in element " << n << ". " << "\n";
        }
    }
    
    elementsFile.close();
} // CBModelLoaderTetgen::LoadElements

void CBModelLoaderTetgen::LoadElementsNeighbors() { // CBElement<T> Class Factory
    throw std::runtime_error(
                             "CBModelLoaderTetgen::LoadElementsNeighbors is no longer supported for tetgen and handled by CBModel.");
    
    std::string neighborsFilename = Base::parameters_->Get<std::string>("Mesh.Tetgen.Neighbors", "");
    if (neighborsFilename != "") {
        std::ifstream neighborsFile;
        neighborsFile.open(neighborsFilename.c_str());
        
        if (!neighborsFile.good()) {
            neighborsFile.close();
            throw std::runtime_error(
                                     "CBModelLoaderTetgen::LoadElementsNeighbors(): File " + neighborsFilename + " is not accessable !");
        } else {
            TInt         numElements    = 0;
            TInt         currentElement = 0;
            TInt         neighbor;
            unsigned int neighborsPerElement = 0;
            
            neighborsFile >> numElements;
            if (numElements != Base::model_->GetElements().size()) {
                std::string elementsFilename =  Base::parameters_->Get<std::string>("Mesh.Tetgen.Elements");
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadElementsNeighbors(): File " + neighborsFilename + " does not fit to elements file " + elementsFilename +
                                         " !");
            }
            neighborsFile >> neighborsPerElement;
            
            for (TInt i = 0; i < numElements; i++) {
                if (neighborsFile.eof()) {
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadElementsNeighbors(): File " + neighborsFilename +
                                             " is corrupt. Check number of elements !");
                }
                
                neighborsFile >> currentElement;
                std::vector<TInt> currentElementsNeighbors;
                for (unsigned int j = 0; j < neighborsPerElement; j++) {
                    neighborsFile >> neighbor;
                    if (neighbor != -1)
                        currentElementsNeighbors.push_back(neighbor - 1);
                }
            }
        }
        neighborsFile.close();
    }
} // CBModelLoaderTetgen::LoadElementsNeighbors

void CBModelLoaderTetgen::LoadSurfaces(std::string surfacesFilename, std::vector<CBElement *> &newElements) { // CBElement<T> Class Factory
                                                                                                              // Format:
                                                                                                              // <amount of elements> <Nodes per Element> (must be 3)> <Attributes> (min 2, max 3)>
                                                                                                              // <surface element Index> <node 1> <node 2> <node 3> <surface index>
    
    if (surfacesFilename != "") {
        std::ifstream surfacesFile;
        surfacesFile.open(surfacesFilename.c_str());
        
        if (!surfacesFile.good()) {
            surfacesFile.close();
            throw std::runtime_error("CBModelLoaderTetgen::LoadSurfaces(): File " + surfacesFilename +
                                     " is not accessable !");
        } else {
            TInt         surfaceElements        = 0;
            TInt         currentSurfaceElement  = 0;
            unsigned int nodesPerSurfaceElement = 0;
            unsigned int attr                   = 0;
            
            std::string str = "";
            
            while (str.size() < 1)
                getline(surfacesFile, str);
            
            std::stringstream ss(str);
            
            ss >> surfaceElements;
            ss >> nodesPerSurfaceElement;
            ss >> attr;
            
            if (attr < 2) {
                throw std::runtime_error(
                                         "CBModelLoaderTetgen::LoadSurfaces(): Attributes has to be two (surface material & surface index) or three (optional surface traction scaling)");
            }
            
            std::vector<TInt> nodesIndexes;
            
            std::vector<TFloat> attributes;
            attributes.resize(attr);
            
            for (TInt i = 0; i < surfaceElements; i++) {
                if (surfacesFile.eof()) {
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadSurfaces(): File " + surfacesFilename +
                                             " is corrupt. Check number of elements !");
                }
                
                std::string str;
                
                while (str.size() < 1)
                    getline(surfacesFile, str);
                
                std::stringstream ss(str);
                std::vector<TFloat> v;
                while (ss) {
                    TFloat t = 0;
                    ss >> t;
                    v.push_back(t);
                }
                
                
                nodesPerSurfaceElement = v.size() - attr - 2;
                nodesIndexes.resize(nodesPerSurfaceElement);
                if ((nodesPerSurfaceElement != 3) && (nodesPerSurfaceElement != 6)) {
                    throw std::runtime_error(
                                             "CBModelLoaderTetgen::LoadSurfaces(): Only 3 or 6 nodes per surfaces are supported yet");
                }
                
                currentSurfaceElement = static_cast<TInt>(v.at(0));
                
                for (unsigned int j = 0; j < nodesPerSurfaceElement; j++)
                    nodesIndexes.at(j) = static_cast<TInt>(v.at(j+1));
                
                
                for (unsigned int j = 0; j < attr; j++)
                    attributes.at(j) =  v.at(j+1+nodesPerSurfaceElement);
                
                CBElement *element;
                std::string surfaceType;
                if (Base::parameters_->IsAvailable("Mesh.Surfaces.Surface_" +
                                                   std::to_string(static_cast<TInt>(attributes.at(1))) + ".Type")) {
                    surfaceType =
                    Base::parameters_->Get<std::string>("Mesh.Surfaces.Surface_" +
                                                        std::to_string(static_cast<TInt>(attributes.at(1))) + ".Type");
                } else {
                    throw std::runtime_error("CBModelLoaderTetgen::LoadSurfaces(): [Mesh.Surfaces.Surface_" +
                                             std::to_string(static_cast<TInt>(attributes.at(
                                                                                            1))) +
                                             ".Type] does not exists in the parameters file");
                }
                
                
                element = Base::model_->GetElementFactory()->New(surfaceType);
                if (element->GetNumberOfNodesIndices() == nodesPerSurfaceElement) {
                    element->SetIndex(newElements.size());
                    element->SetMaterialIndex(static_cast<TInt>(attributes.at(0)));
                    dynamic_cast<CBElementSurface *>(element)->SetSurfaceIndex(static_cast<TInt>(attributes.at(1)));
                    dynamic_cast<CBElementSurface *>(element)->SetSurfaceElementIndex(currentSurfaceElement);
                    if (attr == 3) {
                        dynamic_cast<CBElementSurface *>(element)->SetSurfaceTractionScaling(attributes.at(2));
                    }
                    for (unsigned int j = 0; j < nodesPerSurfaceElement; j++)
                        element->SetNodeIndex(j, nodesIndexes.at(j) - 1);
                    element->SetParameters(parameters_);
                    newElements.push_back(element);
                } else if ((element->GetNumberOfNodesIndices() == 3) && (nodesPerSurfaceElement == 6)) {
                    for (int t = 0; t < 4; t++) {
                        int triangleNodes[3];
                        if (t == 0) {
                            triangleNodes[0] = 0;
                            triangleNodes[1] = 3;
                            triangleNodes[2] = 5;
                        } else if (t == 1) {
                            triangleNodes[0] = 3;
                            triangleNodes[1] = 1;
                            triangleNodes[2] = 4;
                        } else if (t == 2) {
                            triangleNodes[0] = 3;
                            triangleNodes[1] = 4;
                            triangleNodes[2] = 5;
                        } else if (t == 3) {
                            triangleNodes[0] = 5;
                            triangleNodes[1] = 4;
                            triangleNodes[2] = 2;
                        }
                        
                        element->SetIndex(newElements.size());
                        element->SetMaterialIndex(attributes.at(0));
                        dynamic_cast<CBElementSurface *>(element)->SetSurfaceIndex(attributes.at(1));
                        dynamic_cast<CBElementSurface *>(element)->SetSurfaceElementIndex(currentSurfaceElement);
                        element->SetNodeIndex(0, nodesIndexes.at(triangleNodes[0]) - 1);
                        element->SetNodeIndex(1, nodesIndexes.at(triangleNodes[1]) - 1);
                        element->SetNodeIndex(2, nodesIndexes.at(triangleNodes[2]) - 1);
                        newElements.push_back(element);
                    }
                } else {
                    throw std::runtime_error("CBModelLoaderTetgen::LoadSurfaces(): [Mesh.Surfaces.Surface_" +
                                             std::to_string(attributes.at(
                                                                          1)) +
                                             ".Type] does not fit to the amount of node indices in the surface file");
                }
            }
        }
        surfacesFile.close();
    }
} // CBModelLoaderTetgen::LoadSurfaces
