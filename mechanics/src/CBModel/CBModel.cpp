/*
 *  CBModel.cpp
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 13.09.11.
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

extern "C" double dgesvd_(const char *, const char *, int *, int *, double *, int *, double *, double *, int *,
                          double *, int *, double *, int *, int *);

// #endif

#include <algorithm>
#include <unordered_set>

#include "CBModel.h"
#include "CBElement.h"
#include "CBElementSolid.h"
#include "CBElementSurface.h"
#include "DCCtrl.h"

void CBModel::InitNodes(TInt numNodes) {
    nodes_.resize(numNodes);
}

void CBModel::InitMapping() {
    forwardMapping_.resize(nodes_.size());
    backwardMapping_.resize(nodes_.size());
    
    for (int i = 0; i < nodes_.size(); i++) {
        forwardMapping_.at(i) = i;
        backwardMapping_.at(i) = i;
    }
}

void CBModel::CalcNonZeroes() {
    nodeNeighborsNnz_.assign(nodes_.size()*3, 1);
    std::map<TInt, std::set<TInt>> SurIndexMap;
    for (auto &ele : elements_) {
        CBElementSurface *sele = dynamic_cast<CBElementSurface *>(ele);
        
        for (unsigned int i = 0; i < ele->GetNumberOfNodesIndices(); ++i) {
            TInt n = ele->GetNodeIndex(i);
            nodeNeighborsNnz_.at(n*3) += ele->GetNumberOfNodesIndices() - 1;
            nodeNeighborsNnz_.at(n*3+1) += ele->GetNumberOfNodesIndices() - 1;
            nodeNeighborsNnz_.at(n*3+2) += ele->GetNumberOfNodesIndices() - 1;
            if (sele) {
                if (sele->GetType() == "CAVITY") {
                    SurIndexMap[sele->GetSurfaceIndex()].insert(n);
                }
                if (sele->GetType() == "CONTACT_MASTER") {
                    SurIndexMap[sele->GetSurfaceIndex()].insert(n);
                }
                if (sele->GetType() == "CONTACT_SLAVE") {
                    SurIndexMap[sele->GetSurfaceIndex()].insert(n);
                }
                if (sele->GetType() == "CONTACT_ROBIN") {
                    SurIndexMap[sele->GetSurfaceIndex()].insert(n);
                }
            }
        }
    }
    for (auto SurIndex : SurIndexMap) {
        TInt NumNodes = (int)SurIndex.second.size();
        for (auto NodeIndex : SurIndex.second) {
            nodeNeighborsNnz_.at(NodeIndex*3) += NumNodes;
            nodeNeighborsNnz_.at(NodeIndex*3+1) += NumNodes;
            nodeNeighborsNnz_.at(NodeIndex*3+2) += NumNodes;
        }
    }
    
    for (int i = 0; i < nodeNeighborsNnz_.size(); i++) {
        nodeNeighborsNnz_[i] *= 3;
        nodeNeighborsNnz_[i] = nodeNeighborsNnz_[i] > maxNnz_ ? maxNnz_ : nodeNeighborsNnz_[i];
    }
} // CBModel::CalcNonZeroes

void CBModel::InitNodesComponentsBoundaryConditions() {
    nodesComponentsBoundaryConditions_.resize(nodes_.size());
}

void CBModel::ExtractSolidElements() {
    solidElements_.clear();
    for (auto &e : elements_) {
        // if(e->IsProperty(ElementProperty::Solid)) // [ss029] This does not work yet.
        if (dynamic_cast<CBElementSolid *>(e) != 0)
            solidElements_.push_back(static_cast<CBElementSolid *>(e));
    }
}

void CBModel::ExtractSurfaceElements() {
    surfaceElements_.clear();
    for (auto &e : elements_) {
        // if(e->IsProperty(ElementProperty::Surface)) // [ss029] This does not work yet.
        if (dynamic_cast<CBElementSurface *>(e) != 0)
            surfaceElements_.push_back(static_cast<CBElementSurface *>(e));
    }
}

void CBModel::UpdateNode(TInt pos, const Vector3<TFloat> &node) {
    nodes_[pos] = node;
}

void CBModel::CheckElementsIndexes() {
    for (TInt i = 0; i < elements_.size(); i++) {
        if (elements_.at(i)->GetIndex() != i) {
            throw std::runtime_error(
                                     "CBModel::CheckElementsIndexes(): Index of element does not fit to index of elements container !!!");
        }
    }
}

void CBModel::PrintNodes() {
    for (TInt i = 0; i < nodes_.size(); i++) {
        std::cout << "Node " << i << ":\t\t";
        nodes_[i].Print();
    }
}

void CBModel::PrintElements() {
    for (TInt i = 0; i < elements_.size(); i++) {
        std::cout << "Element: " << i << ":\t\t";
        elements_[i]->PrintNodesIndices();
        std::cout << "\n";
    }
}

void CBModel::DetermineNeighbors() {
    DetermineNodeSolidElements();
    DetermineNodeSurfaceElements();
    DetermineSolidElementNeighbors();
    DetermineSurfaceElementNeighbors();
}

void CBModel::DetermineNodeSolidElements() {
    // nodeSolidElements_ contains for each node all solid elements that contain this node
    nodeSolidElements_.resize(nodes_.size());
    
    for (auto &e : solidElements_) {
        for (TInt i = 0; i < e->GetNumberOfNodesIndices(); i++) {
            TInt nodeId = e->GetNodeIndex(i);
            nodeSolidElements_.at(nodeId).insert(e->GetIndex());
        }
    }
}

void CBModel::DetermineNodeSurfaceElements() {
    // nodeSurfaceElements_ contains for each node all surface elements that contain this node
    nodeSurfaceElements_.resize(nodes_.size());
    
    for (auto &e : surfaceElements_) {
        for (TInt i = 0; i < e->GetNumberOfNodesIndices(); i++) {
            TInt nodeId = e->GetNodeIndex(i);
            nodeSurfaceElements_.at(nodeId).insert(e->GetIndex());
        }
    }
}

/// For each solid element, determine all neighboring solid elements
void CBModel::DetermineSolidElementNeighbors() {
    // solidElementNeighborsCommonNode_ contains for each solid element all solid elements that have 1 node in common
    // solidElementNeighborsCommonFace_ contains for each solid element all solid elements that have 1 face in common
    
    if (nodeSolidElements_.size() == 0) {
        throw std::runtime_error(
                                 "CBModel::DetermineSolidElementNeighbors: DetermineNodeSolidElements needs to be called first.");
    }
    
    for (auto &e : solidElements_) {
        solidElementNeighborsCommonNode_.insert(std::pair<TInt, std::unordered_set<TInt>>(e->GetIndex(),
                                                                                          std::unordered_set<TInt>()));
        solidElementNeighborsCommonFace_.insert(std::pair<TInt, std::unordered_set<TInt>>(e->GetIndex(),
                                                                                          std::unordered_set<TInt>()));
    }
    
    for (TInt nodeId = 0; nodeId < nodes_.size(); nodeId++) {
        std::unordered_set<TInt> eleCandidates = nodeSolidElements_.at(nodeId);
        
        // compare two elements each
        // this is fast, eleCandidates.size() does not depend on mesh size
        for (auto ii : eleCandidates) {
            CBElement *e1 = elements_.at(ii);
            std::vector<int> v1; // point ids of e1
            for (int i = 0; i < e1->GetNumberOfNodesIndices(); i++)
                v1.push_back(e1->GetNodeIndex(i));
            
            for (auto jj : eleCandidates) {
                if (ii == jj)
                    continue;
                
                CBElement *e2 = elements_.at(jj);
                std::vector<int> v2; // point ids of e2
                for (int j = 0; j < e2->GetNumberOfNodesIndices(); j++)
                    v2.push_back(e2->GetNodeIndex(j));
                
                // compare number of common nodes per element
                // compute set intersection
                std::sort(v1.begin(), v1.end());
                std::sort(v2.begin(), v2.end());
                std::vector<int> v_intersection;
                std::set_intersection(v1.begin(), v1.end(),
                                      v2.begin(), v2.end(),
                                      std::back_inserter(v_intersection));
                
                // version 1: All elements that have 1 node in common are neighbors.
                // [this is used currently nowhere, but might regularize better than version 2]
                solidElementNeighborsCommonNode_.at(ii).insert(jj);
                solidElementNeighborsCommonNode_.at(jj).insert(ii);
                
                // version 2: All elements that have 1 face in common are neighbors.
                // T10 elements: 6 neighbors in common, T4 elements: 3 neighbors in common
                // [this is used by tf120 in a laplace connectivity matrix for regularization of the inverse tension estimatior CBSolverActiveStressEstimator]
                if (((v_intersection.size() >= 6) && (e1->GetNumberOfNodesIndices() == 10) ) ||
                    ((v_intersection.size() >= 3) && (e1->GetNumberOfNodesIndices() == 4))) { // there should be a better way to distinguish between T10 and T4
                    solidElementNeighborsCommonFace_.at(ii).insert(jj);
                    solidElementNeighborsCommonFace_.at(jj).insert(ii);
                }
            }
        }
    }
    
    // PrintSolidElementNeighbors();
} // CBModel::DetermineSolidElementNeighbors

/// For each surface element, determine all neighboring surface elements
void CBModel::DetermineSurfaceElementNeighbors() {
    // surfaceElementNeighborsCommonNode_ contains for each surface element all surface elements that have 1 node in common
    // surfaceElementNeighborsCommonEdge_ contains for each surface element all surface elements that have 1 edge in common
    
    if (nodeSurfaceElements_.size() == 0) {
        throw std::runtime_error(
                                 "CBModel::DetermineSurfaceElementNeighbors: DetermineNodeSurfaceElements needs to be called first.");
    }
    
    for (auto &e : surfaceElements_) {
        surfaceElementNeighborsCommonNode_.insert(std::pair<TInt, std::unordered_set<TInt>>(e->GetIndex(),
                                                                                            std::unordered_set<TInt>()));
        surfaceElementNeighborsCommonEdge_.insert(std::pair<TInt, std::unordered_set<TInt>>(e->GetIndex(),
                                                                                            std::unordered_set<TInt>()));
    }
    
    for (TInt nodeId = 0; nodeId < nodes_.size(); nodeId++) {
        std::unordered_set<TInt> eleCandidates = nodeSurfaceElements_.at(nodeId);
        
        // compare two elements each
        // this is fast, eleCandidates.size() does not depend on mesh size
        for (auto ii : eleCandidates) {
            CBElement *e1 = elements_.at(ii);
            std::vector<int> v1; // point ids of e1
            for (int i = 0; i < e1->GetNumberOfNodesIndices(); i++)
                v1.push_back(e1->GetNodeIndex(i));
            
            for (auto jj : eleCandidates) {
                if (ii == jj)
                    continue;
                
                CBElement *e2 = elements_.at(jj);
                std::vector<int> v2; // point ids of e2
                for (int j = 0; j < e2->GetNumberOfNodesIndices(); j++)
                    v2.push_back(e2->GetNodeIndex(j));
                
                // compare number of common nodes per element
                // compute set intersection
                std::sort(v1.begin(), v1.end());
                std::sort(v2.begin(), v2.end());
                std::vector<int> v_intersection;
                std::set_intersection(v1.begin(), v1.end(),
                                      v2.begin(), v2.end(),
                                      std::back_inserter(v_intersection));
                
                // version 1: All elements that have 1 node in common are neighbors.
                surfaceElementNeighborsCommonNode_.at(ii).insert(jj);
                surfaceElementNeighborsCommonNode_.at(jj).insert(ii);
                
                // version 2: All elements that have 1 edge in common are neighbors.
                // T6 elements: 3 neighbors in common, T3 elements: 2 neighbors in common
                if (((v_intersection.size() >= 3) && (e1->GetNumberOfNodesIndices() == 6)) ||
                    ((v_intersection.size() >= 2) && (e1->GetNumberOfNodesIndices() == 3))) { // there should be a better way to distinguish between T6 and T3
                    surfaceElementNeighborsCommonEdge_.at(ii).insert(jj);
                    surfaceElementNeighborsCommonEdge_.at(jj).insert(ii);
                }
            }
        }
    }
    
    // PrintSurfaceElementNeighbors();
} // CBModel::DetermineSurfaceElementNeighbors

std::vector<std::unordered_set<TInt>> CBModel::GetNodeNeighbors(TInt material) {
    std::vector<std::unordered_set<TInt>> nodeNeighbors(nodes_.size());
    
    for (auto &e : solidElements_) {
        if ((material != -1) && (e->GetMaterialIndex() != material) )
            continue;
        
        for (TInt i = 0; i < e->GetNumberOfNodesIndices(); i++) {
            TInt nodeId = e->GetNodeIndex(i);
            for (TInt j = 0; j < e->GetNumberOfNodesIndices(); j++) {
                if (i == j)
                    continue;
                
                TInt neighId = e->GetNodeIndex(j);
                nodeNeighbors.at(nodeId).insert(neighId);
            }
        }
    }
    return nodeNeighbors;
}

void CBModel::PrintSolidElementNeighbors() {
    // solid elements only
    std::vector<TInt> mapping;
    
    for (auto &e : elements_) {
        if (dynamic_cast<CBElementSolid *>(e) != 0)
            mapping.push_back(e->GetIndex());
    }
    TInt numSolidElements = mapping.size();
    
    for (int i = 0; i < numSolidElements; i++) {
        std::cout << "Neighbors of solid element " << i << " with common face: ";
        std::unordered_set<TInt> neighbors = solidElementNeighborsCommonFace_.at(elements_.at(mapping[i])->GetIndex());
        for (auto n : neighbors)
            std::cout << " " << n;
        std::cout << std::endl;
    }
    for (int i = 0; i < numSolidElements; i++) {
        std::cout << "Neighbors of solid element " << i << " with common node: ";
        std::unordered_set<TInt> neighbors = solidElementNeighborsCommonNode_.at(elements_.at(mapping[i])->GetIndex());
        for (auto n : neighbors)
            std::cout << " " << n;
        std::cout << std::endl;
    }
}

void CBModel::PrintSurfaceElementNeighbors() {
    // surface elements only
    std::vector<TInt> mapping;
    
    for (auto &e : elements_) {
        if (dynamic_cast<CBElementSurface *>(e) != 0)
            mapping.push_back(e->GetIndex());
    }
    TInt numSurfaceElements = mapping.size();
    
    for (int i = 0; i < numSurfaceElements; i++) {
        std::cout << "Neighbors of surface element " << elements_.at(mapping[i])->GetIndex() << " with common edge: ";
        std::unordered_set<TInt> neighbors = surfaceElementNeighborsCommonEdge_.at(elements_.at(mapping[i])->GetIndex());
        for (auto n : neighbors)
            std::cout << " " << n;
        std::cout << std::endl;
    }
    for (int i = 0; i < numSurfaceElements; i++) {
        std::cout << "Neighbors of surface element " << elements_.at(mapping[i])->GetIndex() << " with common node: ";
        std::unordered_set<TInt> neighbors = surfaceElementNeighborsCommonNode_.at(elements_.at(mapping[i])->GetIndex());
        for (auto n : neighbors)
            std::cout << " " << n;
        std::cout << std::endl;
    }
}

void CBModel::PrintMaterials() {
    for (TInt i = 0; i < elements_.size(); i++) {
        std::cout << "Material of element : " << i << ":\t\t" << (elements_[i]->GetMaterialIndex()) << "\n";
    }
}

void CBModel::SaveNodes(std::string Filename) {
    std::ofstream file;
    
    file.open(Filename.c_str(), std::ios::out | std::ios::binary);
    TFloat *data = new TFloat[3 * nodes_.size()];
    TFloat *p    = data;
    for (std::vector<Vector3<TFloat>>::iterator it = nodes_.begin(); it != nodes_.end(); it++) {
        TFloat *c = it->GetArray();
        memcpy(p, c, 3 * sizeof(TFloat));
        p += 3;
    }
    file.write((char *)data, 3 * sizeof(TFloat) * nodes_.size());
    file.close();
}

void CBModel::SetGlobalData(std::string str, TFloat val) {
    auto it = globalData_.find(str);
    
    if (it == globalData_.end()) {
        globalData_.insert({str, val});
    } else {
        globalData_[str] = val;
    }
}

void CBModel::ApplyDomainDecomposition(std::vector<std::pair<TInt, TInt>> &mapping) {
    TInt numDomains = 0;
    
    if (mapping.size() != 0) {
        std::vector<Vector3<TFloat>> reorderdNodes;
        std::vector<TInt>             reorderdNodesComponentsBoundaryConditions;
        
        reorderdNodes = nodes_;
        reorderdNodesComponentsBoundaryConditions = nodesComponentsBoundaryConditions_;
        
        std::vector<TInt> tmpMapping = backwardMapping_;
        
        for (TInt i = 0; i < nodes_.size(); i++) {
            reorderdNodes.at(mapping.at(i).first) = nodes_.at(i);
            reorderdNodesComponentsBoundaryConditions.at(mapping.at(i).first) = nodesComponentsBoundaryConditions_.at(i);
            
            forwardMapping_.at(tmpMapping.at(mapping.at(i).first)) = i; // using tmpMapping to map from current indices to original to update forwardMapping
            backwardMapping_.at(i) = tmpMapping.at(mapping.at(i).first); // and vice versa
            
            nodesScalarData.SetData("DomainDecomposition", mapping.at(i).first, mapping.at(i).second);
            if (mapping.at(i).second > numDomains)
                numDomains = mapping.at(i).second;
        }
        
        nodes_ = reorderdNodes;
        nodesComponentsBoundaryConditions_ = reorderdNodesComponentsBoundaryConditions;
        
        for (std::vector<CBElement *>::iterator it = elements_.begin(); it != elements_.end(); it++) {
            for (unsigned int j = 0; j < (*it)->GetNumberOfNodesIndices(); j++) {
                TInt n = (*it)->GetNodeIndex(j);
                (*it)->SetNodeIndex(j, mapping.at(n).first);
            }
        }
    }
    
    numDomains++; // Number of Domains is the highest process id + 1 since first process has id zero
    
    nodesRanges_.resize(numDomains+1);
    
    nodesRanges_[0] = 0;
    
    for (unsigned int i = 1; i < numDomains; i++) {
        TInt m = 0;
        for (int j = 0; j < mapping.size(); j++) {
            if (mapping.at(j).second == i)
                m++;
        }
        if (m == 0) {
            throw std::runtime_error(
                                     "CBSolver::DetermineNodesRanges(): Nodes cannot be distributed efficiently, reduce number of processes");
        } else {
            nodesRanges_[i] = nodesRanges_[i - 1] + m;
        }
    }
    nodesRanges_[numDomains] = nodes_.size();
} // CBModel::ApplyDomainDecomposition

void CBModel::ApplySpatialSortPCA() {
    int             mDim = 3, nDim = nodes_.size();
    int             lda  = mDim, ldu = mDim, ldvt = nDim, info, lwork;
    double          wkopt;
    double *work;
    double *sVec = new double[nDim];
    double *uVec = new double[lda * mDim];
    
    //    double* vVecTransposed = new double[ldvt*nDim];
    double *a = new double[lda * nDim];
    Vector3<TFloat> mean(0, 0, 0);
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        mean += nodes_.at(i);
    }
    mean /= nodes_.size();
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        Vector3<TFloat> node = nodes_.at(i);
        node -= mean;
        TFloat *n = node.GetArray();
        for (int j = 0; j < 3; j++) {
            a[3 * i + j] = n[j];
        }
    }
    
    lwork = -1;
    dgesvd_("S", "N", &mDim, &nDim, a, &lda, sVec, uVec, &ldu, 0, &ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work  = new double[lwork];
    dgesvd_("S", "N", &mDim, &nDim, a, &lda, sVec, uVec, &ldu, 0, &ldvt, work, &lwork, &info);
    if (info > 0) {
        throw std::runtime_error("CBModel::ApplySpatialSortPCA(): The algorithm computing SVD failed to converge");
    }
    
    double *score = new double[nodes_.size()];
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        double          s    = 0;
        Vector3<TFloat> node = nodes_.at(i);
        node -= mean;
        
        TFloat *n = node.GetArray();
        
        for (int j = 0; j < 3; j++) {
            s += n[j] * uVec[j];
        }
        score[i] = s;
    }
    
    std::vector<std::pair<TInt, double>> nodesIndexes;
    nodesIndexes.reserve(nodes_.size());
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        std::pair<TInt, double> p;
        p.first  = i;
        p.second = score[i];
        nodesIndexes.push_back(p);
    }
    
    std::sort(nodesIndexes.begin(), nodesIndexes.end(), [](std::pair<TInt, double> a, std::pair<TInt, double> b) {
        return a.second < b.second;
    });
    
    TInt *mapping = new TInt[nodes_.size()];
    
    std::vector<Vector3<TFloat>> sortedNodes;
    std::vector<TInt>             sortedNodesComponentsBoundaryConditions;
    
    sortedNodes.reserve(nodes_.size());
    sortedNodesComponentsBoundaryConditions.reserve(nodes_.size());
    
    std::vector<TInt> tmpMapping = backwardMapping_;
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        sortedNodes.push_back(nodes_.at(nodesIndexes.at(i).first));
        
        forwardMapping_.at(tmpMapping.at(nodesIndexes.at(i).first)) = i; // using tmpMapping to map from current indices to original to update forwardMapping
        backwardMapping_.at(i) = tmpMapping.at(nodesIndexes.at(i).first); // and vice versa
        
        sortedNodesComponentsBoundaryConditions.push_back(nodesComponentsBoundaryConditions_.at(nodesIndexes.at(i).first));
        mapping[nodesIndexes.at(i).first] = i;
    }
    
    for (TInt i = 0; i < nodes_.size(); i++) {
        nodes_.at(i)                   = sortedNodes.at(i);
        nodesComponentsBoundaryConditions_.at(i) = sortedNodesComponentsBoundaryConditions.at(i);
    }
    
    for (std::vector<CBElement *>::iterator it = elements_.begin(); it != elements_.end(); it++) {
        for (unsigned int j = 0; j < (*it)->GetNumberOfNodesIndices(); j++) {
            TInt n = (*it)->GetNodeIndex(j);
            (*it)->SetNodeIndex(j, mapping[n]);
        }
    }
    
    delete[] work;
    delete[] score;
    delete[] sVec;
    delete[] uVec;
    
    //    delete vVecTransposed;
    delete[] a;
} // CBModel::ApplySpatialSortPCA

