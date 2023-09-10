/*
 *  CBModel.h
 *  CardioMechanics
 *
 *  Created by Thomas Fritz on 23.02.10.
 *  Copyright 2010 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#ifndef CB_MODEL_H
#define CB_MODEL_H

#include <vector>
#include <unordered_set>
#include <mutex>


#include "Matrix3.h"

#include "CBModelExporter.h"
#include "CBElementFactory.h"

#include "CBElement.h"
#include "CBElementSolid.h"
#include "CBElementSurface.h"

#include "DCType.h"
using namespace math_pack;

class CBModel {
protected:
    using SizeFunction = std::function<size_t()>;
    
    template<class T> class DataMap {
    private:
        using MapType = std::map<std::string, std::vector<T>>;
        
        SizeFunction sizeFunction;
        MapType data;
        
    public:
        friend class CBModel;
        DataMap(SizeFunction _sizeFunction) : sizeFunction(_sizeFunction) {
        }
        
        void SetData(std::string key, TInt index, T value) {
            auto& list = data[key];
            if (list.size() == 0)
                list.resize(sizeFunction(), T(0));
            list[index] = value;
        }
        
        T GetData(std::string key, TInt index) const {
            auto it = data.find(key);
            if(it != data.end()) {
                return(it->second.at(index));
            } else
                throw std::runtime_error(" CBModel::GetData(std::string name,TInt index): Data " + key + " at index " + std::to_string(index) + " does not exist !!");
        }
        
        std::vector<T>& GetData(std::string key) {
            auto it = data.find(key);
            if(it != data.end()) {
                return(it->second);
            } else
                throw std::runtime_error(" CBModel::GetData(std::string name,TInt index): Data " + key + " does not exist !!");
        }
        
        const MapType& GetData() const {
            return data;
        }
        
        void SetDataMap(std::string key, std::vector<T> value) {
            data[key] = value;
        }
        
        void InitDataMap(std::string key, T value) {
            data[key] = std::vector<T>(sizeFunction(), value);
        }
        
        bool IsData(std::string key) const {
            return data.find(key) != data.end();
        }
        
        void PrepareDataMap(std::string key, T value) {
            if (!IsData(key)) {
                InitDataMap(key, value);
            }
        }
        
        void Clear() {
            data.clear();
        }
        
    protected:
        MapType& GetDataReadWrite() {
            return data;
        }
    };
    
public:
    CBModel() :
    nodesScalarData([this](){return nodes_.size();}),
    elementsScalarData([this](){return elements_.size();}),
    elementsVectorData([this](){return elements_.size();}),
    nodesVectorData([this](){return nodes_.size();}),
    index_(0),
    exporter_(0) {}
    
    ~CBModel(){}
    
    void InitNodes(TInt numNodes);
    void InitMapping();
    void CalcNonZeroes();
    void InitNodesComponentsBoundaryConditions();
    void InitElements(TInt numElements);
    void ExtractSolidElements();
    void ExtractSurfaceElements();
    
    void SetExporter(CBModelExporter* exporter){exporter_ = exporter;}
    CBModelExporter* GetExporter(){return(exporter_);}
    bool HasExporter() { return (exporter_ != 0); }
    
    void UpdateNode(TInt pos, const Vector3<TFloat>& node);
    
    void CheckElementsIndexes();
    
    std::vector<PetscInt> GetNodeNeighborsForNnz() {return nodeNeighborsNnz_;}
    
    void DetermineNeighbors();
    void DetermineNodeSolidElements();
    void DetermineNodeSurfaceElements();
    void DetermineSolidElementNeighbors();
    void DetermineSurfaceElementNeighbors();
    /// Returns for the node 'pos' all solid elements that contain this node
    std::unordered_set<TInt> GetNodeSolidElements(TInt pos){return(nodeSolidElements_.at(pos));}
    /// Returns for the node 'pos' all surface elements that contain this node
    std::unordered_set<TInt> GetNodeSurfaceElements(TInt pos){return(nodeSurfaceElements_.at(pos));}
    /// Returns for the solid element 'pos' all solid elements that have 1 node in common
    std::unordered_set<TInt> GetSolidElementNeighborsCommonNode(TInt pos){return(solidElementNeighborsCommonNode_.at(pos));}
    /// Returns for the solid element 'pos' all solid elements that have 1 face in common
    std::unordered_set<TInt> GetSolidElementNeighborsCommonFace(TInt pos){return(solidElementNeighborsCommonFace_.at(pos));}
    /// Returns for the solid element 'pos' all solid elements that have 1 node in common
    std::unordered_set<TInt> GetSurfaceElementNeighborsCommonNode(TInt pos){return(surfaceElementNeighborsCommonNode_.at(pos));}
    /// Returns for the solid element 'pos' all solid elements that have 1 face in common
    std::unordered_set<TInt> GetSurfaceElementNeighborsCommonEdge(TInt pos){return(surfaceElementNeighborsCommonEdge_.at(pos));}
    /// Returns for all nodes all neighboring nodes that have a solid element with material number 'material' in common (use material = -1 for all materials)
    std::vector<std::unordered_set<TInt> > GetNodeNeighbors(TInt material);
    
    void SetCurrentTime(TFloat currentTime){currentTime_ = currentTime;}
    TFloat GetCurrentTime(){return(currentTime_);}
    
    void SetIndex(TInt index){index_ = index;}
    TInt GetIndex(){return(index_);}
    
    void SetMaxNnz(TInt nnz){maxNnz_ = nnz;}
    TInt GetMaxNnz(){return(maxNnz_);}
    
    void SetGlobalData(std::string str, TFloat val);
    
    void SetNodes(std::vector<Vector3<TFloat>>& newNodes) {nodes_ = newNodes;}
    void SetNodes(std::vector<Vector3<TFloat>>&& newNodes) {nodes_ = newNodes;}
    const std::vector<Vector3<TFloat>>& GetNodes() {return nodes_;}
    
    const std::vector<Vector3<TFloat>>& GetRefNodes() {return refNodes_;}
    
    void SetBoundaryConditions(std::vector<TInt>&& newBC) {nodesComponentsBoundaryConditions_ = newBC;}
    const std::vector<TInt>& GetNodeBoundaryConditions(){return nodesComponentsBoundaryConditions_;}
    
    void SetElements(std::vector<CBElement*>&& newElements) { elements_ = newElements; }
    void SetElements(std::vector<CBElement*>& newElements) { elements_ = newElements; }
    const std::vector<CBElement*>& GetElements() {return elements_;}
    const std::vector<CBElementSolid*>& GetSolidElements() {return solidElements_;}
    const std::vector<CBElementSurface*>& GetSurfaceElements() {return surfaceElements_;}
    
    const std::map<std::string, TFloat>& GetGlobalData() { return globalData_; }
    
    CBElementFactory* GetElementFactory(){return(&elementFactory_);}
    
    void SaveNodes(std::string Filename);
    void LoadNodes();
    
    void LoadElements();
    void SaveElements();
    
    void ApplySpatialSortPCA();
    void ApplyDomainDecomposition(std::vector<std::pair<int,int>>& mapping);
    
    TInt GetForwardMapping(TInt i){return(forwardMapping_.at(i));}
    TInt GetBackwardMapping(TInt i){return(backwardMapping_.at(i));}
    
    std::vector<TInt> GetNodesRanges(){return(nodesRanges_);}
    
    //* Lock CBModel for exclusive use
    std::mutex& GetLock(){return mutex_;}
    
    // Print data to stdout - Use only for debugging
    void PrintNodes();
    
    void PrintElements();
    void PrintSolidElementNeighbors();
    void PrintSurfaceElementNeighbors();
    void PrintElement(TInt i){elements_.at(i)->PrintNodesIndices();}
    
    void PrintMaterials();
    void PrintMaterial(TInt i){std::cout << "Material of element : " << i << ":\t\t" << (elements_.at(i)->GetMaterialIndex()) << "\n";}
    
    void ResetRefNodes(){ refNodes_ = nodes_;}
    
    void SetVelocity(std::vector<Vector3<TFloat>> vel){   velocity_      = vel; IsVelocity_     = true;}
    void SetAcceleration(std::vector<Vector3<TFloat>> acc){acceleration_ = acc; IsAcceleration_ = true;}
    void SetActiveStress(std::vector<TFloat> as){ ActiveStress_ = as;  IsActiveStress_ = true;}
    bool IsVelocity(){return IsVelocity_;}
    bool IsAcceleration(){return IsAcceleration_;}
    bool IsActiveStress(){return IsActiveStress_;}
    std::vector<Vector3<TFloat>> GetVelocity(){if(!IsVelocity_){throw std::runtime_error("CBModel: Velocity not set"); } return velocity_;}
    std::vector<Vector3<TFloat>> GetAcceleration(){if(!IsAcceleration_){throw std::runtime_error("CBModel: Acceleration not set"); } return acceleration_;}
    std::vector<TFloat> GetActiveStress(){if(!IsActiveStress_){throw std::runtime_error("CBModel: ActiveStress not set"); } return ActiveStress_;}
    
    DataMap<TFloat> nodesScalarData;
    DataMap<TFloat> elementsScalarData;
    DataMap<Vector3<TFloat>> elementsVectorData;
    DataMap<Vector3<TFloat>> nodesVectorData;
    
protected:
    TInt       index_;
    std::mutex mutex_;
    
    CBElementFactory elementFactory_;
    std::vector<Vector3<TFloat>> nodes_;
    std::vector<Vector3<TFloat>> refNodes_;
    std::vector<TInt> nodesComponentsBoundaryConditions_;
    std::vector<CBElement*> elements_;
    std::vector<CBElementSolid*> solidElements_;
    std::vector<CBElementSurface*> surfaceElements_; // [ep901] Implement surface list!
    
    std::vector<std::unordered_set<TInt> > nodeSolidElements_;
    std::vector<std::unordered_set<TInt> > nodeSurfaceElements_;
    std::map<TInt, std::unordered_set<TInt> > solidElementNeighborsCommonNode_;
    std::map<TInt, std::unordered_set<TInt> > solidElementNeighborsCommonFace_;
    std::map<TInt, std::unordered_set<TInt> > surfaceElementNeighborsCommonNode_;
    std::map<TInt, std::unordered_set<TInt> > surfaceElementNeighborsCommonEdge_;
    
    std::map<std::string, TFloat> globalData_;
    
    TFloat currentTime_;
    CBModelExporter* exporter_;
    
    std::vector<TInt> nodesRanges_;
    std::vector<TInt> forwardMapping_;
    std::vector<TInt> backwardMapping_;
    std::vector<PetscInt> nodeNeighborsNnz_;
    TInt maxNnz_ = -1;
    
    std::vector<Vector3<TFloat>> velocity_;
    std::vector<Vector3<TFloat>> acceleration_;
    std::vector<TFloat> ActiveStress_;
    
    bool IsVelocity_ = false;
    bool IsAcceleration_ = false;
    bool IsActiveStress_ = false;
    
private:
    CBModel(const CBModel&);
    void operator=(const CBModel&);
    
};
#endif
