/*
 * File: CBModelExporter.h
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


#ifndef CB_MODEL_EXPORTER
#define CB_MODEL_EXPORTER

#include <string>
#include <thread>
#include "ParameterMap.h"
#include "DCType.h"
#include "filesystem.h"
class CBModel;


class CBModelExporter {
public:
    CBModelExporter(CBModel *model, ParameterMap *parameters);
    virtual ~CBModelExporter();
    void Init();
    void DeInit();
    
    void ResetCounter() {counter_ = 0;}
    
    void SetFailedSimCounter() {counter_ = 99999999;}
    
    // Exporter Options:
    TFloat GetExportTimeStep() {return exportTimeStep_;}
    
    void SetExportTimeStep(TFloat exportTimeStep) { exportTimeStep_ = exportTimeStep; }
    
    std::string GetExportDir() { return exportGlobalDataFileDir_; }
    
    std::string GetExportDirPrefix() { return exportFileDir_ + "/" +exportFilenamePrefix_;}
    
    bool GetExportOption(std::string str, bool defaultValue = false);
    void WriteToFile();
    
protected:
    virtual bool ExportModel() = 0;
    void ExportGlobalData();
    CBModel *model_;
    ParameterMap *parameters_;
    TFloat        exportTimeStep_;
    TInt          counter_;
    std::string   exportFileDir_;
    std::string   exportGlobalDataFileDir_;
    std::string   exportFilenamePrefix_;
    std::thread *writeToFileThread_;
    
private:
    CBModelExporter() {}
    
    void WriteToFileThreadFunction();
    bool isInitialized_ = false;
};
#endif // ifndef CB_MODEL_EXPORTER
