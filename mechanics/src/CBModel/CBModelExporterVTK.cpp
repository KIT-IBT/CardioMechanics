/*
 * File: CBModelExporterVTK.cpp
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



#include "filesystem.h"
#include "DCCtrl.h"

#include "CBModelExporterVTK.h"
#include "CBElementSolidT4.h"
#include "CBElementSolidT10.h"
#include "CBElementSurfaceT3.h"
#include "CBElementSurfaceT6.h"
#include "CBModel.h"


using namespace math_pack;

void CBModelExporterVTK::InitPvdFile() {
    std::string timeStepsDir = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_ + "_vtu";
    
    if (!frizzle::filesystem::CreateDirectory(exportFileDir_)) {
        throw(std::string("void CBModelExporterVTK::InitPvdFile(): Path: " +exportFileDir_ +
                          " exists but is not a directory"));
    }
    if (!frizzle::filesystem::CreateDirectory(timeStepsDir)) {
        throw std::runtime_error(
                                 "void CBModelExporterVTK::InitPvdFile(): Path: " +timeStepsDir + " exists but is not a directory");
    }
    
    std::string   pvdFilename = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_ + ".pvd";
    std::ofstream pvdFile(pvdFilename.c_str());
    
    pvdFile <<
    "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    pvdFile << "\t<Collection>\n";
    pvdFile << "\t</Collection>\n";
    pvdFile << "</VTKFile>\n";
    
    pvdFile.close();
}

bool CBModelExporterVTK::ExportModel() {
    std::string   pvdFilename = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_ + ".pvd";
    std::ifstream pvdFileIn(pvdFilename.c_str());
    
    std::string str;
    std::vector<std::string> fileContent;
    
    while (getline(pvdFileIn, str)) {
        fileContent.push_back(str);
    }
    pvdFileIn.close();
    
    std::string timeStepFile = Base::exportFilenamePrefix_ + "." + std::to_string(Base::counter_) + ".vtu";
    std::string timeStepsDir = Base::exportFilenamePrefix_ + "_vtu";
    
    std::ofstream pvdFileOut(pvdFilename.c_str());
    for (std::vector<std::string>::iterator it = fileContent.begin(); it != fileContent.end(); it++) {
        if (it->find(std::string("</Collection>")) != std::string::npos) {
            pvdFileOut << "\t\t<DataSet timestep=\"" << Base::model_->GetCurrentTime() <<
            "\" part=\"0\" group=\"\" file=\"" << timeStepsDir << "/" << timeStepFile << "\"/>\n";
        }
        pvdFileOut << (*it) << std::endl;
    }
    pvdFileOut.close();
    
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> meshNodes = vtkSmartPointer<vtkPoints>::New();
    meshNodes->Allocate(Base::model_->GetNodes().size());
    
    vtkSmartPointer<vtkIntArray> originalPointIds = vtkSmartPointer<vtkIntArray>::New();
    originalPointIds->SetNumberOfComponents(1);
    originalPointIds->Allocate(Base::model_->GetNodes().size());
    originalPointIds->SetName("PointID");
    
    vtkSmartPointer<vtkIntArray> originalCellIds = vtkSmartPointer<vtkIntArray>::New();
    originalCellIds->SetNumberOfComponents(1);
    originalCellIds->Allocate(Base::model_->GetElements().size());
    originalCellIds->SetName("CellID");
    
    vtkSmartPointer<vtkIntArray> surfaceElementIDs = vtkSmartPointer<vtkIntArray>::New();
    surfaceElementIDs->SetNumberOfComponents(1);
    surfaceElementIDs->Allocate(Base::model_->GetElements().size());
    surfaceElementIDs->SetName("SurfaceElementID");
    
    vtkSmartPointer<vtkIntArray> surfaceIDs = vtkSmartPointer<vtkIntArray>::New();
    surfaceIDs->SetNumberOfComponents(1);
    surfaceIDs->Allocate(Base::model_->GetElements().size());
    surfaceIDs->SetName("SurfaceID");
    
    vtkSmartPointer<vtkDoubleArray> surfaceNormals = vtkSmartPointer<vtkDoubleArray>::New();
    surfaceNormals->SetNumberOfComponents(3);
    surfaceNormals->Allocate(Base::model_->GetElements().size());
    surfaceNormals->SetName("SurfaceNormal");
    
    for (TInt i = 0; i < Base::model_->GetNodes().size(); i++) {
        Vector3<TFloat> p = Base::model_->GetNodes().at(model_->GetForwardMapping(i));
        meshNodes->InsertNextPoint(p(0), p(1), p(2));
        originalPointIds->InsertNextValue(i);
    }
    mesh->GetPointData()->AddArray(originalPointIds);
    
    mesh->SetPoints(meshNodes);
    
    for (TInt i = 0; i < Base::model_->GetElements().size(); i++) {
        CBElement *element = Base::model_->GetElements().at(i);
        if (dynamic_cast<CBElementSolid *>(element) != 0) {
            if (dynamic_cast<CBElementSolidT4 *>(element) != 0) {
                vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
                for (int j = 0; j < 4; j++)
                    tetra->GetPointIds()->SetId(j, Base::model_->GetBackwardMapping(element->GetNodeIndex(j)));
                mesh->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
            } else if (dynamic_cast<CBElementSolidT10 *>(element) != 0) {
                vtkSmartPointer<vtkQuadraticTetra> tetra = vtkSmartPointer<vtkQuadraticTetra>::New();
                for (int j = 0; j < 10; j++)
                    tetra->GetPointIds()->SetId(j, Base::model_->GetBackwardMapping(element->GetNodeIndex(j)));
                mesh->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
            }
            
            if (Base::GetExportOption("SurfaceElementIds", false))
                surfaceElementIDs->InsertNextValue(-2);
            if (Base::GetExportOption("SurfaceIds", false))
                surfaceIDs->InsertNextValue(-2);
            if (Base::GetExportOption("SurfaceNormals", false)) {
                double normal[3] = {0, 0, 0};
                surfaceNormals->InsertNextTuple(normal);
            }
        } else if (dynamic_cast<CBElementSurface *>(element) != 0) {
            CBElementSurface *surface = dynamic_cast<CBElementSurface *>(element);
            
            if (dynamic_cast<CBElementSurfaceT3 *>(element) != 0) {
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, Base::model_->GetBackwardMapping(element->GetNodeIndex(0)));
                triangle->GetPointIds()->SetId(1, Base::model_->GetBackwardMapping(element->GetNodeIndex(1)));
                triangle->GetPointIds()->SetId(2, Base::model_->GetBackwardMapping(element->GetNodeIndex(2)));
                mesh->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
            } else if (dynamic_cast<CBElementSurfaceT6 *>(element) != 0) {
                vtkSmartPointer<vtkQuadraticTriangle> triangle = vtkSmartPointer<vtkQuadraticTriangle>::New();
                triangle->GetPointIds()->SetId(0, Base::model_->GetBackwardMapping(element->GetNodeIndex(0)));
                triangle->GetPointIds()->SetId(1, Base::model_->GetBackwardMapping(element->GetNodeIndex(1)));
                triangle->GetPointIds()->SetId(2, Base::model_->GetBackwardMapping(element->GetNodeIndex(2)));
                triangle->GetPointIds()->SetId(3, Base::model_->GetBackwardMapping(element->GetNodeIndex(3)));
                triangle->GetPointIds()->SetId(4, Base::model_->GetBackwardMapping(element->GetNodeIndex(4)));
                triangle->GetPointIds()->SetId(5, Base::model_->GetBackwardMapping(element->GetNodeIndex(5)));
                mesh->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
            } else {
                throw std::runtime_error(
                                         "CBModelExporterVTK::WriteToFile(): Only surface elements with 3 or 6 vertices are supported");
            }
            
            if (Base::GetExportOption("SurfaceElementIds", false))
                surfaceElementIDs->InsertNextValue(surface->GetSurfaceElementIndex());
            if (Base::GetExportOption("SurfaceIds", false))
                surfaceIDs->InsertNextValue(surface->GetSurfaceIndex());
            if (Base::GetExportOption("SurfaceNormals", false)) {
                Vector3<TFloat> p1 = Base::model_->GetNodes().at(element->GetNodeIndex(0));
                Vector3<TFloat> p2 = Base::model_->GetNodes().at(element->GetNodeIndex(1));
                Vector3<TFloat> p3 = Base::model_->GetNodes().at(element->GetNodeIndex(2));
                Vector3<TFloat> p4 = CrossProduct(p1-p2, p1-p3);
                surfaceNormals->InsertNextTuple(p4.GetArray());
            }
        } else {
            model_->GetLock().unlock();
            
            // Base::model_->Unlock(); // [ep901] Why is here an unlock. There is no corresponding lock!
            throw std::runtime_error(
                                     "CBModelExporterVTK::WriteToFile(): Element Type: " + element->GetType() + " not implemented !");
            return false;
        }
        
        originalCellIds->InsertNextValue(i);
    }
    mesh->GetCellData()->AddArray(originalCellIds);
    
    if (Base::GetExportOption("SurfaceElementIds", false))
        mesh->GetCellData()->AddArray(surfaceElementIDs);
    
    if (Base::GetExportOption("SurfaceIds", false))
        mesh->GetCellData()->AddArray(surfaceIDs);
    
    if (Base::GetExportOption("Material", true)) {
        vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
        mat->SetNumberOfComponents(1);
        mat->Allocate(Base::model_->GetElements().size());
        mat->SetName("Material");
        
        for (TInt i = 0; i < Base::model_->GetElements().size(); i++) {
            CBElement *element = Base::model_->GetElements().at(i);
            mat->InsertNextValue(element->GetMaterialIndex());
        }
        mesh->GetCellData()->AddArray(mat);
    }
    
    if (Base::GetExportOption("Fixation", true)) {
        vtkSmartPointer<vtkIntArray> fix = vtkSmartPointer<vtkIntArray>::New();
        fix->SetNumberOfComponents(1);
        fix->Allocate(Base::model_->GetNodes().size());
        fix->SetName("Fixation");
        
        for (TInt i = 0; i < Base::model_->GetNodes().size(); i++) {
            fix->InsertNextValue(Base::model_->GetNodeBoundaryConditions().at(Base::model_->GetForwardMapping(i)));
        }
        mesh->GetPointData()->AddArray(fix);
    }
    
    // Export of Fiber, Sheet, and SheetNormal has been moved to CBSolver::Export.
    // This section is now deprecated.
    //    if(Base::GetExportOption("Fiber", false))
    //    {
    //        vtkSmartPointer<vtkDoubleArray> fiber = vtkSmartPointer<vtkDoubleArray>::New();
    //        fiber->SetNumberOfComponents(3);
    //        fiber->Allocate(3 * Base::model_->GetElements().size());
    //        fiber->SetName("FiberUnrotated");
    //        for(TInt i = 0; i < Base::model_->GetElements().size(); i++)
    //        {
    //            CBElement* element = Base::model_->GetElements().at(i);
    //            Vector3<TFloat>          b       = element->GetBasis()->GetCol(0);
    //            double                   tuple[3];
    //            tuple[0] = b(0);
    //            tuple[1] = b(1);
    //            tuple[2] = b(2);
    //            fiber->InsertNextTuple(tuple);
    //        }
    //        mesh->GetCellData()->AddArray(fiber);
    //    }
    //
    //    if(Base::GetExportOption("Sheet", false))
    //    {
    //        vtkSmartPointer<vtkDoubleArray> sheet = vtkSmartPointer<vtkDoubleArray>::New();
    //        sheet->SetNumberOfComponents(3);
    //        sheet->Allocate(3 * Base::model_->GetElements().size());
    //        sheet->SetName("SheetUnrotated");
    //        for(TInt i = 0; i < Base::model_->GetElements().size(); i++)
    //        {
    //            CBElement* element = Base::model_->GetElements().at(i);
    //            Vector3<TFloat>          b       = element->GetBasis()->GetCol(1);
    //            double                   tuple[3];
    //            tuple[0] = b(0);
    //            tuple[1] = b(1);
    //            tuple[2] = b(2);
    //            sheet->InsertNextTuple(tuple);
    //        }
    //        mesh->GetCellData()->AddArray(sheet);
    //    }
    //
    //    if(Base::GetExportOption("SheetNormal", false))
    //    {
    //        vtkSmartPointer<vtkDoubleArray> sheetnormal = vtkSmartPointer<vtkDoubleArray>::New();
    //        sheetnormal->SetNumberOfComponents(3);
    //        sheetnormal->Allocate(3 * Base::model_->GetElements().size());
    //        sheetnormal->SetName("SheetNormalUnrotated");
    //        for(TInt i = 0; i < Base::model_->GetElements().size(); i++)
    //        {
    //            CBElement* element = Base::model_->GetElements().at(i);
    //            Vector3<TFloat>          b       = element->GetBasis()->GetCol(2);
    //            double                   tuple[3];
    //            tuple[0] = b(0);
    //            tuple[1] = b(1);
    //            tuple[2] = b(2);
    //            sheetnormal->InsertNextTuple(tuple);
    //        }
    //        mesh->GetCellData()->AddArray(sheetnormal);
    //    }
    
    if (Base::GetExportOption("FiberAllQuadraturePoints", false)) {
        // increase this when more quadraturepoints are needed
        int nQ = 5;
        vtkSmartPointer<vtkDoubleArray> Fiber[5];
        vtkSmartPointer<vtkDoubleArray> Sheet[5];
        vtkSmartPointer<vtkDoubleArray> Normal[5];
        for (int i = 0; i < nQ; i++) {
            Fiber[i] = vtkSmartPointer<vtkDoubleArray>::New();
            Fiber[i]->SetNumberOfComponents(3);
            Fiber[i]->Allocate(3*Base::model_->GetElements().size());
            Fiber[i]->SetName(("Fiber" + std::to_string(i)).c_str());
            
            Sheet[i] = vtkSmartPointer<vtkDoubleArray>::New();
            Sheet[i]->SetNumberOfComponents(3);
            Sheet[i]->Allocate(3*Base::model_->GetElements().size());
            Sheet[i]->SetName(("Sheet" + std::to_string(i)).c_str());
            
            Normal[i] = vtkSmartPointer<vtkDoubleArray>::New();
            Normal[i]->SetNumberOfComponents(3);
            Normal[i]->Allocate(3*Base::model_->GetElements().size());
            Normal[i]->SetName(("Normal" + std::to_string(i)).c_str());
        }
        for (TInt i = 0; i < Base::model_->GetElements().size(); i++) {
            CBElement *element = Base::model_->GetElements().at(i);
            for (int j = 0; j < nQ; j++) {
                Vector3<TFloat> f = Vector3<TFloat>(0, 0, 0);
                Vector3<TFloat> s = Vector3<TFloat>(0, 0, 0);
                Vector3<TFloat> n = Vector3<TFloat>(0, 0, 0);
                if ((dynamic_cast<CBElementSolid *>(element) != 0) && (j < element->GetNumberOfQuadraturePoints())) {
                    f = element->GetBasisAtQuadraturePoint(j)->GetCol(0);
                    s = element->GetBasisAtQuadraturePoint(j)->GetCol(1);
                    n = element->GetBasisAtQuadraturePoint(j)->GetCol(2);
                }
                
                double t1[3];
                double t2[3];
                double t3[3];
                
                t1[0] = f(0);
                t1[1] = f(1);
                t1[2] = f(2);
                t2[0] = s(0);
                t2[1] = s(1);
                t2[2] = s(2);
                t3[0] = n(0);
                t3[1] = n(1);
                t3[2] = n(2);
                
                Fiber[j]->InsertNextTuple(t1);
                Sheet[j]->InsertNextTuple(t2);
                Normal[j]->InsertNextTuple(t3);
            }
        }
        for (int i = 0; i < nQ; i++) {
            mesh->GetCellData()->AddArray(Fiber[i]);
            mesh->GetCellData()->AddArray(Sheet[i]);
            mesh->GetCellData()->AddArray(Normal[i]);
        }
    }
    
    for (auto &it : model_->nodesVectorData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->Allocate(3 * Base::model_->GetNodes().size());
        data->SetName(it.first.c_str());
        for (int i = 0; i < it.second.size(); i++) {
            Vector3<TFloat> b = (it.second.at(Base::model_->GetForwardMapping(i)));
            double          tuple[3];
            
            tuple[0] = b(0);
            tuple[1] = b(1);
            tuple[2] = b(2);
            data->InsertNextTuple(tuple);
        }
        mesh->GetPointData()->AddArray(data);
    }
    
    
    for (auto &it : model_->elementsVectorData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->Allocate(3 * Base::model_->GetElements().size());
        data->SetName(it.first.c_str());
        for (auto &it2 : it.second) {
            Vector3<TFloat> b = it2;
            double          tuple[3];
            tuple[0] = b(0);
            tuple[1] = b(1);
            tuple[2] = b(2);
            data->InsertNextTuple(tuple);
        }
        mesh->GetCellData()->AddArray(data);
    }
    
    
    for (auto &it : model_->nodesScalarData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->Allocate(Base::model_->GetNodes().size());
        data->SetName(it.first.c_str());
        for (int i = 0; i < it.second.size(); i++) {
            data->InsertNextValue(it.second.at(Base::model_->GetForwardMapping(i)));
        }
        mesh->GetPointData()->AddArray(data);
    }
    
    
    for (auto &it : model_->elementsScalarData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->Allocate(Base::model_->GetElements().size());
        data->SetName(it.first.c_str());
        for (auto &it2 : it.second) {
            data->InsertNextValue(double(it2));
        }
        mesh->GetCellData()->AddArray(data);
    }
    
    
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(std::string(Base::exportFileDir_ + "/" + timeStepsDir + "/" + timeStepFile).c_str());
#if VTK_MAJOR_VERSION > 5
    writer->SetInputData(mesh);
#else // if VTK_MAJOR_VERSION > 5
    writer->SetInput(mesh);
#endif // if VTK_MAJOR_VERSION > 5
    writer->Write();
    
    ExportTetgenNode();
    ExportTetgenBases();
    
    Base::counter_++;
    return true;
} // CBModelExporterVTK::ExportModel

void CBModelExporterVTK::InitTetgen() {
    std::string nodeDir = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_ + "_node/";
    
    if (!frizzle::filesystem::CreateDirectory(nodeDir))
        throw(std::string("void CBModelExporterVTK::InitTetgen(): Path: " + nodeDir + " exists but is not a directory"));
    
    std::string basesDir = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_ + "_bases/";
    if (!frizzle::filesystem::CreateDirectory(basesDir))
        throw(std::string("void CBModelExporterVTK::InitTetgen(): Path: " + basesDir + " exists but is not a directory"));
}

void CBModelExporterVTK::ExportTetgenNode() {
    std::string nodeFilename = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_;
    
    nodeFilename += "_node/" + Base::exportFilenamePrefix_ + "." + std::to_string(Base::counter_) + ".node";
    
    std::ofstream nodeFile;
    nodeFile.open(nodeFilename);
    if (!nodeFile.good())
        throw std::runtime_error("CBModelExporterTetgen::ExportModel: Couldn't create " + nodeFilename + ".");
    
    size_t numNodes = Base::model_->GetNodes().size();
    nodeFile << numNodes << " 3 1 0" << std::endl;
    
    for (TInt i = 0; i < numNodes; i++) {
        Vector3<TFloat> p = Base::model_->GetNodes().at(model_->GetForwardMapping(i));
        p *= 1e3; // Tetgen is in mm
        
        // Handling of boundary conditions should be extended to:
        // 001 (1): x fixed; 010 (2): y fixed; 100 (4): z fixed; 111 (7): x,y,z fixed
        TInt b = 0;
        if (Base::model_->GetNodeBoundaryConditions().at(Base::model_->GetForwardMapping(i)))
            b = 7;
        
        nodeFile << i+1 << " " << p(0) << " " << p(1) << " " << p(2) << " " << b << std::endl;
    }
    nodeFile.close();
} // CBModelExporterVTK::ExportTetgenNode

void CBModelExporterVTK::ExportTetgenBases() {
    if (model_->GetExporter()->GetExportOption("TetgenBases", false)) {
        // get bases from model and write to a .bases file
        // note: From the solver, the function ExportBasesAtQuadraturePointsToModel() needs to be called first!
        
        std::string basesFilename = Base::exportFileDir_ + "/" + Base::exportFilenamePrefix_;
        basesFilename += "_bases/" + Base::exportFilenamePrefix_ + "." + std::to_string(Base::counter_) + ".bases";
        
        DCCtrl::debug << "Export bases to: " << basesFilename << "\n";
        
        std::ofstream basesFile;
        basesFile.open(basesFilename);
        if (!basesFile.good())
            throw std::runtime_error("CBModelExporterTetgen::ExportModel: Couldn't create " + basesFilename + ".");
        
        TInt nQ = Base::model_->GetSolidElements().at(0)->GetNumberOfQuadraturePoints();
        for (auto &e : Base::model_->GetSolidElements()) {
            if (e->GetNumberOfQuadraturePoints() != nQ)
                throw std::runtime_error("ERROR: Tetgen export is currently only implemented for all T4 or all T10 elements.");
        }
        
        size_t numElements = Base::model_->GetSolidElements().size();
        basesFile << numElements << " " << nQ << std::endl; // HEADER
        
        for (auto &e : Base::model_->GetSolidElements()) {
            basesFile << e->GetIndex()+1;
            for (TInt i = 0; i < nQ; i++) {
                basesFile << " ";
                
                // note: bases can not be taken from the elements, since these are not rotated with the current deformation tensor
                // reading the current deformation tensor requires non-trivial communication between the processes
                std::string qString = std::to_string(i+1);
                Vector3<TFloat> f  = model_->elementsVectorData.GetData("TetgenBasesFiberAtQuadraturePoint" + qString,
                                                                        e->GetIndex());
                Vector3<TFloat> s  = model_->elementsVectorData.GetData("TetgenBasesSheetAtQuadraturePoint" + qString,
                                                                        e->GetIndex());
                Vector3<TFloat> sn = model_->elementsVectorData.GetData("TetgenBasesSheetnormalAtQuadraturePoint" + qString,
                                                                        e->GetIndex());
                
                basesFile << " " <<  f(0) << " " <<  f(1) << " " <<  f(2);
                basesFile << " " <<  s(0) << " " <<  s(1) << " " <<  s(2);
                basesFile << " " << sn(0) << " " << sn(1) << " " << sn(2);
            }
            basesFile << std::endl;
        }
        basesFile.close();
    }
} // CBModelExporterVTK::ExportTetgenBases

vtkSmartPointer<vtkUnstructuredGrid> CBModelExporterVTK::GetMesh() {
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> meshNodes = vtkSmartPointer<vtkPoints>::New();
    
    meshNodes->Allocate(Base::model_->GetNodes().size());
    
    vtkSmartPointer<vtkIntArray> originalPointIds = vtkSmartPointer<vtkIntArray>::New();
    originalPointIds->SetNumberOfComponents(1);
    originalPointIds->Allocate(Base::model_->GetNodes().size());
    originalPointIds->SetName("PointID");
    
    vtkSmartPointer<vtkIntArray> originalCellIds = vtkSmartPointer<vtkIntArray>::New();
    originalCellIds->SetNumberOfComponents(1);
    originalCellIds->Allocate(Base::model_->GetElements().size());
    originalCellIds->SetName("CellID");
    
    vtkSmartPointer<vtkIntArray> surfaceElementIDs = vtkSmartPointer<vtkIntArray>::New();
    surfaceElementIDs->SetNumberOfComponents(1);
    surfaceElementIDs->Allocate(Base::model_->GetElements().size());
    surfaceElementIDs->SetName("SurfaceElementID");
    
    
    vtkSmartPointer<vtkIntArray> surfaceIDs = vtkSmartPointer<vtkIntArray>::New();
    surfaceIDs->SetNumberOfComponents(1);
    surfaceIDs->Allocate(Base::model_->GetElements().size());
    surfaceIDs->SetName("SurfaceID");
    
    vtkSmartPointer<vtkDoubleArray> surfaceNormals = vtkSmartPointer<vtkDoubleArray>::New();
    surfaceNormals->SetNumberOfComponents(3);
    surfaceNormals->Allocate(Base::model_->GetElements().size());
    surfaceNormals->SetName("SurfaceNormal");
    
    for (TInt i = 0; i < Base::model_->GetNodes().size(); i++) {
        Vector3<TFloat> p = Base::model_->GetNodes().at(model_->GetForwardMapping(i));
        meshNodes->InsertNextPoint(p(0), p(1), p(2));
        originalPointIds->InsertNextValue(i);
    }
    mesh->GetPointData()->AddArray(originalPointIds);
    
    mesh->SetPoints(meshNodes);
    
    for (TInt i = 0; i < Base::model_->GetElements().size(); i++) {
        CBElement *element = Base::model_->GetElements().at(i);
        
        if (dynamic_cast<CBElementSolid *>(element) != 0) {
            if (dynamic_cast<CBElementSolidT4 *>(element) != 0) {
                vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
                for (int j = 0; j < 4; j++)
                    tetra->GetPointIds()->SetId(j, Base::model_->GetBackwardMapping(element->GetNodeIndex(j)));
                mesh->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
            } else if (dynamic_cast<CBElementSolidT10 *>(element) != 0) {
                vtkSmartPointer<vtkQuadraticTetra> tetra = vtkSmartPointer<vtkQuadraticTetra>::New();
                for (int j = 0; j < 10; j++)
                    tetra->GetPointIds()->SetId(j, Base::model_->GetBackwardMapping(element->GetNodeIndex(j)));
                mesh->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
            }
            
            if (Base::GetExportOption("SurfaceIds", false))
                surfaceIDs->InsertNextValue(0);
            if (Base::GetExportOption("SurfaceElementIds", false))
                surfaceElementIDs->InsertNextValue(0);
        } else if (dynamic_cast<CBElementSurface *>(element) != 0) {
            CBElementSurface *surface = dynamic_cast<CBElementSurface *>(element);
            
            if (dynamic_cast<CBElementSurfaceT3 *>(element) != 0) {
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, Base::model_->GetBackwardMapping(element->GetNodeIndex(0)));
                triangle->GetPointIds()->SetId(1, Base::model_->GetBackwardMapping(element->GetNodeIndex(1)));
                triangle->GetPointIds()->SetId(2, Base::model_->GetBackwardMapping(element->GetNodeIndex(2)));
                mesh->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
            } else if (dynamic_cast<CBElementSurfaceT6 *>(element) != 0) {
                vtkSmartPointer<vtkQuadraticTriangle> triangle = vtkSmartPointer<vtkQuadraticTriangle>::New();
                triangle->GetPointIds()->SetId(0, Base::model_->GetBackwardMapping(element->GetNodeIndex(0)));
                triangle->GetPointIds()->SetId(1, Base::model_->GetBackwardMapping(element->GetNodeIndex(1)));
                triangle->GetPointIds()->SetId(2, Base::model_->GetBackwardMapping(element->GetNodeIndex(2)));
                triangle->GetPointIds()->SetId(0, Base::model_->GetBackwardMapping(element->GetNodeIndex(3)));
                triangle->GetPointIds()->SetId(1, Base::model_->GetBackwardMapping(element->GetNodeIndex(4)));
                triangle->GetPointIds()->SetId(2, Base::model_->GetBackwardMapping(element->GetNodeIndex(5)));
                mesh->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
            } else {
                throw std::runtime_error(
                                         "CBModelExporterVTK::WriteToFile(): Only surface elements with 3 or 6 vertices are supported");
            }
            
            if (Base::GetExportOption("SurfaceElementIds", false))
                surfaceElementIDs->InsertNextValue(surface->GetSurfaceElementIndex());
            if (Base::GetExportOption("SurfaceIds", false))
                surfaceIDs->InsertNextValue(surface->GetSurfaceIndex());
            if (Base::GetExportOption("SurfaceNormals", false)) {
                Vector3<TFloat> p1 = Base::model_->GetNodes().at(element->GetNodeIndex(0));
                Vector3<TFloat> p2 = Base::model_->GetNodes().at(element->GetNodeIndex(1));
                Vector3<TFloat> p3 = Base::model_->GetNodes().at(element->GetNodeIndex(2));
                Vector3<TFloat> p4 = CrossProduct(p1-p2, p1-p3);
                surfaceNormals->InsertNextTuple(p4.GetArray());
            }
        } else {
            // [ep901] Why is here an unlock. There is no corresponding lock!
            // Base::model_->Unlock();
            model_->GetLock().unlock();
            throw std::runtime_error(
                                     "CBModelExporterVTK::WriteToFile(): Element Type: " + element->GetType() + " not implemented !");
            return vtkSmartPointer<vtkUnstructuredGrid>(NULL);
        }
        
        originalCellIds->InsertNextValue(i);
    }
    mesh->GetCellData()->AddArray(originalCellIds);
    
    // always export material to mesh
    {
        vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
        mat->SetNumberOfComponents(1);
        mat->Allocate(Base::model_->GetElements().size());
        mat->SetName("Material");
        
        for (TInt i = 0; i < Base::model_->GetElements().size(); i++) {
            CBElement *element = Base::model_->GetElements().at(i);
            mat->InsertNextValue(element->GetMaterialIndex());
        }
        mesh->GetCellData()->AddArray(mat);
    }
    
    // always export fixation to mesh
    {
        vtkSmartPointer<vtkIntArray> fix = vtkSmartPointer<vtkIntArray>::New();
        fix->SetNumberOfComponents(1);
        fix->Allocate(Base::model_->GetNodes().size());
        fix->SetName("Fixation");
        
        for (TInt i = 0; i < Base::model_->GetNodes().size(); i++) {
            fix->InsertNextValue(Base::model_->GetNodeBoundaryConditions().at(Base::model_->GetForwardMapping(i)));
        }
        mesh->GetPointData()->AddArray(fix);
    }
    
    for (auto &it : model_->nodesVectorData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->Allocate(3 * Base::model_->GetNodes().size());
        data->SetName(it.first.c_str());
        for (int i = 0; i < it.second.size(); i++) {
            Vector3<TFloat> b = (it.second.at(Base::model_->GetForwardMapping(i)));
            double          tuple[3];
            
            tuple[0] = b(0);
            tuple[1] = b(1);
            tuple[2] = b(2);
            data->InsertNextTuple(tuple);
        }
        mesh->GetPointData()->AddArray(data);
    }
    
    
    for (auto &it : model_->elementsVectorData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->Allocate(3 * Base::model_->GetElements().size());
        data->SetName(it.first.c_str());
        for (auto &it2 : it.second) {
            Vector3<TFloat> b = it2;
            double          tuple[3];
            tuple[0] = b(0);
            tuple[1] = b(1);
            tuple[2] = b(2);
            data->InsertNextTuple(tuple);
        }
        mesh->GetCellData()->AddArray(data);
    }
    
    
    for (auto &it : model_->nodesScalarData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->Allocate(Base::model_->GetNodes().size());
        data->SetName(it.first.c_str());
        for (int i = 0; i < it.second.size(); i++) {
            data->InsertNextValue(it.second.at(Base::model_->GetForwardMapping(i)));
        }
        mesh->GetPointData()->AddArray(data);
    }
    
    
    for (auto &it : model_->elementsScalarData.GetData()) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->Allocate(Base::model_->GetElements().size());
        data->SetName(it.first.c_str());
        for (auto &it2 : it.second) {
            data->InsertNextValue(double(it2));
        }
        mesh->GetCellData()->AddArray(data);
    }
    
    return mesh;
} // CBModelExporterVTK::GetMesh
