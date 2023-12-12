/*
 * File: acltProjectFile.h
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

#ifndef ACLTPROJECTFILE
#define ACLTPROJECTFILE

class acltProjectEntry {
public:
    string name;
    string explanation;
};

class acltProjectFile {
public:
    acltProjectFile() {init();}
    
    inline void printHelp();
    inline int getIndex(string);
    
private:
    vector<acltProjectEntry> PE;
    inline void init();
};

void acltProjectFile::init() {
    acltProjectEntry lPE;
    
    lPE.name = "Resprefix";       lPE.explanation = "The prefix of the result data";
    PE.push_back(lPE);
    lPE.name        = "Condition";
    lPE.explanation =
    "Name of the condition file (entries: <vector position> <amplitude> <rate> <duration> <temporal offset> II|IE|IF|UI|UE|UF)";
    PE.push_back(lPE);
    lPE.name        = "Sensor";
    lPE.explanation = "Name of the sensor file (entries: <vector position> <filename> <beginsave> <dtsave> <type>)";
    PE.push_back(lPE);
    lPE.name = "CalcLength";      lPE.explanation = "Calculation length [s]";
    PE.push_back(lPE);
    lPE.name = "DTCell";          lPE.explanation = "Time step of cell model calculation [s]";
    PE.push_back(lPE);
    lPE.name = "DTExtra";         lPE.explanation = "Time step of extracellular calculation [s]";
    PE.push_back(lPE);
    lPE.name = "DTIntra";         lPE.explanation = "Time step of intracellular calculation [s]";
    PE.push_back(lPE);
    lPE.name = "BeginSave";       lPE.explanation = "Time at which saving of results sterts [s]";
    PE.push_back(lPE);
    lPE.name = "DTSave";          lPE.explanation = "Time step of saving [s]";
    PE.push_back(lPE);
    lPE.name        = "CellModelFile";
    lPE.explanation = "Name of the cell model file (entries: <tissue class> <elphpy model> <force model>)";
    PE.push_back(lPE);
    lPE.name = "Material";        lPE.explanation = "Name of vector describing tissue classes";
    PE.push_back(lPE);
    lPE.name = "MaterialFibro";   lPE.explanation = "Name of vector describing fibroblast classes";
    PE.push_back(lPE);
    lPE.name = "MatrixIntra";     lPE.explanation = "Name of intracellular conductivity matrix";
    PE.push_back(lPE);
    lPE.name = "MatrixExtra";     lPE.explanation = "Name of extracellular conductivity matrix";
    PE.push_back(lPE);
    lPE.name        = "MatrixExtraCombined";
    lPE.explanation = "Name of extracellular conductivity matrix that is already combined with Intra (and Fibro)";
    PE.push_back(lPE);
    lPE.name = "MatrixFibro";     lPE.explanation = "Name of fibroblast conductivity matrix";
    PE.push_back(lPE);
    lPE.name = "BetaMyoFib";      lPE.explanation = "Number of Myo-Fibro gap junctions [1/m^3]";
    PE.push_back(lPE);
    lPE.name = "RMyoFib";         lPE.explanation = "Resistor of single myo-fibro gap junction [Ohm]";
    PE.push_back(lPE);
    lPE.name = "VolumeMyo";       lPE.explanation = "Relative amount of myocyte volume";
    PE.push_back(lPE);
    lPE.name = "VolumeExtra";     lPE.explanation = "Relative amount of extracellular volume";
    PE.push_back(lPE);
    lPE.name = "VolumeFibro";     lPE.explanation = "Relative amount of fibroblast volume";
    PE.push_back(lPE);
    lPE.name = "Domain";          lPE.explanation = "Type, either \"Mono\", \"Bi\", or \"Tri\"";
    PE.push_back(lPE);
    lPE.name = "Verbose";         lPE.explanation = "Print verbose information into terminal";
    PE.push_back(lPE);
    lPE.name = "Protocol";        lPE.explanation = "Name of verbose information file";
    PE.push_back(lPE);
    lPE.name        = "HeteroFileIntra";
    lPE.explanation = "Name of the heterogeneous cell model file (entries: <vector file name> <parameter name>)";
    PE.push_back(lPE);
    lPE.name        = "HeteroFileFibro";
    lPE.explanation = "Name of the heterogeneous fibro model file (entries: <vector file name> <parameter name>)";
    PE.push_back(lPE);
    lPE.name = "LoadBackup";      lPE.explanation = "Loads backup from standard file of <file> if declared";
    PE.push_back(lPE);
    lPE.name = "SaveBackup";      lPE.explanation = "Saves backup to standard file of <file> if declared";
    PE.push_back(lPE);
    lPE.name = "DTBackup";        lPE.explanation = "Time step of backup";
    PE.push_back(lPE);
    lPE.name        = "Results";
    lPE.explanation =
    "List of variables to be saved (e.g., Vm, Ve, AT, m, Na_i, I_Na, and any variable of the cell model)";
    PE.push_back(lPE);
    lPE.name        = "Implicit";
    lPE.explanation = "Number of Jacobi iterations for implicit calculation (!Only for Monodomain)";  PE.push_back(lPE);
    lPE.name        = "Force";           lPE.explanation = "Also calculate and output force";
    PE.push_back(lPE);
    
    // lPE.name="ResultsFibro";    lPE.explanation="List of variables to be saved (i.e. Vf, If, and any variable of the
    // fibro model)";PE.push_back(lPE);
    lPE.name = "MassMatrixIntra"; lPE.explanation = "Name of intracellular mass matrix";
    PE.push_back(lPE);
    lPE.name        = "ThetaIntra";
    lPE.explanation =
    "Theta for intra PDE solver scheme. [0.0..1.0] 0=explicit, 1=implicit (default: 0 (FD) or 0.5 (FE))";
    PE.push_back(lPE);
    lPE.name        = "MembraneCapacitance";
    lPE.explanation = "Membrane capacitance per unit area, number or vector (F/m^2, default: cell model)";
    PE.push_back(lPE);
    lPE.name        = "SurfaceToVolume";
    lPE.explanation = "Myocyte surface to volume ratio, number or vector (1/m, default: cell model)"; PE.push_back(lPE);
    lPE.name        = "Gauss";           lPE.explanation = "Matrices for gauss interpolation and integration";
    PE.push_back(lPE);
    lPE.name = "ActivationThreshold"; lPE.explanation = "TMV threshold for activation time (default: 0V).";
    PE.push_back(lPE);
    lPE.name        = "CurrentScheme";
    lPE.explanation = "Godunow, Strang, SVI, ICI, None. SVI and None require Gauss matrices."; PE.push_back(lPE);
    lPE.name        = "IntraIndexSet";
    lPE.explanation = "PETSc IndexSet of the intracellular domain (if smaller than the xtracellular domain)";
    PE.push_back(lPE);
    lPE.name = "Compress";        lPE.explanation = "Compress output vectors using deflate/gzip."; PE.push_back(lPE);
    lPE.name = "NoExport";        lPE.explanation = "Deactivate Export"; PE.push_back(lPE);
}  // acltProjectFile::init

void acltProjectFile::printHelp() {
    vector<acltProjectEntry>::iterator iPE, bPE = PE.begin(), ePE = PE.end();
    for (iPE = bPE; iPE != ePE; iPE++) {
        cerr << "\t"<< left << setw(15) << (*iPE).name << right << (*iPE).explanation << endl;
    }
}

int acltProjectFile::getIndex(string entry) {
    vector<acltProjectEntry>::iterator iPE, bPE = PE.begin(), ePE = PE.end();
    int i = 0;
    for (iPE = bPE; iPE != ePE; iPE++, i++) {
        if (!entry.compare((*iPE).name))
            return i;
    }
    return -1;
}

#endif  // ifndef ACLTPROJECTFILE
