/*
 * File: ExtractSurfaceNodesFromVTU.cpp
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


// Extract the coordinates of the surface nodes of a deformed mesh (.vtu) for the
// requested surface indices, as listed in a Tetgen .sur file, and write them as a
// plain-text list. Used to build the per-timestep target surfaces consumed by the
// inverse active-stress estimator (see CBPointsCtrl).

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh.vtu> <nodes file> <surface file> '<surface indices>' <output>" << std::endl;
        std::cerr << "Example: mesh.vtu mesh.node mesh.sur '1 2' surface.dat" << std::endl;
        return EXIT_FAILURE;
    }

    // Read the deformed mesh points from the VTU.
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    vtkUnstructuredGrid *grid = reader->GetOutput();
    vtkSmartPointer<vtkPoints> nodes = grid->GetPoints();

    // Surface indices to extract.
    std::vector<int> surfIdxToExtract;
    std::string surfArg(argv[4]);
    std::stringstream stream(surfArg);
    for (int idx; stream >> idx; )
        surfIdxToExtract.push_back(idx);

    std::ifstream surfaceFile(argv[3]);
    std::string str;

    // Header: number of triangles, dimension, number of surface attributes.
    getline(surfaceFile, str);

    // Collect the (1-based) node indices of every triangle whose surface index matches.
    std::set<int> nodesIndices;
    while (getline(surfaceFile, str)) {
        std::stringstream ss(str);
        std::vector<int> currentLine;
        for (int v; ss >> v; )
            currentLine.push_back(v);
        if (currentLine.empty())
            continue;

        int currentSurfIdx = currentLine.back();
        if (std::find(surfIdxToExtract.begin(), surfIdxToExtract.end(), currentSurfIdx) != surfIdxToExtract.end()) {
            // Column 0 is the triangle id; the next three are the edge nodes.
            nodesIndices.insert(currentLine[1]);
            nodesIndices.insert(currentLine[2]);
            nodesIndices.insert(currentLine[3]);
        }
    }
    surfaceFile.close();

    std::ofstream outputFile(argv[5], std::ofstream::out);
    if (!outputFile.good())
        throw std::runtime_error("Cannot write to output file");

    outputFile << std::to_string(3 * nodesIndices.size()) << "\n";
    for (int i : nodesIndices) {
        double *p = nodes->GetPoint(i - 1);
        outputFile << std::to_string(p[0]) << " "
                   << std::to_string(p[1]) << " "
                   << std::to_string(p[2]) << "\n";
    }
    outputFile.close();

    return 0;
}
