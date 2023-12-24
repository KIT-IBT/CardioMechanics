/*
 * File: VTPExtractSurfaceFromMesh.cpp
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


#include<kaLattice.h>
#include<kaPoint.h>
#include<kaMatrixN.h>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<set>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkTetra.h>
#include <vtkQuadraticTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLReader.h>
#include <vtkFillHolesFilter.h>
#include<vtkPolyDataNormals.h>
#include<vtkFeatureEdges.h>
#include<vtkCleanPolyData.h>
#include<vtkPolyDataConnectivityFilter.h>
#include<vtkAppendPolyData.h>
#include<vtkStripper.h>
#include<vtkTriangleFilter.h>
#include<vtkPointLocator.h>
#include<vtkSelectEnclosedPoints.h>
#include<vtkIntArray.h>
#include<vtkMath.h>
#include<vtkLine.h>
#include<vtkMassProperties.h>


long int FindCorrespondingVoxel(kaPoint<double> a, kaLattice<unsigned char>* b)
{
	
	double c;
	kaPoint<int> d;
	
	
	kaMatrixN<double,3> inv;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			inv.a(i,j) = b->m.a(i,j);
	
	inv.Inverse();
	
	a.x -= b->m.a(0,3);
	a.y -= b->m.a(1,3);
	a.z -= b->m.a(2,3);
	
	c = (a.x * inv.a(0,0) + a.y * inv.a(0,1) + a.z * inv.a(0,2));
	if((c-int(c)) >= 0.5 && (c-int(c))< 1)c += 1.0;
	d.x = int(c);
	
	c = (a.x * inv.a(1,0) + a.y * inv.a(1,1) + a.z * inv.a(1,2));
	if((c-int(c)) >= 0.5 && (c-int(c))< 1)c += 1.0;
	d.y = int(c);
	
	c = (a.x * inv.a(2,0) + a.y * inv.a(2,1) + a.z * inv.a(2,2));
	if((c-int(c)) >= 0.5 && (c-int(c))< 1)c += 1.0;
	d.z = int(c);
	
	long int e = d.z*b->xyLattice + d.y*b->xLattice+d.x;
	return e;
}

void LoadNodes(std::ifstream& nodesFile,std::vector<kaPoint<double> >& nodes,double unit=1)
{
	long int numNodes=0;
	long int currentNode=0;
	int dim=0;
	int attr=0;
	int boundaryMarker=0;
	
	nodesFile >> numNodes;
	nodesFile >> dim;		
	nodesFile >> attr;
	nodesFile >> boundaryMarker;
	
	if(attr > 1 || boundaryMarker != 0)
		throw(std::string("More than one attribute or Boundary Marker not yet implemented !"));
	if(dim != 3)
		throw std::runtime_error("Only for 3D model !");
	
	for(long int i=0; i < numNodes; i++)
	{
		if(nodesFile.eof())
			throw(std::string("Nodes file is corrupt. Check number of nodes !"));
		double x,y,z;
		double dummy;
		nodesFile >> currentNode;
		nodesFile >> x;
		nodesFile >> y;
		nodesFile >> z;
		if(attr == 1)
			nodesFile >> dummy;
		nodes.push_back(kaPoint<double>(x*unit,y*unit,z*unit));
	}
}

void LoadElements(std::ifstream& elementsFile,std::vector<std::vector<long int> >& elements)
{
	int numElements=0;
	int nodesPerElement=0;
	long int currentElement=0;
	int attr=0;
	
	elementsFile >> numElements;
	elementsFile >> nodesPerElement;		
	elementsFile >> attr;
	
	if(nodesPerElement != 4 && nodesPerElement != 10)
		throw(std::string("Unkown Volume Element"));
    
	
    if(attr < 1)
        throw(std::string("Material information is needed"));
	
	
	for(long int i=0; i < numElements; i++)
	{
		
		// Load Elements
		if(elementsFile.eof())
			throw(std::string("Elements file is corrupt. Check number of nodes !"));
		
		
		elementsFile >> currentElement;
		std::vector<long int>ele;
		for(int j=0;j<nodesPerElement;j++)
		{
			long int n;
			elementsFile >> n;
			ele.push_back(n-1);
		}
        int m = 0;

        elementsFile >> m; // Materials
        ele.push_back(m);
		elements.push_back(ele);
        
		for(int j=0;j<attr-1;j++) // Ignore remaining attributes;
		{
			long int nirvana;
			elementsFile >> nirvana;
		}
	}
}

void GenerateUnstructuredGrid(std::vector<std::vector<long int> >& elements,std::vector<kaPoint<double> >& nodes, vtkSmartPointer<vtkUnstructuredGrid>& mesh)
{
	vtkSmartPointer<vtkPoints> meshNodes = vtkSmartPointer<vtkPoints>::New();
	meshNodes->Allocate(nodes.size());
	mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
    mat->SetNumberOfComponents(1);
    mat->SetName("Material");
    
	for(std::vector<kaPoint<double> >::iterator it = nodes.begin(); it!= nodes.end();it++)
	{
		meshNodes->InsertNextPoint(it->x,it->y,it->z);
	}
	mesh->SetPoints(meshNodes);
	for(std::vector<std::vector<long int> >::iterator it = elements.begin(); it!= elements.end();it++)
	{
		if(it->size()==5)
        {
            vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
            for(int i=0; i < 4; i++)
                tetra->GetPointIds()->SetId(i,it->at(i));
            mat->InsertNextValue(it->at(4));
            mesh->InsertNextCell(tetra->GetCellType(),tetra->GetPointIds());
            
        }
        else if(it->size()==11)
        {
            vtkSmartPointer<vtkQuadraticTetra> tetra = vtkSmartPointer<vtkQuadraticTetra>::New();
            for(int i=0; i < 10; i++)
                tetra->GetPointIds()->SetId(i,it->at(i));
            mat->InsertNextValue(it->at(10));
            mesh->InsertNextCell(tetra->GetCellType(),tetra->GetPointIds());
        }
        else
        {
            std::cout << "Only 4 or 10 nodes elements are allowed !";
            exit(0);
        }
	}
	mesh->GetCellData()->AddArray(mat);
}

void ExtractCavity(vtkPolyData* mesh, kaLattice<unsigned char>* lat,int value, int ele,double radius,vtkSmartPointer<vtkPolyData>& cavityMesh,const std::vector<int>& wantedMaterials)
{
	std::vector<std::vector<long int> > cavityElements;
	vtkSmartPointer<vtkPoints> nodes = mesh->GetPoints();

	long int numNodes = mesh->GetNumberOfPoints();
    long int numCells = mesh->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
	cavityMesh = vtkSmartPointer<vtkPolyData>::New();

	std::vector<long int> nodesMapping;
	nodesMapping.resize(numNodes,0);
	
    vtkIntArray* mat= vtkIntArray::SafeDownCast((mesh->GetCellData()->GetArray("Material")));
    
	for(int i=0;i<numCells;i++)
	{
		mesh->GetCellPoints(i,pt);
        if(wantedMaterials.size()!=0)
        {
            
            int m = mat->GetValue(i);
            if(find(wantedMaterials.begin(),wantedMaterials.end(),m)==wantedMaterials.end())
                continue;
        }
		bool isCavityNode[3] = {false,false,false};
		for(int j=0;j<3;j++)
		{
			double* p = nodes->GetPoint(pt->GetId(j));
			for(int k=-1; k<2; k++)
				for(int l=-1;l<2; l++)
					for(int m=-1;m<2; m++)
					{			
						kaPoint<double>a(p[0]+k*radius,p[1]+l*radius,p[2]+m*radius);
						long int voxel = FindCorrespondingVoxel(a, lat);
						if(voxel > 0 && voxel < lat->xyzLattice)
							if(lat->lat[voxel] == value)
								isCavityNode[j]=true;
					}
		}
		if(isCavityNode[0] && isCavityNode[1] && isCavityNode[2])
		{
			std::vector<long int> nodesIndexes;
			for(int j=2; j >= 0; j--)
			{
				nodesIndexes.push_back(pt->GetId(j));
			}
			cavityElements.push_back(nodesIndexes);
		}
		
	}
	
	vtkSmartPointer<vtkPoints> cavityPoints = vtkSmartPointer<vtkPoints>::New();
	std::set<long int> processedNodes;
	for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
	{
		for(std::vector<long int>::iterator it2 = it->begin(); it2 != it->end(); it2++)
		{
			pair <std::set<long int>::iterator, bool> ret =  processedNodes.insert(*it2);
			if(ret.second)
			{
				double* p = nodes->GetPoint(*it2);
				long int i = cavityPoints->InsertNextPoint(p[0],p[1],p[2]);
				nodesMapping[*it2] = i;
				*it2 = i;
			}
			else
			{
				*it2 = nodesMapping[*it2];
			}

		}
	}
	cavityMesh->SetPoints(cavityPoints);

	vtkSmartPointer<vtkCellArray> cavityTriangles = vtkSmartPointer<vtkCellArray>::New();
	
	for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
	{
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		triangle->GetPointIds()->SetId(0,it->at(0));
		triangle->GetPointIds()->SetId(1,it->at(1)); 
		triangle->GetPointIds()->SetId(2,it->at(2)); 
		cavityTriangles->InsertNextCell(triangle);
	}	
	cavityMesh->SetPolys(cavityTriangles);
}

void RemoveOneEdgeBoundaryCells(vtkSmartPointer<vtkPolyData>& mesh)
{
	std::cout << "Smoothing Boundary" << std::endl;
	vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
	
	vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> ne = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> kn = vtkSmartPointer<vtkIdList>::New();
	ne->SetNumberOfIds(1);
	kn->SetNumberOfIds(2);
	
	std::vector<std::vector<long int> > cavityElements;
	bool goAhead = true;
	while(goAhead)
	{
		goAhead = false;
		cavityElements.clear();
		long int numNodes = mesh->GetNumberOfPoints();

		std::vector<long int> nodesMapping;
		nodesMapping.resize(numNodes,0);

		vtkSmartPointer<vtkPoints> nodes = mesh->GetPoints();
		
		long int numCells = mesh->GetNumberOfCells();
		for(int i=0;i<numCells;i++)
		{
			mesh->GetCellPoints(i,pt);
			int boundaries = 0;
			for(int j=0;j<3;j++)
			{
				kn->SetId(0,pt->GetId((j==0)?0:2));
				kn->SetId(1,pt->GetId((j==1)?0:1));
				mesh->GetCellNeighbors(i, kn, ne);
				if( ne->GetNumberOfIds() == 0)
				{
					boundaries++;
					if(boundaries == 2)
						goAhead = true;
				}
			}
			if(boundaries < 2)
			{
				std::vector<long int> nodesIndexes;
				for(int j=0; j < 3; j++)
				{
					nodesIndexes.push_back(pt->GetId(j));
				}
				cavityElements.push_back(nodesIndexes);
			}
			
		}
		vtkSmartPointer<vtkPoints> cavityPoints = vtkSmartPointer<vtkPoints>::New();
		std::set<long int> processedNodes;
		for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
		{
			for(std::vector<long int>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			{
				pair <std::set<long int>::iterator, bool> ret =  processedNodes.insert(*it2);
				if(ret.second)
				{
					double* p = nodes->GetPoint(*it2);
					long int i = cavityPoints->InsertNextPoint(p[0],p[1],p[2]);
					nodesMapping[*it2] = i;
					*it2 = i;
				}
				else
				{
					*it2 = nodesMapping[*it2];
				}
				
			}
		}
		temp->SetPoints(cavityPoints);
		vtkSmartPointer<vtkCellArray> cavityTriangles = vtkSmartPointer<vtkCellArray>::New();
		
		for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
		{
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId(0,it->at(0));
			triangle->GetPointIds()->SetId(1,it->at(1)); 
			triangle->GetPointIds()->SetId(2,it->at(2)); 
			cavityTriangles->InsertNextCell(triangle);
		}	
		temp->SetPolys(cavityTriangles);
		mesh->DeepCopy(temp);
	}
}

void ShaveBoundary(vtkSmartPointer<vtkPolyData>& mesh, int passes)
{
	std::cout << "Shaving Boundary with " << passes << " passes" << std::endl;
	vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> ne = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> kn = vtkSmartPointer<vtkIdList>::New();
	ne->SetNumberOfIds(1);
	kn->SetNumberOfIds(2);
	std::vector<std::vector<long int> > cavityElements;
	
	for(int i=0; i < passes; i++)
	{
		cavityElements.clear();
		long int numNodes = mesh->GetNumberOfPoints();
		vtkSmartPointer<vtkPoints> nodes = mesh->GetPoints();
		
		std::vector<long int> nodesMapping;
		nodesMapping.resize(numNodes,0);
		
		
		
		long int numCells = mesh->GetNumberOfCells();
		for(int i=0;i<numCells;i++)
		{
			mesh->GetCellPoints(i,pt);
			bool isBoundaryElement = false;
			for(int j=0;j<3;j++)
			{
				kn->SetId(0,pt->GetId((j==0)?0:2));
				kn->SetId(1,pt->GetId((j==1)?0:1));
				mesh->GetCellNeighbors(i, kn, ne);
				if( ne->GetNumberOfIds() == 0)
				{
					isBoundaryElement=true;
				}
			}
			if(!isBoundaryElement)
			{
				std::vector<long int> nodesIndexes;
				for(int j=0; j < 3; j++)
				{
					nodesIndexes.push_back(pt->GetId(j));
				}
				cavityElements.push_back(nodesIndexes);
			}
			
		}
		vtkSmartPointer<vtkPoints> cavityPoints = vtkSmartPointer<vtkPoints>::New();
		std::set<long int> processedNodes;
		for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
		{
			for(std::vector<long int>::iterator it2 = it->begin(); it2 != it->end(); it2++)
			{
				pair <std::set<long int>::iterator, bool> ret =  processedNodes.insert(*it2);
				if(ret.second)
				{
					double* p = nodes->GetPoint(*it2);
					long int i = cavityPoints->InsertNextPoint(p[0],p[1],p[2]);
					nodesMapping[*it2] = i;
					*it2 = i;
				}
				else
				{
					*it2 = nodesMapping[*it2];
				}
				
			}
		}

		temp->SetPoints(cavityPoints);
		vtkSmartPointer<vtkCellArray> cavityTriangles = vtkSmartPointer<vtkCellArray>::New();
		
		for(std::vector<std::vector<long int> >::iterator it=cavityElements.begin(); it != cavityElements.end(); it++)
		{
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId(0,it->at(0));
			triangle->GetPointIds()->SetId(1,it->at(1)); 
			triangle->GetPointIds()->SetId(2,it->at(2)); 
			cavityTriangles->InsertNextCell(triangle);
		}	
		temp->SetPolys(cavityTriangles);
		mesh->DeepCopy(temp);
	}
}

void FillHoles(vtkSmartPointer<vtkPolyData>& mesh, double holesSize)
{
	std::cout << "Filling holes" << std::endl;
	vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter = vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHolesFilter->SetHoleSize(holesSize);	
	
	vtkSmartPointer<vtkFeatureEdges> boundaryEdges = vtkSmartPointer<vtkFeatureEdges>::New(); 
	boundaryEdges->FeatureEdgesOff(); 
	boundaryEdges->BoundaryEdgesOn();
	boundaryEdges->NonManifoldEdgesOn();
#if VTK_MAJOR_VERSION > 5
	boundaryEdges->SetInputData(mesh);
#else
	boundaryEdges->SetInput(mesh);
#endif
	boundaryEdges->Update();

	int numHoles = boundaryEdges->GetOutput()->GetNumberOfCells();
	int passes = 1;
	while(numHoles != 0)
	{
#if VTK_MAJOR_VERSION > 5
		fillHolesFilter->SetInputData(mesh);
#else
		fillHolesFilter->SetInputConnection(mesh->GetProducerPort());
#endif
		fillHolesFilter->Update();
		
		temp->DeepCopy(fillHolesFilter->GetOutput());
		mesh->DeepCopy(temp);
		
#if VTK_MAJOR_VERSION > 5
		boundaryEdges->SetInputData(mesh); 
#else
		boundaryEdges->SetInput(mesh); 
#endif
		boundaryEdges->Update();
		int i = numHoles;
		numHoles = (boundaryEdges->GetOutput()->GetNumberOfCells());
		if(i == numHoles)
		{
			std::cout << "Filling the holes failed" << std::endl;
			break;
		}
		passes++;
		std::cout << "Pass: " << passes << " Number of holes left: " << numHoles << std::endl;
	}
}

void CheckNormalOrientation(vtkSmartPointer<vtkPolyData>& mesh)
{
	vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyDataNormals> normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
    
    
    
	normalFilter->ConsistencyOn();
#if VTK_MAJOR_VERSION > 5
	normalFilter->SetInputData(mesh);
#else
	normalFilter->SetInput(mesh);
#endif
	normalFilter->AutoOrientNormalsOn();
	normalFilter->ComputeCellNormalsOn();
	normalFilter->Update();
	temp->DeepCopy(normalFilter->GetOutput());
	mesh->DeepCopy(temp);
	
	long int numElements = mesh->GetNumberOfCells();
	vtkSmartPointer<vtkPoints> nodes = mesh->GetPoints();

	vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> orientation = vtkSmartPointer<vtkDoubleArray>::New();

	orientation->SetNumberOfComponents(3);
	orientation->Allocate(3*numElements);
	orientation->SetName("Check Orientation");
	
	for(long int i=0; i<numElements; i++)
	{
		vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(i,pt);

		double p0[3]; 
		nodes->GetPoint(pt->GetId(0),p0);
		double p1[3]; 
		nodes->GetPoint(pt->GetId(1),p1);		
		double p2[3];
		nodes->GetPoint(pt->GetId(2),p2);		

		
		double p10[3] = {p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]};
		double p20[3] = {p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2]};
		double p3[3] = {0,0,0};
		vtkMath::Cross(p10,p20,p3);
		orientation->InsertNextTuple(p3);
	}
	mesh->GetCellData()->AddArray(orientation);	
}

void DetermineNodeMapping(vtkSmartPointer<vtkUnstructuredGrid> grid, vtkSmartPointer<vtkPolyData>& poly, std::vector<long>& nodesMapping)
{
	nodesMapping.clear();
	
	long int numPolyPoints = poly->GetNumberOfPoints();
	vtkSmartPointer<vtkPoints> polyPoints = poly->GetPoints();

	vtkSmartPointer<vtkPointLocator> locator = 	vtkSmartPointer<vtkPointLocator>::New();
	locator->SetDataSet( grid );
	locator->Update();    
	
	for(long int i=0; i < numPolyPoints; i++)
	{
        if(i%100 == 0)
            std::cout << i <<  " of " << numPolyPoints << "\n";
		double* p = polyPoints->GetPoint(i);
		long int index = locator->FindClosestPoint(p);
		nodesMapping.push_back(index);
	}	
}

int main(int argc, char* argv[])
{
	
	if(argc < 4)
	{
		std::cout << "VTPExtractSurfaceFromMesh <nodes_file.node> <elements_file.ele> <vtp.file> <element_index> <surface_index> <output_prefix> " << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "\t-free_surface_element_index <val>" << std::endl;
		std::cout << "\t-unit <val>" << std::endl;
		std::cout << "\t-append" << std::endl;
		std::cout << "\t-first_index <val>" << std::endl;
		std::cout << "\t-export_nodes_indices <exp_nodes_file.node>" << std::endl;
		std::cout << "\t-export_stand_alone <standalone_prefix>" << std::endl;
		exit(0);
	}
	try 
	{
		
		bool append = false;
		int firstIndex = -1;
        bool standAlone = false;
        bool exportNodesIndices = false;
        double unit=1;
        int freeSurfaceElementIndex = -1;
        std::string exportNodesIndicesFilename;
        std::string standAloneFilenamePrefix;
        vector<int> materials;
        
		for(int i=1; i<argc; i++) 
		{
			
            if(strcmp(argv[i], "-export_stand_alone")==0 && i<argc-1)
			{
                standAloneFilenamePrefix= std::string(argv[++i]);
				standAlone=true;
			}

			if(strcmp(argv[i], "-append")==0)
			{
				append=true;
			}
     
            
            if(strcmp(argv[i], "-free_surface_element_index")==0 && i<argc-1)
			{
                freeSurfaceElementIndex = atoi(argv[++i]);
			}
            
            if(strcmp(argv[i], "-unit")==0 && i<argc-1)
			{
                unit = atof(argv[++i]);
			}
            
            if(strcmp(argv[i], "-export_nodes_indices")==0 && i<argc-1)
			{
                exportNodesIndicesFilename = std::string(argv[++i]);
                exportNodesIndices = true;
			}
            
			if(strcmp(argv[i], "-first_index")==0 && i<argc-1)
			{
				int param =atoi(argv[++i]); 
				{
					if(param < 0)
					{
						std::cout << "Index of first surface element must be >= 0" << std::endl;
						exit(0);
					}
					else
					{
						firstIndex = param;
					}

				}
			}
		}
        
				
        std::ifstream nodesFile(argv[1]);
		if(!nodesFile.good())
			throw(std::string("Cannot open nodes file\n"));
		std::ifstream elementsFile(argv[2]);
		if(!elementsFile.good())
			throw(std::string("Cannot open elements file\n"));
		
		//Load Nodes to vector
		
		std::vector<kaPoint<double> > nodes;
		LoadNodes(nodesFile,nodes,unit);
		nodesFile.close();
		
		std::vector<std::vector<long int> > elements;
		LoadElements(elementsFile,elements);
        
		vtkSmartPointer<vtkUnstructuredGrid> mesh;
		GenerateUnstructuredGrid(elements,nodes,mesh);
        
        string Filename = argv[3];
        int SuffixPos = Filename.find_last_of(".");
        string Suffix = Filename.substr(SuffixPos+1,Filename.size()-SuffixPos);
        vtkSmartPointer<vtkPolyData> cavityMesh = vtkSmartPointer<vtkPolyData>::New();
        if(Suffix == "vtp"){
            std::cout << "Reading .vtp\n";
            vtkSmartPointer<vtkXMLPolyDataReader> vtpreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            vtpreader->SetFileName(argv[3]);
            vtpreader->Update();
            cavityMesh = vtpreader->GetOutput();
        }else if(Suffix == "stl"){
            std::cout << "Reading .stl\n";
            vtkSmartPointer<vtkSTLReader> stlreader = vtkSmartPointer<vtkSTLReader>::New();
            stlreader->SetFileName(argv[3]);
            stlreader->Update();
            cavityMesh = stlreader->GetOutput();
        }
        
        std::cout << cavityMesh->GetNumberOfCells() << endl;
        
		CheckNormalOrientation(cavityMesh);

		vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
#if VTK_MAJOR_VERSION > 5
		massProperties->SetInputData(cavityMesh);
#else
		massProperties->SetInput(cavityMesh);
#endif
		massProperties->Update();

		std::cout<< "Cavity Volume: " << massProperties->GetVolume() << std::endl;
		std::cout<< "Cavity Surface Area: " << massProperties->GetSurfaceArea() << std::endl;
		
		std::vector<long> nodesMapping;
        
        // Either find the surface nodes in the existing mesh or create a stand alone surface mesh with triangles and points, in that case the nodes indices of the surface mesh
        // don't need to be mapped
        
        DetermineNodeMapping(mesh, cavityMesh, nodesMapping);
        std::vector<bool> freeSurface;

        vtkSmartPointer<vtkIntArray> mat = vtkSmartPointer<vtkIntArray>::New();
        mat->SetNumberOfComponents(1);
        mat->SetName("Material");
       
        if(freeSurfaceElementIndex != -1)
            for(int i=0; i < cavityMesh->GetNumberOfCells(); i++)
            {
                vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
                cavityMesh->GetCellPoints(i,pt);
                
                vtkIdType p[3] = {nodesMapping.at(pt->GetId(0)),nodesMapping.at(pt->GetId(1)),nodesMapping.at(pt->GetId(2))};
                bool hasNeighbor = false;
                for(int j=0; j<mesh->GetNumberOfCells();j++)
                {
                    vtkSmartPointer<vtkIdList> pt2 = vtkSmartPointer<vtkIdList>::New();
                    mesh->GetCellPoints(j,pt2);
                    vtkIdType p2[4] = {pt2->GetId(0),pt2->GetId(1),pt2->GetId(2),pt2->GetId(3)};
                    int n=0;
                    for(int k=0;k < 3; k++)
                        for(int l=0; l <4; l++)
                            if(p2[l]==p[k])
                                n++;
                    if(n==3)
                        hasNeighbor = true;
                }
                if(hasNeighbor)
                    freeSurface.push_back(false);
                else
                    freeSurface.push_back(true);
                
            }
        else
            for(int i=0; i < cavityMesh->GetNumberOfCells(); i++)
                    freeSurface.push_back(false);
           
		std::string surfaceFile = std::string(argv[6]) + std::string(".sur");
		
		std::ofstream file;
		
		int amountOfElements = 0;
		int lastIndex = 0;
		
		long int numElements = cavityMesh->GetNumberOfCells();
		vtkSmartPointer<vtkIdList> pt = vtkSmartPointer<vtkIdList>::New();
        
		if(append)
		{
            if(standAlone)
                throw std::string("Appending isn't implemented for standalone surface meshes");
            
			std::ifstream inFile(surfaceFile.c_str());
			std::string str;
			getline(inFile, str);
			std::stringstream ss(str);
			ss >> amountOfElements;
			std::vector<std::string> lines;
			while(getline(inFile,str))
			{
				lines.push_back(str);
				std::stringstream ss(str);
				ss >> lastIndex;
			}
			inFile.close();
			file.open((surfaceFile).c_str());
			
			std::cout << "Amount of Elements: " << amountOfElements << " Last Index: " << lastIndex << std::endl;
			file << numElements+amountOfElements << " 3 2" << std::endl;
			for(std::vector<std::string>::iterator it=lines.begin();it!=lines.end();it++)
			{
				file <<  (*it) << std::endl;
			}
		}
		else
		{
			file.open(surfaceFile.c_str());
			file << numElements << " 3 2" << std::endl;
		}

        int elementIndex = atoi(argv[4]);
		
        if(freeSurfaceElementIndex == -1)
            freeSurfaceElementIndex = elementIndex;
        
		for(long int i=0; i<numElements; i++)
		{
			cavityMesh->GetCellPoints(i,pt);
			file <<  lastIndex + i+1 << " " << nodesMapping.at(pt->GetId(0))+1 << " " << nodesMapping.at(pt->GetId(1))+1 << " " << nodesMapping.at(pt->GetId(2))+1 << " " <<  (freeSurface[i]?freeSurfaceElementIndex:elementIndex) << " " << atoi(argv[5]) << std::endl;
		}
		file.close();
		
        if(standAlone)
        {
            std::ofstream nodesFile(std::string(standAloneFilenamePrefix+".nodes").c_str());
            std::ofstream surFile(std::string(standAloneFilenamePrefix+".sur").c_str());
            
            nodesFile << cavityMesh->GetNumberOfPoints() << " 3 1 0\n";
            for(int i=0; i < cavityMesh->GetNumberOfPoints(); i++)
            {
                nodesFile << i+1 << " " << cavityMesh->GetPoint(i)[0] << " " << cavityMesh->GetPoint(i)[1] << " " << cavityMesh->GetPoint(i)[2] << " 7 \n";
            }
            
            surFile << numElements << " 3 2" << std::endl;
            
            for(long int i=0; i<numElements; i++)
            {
                cavityMesh->GetCellPoints(i,pt);
                surFile <<  i+1 << " " << pt->GetId(0)+1 << " " << pt->GetId(1)+1 << " " << pt->GetId(2)+1 << " " << atoi(argv[2]) << " " << atoi(argv[3]) << std::endl;
            }
            
            nodesFile.close();
            surFile.close();
        }
        
        for(int i=0; i < numElements; i++)
            mat->InsertNextValue((freeSurface[i]?freeSurfaceElementIndex:elementIndex));

        cavityMesh->GetCellData()->AddArray(mat);
        if(exportNodesIndices)
        {
            std::ofstream nodesIndices(std::string(exportNodesIndicesFilename).c_str());
            for(auto i: nodesMapping)
                nodesIndices << i+1 << "\n";
            nodesIndices.close();
        }
		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		std::string vtpFile = std::string(argv[6]) + std::string(".vtp");
		writer->SetFileName(vtpFile.c_str());
#if VTK_MAJOR_VERSION > 5
                writer->SetInputData(cavityMesh);
#else
		writer->SetInput(cavityMesh);
#endif
		writer->Write();
		
		

		
		elementsFile.close();
		
	}	
	
	catch (std::string e) 
	{
		std::cout << e;
		exit(1);
	}
	catch(...) 
	{
        cerr << "Exception of unknown type!\n";
		exit(1);
    }
	

}

