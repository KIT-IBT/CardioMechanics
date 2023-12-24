/*
 * File: ConvertT4toT10.cpp
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

#include<iostream>
#include<vector>
#include<fstream>
#include<map>
#include<stdexcept>
#include<cmath>



class Node
{
public:	
	void clear(){x = 0; y = 0; z = 0;attributes.clear();}
	double x;
	double y;
	double z;
	std::vector<long int> attributes;
};

class T4Element
{
public:
	void clear(){for(int i=0;i<4;i++)nodesIndices[i]=0;attributes.clear();}
	double nodesIndices[4];
	std::vector<long int> attributes;
};

class T10Element 
{
public:
	void clear(){for(int i=0;i<10;i++)nodesIndices[i]=0;attributes.clear();}
	double nodesIndices[10];
	std::vector<long int> attributes;
};

class T3Triangle
{
public:
	void clear(){for(int i=0;i<3;i++)nodesIndices[i]=0;attributes.clear();}
	double nodesIndices[3];
	std::vector<long int> attributes;
};

class T6Triangle
{
public:
	void clear(){for(int i=0;i<6;i++)nodesIndices[i]=0;attributes.clear();}
	double nodesIndices[6];
	std::vector<long int> attributes;
};

bool IsSurfaceOf(T3Triangle& t3, T10Element& t10)
{
    int numFound = 0;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<10; j++)
        {
            if (t3.nodesIndices[i] == t10.nodesIndices[j])
                numFound++;
        }
    }
    if (numFound == 3)
        return true;
    return false;
};

void LoadNodes(std::ifstream& nodesFile,std::vector<Node>& nodes)
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
		
		Node node;
		
		nodesFile >> currentNode;
		nodesFile >> node.x;
		nodesFile >> node.y;
		nodesFile >> node.z;
		if(attr == 1)
		{
			long int a;
			nodesFile >> a;
			node.attributes.push_back(a);
		}
		nodes.push_back(node);
	}
}

void LoadElements(std::ifstream& elementsFile,std::vector<T4Element>& elements)
{
	int numElements=0;
	int nodesPerElement=0;
	long int currentElement=0;
	int attr=0;
	
	elementsFile >> numElements;
	elementsFile >> nodesPerElement;		
	elementsFile >> attr;
	
	if(nodesPerElement != 4)
		throw(std::string("Just for 4-node tetrahedrons"));
	
	//	if(attr != 0)
	//	throw(std::string("Attributes or Boundary Marker not yet implemented !"));
	
	
	for(long int i=0; i < numElements; i++)
	{
		
		// Load Elements
		if(elementsFile.eof())
			throw(std::string("Elements file is corrupt. Check number of nodes !"));
		
		
		elementsFile >> currentElement;
		T4Element ele;

		for(int j=0;j<nodesPerElement;j++)
		{
			long int n;
			elementsFile >> n;
			ele.nodesIndices[j] = n-1;
		}

		for(int j=0;j<attr;j++)
		{
			long int a;
			elementsFile >> a;
			ele.attributes.push_back(a);
		}
        
		elements.push_back(ele);
	}
}

void LoadTriangles(std::ifstream& elementsFile, std::vector<T3Triangle>& elements)
{
	int numElements=0;
	int nodesPerElement=0;
	long int currentElement=0;
	int attr=0;
	
	elementsFile >> numElements;
	elementsFile >> nodesPerElement;		
	elementsFile >> attr;
	
	if(nodesPerElement != 3)
		throw(std::string("Just for 3-node triangles"));
	
	//	if(attr != 0)
	//	throw(std::string("Attributes or Boundary Marker not yet implemented !"));
	
	
	for(long int i=0; i < numElements; i++)
	{
		
		// Load Elements
		if(elementsFile.eof())
			throw(std::string("Surface file is corrupt. Check number of nodes !"));
		
		
		elementsFile >> currentElement;
		T3Triangle tri;

		for(int j=0;j<nodesPerElement;j++)
		{
			long int n;
			elementsFile >> n;
			tri.nodesIndices[j] = n-1;
		}

		for(int j=0;j<attr;j++)
		{
			long int a;
			elementsFile >> a;
			tri.attributes.push_back(a);
		}
        
		elements.push_back(tri);
	}
}

void CheckNodeSorting(std::vector<Node>& nodes, std::vector<T4Element>& t4)
{
    for(auto& e:t4)
    {
        
        double nodesCoords[12] = {
                                    nodes.at(e.nodesIndices[0]).x,nodes.at(e.nodesIndices[0]).y,nodes.at(e.nodesIndices[0]).z,
                                    nodes.at(e.nodesIndices[1]).x,nodes.at(e.nodesIndices[1]).y,nodes.at(e.nodesIndices[1]).z,
                                    nodes.at(e.nodesIndices[2]).x,nodes.at(e.nodesIndices[2]).y,nodes.at(e.nodesIndices[2]).z,
                                    nodes.at(e.nodesIndices[3]).x,nodes.at(e.nodesIndices[3]).y,nodes.at(e.nodesIndices[3]).z
                                 };
       
        
        double z41 = (nodesCoords[11] - nodesCoords[2]);
        double z31 = (nodesCoords[ 8] - nodesCoords[2]);
        double z21 = (nodesCoords[ 5] - nodesCoords[2]);
        double y41 = (nodesCoords[10] - nodesCoords[1]);
        double y31 = (nodesCoords[ 7] - nodesCoords[1]);
        double y21 = (nodesCoords[ 4] - nodesCoords[1]);
        double x41 = (nodesCoords[ 9] - nodesCoords[0]);
        double x31 = (nodesCoords[ 6] - nodesCoords[0]);
        double x21 = (nodesCoords[ 3] - nodesCoords[0]);
        
        double v = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);
        if(v < 0)
        {
            int a = e.nodesIndices[1];
            e.nodesIndices[1] = e.nodesIndices[2];
            e.nodesIndices[2] = a;
        }

    }
}

void ConvertT4toT10(std::vector<Node>& nodes, std::vector<T4Element>& t4, std::vector<T10Element>& t10)
{
	std::map<std::pair<long int,long int>,long int> m;
	int cEle = 0;
	for(std::vector<T4Element>::iterator it=t4.begin(); it!=t4.end(); it++)
	{
		cEle++;
		std::cout << cEle << " of " << t4.size() << std::endl;
		T10Element e;
		e.clear();
		for(int i=0; i < 4; i++)
			e.nodesIndices[i] = it->nodesIndices[i];
		int nodesIndex = 4;
		
		for(int i=0; i < 6;i++)
		{
			Node n;
			n.clear();
			
			int j=-1;
			int k=-1;
			switch (i)
			{
				case 0:
					j=0;k=1;
					break;
				case 1:
					j=2;k=1;
					break;
				case 2:
					j=0;k=2;
					break;
				case 3:
					j=0;k=3;
					break;
				case 4:
					j=3;k=1;
					break;
				case 5:
					j=3;k=2;
					break;
				default:
					break;
			}
			{
				int smallerIndex,biggerIndex;
				if(it->nodesIndices[j] < it->nodesIndices[k])
				{
					smallerIndex = it->nodesIndices[j];
					biggerIndex = it->nodesIndices[k];
				}
				else
				{
					smallerIndex = it->nodesIndices[k];
					biggerIndex = it->nodesIndices[j];
				}

				std::map<std::pair<long int,long int>,long int>::iterator it2 = m.find(std::pair<long int,long int>(smallerIndex,biggerIndex));
				if(it2 == m.end())
				{													
					n.x = ( nodes.at(it->nodesIndices[j]).x + nodes.at(it->nodesIndices[k]).x ) / 2;
					n.y = ( nodes.at(it->nodesIndices[j]).y + nodes.at(it->nodesIndices[k]).y ) / 2;
					n.z = ( nodes.at(it->nodesIndices[j]).z + nodes.at(it->nodesIndices[k]).z ) / 2;
				
					if(nodes.at(0).attributes.size() != 0)
					{
						
						std::vector<long int>::iterator a = nodes.at(it->nodesIndices[j]).attributes.begin();
						std::vector<long int>::iterator b = nodes.at(it->nodesIndices[k]).attributes.begin();
						
						while (a != nodes.at(it->nodesIndices[j]).attributes.end())
						{
							if(*a == *b)
								n.attributes.push_back(*a);
							else
								n.attributes.push_back(0);
							a++;
							b++;
						}
					}
					nodes.push_back(n);
					e.nodesIndices[nodesIndex] = nodes.size()-1;
					m.insert(std::pair< std::pair<long int,long int>,long int >( std::pair<long int,long int>(smallerIndex,biggerIndex),e.nodesIndices[nodesIndex]));
				}
				else
				{
					e.nodesIndices[nodesIndex] = it2->second;
				}
				nodesIndex++;
			}
		}

		e.attributes = it->attributes;
		t10.push_back(e);
	}
}


double norm2(Node p, Node q) {
    double norm2 = 0;
    std::vector<double> p1;
    p1.push_back(p.x);
    p1.push_back(p.y);
    p1.push_back(p.z);
    std::vector<double> p2;
    p2.push_back(q.x);
    p2.push_back(q.y);
    p2.push_back(q.z);
    for (int i=0; i<3; i++)
        norm2 += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return sqrt(norm2);
}


int FindNodeIndex(std::vector<Node>& nodes, Node& p, Node& q, T10Element& t10)
{
    Node r;
    r.x = (p.x + q.x)/2.;
    r.y = (p.y + q.y)/2.;
    r.z = (p.z + q.z)/2.;
    double dist2 = norm2(r, nodes[t10.nodesIndices[0]]);
    double dist2Min = dist2;
    int idxMin = 0;
    for (int i=0; i<10; i++)
    {
        dist2 = norm2(r, nodes[t10.nodesIndices[i]]);
        if (dist2 < dist2Min)
        {
            dist2Min = dist2;
            idxMin = t10.nodesIndices[i];
        }
    }
    return idxMin;
}


void ConvertT3toT6(std::vector<Node>& nodes, std::vector<T3Triangle>& t3, std::vector<T10Element>& t10, std::vector<T6Triangle>& t6)
{

    T6Triangle e;
    std::vector<double> p(3,0);
    int n4 = -1;
    int n5 = -1;
    int n6 = -1;
    for (int i=0; i<t3.size(); i++)
        for (int j=0; j<t10.size(); j++)
            // find corresponding t10 element
            if (IsSurfaceOf(t3[i], t10[j])) {
                e.clear();
                // find indices of intermediate points
                n4 = FindNodeIndex(nodes, nodes[t3[i].nodesIndices[0]], nodes[t3[i].nodesIndices[1]], t10[j]);
                n5 = FindNodeIndex(nodes, nodes[t3[i].nodesIndices[1]], nodes[t3[i].nodesIndices[2]], t10[j]);
                n6 = FindNodeIndex(nodes, nodes[t3[i].nodesIndices[2]], nodes[t3[i].nodesIndices[0]], t10[j]);
                e.nodesIndices[0] = t3[i].nodesIndices[0];
                e.nodesIndices[1] = t3[i].nodesIndices[1];
                e.nodesIndices[2] = t3[i].nodesIndices[2];
                e.nodesIndices[3] = n4;
                e.nodesIndices[4] = n5;
                e.nodesIndices[5] = n6;
                e.attributes = t3[i].attributes;
                t6.push_back(e);
            }
}


int main(int argc, char* argv[])
{

	if(argc < 4)
	{
		std::cout << "ConvertT4toT10 <Nodes_File.node> <Elements_File.ele> <Output>" << std::endl;
		std::cout << "ConvertT4toT10 <Nodes_File.node> <Elements_File.ele> <surfaceT3.sur> <Output>" << std::endl;
		exit(0);
	}
	
    bool useSur = false;
    if (argc==4) {
        useSur = false;
    }
    if (argc==5) {
        useSur = true;
    }

	try 
	{

        char* nodesInFname = NULL;
        char* eleInFname   = NULL;
        char* surInFname   = NULL;
        char* vtuOutFname  = NULL;
        char nullString[1] = {'\0'};
        nodesInFname = argv[1];
        eleInFname = argv[2];
        if (argc==4)
        {
            surInFname = nullString;
            vtuOutFname = argv[3];
        }
        if (argc==5)
        {
            surInFname = argv[3];
            vtuOutFname = argv[4];
        }
        std::cout << "ou0 " << argc << std::endl;
        std::cout << "ou1 " << nodesInFname << std::endl;
        std::cout << "ou2 " << eleInFname << std::endl;
        std::cout << "ou3 " << surInFname << std::endl;
        std::cout << "ou4 " << vtuOutFname << std::endl;

		std::ifstream nodesFile(nodesInFname);
		if(!nodesFile.good())
			throw(std::string("Cannot open nodes file\n"));
		std::ifstream elementsFile(eleInFname);
		if(!elementsFile.good())
			throw(std::string("Cannot open elements file\n"));
        std::ifstream surFile(surInFname);
        if (useSur)
        {
            if(!surFile.good())
                throw(std::string("Cannot open surface file\n"));
        }
	
		std::vector<Node> nodes;
		LoadNodes(nodesFile,nodes);				
		nodesFile.close();
		
		std::vector<T4Element> t4;
		LoadElements(elementsFile,t4);
		elementsFile.close();

        std::vector<T3Triangle> t3;
        if (useSur)
        {
            LoadTriangles(surFile,t3);
            surFile.close();
        }

		std::vector<T10Element> t10;
        CheckNodeSorting(nodes,t4);
		std::vector<T6Triangle> t6;
        
		ConvertT4toT10(nodes,t4,t10);

        if (useSur)
            ConvertT3toT6(nodes,t3,t10,t6);
		
		std::ofstream eleOutput(std::string(vtuOutFname + std::string(".ele")).c_str(),std::ios::out | std::ios::binary);
		if(!eleOutput.good())
			throw(std::string("Cannot write to output file\n"));

		std::ofstream nodesOutput(std::string(vtuOutFname + std::string(".node")).c_str(),std::ios::out | std::ios::binary);
		if(!nodesOutput.good())
			throw(std::string("Cannot write to output file\n"));

		std::ofstream surOutput;
        if (useSur)
        {
            surOutput.open(std::string(vtuOutFname + std::string(".sur")).c_str(),std::ios::out | std::ios::binary);
            if(!surOutput.good())
                throw(std::string("Cannot write to output file\n"));
        }

        // output nodes

		nodesOutput << nodes.size() << " 3 " << nodes.at(0).attributes.size() << " 0 " << std::endl;
		for(int i=0; i < nodes.size(); i++)
		{
			nodesOutput << i+1 << " " << nodes.at(i).x << " " << nodes.at(i).y << " " << nodes.at(i).z;
			for(int j=0; j < nodes.at(i).attributes.size(); j++)
				nodesOutput << " " << nodes.at(i).attributes.at(j);
			nodesOutput << std::endl;
		}
		nodesOutput.close();

        // output elements
		
		eleOutput << t10.size() << " 10 " << t10.at(0).attributes.size() << std::endl;
		for(int i=0; i < t10.size(); i++)
		{
			eleOutput << i+1;
			for(int j=0; j < 10; j++)
			{
				eleOutput << " " << t10.at(i).nodesIndices[j]+1;
			}
			for(int j=0; j < t10.at(i).attributes.size(); j++)
				eleOutput << " " << t10.at(i).attributes.at(j);
			eleOutput << std::endl;
		}
		eleOutput.close();

        // output surface

        if (useSur)
        {
            surOutput << t6.size() << " 6 " << t6.at(0).attributes.size() << std::endl;
            for(int i=0; i < t6.size(); i++)
            {
                surOutput << i+1;
                for(int j=0; j < 6; j++)
                {
                    surOutput << " " << t6.at(i).nodesIndices[j]+1;
                }
                for(int j=0; j < t6.at(i).attributes.size(); j++)
                    surOutput << " " << t6.at(i).attributes.at(j);
                surOutput << std::endl;
            }
            surOutput.close();
        }
		

	}	
	
	catch (std::string e) 
	{
		std::cout << e;
		exit(1);
	}
	catch(...) 
	{
        std::cerr << "Exception of unknown type!\n";
		exit(1);
    }
	
}
