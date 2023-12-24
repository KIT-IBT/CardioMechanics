/*
 * File: FixNodes.cpp
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
#include<fstream>
#include<sstream>
#include<string>

using namespace std;

int main(int argc, char* argv[])
{   
    if(argc < 4)
    {
        cout << "FixNodes <input.node> <output.node> <fixation.list> \n\t-val <val> (optional) \n\t-keep_original (optional) Modifiy only nodes from list, if not set, all other nodes are set to 0 (unfixed)" << endl;
        
        exit(0);
    }
    
	ifstream inFile(argv[1]);
	ofstream outFile(argv[2]);
	ifstream fixFile(argv[3]);

    
    int val = 7;
    bool keepOriginal = false;
    
    for(int i=0; i < argc; i++)
    {
        if(string(argv[i])=="-val" && i+1 < argc)
            val = atoi(argv[++i]);
        if(string(argv[i])=="-keep_original")
            keepOriginal = true;
    }

    if(!inFile.good() || !fixFile.good() )
    {
        cout << "Check the file names, at least one is wrong, but i won't tell you which it is, arrr arr arr!" << endl;
        exit(0);
    }

	
	int numNodes,dim,attr,nirvana;
	inFile >> numNodes;
	inFile >> dim;
	inFile >> attr;
	inFile >> nirvana;
	
	bool* fixedNodes = new bool[numNodes];
	for(int i=0; i < numNodes; i++)
		fixedNodes[i] = false;

	string str;

	while (getline(fixFile, str))
	{
		stringstream ss(str);
		int i;
		ss >> i;
        /*if(!ss.eof())
        {
            cout << "fixation.list is corrupt" << endl;
            exit(0);
        }*/
		if(i+1 <= numNodes)
			fixedNodes[i+1] = true;
		else
		{
			cout << "Node " << i << " is out of range " << endl;
			exit(0);
		}
	}
	
	outFile << numNodes << " " << dim << " 1 0" << endl;
	

	getline(inFile, str); // ignore first line
	
	while (getline(inFile, str))
	{
		
		stringstream ss(str);
		int i;
		double x,y,z;
		ss >> i;
		ss >> x;
		ss >> y;
		ss >> z;
       
		outFile << i << " " << x << " " << y << " " << z << " ";
		if(attr == 1)
		{
			int a;
			ss >> a;
			if(fixedNodes[i])
				outFile << val << endl;
			else
                if(keepOriginal)
                    outFile << a << endl;
                else
                    outFile << 0 << endl;
		}
		else
		{
			if(fixedNodes[i])
				outFile << val << endl;
			else
				outFile << 0 << endl;
		}
	}
	inFile.close();
	fixFile.close();
	outFile.close();
}
