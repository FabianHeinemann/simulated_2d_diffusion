#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "CSimulatedDiffusion.h"

int main(int argc, const char* argv[])
{
	string entryStr;
	CSimulatedDiffusion* simulation;

	if (argc >= 2)
	{		
		std::ifstream fin(argv[1]);
		YAML::Parser parser(fin);
		YAML::Node	descriptionNode;

		// Go through all simulations defined in the yaml file
		while(parser.GetNextDocument(descriptionNode))
		{					
			simulation = new CSimulatedDiffusion;
			// Load the voronoi mesh
			descriptionNode["voronoimesh"] >> entryStr;
			if (strcmp(entryStr.c_str(), "none") != 0)
			{
				simulation->loadFile(entryStr);
			}
			else
			{
				cout << "No mesh" << endl;
			}
			cout << endl;

			// Name of the result file
			descriptionNode["resultfile"] >> entryStr;			
			simulation->outputFileStr = entryStr;

			// pjump
			descriptionNode["pjump"] >> entryStr;
			simulation->pjump = atof(entryStr.c_str());

			// Time
			descriptionNode["Tmax"] >> entryStr;
			simulation->Tmax = atof(entryStr.c_str());

			// Write simulation parameters
			cout << "Output file: \t" << simulation->outputFileStr << endl;
			cout << "Pjump: \t\t" << simulation->pjump << endl;
			cout << "Tmax: \t\t" << simulation->Tmax << " s" << endl;
			
			// Start the simulation
			cout << endl;
			cout << "Simulating..." << endl << endl;
			simulation->start();
			
			cout << "Done!" << endl << endl;
			cout << "-------------------------------" << endl << endl;

			// Delete the simulation to free memory
			delete simulation;
		}	

		cout << "All simulations finished." << endl << endl;
		cin >> argc;
	}
	return 0;
}

