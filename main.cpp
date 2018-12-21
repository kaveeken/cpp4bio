/*****************************************
              main.cpp

  Main program for simulating collective antibiotic resistance.
  This program simulates growth and diffusion of bacteria in a 2d system,
  when under influence of an antibiotic and a growth-limiting resource.
  At set intervals, .csv files containing bacterial densities are written.
  
  At runtime, this program requires a config.txt file. 
  See readConfig.cpp for details.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <iostream>
#include <fstream>
#include <vector> 
#include <string>
#include <exception>
#include <cstdlib>
#include <cmath>

#include "project.h"

// these are read by readConfig(), maybe there is a better way to handle them so we can const them
int    iN      = 0;    // matrix size
double tEnd    = 0.0;  // simulation length
double dtSav   = 0.0;  // interval at which is written to file
double diffXpH = 0.0;  // diffusion speed of bacteria
double diffYpH = 0.0;  // diffusion speed of Cm
double diffZpH = 0.0;  // diffusion speed of resource
int    iPopCol = 0;    // amount of populated columns

// these are extern in header file
const int nvar 	= 6; // # of variables

// integration parameters
const double dt0 	= 0.001;	// starting step size


int main()
{
	try {
		// get configuration variables
		readConfig(); 

		// things to track:
		int iFile = 0;     // iterator for generating filenames
		double maxS = 0.0; // track maximum values
		double maxR = 0.0; // for downstream normalization purposes

		// make the initial system
		std::vector<std::vector<std::vector<double> > > matX = buildSystem();

		// ********** simulation *************
		// start numerically integrating
		int nOK = 0, nStep = 0; // number of OK steps, number of total steps
		double dtMin = dt0, dtMax = kdMinH;
		for(double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep) {
			// we should maybe get out of this loop when system has reached a steady state
			// we would need to statistically compare the current system with the previous one
			// which would be annoying
			// the simulation seems to slow down on (some) steady state systems, probably due to densities summing to almost 1
		 	// this summing to 1 problem is sort of handled in matrixDPAS()
			
			// run and evaluate numerical integrator
			if(matrixDPAS(t,matX,dt)) { 
				++nOK; // step OK 
				diffundMatrix(matX, dt); // diffusion step
				inflow(matX);
				//std::cout << t << std::endl;
			} else if(dt < dtMin)
				dtMin = dt;
			else if(dt > dtMax)
				dtMax = dt;

			// sequence of padded number codes to track output
			// right now this limits the amount of writes we can do
			// can be changed to 01 to 09, followed by just using the value of i up to 99
			// if we adjust below logic
			std::string numbers = "01020304050607080910111213141516171819202122";
			// this is made every loop, better to pull it outside of loop

			// ******* writing to file ********** 
			if(t > tsav){ 
				std::cout << t << std::endl; // this is the only hint to see if the program is running correctly

				int index = (iFile) * 2;

				// generate filename for susceptible strain iteratively
				std::string strFileNameS = "data_" 
					     + numbers.substr(index,2) + "_S.csv";
				std::string strErrorS = "unable to open " 
							  + strFileNameS + "\n";

				// open file based on above filename
				std::ofstream fileStreamS(strFileNameS);
				if(!fileStreamS.is_open())
					throw std::runtime_error(strErrorS);

				// start writing to file
				for(int i = 0; i < iN; ++i) {
					for(int j = 0; j < iN; ++j) {
						if(matX[i][j][0] > maxS)
							maxS = matX[i][j][0];
						fileStreamS 
							<< matX[i][j][0] << ','; // leaves trailing , which is sometimes annoying
					}
					fileStreamS << std::endl;
				}
				fileStreamS.close();


				// same as above, but for resistant strain
				std::string strFileNameR = "data_" 
					     + numbers.substr(index,2) + "_R.csv";
				std::string strErrorR = "unable to open " 
							  + strFileNameR + "\n";

				std::ofstream fileStreamR(strFileNameR);
				if(!fileStreamR.is_open())
					throw std::runtime_error(strErrorR);

				for(int i = 0; i < iN; ++i) {
					for(int j = 0; j < iN; ++j) {
						if(matX[i][j][1] > maxR)
							maxR = matX[i][j][1];
						fileStreamR 
						        << matX[i][j][1] << ',';
					}
					fileStreamR << std::endl;
				}
				fileStreamR.close();

				tsav += dtSav; // set time of next save

				++iFile;
			}
		}
		// report integration statistics
		std::cout << "steps: " << nStep << std::endl
			<< "proportion bad steps: " << 1.0 - nOK * 1.0 / nStep 
						    << std::endl
			<< "avg step size: " << tEnd / nStep << std::endl
			<< "min step size: " << dtMin << std::endl // this reports the wrong value
			<< "max step size: " << dtMax << std::endl;

		double maxMax = maxR; // another variable for normalization purposes
		if(maxS > maxR)
			maxMax = maxS;
	

		// write a file containing some relevant information
		std::ofstream fileStream("parameters.txt");
		if(!fileStream.is_open())
			throw std::runtime_error("unable to open parameters.txt\n");
		fileStream <<  "// Integration parameters\n"
			<< "Matrix size iN: \t\t\t" << iN << std::endl
			<< "Number of variables nvar: \t\t" << nvar << std::endl
			<< "Error tolerance: \t\t\t" << tolerance << std::endl	
                        << "Max step size decrease: \t\t" << kdShrinkMax << std::endl
                        << "Max step size increase: \t\t" << kdGrowMax << std::endl
                        << "Step size control safety factor\t\t" << kdSafety    << std::endl
                        << "Minimum step size: \t\t\t" << kdMinH << std::endl
                               	
                        << "Minimum population density: \t\t" << lowBound 	 << std::endl
			<< "// Model Parameters\n"
			<< "Max growth rate r: \t\t\t" << r << std::endl
			<< "Relative growth rate eta: \t\t\t" << eta << std::endl
			<< "Resource consumption rate c: \t\t\t" << c << std::endl
			<< "Growth function half-saturation rate kz: \t\t" << kz << std::endl
			<< "Growth function inhibitory constant hy: \t\t" << hy << std::endl
			<< "Passive transport of Cm between cell compartments p: \t" << p << std::endl
			<< "Maximum rate d: \t\t\t" << d << std::endl
			<< "Half-saturation constant ky: \t\t" << ky << std::endl

			<< "// Initial Conditions\n" // this one is tricky because we would like to qualify the initial density distribution

			<< "// Other\n"
			<< "Bacteria diffusion per time: \t\t" << diffXpH << std::endl
			<< "Cm diffusion per time: \t\t\t" << diffYpH << std::endl
			<< "Resource diffusion per time: \t\t" << diffZpH << std::endl
			<< "Maximum S density: \t\t\t" << maxS << std::endl
			<< "Maximum R density: \t\t\t" << maxR << std::endl
			<< "Number of files written (before extinction):\t" << iFile << std::endl
			<< "Simulation time: \t\t\t" << tEnd << std::endl

			<< "// Integration Statistics\n" 
		        << "steps: " << nStep << std::endl
			<< "proportion bad steps: " << 1.0 - nOK * 1.0 / nStep 
						    << std::endl
			<< "avg step size: " << tEnd / nStep << std::endl
			<< "min step size: " << dtMin << std::endl
			<< "max step size: " << dtMax << std::endl
			<< "// For reading in bash\n" 
			<< "sMax = " << maxS << std::endl 
			<< "rMax = " << maxR << std::endl
			<< "maxMax = " << maxMax << std::endl; 
			
	}
	// every error is handled the same
	catch(std::exception &error) {
		std::cerr << "error: " << error.what();
		exit(EXIT_FAILURE);
	}
	return 0;
}
