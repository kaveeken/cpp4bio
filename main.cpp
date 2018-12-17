#include <iostream>
#include <fstream>
#include <vector> 
#include <string>
#include <exception>
#include <cstdlib>
#include <cmath>

#include "project.h"

// these are read by readConfig(), maybe there is a better way to handle them
int    iN      = 0;    // matrix size
double tEnd    = 0.0;  // simulation length
double dtSav   = 0.0;  // interval at which is written to file
double diffXpH = 0.0;  // diffusion speed of bacteria
double diffYpH = 0.0;  // diffusion speed of Cm
double diffZpH = 0.0;  // diffusion speed of resource

// these are extern in header file
const int nvar 	= 6; // # of variables

// integration parameters
const double dt0 	= 0.001;	// starting step size


int main()
{
	try {
		readConfig(); // get configuration variables
		// things to track
		int iFile = 0;     // iterator for generating filenames
		double maxS = 0.0; // track maximum value
		double maxR = 0.0; // for normalization purposes

		//*********** initial conditions ********
		// make vector containing
		// Xs, Xr, z, ym, ys, yr
		std::vector<double> x(nvar);
		// set initial conditions
		x[0] = 0.05;
		x[1] = 0.0064; // results very dependent on this initial condition
		x[2] = 0.95;
		x[3] = 0.95;
		x[4] = 0.90;
		x[5] = 0.8;

		std::vector<double> xEmpty = x; // small values break diffundMatrix()
		xEmpty[0] = 0.0;
		xEmpty[1] = 0.0;

		std::vector<std::vector<double> > rowX(iN,xEmpty); // row vector
		std::vector<std::vector<std::vector<double > > > matX(iN,rowX); // matrix


		for(int i = 0; i < 9; ++i)
			for(int j = 0; j < 12; ++j)
			       matX[i][j] = x;	

		for(int j = 5; j < 8; ++j){
			matX[9][j][1] = 0.2;
			matX[10][j][1] = 0.2;
			matX[11][j][1] = 0.2;
		}


		// ********** simulation *************
		// start numerically integrating
		int nOK = 0, nStep = 0; // number of OK steps, number of total steps
		double dtMin = dt0, dtMax = kdMinH;
		for(double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep) {
			/*
			std::cout << matX[0][10][0] << ',' << matX[0][10][1] << ','
				<< matX[0][10][3]  << ',' << matX[0][10][2] << ',' << matX[0][10][4] << std::endl
				<< matX[1][10][0] << '\t' << matX[1][10][4] << std::endl
				<< matX[2][10][0] << '\t' << matX[2][10][4] << std::endl
				<< matX[3][10][0] << '\t' << matX[2][10][4] << std::endl;
				*/
			//std::cout << dt << std::endl;
			//std::cout << matX[11][0][0] << ' ' << matX[11][1][0] << std::endl
			//	<< matX[10][0][0] << ' ' <<  matX[10][1][0] << std::endl;

			// run and evaluate numerical integrator
			if(matrixDPAS(t,matX,dt)) { // this should encompass more of the below maybe
				++nOK; // step OK 
				diffundMatrix(matX, dt); // diffusion step
				inflow(matX);
//				std::cout << t << std::endl;
			} else 

			if(dt < dtMin) 
				dtMin = dt;
			else if(dt > dtMax)
				dtMax = dt;

			// sequence of number codes to track output
			// this can just be 01 to 09, followed by just the value of i up to 99
			std::string numbers = "01020304050607080910111213141516171819202122";
			if(t > tsav){ // this maybe happens on bad integration step if it is the first step
				std::cout << t << std::endl;
				// make filename
				int index = (iFile) * 2;

				std::string strFileNameS = "data_S_" 
					     + numbers.substr(index,2) + ".csv";
				std::string strErrorS = "unable to open " 
							  + strFileNameS + "\n";

				// open file
				std::ofstream fileStreamS(strFileNameS);
				if(!fileStreamS.is_open())
					throw std::runtime_error(strErrorS);

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


				std::string strFileNameR = "data_R_" 
					     + numbers.substr(index,2) + ".csv";
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

				tsav += dtSav;

				++iFile;
			}
		}
		// report integration statistics
		std::cout << "steps: " << nStep << std::endl
			<< "proportion bad steps: " << 1.0 - nOK * 1.0 / nStep 
						    << std::endl
			<< "avg step size: " << tEnd / nStep << std::endl
			<< "min step size: " << dtMin << std::endl
			<< "max step size: " << dtMax << std::endl;

		// write a file containing the relevant information
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

			<< "// Initial Conditions\n"

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
			<< "max step size: " << dtMax << std::endl;
	}
	catch(std::exception &error) {
		std::cerr << "error: " << error.what();
		exit(EXIT_FAILURE);
	}
	return 0;
}
