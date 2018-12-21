/*****************************************
              buildSystem.cpp

  Function to handle initialization of a 2d spatial system.
  Used in simulating collective abtibiotic resistance.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <vector>
#include <iostream>

#include "project.h"

// this function initializes the system as a matrix with certain conditions and a local population of bacteria
// some of the values here should probably be read from config file
std::vector<std::vector<std::vector<double> > > buildSystem()
{
	// *** single-cell starting conditions ***
	// populated cell conditions
	std::vector<double> vecPop(nvar);
	vecPop[0] = 0.005; 	// density of susceptible bacteria
	vecPop[1] = 0.1;	// density of resistant bacteria	
	vecPop[2] = 0.95;	// concentration of growth-limiting resource
	vecPop[3] = 1.3;	// concentration of antibiotic
	vecPop[4] = 1.3;	// internal concentration of antibiotic for susceptible strain
	vecPop[5] = 0.8;	// internal concentration of antibiotic for resistant strain
	if(vecPop[0] + vecPop[1] > 1)
		throw std::runtime_error("Initial bacteria density > 1\n");
	// conditions w/o bacteria
	std::vector<double> vecEmp(nvar);
	vecEmp[0] = 0.0;
	vecEmp[1] = 0.0; 
	vecEmp[2] = 0.95;
	vecEmp[3] = 1.3;
	vecEmp[4] = 0.0; 
	vecEmp[5] = 0.0;

	// *** building a(n empty) system ***
	std::vector<std::vector<double> > rowX(iN,vecEmp); // row vector
	std::vector<std::vector<std::vector<double > > > matX(iN,rowX); // matrix

	// *** populating some place ***
	for(int i = iN / 2 - 1; i < iN / 2 + 1; ++i) 
		for(int j = 0; j < iPopCol; ++j) 
			matX[i][j] = vecPop;

	return matX;
}
