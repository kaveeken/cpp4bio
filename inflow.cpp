/*****************************************
              inflow.cpp

  Function to handle inflow of resource and antibiotic into a system.
  Used in simulating collective abtibiotic resistance.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <vector>
#include <exception>
#include <cstdlib>

#include "project.h"

// inflow along upper i border of the matrix  is done by keeping the concentrations at a fixed value,
// as if the compartment can freely take chemicals from a neighbouring cappillary,
// which is assumed to have constant amounts of these.
void inflow(std::vector<std::vector<std::vector<double > > > &matrix)
{
	for(int j = 0; j < iN; ++j) {
		matrix[0][j][2] = 0.95;  // inflow
		matrix[0][j][3] = 1.3;
	       				// these should really be read from config
	}
}

