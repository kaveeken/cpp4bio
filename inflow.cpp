#include <iostream>
#include <vector>
#include <exception>
#include <cstdlib>

#include "project.h"

// inflow along upper i border
void inflow(std::vector<std::vector<std::vector<double > > > &matrix)
{
	for(int j = 0; j < iN; ++j) {
		matrix[0][j][2] = 0.8;
		matrix[0][j][3] = 0.95;
	}
}

