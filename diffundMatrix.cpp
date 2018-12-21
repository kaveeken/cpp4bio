/*****************************************
              diffundMatrix.cpp

  Function to handle diffusion in a 2d system.
  Used in simulating collective abtibiotic resistance.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <iostream>
#include <vector>
#include <exception>
#include <cstdlib>

#include "project.h"

//function to simulate diffusion between compartments (of Xs, Xr and Ym)
void diffundMatrix(std::vector<std::vector<std::vector<double> > > &matrix, // original matrix
		   const double &h) // step size
{
	const double diffX = h * diffXpH;
	const double diffY = h * diffYpH;
	const double diffZ = h * diffZpH;
	std::vector < std::vector < std::vector<double> > > matrixUpdate 
								       = matrix; // sequentialy updated
	//generate vectors used in the loops
	std::vector<double> concAbove(iN);
	std::vector<double> concBelow(iN);
	for (int i = 0; i < iN; ++i) {
		for (int j = 0; j < iN; ++j) {
			//if-else statements to close matrix above and below
			// closed borders around i-axis
			if (i != 0)
				concAbove = matrix[i - 1][j];
			else
				concAbove = matrix[i][j];
			if (i != iN - 1)
				concBelow = matrix[i + 1][j];
			else
				concBelow = matrix[i][j];
			//but left and right side of matrix are connected
			
			// compartments are evaluated by looking at what they receive from their neighbours
			//update of Xs
			//xflow [0] is above, [1] is below, [2] is left, [3] is right
			std::vector<double> XsFlow(4);
			XsFlow[0] =	diffX *(concAbove[0] - matrix[i][j][0]); // receive diffusion speed times concentration difference
			XsFlow[1] =	diffX *(concBelow[0] - matrix[i][j][0]);
			XsFlow[2] =	diffX *(matrix[i][(j - 1 + iN) % iN][0]  // (j - 1 + iN) % iN tracks if we are at a border
				       			     - matrix[i][j][0]);
			XsFlow[3] =	diffX * (matrix[i][(j + 1) % iN][0] 
					                     - matrix[i][j][0]);

			matrixUpdate[i][j][0] += XsFlow[0] + XsFlow[1] 
					       + XsFlow[2] + XsFlow[3]; 
			//update of Xr
			// same as above
			//xflow [0] is above, [1] is below, [2] is left, [3] is right
			std::vector<double> XrFlow(4);
			XrFlow[0] =	diffX *(concAbove[1] - matrix[i][j][1]);
			XrFlow[1] =	diffX *(concBelow[1] - matrix[i][j][1]);
			XrFlow[2] =	diffX *(matrix[i][(j - 1 + iN) % iN][1] 
					                     - matrix[i][j][1]);
			XrFlow[3] =	diffX * (matrix[i][(j + 1) % iN][1] 
							     - matrix[i][j][1]);
			// add 
			matrixUpdate[i][j][1] += XrFlow[0] + XrFlow[1]  // inflow the compartment receives from every direction is added
				               + XrFlow[2] + XrFlow[3]; 
			//update of Ym
			matrixUpdate[i][j][3] += diffY * -(4 * matrix[i][j][3] 
				               - concAbove[3] - concBelow[3]
					       - matrix[i][(j - 1 + iN) % iN][3]
					       - matrix[i][(j + 1) % iN][3]);
			//update of Z
			matrixUpdate[i][j][2] += diffZ * -(4 * matrix[i][j][2] 
				               - concAbove[2] - concBelow[2]
					       - matrix[i][(j - 1 + iN) % iN][2]
					       - matrix[i][(j + 1) % iN][2]);
			//update of Ys
			//neglect ys of neighbour if flow is towards neighbour
			for (int k = 0; k < 4; ++k)
				if (XsFlow[k] < 0.0)
					XsFlow[k] = 0.0;
			if(matrixUpdate[i][j][0] > 1.0e-5) {
				matrixUpdate[i][j][4] =  // the new concentration is determined using weighted average
					( matrix[i][j][0] * matrix[i][j][4]  
					+ XsFlow[0] * concAbove[0] 
					+ XsFlow[1] * concBelow[0] 
					+ XsFlow[2] * matrix[i][(j - 1 + iN) % iN][4]
					+ XsFlow[3] * matrix[i][(j + 1) % iN][4] ) 
					/ (matrixUpdate[i][j][0]);
			} else{ // population is dead to us (below accepted level)
				matrixUpdate[i][j][4] = 0.0;
				matrixUpdate[i][j][0] = 0.0;
			}

			//update of Yr
			// same as above
			for (int k = 0; k < 4; ++k)
				if (XrFlow[k] < 0)
					XrFlow[k] = 0;
			if(matrixUpdate[i][j][1] > 1.0e-5) {
				matrixUpdate[i][j][5] = 
					( matrix[i][j][1] * matrix[i][j][5] 
					+ XrFlow[0] * concAbove[1]
					+ XrFlow[1] * concBelow[1] 
					+ XrFlow[2] * matrix[i][(j - 1 + iN) % iN][5]
					+ XrFlow[3] * matrix[i][(j + 1) % iN][5] ) 
					/ (matrixUpdate[i][j][1]);
			} else {
				matrixUpdate[i][j][5] = 0.0;
				matrixUpdate[i][j][1] = 0.0;
			}
		}
	}
	matrix = matrixUpdate;	//update original matrix
}

