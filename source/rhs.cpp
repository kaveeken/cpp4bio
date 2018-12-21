/*****************************************
              buildSystem.cpp

  Function to compute derivatives of a system of differential equations.
  Used in simulating collective abtibiotic resistance.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <ctgmath>

#include "project.h"

// model specification 
// maybe put these in a file
const double r 	 = 20.0;		//max growth rate
const double eta = 0.9;			//relative growth rate
const double c 	 = 1.0;			//resource consumation rate
const double kz  = 4.0;			//half saturation constant of the growth function
const double hy  = 0.3125;		//inhibitory constant of the growth function (div by y0)
const double p 	 = 50.0;		//pasive transport of Cm between cell compartments
const double d 	 = 37.5;		//maximum rate and...
const double ky  = 3.125;		//...half saturation constant...

// right-hand sides of a system of ODEs
void rhs(const double &t, 		// current time
	 const std::vector<double> &x, 	// current variables
	 std::vector<double> &dxdt) 	// (empty) vector that recieves values of derivatives at time t
{
	//calculate derivatives of Xs and Ys
	if (x[0] != 0.0) {
		dxdt[0] = r * x[2] / (kz + x[2]) * hy / (hy + x[4]) * x[0] 
							            - x[0];
		dxdt[4] = p * (x[3] - x[4]) /*- x[4] */
					    - x[4] * dxdt[0] / x[0];
	}
	else {		//neglect population if it is extinct
		dxdt[0] = 0.0;
		dxdt[4] = 0.0;
	}
	//calculate derivatives of Xr and Yr
	if (x[1] != 0.0) {
		dxdt[1] = eta * r * x[2] / (kz + x[2]) 
				  * hy / (hy + x[5]) * x[1] - x[1];
		dxdt[5] = p * (x[3] - x[5]) - d * x[5] / (ky + x[5]) 
					    /*- x[5] */ - x[5] * dxdt[1] / x[1];
	}
	else {		//neglect population if it is extinct
		dxdt[1] = 0.0;
		dxdt[5] = 0.0;
	}

	// derivative of Z
	dxdt[2] = /*(1.0 - x[2])*/ - c * r * x[2] / (kz + x[2])
		* (x[0] * hy / (hy + x[4]) + x[1] * hy / (hy + x[5]));  
	
	// derivative of Ym
	if(x[0] + x[1] < 1.0) {
		dxdt[3] = /*1.0 / (1.0 - x[0] - x[1])*/
			- p * (x[0] * (x[3] - x[4]) + x[1] * (x[3] - x[5]))
			/ (1.0 - x[0] - x[1])
			/*- x[3]*/ - x[3] * (dxdt[0] + dxdt[1]) / (x[0] + x[1] - 1.0);
	}
	else dxdt[3] = 0.0; // this was for debugging 
}
