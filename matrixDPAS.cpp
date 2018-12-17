#include <iostream>
#include <vector>
#include <cmath>
#include <exception>
#include <cstdlib>

#include "project.h"

// integration parameters
const double tolerance 	= 1.0e-7;	// accepted local error

// integration routine
const double kdShrinkMax = 0.1;		// max decrease in step size
const double kdGrowMax   = 5.0; 	// max increase in step size
const double kdSafety    = 0.9; 	// safety factor in adaptive step size control (?)
const double kdMinH      = 1.0e-6; 	// minimum step size
					// shared with extern in header
const double lowBound 	 = 1.0e-5;

// for optimizing the number of vector instantiations, these are taken out of function body 
// and shared in project.h 
// saves ~30% at runtime
std::vector<double> x(nvar);
std::vector<double> xtmp(nvar);  
std::vector<double> dxdtA(nvar);
std::vector<double> dxdtB(nvar);
std::vector<double> dxdtC(nvar);
std::vector<double> dxdtD(nvar);
std::vector<double> dxdtE(nvar);
std::vector<double> dxdtF(nvar);
std::vector<double> dxdtG(nvar);
std::vector<double> xHi(nvar);
std::vector<double> xLo(nvar);



bool matrixDPAS(double &t, // current time
	       	std::vector<std::vector<std::vector<double> > > &matX, // matrix to be updated
		double &h) // time-step size
// for every element in a matrix,
// computes the next step of the Dormand-Price routine for numerical integration by updating the given arguments
// returns false if the error found is too large, in which case it does not update t and x
// can be optimized by retaining one of the solutions to replace one of the following steps solutions
{
	// is a temporary matrix even needed?
	std::vector<std::vector<std::vector<double> > > matTmp = matX;
	double errMax = 0.0;
	int culprit = 0; // variable with biggest error
	int culpRow = 0; // location of biggest error
	int culpCol = 0;
	// double errMax = 0.0;
	double valLo = 0.0;
	double valHi = 0.0;
	for(int row = 0; row < iN; ++row) { // for every row
		for(int col = 0; col < iN; ++col) { // and every column
			x = matX[row][col]; // take values from given matrix
			for(int k = 0; k < 2; ++k)
				if(x[k] < lowBound)
					x[k] = 0.0;
	
			// step 1 a
			rhs(t,x,dxdtA); // determine derivatives for first step

			// step 2 b
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + 0.2 * h * dxdtA[i];
			rhs(t + 0.5 * h ,xtmp ,dxdtB); 

			// step 3 c
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + 0.3 * h * ( 3.0 * dxdtA[i] 
				               		   + 9.0 * dxdtB[i] ) 
									 / 40.0;
			rhs(t + 0.75 * h ,xtmp ,dxdtC);

			// step 4 d
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + 0.8 * h * ( 44.0 * dxdtA[i] 
							  - 138.0 * dxdtB[i] 
						          + 160.0 * dxdtC[i] ) 
					 				 / 45.0;
			rhs(t + h ,xtmp ,dxdtD);

			// step 5 e
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + 8 * h * ( 19372 * dxdtA[i] 
						         - 76080 * dxdtB[i]
							 + 64448 * dxdtC[i] 
							 - 1908 * dxdtD[i] )
							       / (9.0 * 6561.0);
			rhs(t + h ,xtmp ,dxdtE);

			// step 6 f
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + h  
				      *	( 9017.0 / 3168.0  * dxdtA[i]  
				         - 355.0 / 33.0    * dxdtB[i] 
				       + 46732.0 / 5247.0  * dxdtC[i]
				          + 49.0 / 176.0   * dxdtD[i]
					- 5103.0 / 18656.0 * dxdtE[i] );
			rhs(t + h ,xtmp ,dxdtF);
			
			// step 7 g
			for(int i = 0; i < nvar; ++i) 
				xtmp[i] = x[i] + h  
				             * ( 35.0 / 384.0 	* dxdtA[i] 
					      + 500.0 / 1113.0 	* dxdtC[i]
					      + 125.0 / 192.0 	* dxdtD[i]
					     - 2187.0 / 6784.0 	* dxdtE[i]
					       + 11.0 / 84.0 	* dxdtF[i] );
			rhs(t + h ,xtmp ,dxdtG);

			// compute x_lo: 5th order solution
			for(int i = 0; i < nvar; ++i) 
				xLo[i] = x[i] + h 
					  * ( 35.0 / 384.0 	* dxdtA[i] 
					   + 500.0 / 1113.0 	* dxdtC[i] 
					   + 125.0 / 192.0	* dxdtD[i] 
					  - 2187.0 / 6784.0	* dxdtE[i]
					    + 11.0 / 84.0	* dxdtF[i] );
			// compute x_hi: solution for error estimation
			for(int i = 0; i < nvar; ++i)
				xHi[i] = x[i] + h 
				       * ( 5179.0 / 57600.0	* dxdtA[i] 
					 + 7571.0 / 16695.0	* dxdtC[i] 
					  + 393.0 / 640.0	* dxdtD[i]
			                - 92097.0 / 339200.0	* dxdtE[i]
				          + 187.0 / 2100.0 	* dxdtF[i]
				            + 1.0 / 40.0 	* dxdtG[i] );

			//int culprit = 0;
			// double errMax = 0.0;
			//double valLo = 0.0;
			//double valHi = 0.0;
			for(int i = 0; i < nvar; ++i) {
				// error
				double erri = fabs(0.5 * h * (xLo[i] - xHi[i]))
					            		   / tolerance;
				if(erri > errMax){
					errMax = erri;
					culprit = i;
					valLo = xLo[i];
					valHi = xHi[i];
					culpRow = row;
					culpCol = col;
				}
			}
			// cut off dead populations
			if(xLo[0] < lowBound) {
				//xLo[4] = 0.0;
				xLo[0] = 0.0;
			}
			if(xLo[0] >= 1.0){
				std::cout << "asdsad\n";
				xLo[0] = 0.99;
			}
			if(xLo[1] < lowBound) {
				xLo[1] = 0.0;
				//xLo[5] = 0.0;
			}
			if(xLo[1] >= 1.0) {
				xLo[1] = 0.99;
				std::cout << "asdsad\n";
			}
			matTmp[row][col] = xLo;
		}
	}

		// adjust step size
	const double fct = errMax > 0.0 ? kdSafety / pow(errMax, 1.0 / 5.0) 
					: kdGrowMax; 
	if(errMax > 1.0) {
		// reduce size
		if(fct < kdShrinkMax)
			h *= kdShrinkMax;
		else
			h *= fct;
		if(h < kdMinH){
			std::cout << "culprit: " << culprit << ", error: " 
				     			  << errMax << std::endl
				<< "valLo: " << valLo << " valHi: " 
							<< valHi << std::endl
				<< "culpRow: " << culpRow << " culpCol: " 
							<< culpCol << std::endl;
			throw std::runtime_error(
				      "step size underflow in matrixDPAS().\n");
		}
		return false; // step rejected
		
	} else {
//		std::cout << "fct: " << fct << std::endl;
		// update solution
		matX = matTmp;
		t += h;
		// increase step size
		if(fct > kdGrowMax)
			h *= kdGrowMax;
		else 
			h *= fct;
		return true; // step accepted
	}
}
