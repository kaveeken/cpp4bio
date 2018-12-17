#ifndef HEADER_H
#define HEADER_H

#include <vector>

//******** shared variables ********
// system size
extern int iN; 		// assumes square system
const extern int nvar;
extern double tEnd;
extern double dtSav;
// model parameters
const extern double r;		//max growth rate
const extern double eta;	//relative growth rate
const extern double c;		//resource consumation rate
const extern double kz;		//half saturation constant of the growth function
const extern double hy;		//inhibitory constant of the growth function (div by y0)
const extern double p;		//pasive transport of Cm between cell compartments
const extern double d;		//maximum rate and...
const extern double ky;		//...half saturation constant...
// integration parameters
const extern double tolerance;	
const extern double kdShrinkMax;
const extern double kdGrowMax; 
const extern double kdSafety; 
const extern double kdMinH;    
const extern double lowBound; 	
const extern double kdMinH; 
// diffusion
extern double diffXpH;		//diffuse speed of bacteria per time
extern double diffYpH;		//diffuse speed of Cm
extern double diffZpH;		//diffuse speed of Cm

// matrixDPAS vector instantiations
extern std::vector<double> x;
extern std::vector<double> xtmp;  
extern std::vector<double> dxdtA;
extern std::vector<double> dxdtB;
extern std::vector<double> dxdtC;
extern std::vector<double> dxdtD;
extern std::vector<double> dxdtE;
extern std::vector<double> dxdtF;
extern std::vector<double> dxdtE;
extern std::vector<double> xHi;
extern std::vector<double> xLo;

//************ shared functions **********
void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt);

bool matrixDPAS(double &t, std::vector<std::vector<std::vector<double> > > &matX, double &h);

void diffundMatrix(std::vector<std::vector<std::vector<double> > > &matrix, const double &h);

void inflow(std::vector<std::vector<std::vector<double> > > &matrix);

void readConfig();

#endif
