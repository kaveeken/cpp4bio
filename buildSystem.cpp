#include <vector>

#include "project.h"

std::vector<std::vector<std::vector<double> > > buildSystem(const int &which)
{
	// "default" cell conditions
	std::vector<double> vecDef(nvar);
	vecDef[0] = 0.05;
	vecDef[1] = 0.0064;
	vecDef[2] = 0.95;
	vecDef[3] = 0.95;
	vecDef[4] = 0.90;
	vecDef[5] = 0.8;
	// "empty" conditions w/o bacteria
	std::vector<double> vecEmp(nvar);
	vecEmp[0] = 0.0;
	vecEmp[1] = 0.0; 
	vecEmp[2] = 0.95;
	vecEmp[3] = 0.95;
	vecEmp[4] = 0.90;
	vecEmp[5] = 0.8;
	// many resistant bacteria
	std::vector<double> vecRes(nvar);
	vecRes[0] = 0.0;
	vecRes[1] = 0.4; 
	vecRes[2] = 0.95;
	vecRes[3] = 0.95;
	vecRes[4] = 0.90;
	vecRes[5] = 0.8;


