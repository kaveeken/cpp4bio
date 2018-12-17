#include  <iostream>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <vector>

#include "project.h"

/* the configuration file should contain line separated values (doubles or integers).
 * any other text breaks it
 * maybe incorporate a sed command in bash jobscript to convert from a more readable file format
 *
 * current sequence is:
 * iN
 * tEnd
 * dtSav
 * diffXpH
 * diffYpH
 * diffZpH
 */

void readConfig()
{
	// open filestream
	std::ifstream ifs("config.txt");
	if(!ifs.is_open()) 
		throw std::runtime_error("unable to open config.txt");
	
	// read values
	int iLen = 6;
	std::vector<double> vecConf(iLen);
	for(int i = 0; i < iLen;) {
		double x; // >> vecConf[i] doesnt work for some reason
		ifs >> x;
		if(!ifs.fail()) {
			vecConf[i] = x;
			++i;
		}
	}
	
	// instantiate variables
	iN  = int(vecConf[0]); // matrix size
	tEnd    = vecConf[1];  // simulation length
	dtSav   = vecConf[2];  // interval at which is written to file
	diffXpH = vecConf[3];  // diffusion speed of bacteria
	diffYpH = vecConf[4];  // diffusion speed of Cm
	diffZpH = vecConf[5];  // diffusion speed of resource
}

/*
int main()
{
	try {
		readConfig();
	}
	catch(std::exception &error) {
		std::cerr << "error: " << error.what();
		exit(EXIT_FAILURE);
	}

	return 0;
}
*/
