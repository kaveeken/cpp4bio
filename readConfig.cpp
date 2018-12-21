/*****************************************
              readConfig.cpp

  Function that reads values from a file.
  Used in simulating collective abtibiotic resistance.

  Auke van der Meij
  Kris Veeken

  Project C++ for biologists

  21/12/2018

  ***************************************/

#include <fstream>
#include <exception>
#include <cstdlib>
#include <vector>

#include "project.h"

/* this function reads a file named config.txt and uses it to instantiate variables of interest
 *
 * the configuration file should contain line separated values (doubles or integers).
 * space separated values should work as well, but any other character breaks it.
 * maybe incorporate a sed command in bash jobscript to convert from a more readable file format
 *
 * maybe find a way to const these read variables
 *
 * current sequence is: sane(?) values:
 * iN			something like 10 - 100 depending on diffusion speeds
 * tEnd 		around 200? sometimes more
 * dtSav		relative to tEnd. at current, main() supports up to 22 writes			 
 * diffXpH		all simulations are ran at 0.1
 * diffYpH		we used 0.5-1.5
 * diffZpH		we used 2-3 i think
 * "cluster" size 	1 - ? probably relative to iN
 */

void readConfig()
{
	// open filestream
	std::ifstream ifs("config.txt");
	if(!ifs.is_open()) 
		throw std::runtime_error("unable to open config.txt\n");
	
	// read values
	int iLen = 7;
	std::vector<double> vecConf(iLen);
	for(int i = 0; i < iLen;) {
		double x; // >> vecConf[i] doesnt work for some reason
		ifs >> x;
		if(!ifs.fail()) {
			vecConf[i] = x;
			++i;
		}
	}

	if(vecConf[0] < 2)
		throw std::runtime_error("iN (matrix size) too small\n");
	if(vecConf[6] > vecConf[0]) 
		throw std::runtime_error("iPopCol out of bounds for iN\nConsider larger iN or smaller iPopCol\n");
	if(vecConf[1] >= vecConf[2] * 22.0)
		throw std::runtime_error("Not enough writes for these tEnd and dtSav.\nMake sure tEnd < 22 * dtSav\n");
	
	// instantiate variables
	iN  = int(vecConf[0]); // matrix size
	tEnd    = vecConf[1];  // simulation length
	dtSav   = vecConf[2];  // interval at which is written to file
	diffXpH = vecConf[3];  // diffusion speed of bacteria
	diffYpH = vecConf[4];  // diffusion speed of Cm
	diffZpH = vecConf[5];  // diffusion speed of resource
	iPopCol = int(vecConf[6]);  // number of populated columns
}
