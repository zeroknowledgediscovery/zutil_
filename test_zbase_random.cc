/*
 * Test the new function in the Yi branch of zbase 
 * on generating random PFSA/XPFSA.
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>
#include <list> 
#include <stack> 

#include <dirent.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;


int main(int argc, char** argv)
{ 
	size_t num_states = stoi(argv[1]);
	size_t alphabet = stoi(argv[2]);	
	bool pos = stoi(argv[3]);

	bool POS = (pos > 0)? true : false;
	double rate = 0.1;
	unsigned int res = 10;

	connx aut = SCC_UTIL__::get_random_aut(num_states, alphabet, POS, rate, true);
	pitilde pit = SCC_UTIL__::get_random_pit(aut, res);
	pitilde Xpit = SCC_UTIL__::get_uniform(num_states, alphabet + 1);
	
	PFSA G(pit, aut);
	G.set_Xpit(Xpit);
	G.mc_print();
	PFSA H(~G);
	
	cout << num_states << ", " << H.get_aut().size() << endl;
	
	return 0;
}
