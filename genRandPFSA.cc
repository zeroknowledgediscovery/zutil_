#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

#include <dirent.h>
#include <time.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"
#include <boost/filesystem.hpp>

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

void save(PFSA& G, string filename, string type)
{
	connx aut = G.get_aut();
	pitilde pit;
	if (type.compare("PFSA") == 0)
	{
		pit = G.get_pit();	
	}
	else if (type.compare("XPFSA") == 0)
	{
		pit = G.get_Xpit();
	}
	else
	{
		cout << "Incorrect save type!" << endl;
		exit(1);
	}

	ofstream output(filename);
	output << "#CONNX" << endl;
	for (size_t i = 0; i < aut.size(); i++)
	{
		for (symbol s(0); s < aut[i].size(); s++)
		{
			output << aut[i][s];
			(s < aut[i].size() - 1) ? output << " " : output << endl;
		}
	}
	
	output << "\n#PITILDE" << endl;
	for (size_t i = 0; i < pit.size(); i++)
	{
		for (size_t j = 0; j < pit[i].size(); j++)
		{
			output << pit[i][j];
			(j < pit[i].size() - 1) ? output << " " : output << endl;
		}
	}
	output.close();
}

int main(int argc, char** argv)
{	
	size_t alphabet;
	size_t num_states;
	string PFSA_fname;
	bool verbose;

	const string version = "Generate random PFSA";
	const string EMPTY_ARG_MESSAGE = "Exiting. Type -h or --help for usage";
	
	options_description desc("");
  	desc.add_options()
    ("help,h", "Print help message.")
    ("version,V", "Print version number.")
    ("alphabet,a", value<size_t>(&alphabet)->default_value(2), "alphabet size.")
    ("mum_states,n", value<size_t>(&num_states)->default_value(2), "number of PFSA states.")
    ("verbose,v", value<bool>(&verbose)->default_value(false), "1 for print PFSA.")
	("PFSA_fname,p", value<string>(&PFSA_fname)->default_value("PFSA.cfg"), "PFSA filename.");
  	positional_options_description p;
  	variables_map vm;
  	if (argc == 1)
	{
    	cout << EMPTY_ARG_MESSAGE << endl;
	}
  	try
    {
		store(command_line_parser(argc, argv).options(desc).run(), vm);
      	notify(vm);
    } 
  	catch (std::exception &e)
    {
      	cout << endl << e.what() << endl << desc << endl;
    }
  	if (vm.count("help"))
    {
      	cout << desc << endl;
      	return 1;
    }
  	if (vm.count("version"))
    {
      	cout << version << endl; 
      	return 1;
    }

	connx aut = SCC_UTIL__::get_random_aut(num_states, alphabet);
	pitilde pit = SCC_UTIL__::get_random_pit(aut);
	PFSA G = PFSA(pit, aut);
	if (verbose)
		G.describe("PFSA");
	save(G, PFSA_fname, "PFSA");

	return 0;
}
