#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

#include <dirent.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;


string get_name(string path)
{
	size_t start = path.find("PFSA");
	size_t end = path.find(".cfg");

	if ((start != 0) && (path[start - 1] == 'X'))
	{
		start = start - 1;
	}
	size_t len = end - start;

	return path.substr(start, len);
}


int main(int argc, char** argv)
{
  	const string version="TransducerStudy-llk v0 2020 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
  
	string PFSA_filename;
	string XPFSA_filename;
	size_t length;	// sequence length
	size_t num_samples;
	string random_pfsa_parameters;
	string random_xpfsa_parameters;
	string result_folder;


	double min;
	double max;
	size_t num_steps;
	// step_size is calculated as (max - min) / num_steps, however since max is also included	

	double step_size; 
	// Inferred the parameters
	PFSA G;
	PFSA H; 
	size_t num_states;
	size_t alphabet;
	size_t x_num_states;
	size_t x_alphabet;
	connx aut;
	pitilde pit;
	connx x_aut;
	pitilde Xpit;
 	
	options_description desc("### TransducerStudy-llk zed.uchicago.edu 2020 ###\n\
--------------------------\n\
For a given or randomly generated XPFSA transducer, and a given or randomly \n\
generated PFSA, calculate \n\
\n\
\n\
--------------------------\n\
Example (in zutil_ folder): ./bin_practice/test_zbase_llk\
 -p cfgfiles/PFSA_M2_0.cfg\
 -x cfgfiles/XPFSA_2_3_2.cfg\
 -l 1000 -n 1");
  	desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("pfsa,p", value<string>(&PFSA_filename), "PFSA file")
    ("random_pfsa_parameters,q", value<string>(&random_pfsa_parameters), "Parameters for randomly \
generated PFSA. Input in the form \"input_num_states input_alphabet\"")
    ("xpfsa,x", value<string>(&XPFSA_filename), "XPFSA file")
    ("random_xpfsa_parameters,y", value<string>(&random_xpfsa_parameters), "Parameters for randomly \
generated XPFSA. Input in the form \"output_num_states output_alphabet\"")
    ("length,l", value<size_t>(&length)->default_value(1000), "Length of the sequence")
    ("num_samples,n", value<size_t>(&num_samples)->default_value(1), "If number of samples\
is 1, then only used the input PFSA or the randomly generated PFSA, otherwise \
num_samples many PFSA will be generated using the aut of the input PFSA or the\
randomly generated PFSA.")
	("result_folder,r", value<string>(&result_folder), "If no result folder specified\
a folder ./TS_result_yy:mm:dd-HH:MM:SS will be used.");
  	positional_options_description p;
  	variables_map vm;
  	if (argc == 1)
	{
    	cout <<"empty arg, type -h or --help for manual" << endl;
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

  	if ((vm.count("pfsa") == 0) && (vm.count("random_pfsa_parameters") == 0))
	{
    	cout << "Either specify a PFSA file or specify parameters for randomly generated PFSA"<< endl;
		return 1;
	}
	else if ((vm.count("pfsa") == 1))
	{
		if (vm.count("random_pfsa_parameters") == 1)
		{
			cout << "parameter pfsa(p) will override parameter random_pfsa_parameters(q)";
		}
		G = SCC_UTIL__::read_mc(PFSA_filename, "PFSA");
	
		aut = G.get_aut();
		num_states = aut.size();
		alphabet = aut[0].size();
	}
	else
	{
		stringstream ss;
		ss << random_pfsa_parameters;
		string a;
		ss >> a;
		num_states = (size_t)stoi(a);
		ss >> a;
		alphabet = (size_t)stoi(a);
		
		aut = SCC_UTIL__::get_random_aut(num_states, alphabet);
		pit = SCC_UTIL__::get_random_pit(aut);
		G = PFSA(pit, aut);
	}
  	
	if ((vm.count("xpfsa") == 0) && (vm.count("random_xpfsa_parameters") == 0))
	{
    	cout << "Either specify a XPFSA file or specify parameters for randomly generated PFSA" << endl;
		return 0;
	}
	else if ((vm.count("xpfsa") == 1))
	{
		if (vm.count("random_xpfsa_parameters") == 1)
		{
			cout << "parameter xpfsa(x) will override parameter random_xpfsa_parameters(y)";
		}
		H = SCC_UTIL__::read_mc(XPFSA_filename, "XPFSA");
		
		x_aut = H.get_aut();
		x_num_states = x_aut.size();
		x_alphabet = x_aut[0].size();
	}
	else
	{
		stringstream ss;
		ss << random_xpfsa_parameters;
		string a;
		ss >> a;
		x_num_states = (size_t)stoi(a);
		ss >> a;
		x_alphabet = (size_t)stoi(a);
		
		x_aut = SCC_UTIL__::get_random_aut(x_num_states, alphabet);
		pitilde x_pit = SCC_UTIL__::get_uniform(x_num_states, alphabet);
		Xpit = SCC_UTIL__::get_random_pit(x_num_states, x_alphabet);
		
		H = PFSA(x_pit, x_aut); 
		H.set_Xpit(Xpit);
	}

	G.describe("PFSA G");
	H.describe("Transducer");

	
	// string PFSA_file = argv[1];
	// PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");

	// string XPFSA_file = argv[2];	
	// PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "XPFSA");

	// 
	// double min = -2;
	// double max = 2;
	// size_t num_steps = 100;
	// 
	// double range = max - min;
	// double step_size = range / num_steps;
 	// vector<double> scale(num_steps + 1, 0.);
	// 
	// for (size_t i = 0; i <= num_steps; i++)
	// {
	// 	scale[i] = min + i * step_size;
	// }
	// 
	// cout << "\n############################### Test Matrix: ###########################" << endl;
	// string output = "llk_matrices/llk-" + get_name(PFSA_file) + "-" + get_name(XPFSA_file); 
	// test_matrix(G, H, scale, length, output);
	// cout << "############################# Test Matrix END: #########################\n" << endl;
	
  	return 0;
}
