/*!
 This is code is used to calculate the llk matrix so as to 
 study the eigen structure of XPFSA transducer.
 With input PFSA G and XPFSA transducer H, we generated a linear
 space of sample points between \c min, \c max, with number of 
 steps equalling \c num_steps. We then calculate x * G for each 
 sample points and generate a sequence from x * G and transduce
 the sequence with H. We then calculate the llk of y * G generating
 the output for all y in the linear sample space.

 The code output two files one save the llk matrix, and the other,
 the entropy rates of the x * G. The saved results are scaled since
 the first row of the output is the linear space of sample points.  
 */

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

void test_matrix(
	PFSA& G, 
	PFSA& H, 
	vector<double> scale,
	size_t length, 
	string filename)
{
	size_t num_steps = scale.size();

	vector<PFSA> PFSA_vec;
	vector<double> er_vec;
	for (double s : scale)
	{
		PFSA sG(G * s);	
		PFSA_vec.push_back(sG);
		er_vec.push_back(sG.entropy());
	}

	vector<double> llk_row(num_steps, 0);
	vector<vector<double>> llk_matrix(num_steps, llk_row);
	
	for (size_t i = 0; i < num_steps; i++)
	{
		PFSA Comp(PFSA_vec[i]||H);
		for (size_t j = 0; j < num_steps; j++) 
		{
			symbol_list_ sl = PFSA_vec[i].gen_transduced_data(H, length);
			double llk = PFSA_vec[j].log_likelihood(sl);
			llk_matrix[i][j] = llk;
		}		
	}
	// output result.
	// The first row is the scale.
	// The second row are the entropy rates.
	// The remaining is the llk matrix.
	ofstream output(filename);
	for (double s : scale)
	{
		output << s << " ";
	}
	output << endl;
	for(double er : er_vec)
	{
		output << er << " ";
	}
	output << endl;
	for (vector<double> vec : llk_matrix)
	{
		for(double llk : vec)
		{
			output << llk << " ";
		}
		output << endl;
	}
}


int main(int argc, char** argv)
{
//   	const string version="TransducerStudy-llk v0 2020 zed.uchicago.edu";
//   	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
//   
// 	string PFSA_filename;
// 	string XPFSA_filename;
// 	size_t length;	// sequence length
// 	size_t num_samples;
// 	string random_parameters;
//  	
// 	options_description desc("### TransducerStudy-llk zed.uchicago.edu 2020 ###\n\
// --------------------------\n\
// \n\
// \n\
// \n\
// \n\
// --------------------------\n\
// Example (in zutil_ folder): ./bin_practice/test_zbase_llk\
//  -p cfgfiles/PFSA_M2_0.cfg\
//  -x cfgfiles/XPFSA_2_3_2.cfg\
//  -l 1000 -n 1");
//   	desc.add_options()
//     ("help,h", "print help message.")
//     ("version,V", "print version number")
//     ("pfsa,p", value<string>(&PFSA_filename)->("random"), "PFSA file")
//     ("xpfsa,x", value<string>(&XPFSA_filename)->("random"), "XPFSA file")
//     ("length,l", value<size_t>(&length)->default_value(1000), "Length of the sequence")
//     ("random_parameters, r", value<string>(&random_parameters), "Parameters for randomly \
// generated PFSA and XPFSA. Input in the form \"input_num_states input_alphabet output_num_states output_alphabet\"")
//     ("num_samples,n", value<size_t>(&num_samples)->default_value(1), "If number of samples\
// is 1, then only used the input PFSA or the randomly generated PFSA, otherwise \
// num_samples many PFSA will ge generated using the aut of the input PFSA or \
// randomly generated PFSA.");
//   positional_options_description p;
//   variables_map vm;
//   if (argc == 1)
//     cout <<"empty arg, type -h or --help" << endl;
//   try
//     {
//       store(command_line_parser(argc, argv)
// 	    .options(desc)
// 	    .run(), vm);
//       notify(vm);
//     } 
//   catch (std::exception &e)
//     {
//       cout << endl << e.what() 
// 	   << endl << desc << endl;
//     }
//   if (vm.count("help"))
//     {
//       cout << desc << endl;
//       return 1;
//     }
//   if (vm.count("version"))
//     {
//       cout << version << endl; 
//       return 1;
//     }
// 
//   if (DATA_DIR=="row")
//     DATA_DIR="across";
// 
//   if (DATA_DIR=="column")
//     DATA_DIR="up";
	
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
