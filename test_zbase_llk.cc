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
	const string version="Transducer Eigen-Study-llk v0 2020 zed.uchicago.edu";
	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
 
	string PFSA_filename;
	string XPFSA_filename;
	size_t length;
	double min;
	double max;
	size_t num_samples;
 	
	options_description desc(
"\n\
-----------------------### TransducerStudy-llk zed.uchicago.edu 2020 ###---------------------------------\n\
---------------------------------------------------------------------------------------------------------\n\
This is code is used to calculate the llk matrix so as to study the eigen structure of XPFSA transducer.\n\n\
The input includes:\n\
\t1) A PFSA G;\n\
\t2) An XPFSA transducer H;\n\
\t3) A sequence length, length;\n\
\t4) A minimum value(inclusive), min, for the range of scalar, with default -2;\n\
\t5) A maximum value(inclusive), max, for the range of scaler, with default 2;\n\
\t\t (max must be strictly greater than min);\n\
\t6) Number of sampling points, num_samples, between min and max;\n\
\t\t (num_samples cannot be less than 2 with default value 2).\n\n\
The output is a file that has num_samples + 2 lines with:\n\
\ta) the first line being the sample points;\n\
\tb) the second line being the entropy rate of x * G where x is a sample point;\n\
\tc) the num_samples by num_samples matrix of llk of x * G generating H(y * G).\n\n\
Example (in zutil_ folder):\n\
./bin_practice/test_zbase_llk -p cfgfiles/PFSA_M2_0.cfg -x cfgfiles/XPFSA_2_3_2.cfg -l 1000 -n 100\n\
---------------------------------------------------------------------------------------------------------\n"
);
  	positional_options_description p;
	desc.add_options()
    ("help,h", "Print help message.")
    ("version,V", "Print version number.")
    ("pfsa,p", value<string>(&PFSA_filename), "PFSA filename.")
    ("xpfsa,x", value<string>(&XPFSA_filename), "XPFSA filename.")
    ("length,l", value<size_t>(&length)->default_value(1000), "Length of the sequence.")
    ("min,a", value<double>(&min)->default_value(-2), "Minimum value of the range.")
    ("max,b", value<double>(&max)->default_value(2), "Maximum value of the range. Must be strictly bigger than min.")
    ("num_samples,n", value<size_t>(&num_samples)->default_value(2), "Number of sample points.\
 The first sample point is min while the last sample point is max.");
  	variables_map vm;
  	if (argc == 1)
	{
    	cout <<"empty arg, type -h or --help" << endl;
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

	// Other parameter errors
	if (vm.count("num_samples"))
	{
		if (num_samples < 2)
		{
			cout << "num_samples must be at least 2!" << endl;
			return 1;
		}
	}
	if (max <= min)
	{
		cout << "max(b) must be strictly bigger than min(a)!" << endl;
		return 1;
	}

	PFSA G = SCC_UTIL__::read_mc(PFSA_filename, "PFSA");
	PFSA H = SCC_UTIL__::read_mc(XPFSA_filename, "XPFSA");
	H.describe("transducer");

	double range = max - min;
	double step_size = range / (num_samples - 1);
 	vector<double> scale(num_samples, 0.);
	
	for (size_t i = 0; i < num_samples; i++)
	{
		scale[i] = min + i * step_size;
	}
	
	string output = "llk_matrices/llk-" 
		+ get_name(PFSA_filename)
		+ "-" + get_name(XPFSA_filename)
		+ "-" + to_string(length)
		+ "-" + to_string(num_samples); 
	test_matrix(G, H, scale, length, output);
	
  	return 0;
}
