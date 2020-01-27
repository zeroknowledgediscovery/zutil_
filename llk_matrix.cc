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

const double EPS = 1e-8;


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
	string PFSA_file = argv[1];
	PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
	
	string XPFSA_file = argv[2];	
	PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "XPFSA");

	int length = stoi(argv[3]);	// sequence length
	
	double min = -2;
	double max = 2;
	double range = max - min;
	double num_steps = 100;
	double step_size = range / num_steps;
 
	vector<PFSA> PFSA_vec;
	vector<double> er_vec;
	vector<double> scale;
	for (size_t i = 0; i <=num_steps; i++)
	{
		double x = min + i * step_size;
		PFSA xG(G * x);
		
		PFSA_vec.push_back(xG);
		er_vec.push_back(xG.entropy());
		scale.push_back(x);
	}
	
	vector<double> llk_row(num_steps + 1, 0);
	vector<vector<double>> llk_matrix(num_steps + 1, llk_row);
	

	for (size_t i = 0; i <=num_steps; i++)
	{
		for (size_t j = 0; j <= num_steps; j++) 
		{
			symbol_list_ output = PFSA_vec[i].gen_transduced_data(H, length);

			double llk = PFSA_vec[j].log_likelihood(output);		
			llk_matrix[i][j] = llk;
		}
	}

	
	string PFSA_name = get_name(PFSA_file);
	string XPFSA_name = get_name(XPFSA_file);

	string llk_matrix_filename = "llk_matrices/llk_matrix-scaled-" + PFSA_name + "-" + XPFSA_name;
	ofstream llk_matrix_file(llk_matrix_filename);
	for (double s : scale)
	{
		llk_matrix_file << s << " ";
	}
	llk_matrix_file << endl;
	for (vector<double> vec : llk_matrix)
	{
		for(double llk : vec)
		{
			llk_matrix_file << llk << " ";
		}
		llk_matrix_file << endl;
	}

	string er_filename = "llk_matrices/er-scaled-" + PFSA_name;
	ofstream er_file(er_filename);
	for (double s : scale)
	{
		er_file << s << " ";
	}
	er_file << endl;
	for(double er : er_vec)
	{
		er_file << er << " ";
	}
  	return 0;
}
