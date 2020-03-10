#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

#include <dirent.h>

#include <gsl/gsl_blas.h>
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


void test(PFSA& G, PFSA& H, size_t length)
{
	symbol_list_ sl = G.gen_data(length).get_symbol_list();
	
	PFSA Comp(G||H);
	// Comp.describe("G||H");

	double llk_G = G.log_likelihood(sl);
	cout << "llk of G  = " << llk_G << endl;
	
	double llk_HMM = Comp.HMM_log_likelihood(sl);
	cout << "llk of H(G) = " << llk_HMM << endl;
}

void test_stat(PFSA& G, PFSA& H, size_t length, size_t num_runs)
{
	vector<double> llks;
	PFSA Comp(G||H);
	// Comp.describe("G||H");
	
	for (size_t n = 0; n < num_runs; n++)
	{
		symbol_list_ sl = G.gen_data(length).get_symbol_list();
		double llk = Comp.HMM_log_likelihood(sl);
		// cout << n << ": " << llk << endl;
		llks.push_back(llk);
	}

	double mean = gsl_stats_mean(&llks[0], 1, llks.size());
	double std = gsl_stats_sd(&llks[0], 1, llks.size());	

	cout << "mean = " << mean << endl;
	cout << "std = " << std << endl;
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
			symbol_list_ sl = PFSA_vec[j].gen_data(length).get_symbol_list();
			double llk = Comp.HMM_log_likelihood(sl);
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
	string PFSA_file = argv[1];
	PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
	G.describe("XPFSA G");
	
	string XPFSA_file = argv[2];	
	PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "XPFSA");
	H.describe("XPFSA H");

	size_t length = (size_t)stoi(argv[3]);
	size_t num_runs = (size_t)stoi(argv[4]);

	cout << "\n############################### Test: ###########################" << endl;
	test(G, H, length);
	cout << "############################# Test END: #########################\n" << endl;
	

	cout << "\n############################### Test STAT: ###########################" << endl;
	test_stat(G, H, length, num_runs);
	cout << "############################# Test STAT END: #########################\n" << endl;

	double min = -2;
	double max = 2;
	size_t num_steps = 100;
	
	double range = max - min;
	double step_size = range / num_steps;
 	vector<double> scale(num_steps + 1, 0.);
	
	for (size_t i = 0; i <= num_steps; i++)
	{
		scale[i] = min + i * step_size;
	}

	cout << "\n############################### Test Matrix: ###########################" << endl;
	string output = "llk_matrices/HHM_llk-" + get_name(PFSA_file) + "-" + get_name(XPFSA_file); 
	test_matrix(G, H, scale, length, output);
	cout << "############################# Test Matrix END: #########################\n" << endl;
	
  	return 0;
}
