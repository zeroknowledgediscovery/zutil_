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

void test(
	PFSA& G, 
	PFSA& H, 
	size_t length,
	size_t depth, 
	size_t num_runs)
{
	connx aut = G.get_aut();
	size_t alphabet = aut[0].size();


	map<unsigned int, double> DG = G.get_distr(depth);
	map<unsigned int, double> DH = H.get_distr(depth);

	double dist = SCC_UTIL__::get_dist(DG, DH, alphabet, depth);
	cout << depth << "-order distance =" << dist << endl;
	symbol_list_ sl = H.gen_data(length).get_symbol_list();
	
	map<unsigned int, double> DH_sl = SCC_UTIL__::get_distr_sl(sl, alphabet, depth);
	double dist_sl = SCC_UTIL__::get_dist(DG, DH_sl, alphabet, depth);
	cout << depth << "-order empirical distance =" << dist_sl << endl;

	pair<double, double> stat_pair = SCC_UTIL__::get_stat(G, length, depth, num_runs);

	double mean = stat_pair.first;
	double std = stat_pair.second;
	cout << "mean = " << mean << endl;
	cout << "std = " << std << endl;
	
	double zValue = (dist - mean) / std;
	double zValue_sl = (dist_sl - mean) / std;
	
	cout << depth << "-order zvalue =" << zValue << endl;
	cout << depth << "-order empirical zvalue =" << zValue_sl << endl;
}


void test_matrix(
	PFSA& G, 
	PFSA& H, 
	vector<double> scale, 
	size_t length, 
	size_t depth,
	size_t num_runs,
	string filename)
{
	connx aut = G.get_aut();
	size_t alphabet = aut[0].size();
	

	vector<PFSA> PFSA_vec;
	vector<double> mean_vec;
	vector<double> std_vec;
	vector<map<unsigned int, double>> dist_vec;
	for (double s : scale)
	{
		PFSA sG(G * s);	
		PFSA_vec.push_back(sG);
		pair<double, double> stat_pair = SCC_UTIL__::get_stat(sG, length, depth, num_runs);
		mean_vec.push_back(stat_pair.first);
		std_vec.push_back(stat_pair.second);
		dist_vec.push_back(sG.get_distr(depth));
	}

	size_t num_steps = scale.size();
	vector<double> dist_row(num_steps, 0);
	vector<vector<double>> dist_matrix(num_steps, dist_row);
	
	for (size_t i = 0; i < num_steps; i++)
	{
		for (size_t j = 0; j < num_steps; j++) 
		{
			symbol_list_ sl = PFSA_vec[i].gen_transduced_data(H, length);
			map<unsigned int, double> distr = SCC_UTIL__::get_distr_sl(sl, alphabet, depth);
			double dist = SCC_UTIL__::get_dist(dist_vec[j], distr, alphabet, depth);
			dist_matrix[i][j] = dist;
		}		
	}
	// output result.
	// The first row is the scale;
	// The second row is the mean;
	// The third row is the standard deviation;
	// The remaining is the llk matrix.
	ofstream output(filename);
	for (double s : scale)
	{
		output << s << " ";
	}
	output << endl;
	for(double mean : mean_vec)
	{
		output << mean << " ";
	}
	output << endl;
	for(double std : std_vec)
	{
		output << std << " ";
	}
	output << endl;
	for (vector<double> vec : dist_matrix)
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
	size_t depth = (size_t)stoi(argv[4]);
	size_t num_runs = (size_t)stoi(argv[5]);

	// cout << "\n############################### Test: ###########################" << endl;
	// test(G, H, length, depth, num_runs);
	// cout << "############################# Test END: #########################\n" << endl;
	
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
	string output = "llk_matrices/dist-" + get_name(PFSA_file) + "-" + get_name(XPFSA_file); 
	test_matrix(G, H, scale, length, depth, num_runs, output);
	cout << "############################# Test Matrix END: #########################\n" << endl;
	
  	return 0;
}
