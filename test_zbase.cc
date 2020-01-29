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

double compare(
	PFSA& G,
	const pair<double, double>& stat,
	const symbol_list_& sl, 
	unsigned int depth=5)
{
	int alphabet = G.get_aut()[0].size();
	map<unsigned int, double> distr = G.get_distr(depth);
	map<unsigned int, double> distr_sl = SCC_UTIL__::get_distr_sl(sl, alphabet, depth);
	double dist = SCC_UTIL__::get_dist(distr, distr_sl, alphabet, depth);
	
	double mean = stat.first;
	double std = stat.second;

	double zValue = (dist - mean) / std;
	
	return zValue;
}


int main(int argc, char** argv)
{ 
	string PFSA_file = argv[1];
	string XPFSA_file = argv[2];	
	int length = stoi(argv[3]);		// sequence length
	int depth = stoi(argv[4]);		// depth for the norm
	int num_runs = stoi(argv[5]);	// number of runs to the mean and standard deviation

	
	PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
	G.describe("G");
	PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "XPFSA");
	H.describe("H");

	
	Symbolic_string_ input_ = G.gen_data(length);
	Symbolic_string_ output_ = input_ * H;
	symbol_list_ input = input_.get_symbol_list();
	symbol_list_ output = output_.get_symbol_list();
	cout << "Input:\n" << input_ << endl;
	cout << "Output:\n" << output_ << endl;
	
	cout << "entropy rate of G = " << G.entropy() << endl;
	cout << "llk G generating input = " << G.log_likelihood(input) << endl;		
	cout << "llk G generating output = " << G.log_likelihood(output) << endl;		
	 
	pair<double, double> stat = SCC_UTIL__::get_stat(G, length, depth, num_runs);
	cout << "mean = " << stat.first << endl;
	cout << "standard deviation = " << stat.second << endl;
	cout << "z-value of input = " << compare(G, stat, input, depth)  << endl;
	cout << "z-value of output = " << compare(G, stat, output, depth)  << endl;
	
	pair<PFSA, vector<double>> sync_pair = G.synchronous_product(H);
	PFSA Prod = sync_pair.first;
	vector<double> proj_sd = sync_pair.second;

	cout << "Synchronous product:" << endl;
	Prod.mc_print();
	cout << "projective distribution = " << proj_sd << endl;

	PFSA Comp(~Prod);
	Comp.mc_print();
	
	vector<vector<double>> Omega; 
	cout << "Omega of H with FWN Input:" << endl;
	Omega = H.get_Omega();
	for (vector<double> OmegaRow : Omega)
	{
		cout << "\t" << OmegaRow;
	}
	
	cout << "Omega of H with Input G:" << endl;
	Omega = H.get_Omega(G);
	for (vector<double> OmegaRow : Omega)
	{
		cout << "\t" << OmegaRow;
	}
  	return 0;
}
