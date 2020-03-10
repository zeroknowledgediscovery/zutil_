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


void get_llk_matrix(
	PFSA& G, 				// input PFSA
	PFSA& T, 				// XPFSA transducer
	vector<double> scale, 	// scale for x in xG
	size_t length,			// sequence length
	string filename)		// filename to save the result
{
	vector<PFSA> PFSA_vec;
	vector<double> er_vec;
	size_t num_steps = scale.size();
	for (double s : scale)
	{
		PFSA sG(G * s);
		PFSA_vec.push_back(sG);
		er_vec.push_back(sG.entropy());
	}
	
	// Compute the llk matrix
	vector<double> llk_row(num_steps + 1, 0);
	vector<vector<double>> llk_matrix(num_steps + 1, llk_row);
	
	for (size_t i = 0; i <=num_steps; i++)
	{
		for (size_t j = 0; j <= num_steps; j++) 
		{
			symbol_list_ output = PFSA_vec[i].gen_transduced_data(T, length);
			double llk = PFSA_vec[j].log_likelihood(output);		
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
	// int length = stoi(argv[1]);
	// double min = stod(argv[2]);
	// double max = stod(argv[3]);
	// double num_steps = stoi(argv[4]);
	// string filename = argv[5];
	// 
	// double range = max - min;
	// double step_size = range / num_steps;

 	// Get random transducer
 	size_t num_statesT = 4;
	size_t alphabet = 2;
 	connx autT = SCC_UTIL__::get_random_aut(num_statesT, alphabet, true);
	pitilde pitT= SCC_UTIL__::get_uniform(num_statesT, alphabet);
	pitilde XpitT = SCC_UTIL__::get_random_pit(num_statesT, alphabet, true);
	PFSA T(pitT, autT);
	T.set_Xpit(XpitT);
	T.mc_print();

	vector<vector<double>> Omega;

	Omega = T.get_Omega();
	for (vector<double> omega : Omega)
	{
		cout << "\t" << omega << endl;
	}


	// Generate random M2
	connx aut;
	aut[0][symbol(0)] = 0;
	aut[0][symbol(1)] = 1;
	aut[1][symbol(0)] = 0;
	aut[1][symbol(1)] = 1;
	pitilde pit = SCC_UTIL__::get_random_pit(aut);

	PFSA G(pit, aut);
	G.mc_print();
	
	Omega = T.get_Omega(G);
	for (vector<double> omega : Omega)
	{
		cout << "\t" << omega << endl;
	}

  	return 0;
}
