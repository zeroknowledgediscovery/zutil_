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

symbol_list_ transduce(PFSA& G, PFSA& H, int length)
{
	Symbolic_string_ input = G.gen_data(length);
	symbol_list_ input_ = input.get_symbol_list();

	Symbolic_string_ output = input * H;
	symbol_list_ output_ = output.get_symbol_list();
	
	return output_;
}


double entropy(PFSA& G)
{
	vector<double> stationary_dist(G.get_Stationary());
	pitilde pit(G.get_pit());
	
	double ent = 0;
	for (size_t s=0; s < pit.size(); s++)
	{
		ent += stationary_dist[s] * dist_entropy(pit[s]);
	}	
	return ent;
}


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


pair<PFSA, vector<double>> synchronous_product(PFSA& G, PFSA& H)
{
	connx autG = G.get_aut();
	connx autH = H.get_aut();
	pitilde pitG = G.get_pit();
	pitilde XpitH = H.get_Xpit();

	unsigned int alphabet = autG[0].size();
	
	map<unsigned int, map<unsigned int, unsigned int>> state_index;
	unsigned int count = 0;
	for (unsigned int i = 0; i < autG.size(); i++)
	{
		for (unsigned int j = 0; j < autH.size(); j++)
		{
			state_index[i][j] = count++;
		}
	}
	
	connx autProd;
	pitilde pitProd;
	pitilde XpitProd;

	for (unsigned int i = 0; i < autG.size(); i++)
	{
		for (unsigned int j = 0; j < autH.size(); j++)
		{
			unsigned int from_state = state_index[i][j];
			for (symbol symb(0); symb < alphabet; symb++)
			{
				int k = autG[i][symb];
				int l = autH[j][symb];
				int to_state;
				if ((k >= 0) and (l >= 0))
				{
					to_state = state_index[k][l];
				}
				else if ((k < 0) and (l < 0))
				{
					to_state = -1;
				}
				else
				{
					// when one of k or l is -1, the composition is impossible.
					// If composition is impossible, return a trivial pair.
					cout << "Composition is impossible. Return a trivial pair" << endl;
					return make_pair(PFSA (), vector<double>()); 
				}
				autProd[from_state][symb] = to_state;
			}
			pitProd[from_state] = vector<double>(pitG[i]);
			XpitProd[from_state] = vector<double>(XpitH[j]);
		}
	}
	PFSA Prod(pitProd, autProd);
	Prod.set_Xpit(XpitProd);

	vector<double> sd(Prod.get_Stationary());
	vector<double> proj_sd(XpitH.size(), 0.);

	for (unsigned int i = 0; i < autG.size(); i++)
	{
		for (unsigned int j = 0; j < autH.size(); j++)
		{
			proj_sd[j] += sd[state_index[i][j]];
		}
	}
	
	return make_pair(Prod, proj_sd);	
}


PFSA get_FWN(
	size_t num_states=1, 
	size_t alphabet=2)
{
	connx aut;
	pitilde pit;
	pitilde Xpit;
	
	vector<double> pitrow(alphabet, 1. / (double)alphabet);

	for (size_t i = 0; i < num_states; i++)
	{
		for (symbol symb(0); symb < alphabet; symb++)
		{
			aut[i][symb] = 0;
			pit[i] = pitrow;
			Xpit[i] = pitrow; 
		}
	}

	PFSA FWN(pit, aut);
	FWN.set_Xpit(Xpit);

	return FWN;
}

vector<vector<double>> get_Omega(PFSA& H, PFSA G=PFSA())
{
	connx aut = H.get_aut();
	pitilde Xpit = H.get_Xpit();

	size_t num_states = aut.size();
	size_t alphabet_i = aut[0].size();	// input alphabet size;
	size_t alphabet_o = Xpit[0].size();	// output alphabet size;
	
	// If no input PFSA specified, use a flat white noise;
	if (G.empty())
	{
		G = get_FWN(1, alphabet_i); // flat white noise with 1 state;
	}
	
	// Get projective probability;
	vector<double> proj_sd = synchronous_product(G, H).second;

	// Buffing Omega;
	vector<double> OmegaRow(alphabet_o, 0);
	vector<vector<double>> Omega(alphabet_i, OmegaRow);

	// Evaluating Omega;
	for (symbol symb_i(0); symb_i < alphabet_i; symb_i++)
	{
		for (symbol symb_o(0); symb_o < alphabet_o; symb_o++)
		{
			for (size_t cur_state = 0; cur_state < num_states; cur_state++)
			{
				int nxt_state = aut[cur_state][symb_i];
				Omega[symb_i][symb_o] += proj_sd[cur_state] * (Xpit[nxt_state][symb_o] - Xpit[cur_state][symb_o]);
			}
		}
	}
	
	return Omega;
}


void describe(PFSA& G, string name="PFSA")
{
	connx aut = G.get_aut();
	pitilde pit = G.get_pit();
	pitilde Xpit = G.get_Xpit();

	size_t num_states = aut.size();
	size_t alphabet = aut[0].size();


	cout << "\n############## DESCRIPTION of " << name << " START ##############" << endl;
	
	cout << "number of states = " << num_states << endl;
	cout << "alphabet size = " << alphabet << endl;

	cout << "\naut:" << endl;
	for (size_t i = 0; i < num_states; i++)
	{
		cout << "\t";
		for (symbol symb(0); symb < alphabet; symb++)
		{
			cout << aut[i][symb] << " ";
		} 
		cout << endl;
	}
	
	cout << "\npitilde:" << endl;
	for (size_t i = 0; i < num_states; i++)
	{
		cout << "\t" << pit[i] << endl;
	}
	
	if (Xpit.size() > 0)
	{
		cout << "\noutput probability:" << endl;
		for (size_t i = 0; i < num_states; i++)
		{
			cout << "\t" << Xpit[i] << endl;
		}
	}
	
	cout << "############### DESCRIPTION of " << name << " END ###############\n" << endl;
}


double HMM_log_likelihood(PFSA& H, symbol_list_ sl)
{
	// The PFSA H has to have Xpit
	connx aut = H.get_aut();
	pitilde pit = H.get_pit();
	pitilde Xpit = H.get_Xpit();
	size_t num_states = aut.size();
	
	// size_t alphabet_i = aut[0].size();
	// size_t alphabet_o = Xpit[0].size();

	vector<double> distr = H.get_Stationary();
	vector<vector<double>> PI = H.get_PI();

	double llk = 0.;
	for (symbol symb : sl)
	{
		double sum = 0.;
		vector<double> new_distr(num_states, 0.);
		for (size_t i = 0; i < num_states; i++)
		{
			double tmp = distr[i] * Xpit[i][symb];
			sum += tmp;
			for (size_t j = 0; j < num_states; j++)
			{
				new_distr[j] += tmp * PI[i][j];
			}
		}
		for (size_t i = 0; i < num_states; i++)
		{
			new_distr[i] /= sum;
		}
		llk -= log2(sum);
		distr = new_distr;
	}
	return llk / sl.size(); 
}

void test(PFSA& G, PFSA& H, size_t length)
{
	symbol_list_ sl = G.gen_data(length).get_symbol_list();
	
	PFSA Comp(G||H);
	describe(Comp, "Composition of G and H using the library");


	double llk_G = G.log_likelihood(sl);
	cout << "llk of G  = " << llk_G << endl;
	
	double llk_HMM = HMM_log_likelihood(Comp, sl);
	cout << "llk of H(G) = " << llk_HMM << endl;
}

void test_stat(PFSA& G, PFSA& H, size_t num_runs, size_t length)
{
	vector<double> llks;
	PFSA Comp(G||H);
	for (size_t n = 0; n < num_runs; n++)
	{
		symbol_list_ sl = G.gen_data(length).get_symbol_list();
		double llk = HMM_log_likelihood(Comp, sl);
		cout << n << ": " << llk << endl;
		llks.push_back(llk);
	}

	double mean = gsl_stats_mean(&llks[0], 1, llks.size());
	double std = gsl_stats_sd(&llks[0], 1, llks.size());	

	cout << "mean = " << mean << endl;
	cout << "std = " << std << endl;
}


int main(int argc, char** argv)
{ 
	string PFSA_file = argv[1];
	PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
	describe(G, "XPFSA G");
	
	string XPFSA_file = argv[2];	
	PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "XPFSA");
	describe(H, "XPFSA H");

	size_t length = (size_t)stoi(argv[3]);


	double min = -2;
	double max = 2;
	double range = max - min;
	double num_steps = 100;
	double step_size = range / num_steps;
 
	vector<double> scale;
	vector<PFSA> PFSA_vec;
	vector<double> er_vec;
	for (size_t i = 0; i <=num_steps; i++)
	{
		double x = min + i * step_size;
		scale.push_back(x);
		PFSA xG(G * x);	
		PFSA_vec.push_back(xG);
		er_vec.push_back(entropy(xG));
	}

	vector<double> llk_row(num_steps + 1, 0);
	vector<vector<double>> llk_matrix(num_steps + 1, llk_row);
	
	for (size_t i = 0; i <=num_steps; i++)
	{
		PFSA Comp(PFSA_vec[i]||H);
		for (size_t j = 0; j <= num_steps; j++) 
		{
			symbol_list_ sl = PFSA_vec[j].gen_data(length).get_symbol_list();
		}		
	}
	
	// size_t num_runs = (size_t)stoi(argv[4]);


	bool DEBUG = false;
	if (DEBUG)
	{
		PFSA Comp(G||H);
		describe(Comp, "Composition of G and H using the library");
		cout << Comp.get_Stationary() << endl;
	}

  	return 0;
}
