#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>

#include <dirent.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;


const bool DEBUG=false;

vector<PFSA> read_PFSA(string PFSA_dir)
{
	DIR *pDir;
	pDir = opendir(PFSA_dir.c_str());

	if (pDir == NULL)
	{
		cout << "Cannot open " << PFSA_dir << "!";
		exit(1);
	}

	vector<PFSA> PFSA_vec;
	struct dirent *pDirent;

	while ((pDirent = readdir(pDir)) != NULL)
	{
		string filename(pDirent->d_name);
		if (filename.compare(".") * filename.compare("..") != 0)
		{
			string PFSA_file = PFSA_dir + "/" + filename;
			PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
			PFSA_vec.push_back(G);
			
			if (DEBUG)
			{	
				cout << PFSA_file << endl;
  				vector<double> stationary_dist(G.get_Stationary());
      			cout << "Stationary Distribution: " << stationary_dist << endl;
      			
				connx delta(G.get_aut());
				unsigned int num_states = delta.size();
				unsigned int alphabet_size  = delta[1].size();

				cout << "Transition Matrix: " << endl;
      			for (size_t i = 0; i < num_states; i++)
				{
					cout << "\t";
					for (symbol s(0); s < alphabet_size; s++)
					{
						cout << delta[i][s] << " ";
					}
					cout << endl;
				}

  				vector<vector<double>>  pitilde(G.get_PI());
				cout << "Transition Probability Matrix: " << endl;
      			for (size_t i = 0; i < num_states; i++)
				{
					cout << "\t" << pitilde[i] << endl;
				}	
				cout << endl;
			}
		}
	}

	closedir(pDir);

	return PFSA_vec;
}

double get_log_likelihood(PFSA& G, 
						const vector<double>& stationary, 
						const symbol_list_& seq)
{
  	if(seq.empty())
	{
    	return 0.0;
	}
  
  	double llk = 0.0;
  	size_t num_states = G.get_aut().size();
 
 	vector<double> curr_distr(stationary);

  	for(symbol s : seq)
    {
      	double pr = 0.0;
      	for(unsigned int st = 0; st < num_states; st++)
		{
			pr += curr_distr[st] * G.get_pit()[st][s];
		}

      	llk += log2(pr);

      	vector <double> nxt_distr(num_states, 0.0);
	  	double sum = 0.0;

      	for (size_t st = 0; st < num_states; st++)
		{
			int nxt_state = G.get_aut()[st][s];
			
			if (nxt_state != -1)
			{
				double temp = curr_distr[st] * G.get_pit()[st][s];
				nxt_distr[nxt_state] += temp; 
				sum += temp;
			}
		}
      
      	for (size_t st = 0; st < num_states; st++)
		{
			curr_distr[st] = nxt_distr[st] / sum;	
		}
	}
	return -llk / seq.size();
}


double total_variation_distance(const vector<double>& vec1, const vector<double>& vec2)
{
	double distance = 0;
	for (size_t i = 0; i < vec1.size(); i++)
	{
		distance += abs(vec1[i] - vec2[i]);
	}
	return distance;
}


vector<vector<double>> get_llk_distance_matrix(vector<PFSA>& PFSA_vec, const vector<symbol_list_>& seqs)
{
	vector<vector<double>> stationary_distrs;
	for (PFSA G: PFSA_vec)
	{
		stationary_distrs.push_back(G.get_Stationary());
	}

	vector<vector<double>> llk_matrix;
	for (size_t i = 0; i < seqs.size(); i++)
	{
		vector<double> llk_vec;
		for (size_t j = 0; j < PFSA_vec.size(); j++)
		{
			llk_vec.push_back(get_log_likelihood(PFSA_vec[j], stationary_distrs[j], seqs[i]));
		}
		llk_matrix.push_back(llk_vec);
	}
	
	vector<vector<double>> distance_matrix;
	for (size_t i = 0; i < seqs.size(); i++)
	{
		vector<double> distance_vec;
		for(size_t j = 0; j < seqs.size(); j++)
		{
			distance_vec.push_back(total_variation_distance(llk_matrix[i], llk_matrix[j]));
		}
		distance_matrix.push_back(distance_vec);
	}
	return distance_matrix;
}



int main(int argc, char** argv)
{ 
	string PFSA_name = argv[1];
	string XPFSA_name = argv[2];
	int length = stoi(argv[3]);

	string PFSA_file = "/home/yhuang10/D3M/zutil_/cfgfiles/" + PFSA_name + ".cfg" ;
	PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
	
	string XPFSA_file = "/home/yhuang10/D3M/zutil_/cfgfiles/" + XPFSA_name + ".cfg" ;
	PFSA H = SCC_UTIL__::read_mc(XPFSA_file, "PFSA");

	Symbolic_string_ input = G.gen_data(length);
	cout << "Input:" << endl;
	cout << input << endl;

	Symbolic_string_ output = input * H;

	cout << "Output:" << endl;
	cout << output << endl;

	bool DEBUG = true;

	if (DEBUG)
	{	
		cout << XPFSA_file << endl;
		vector<double> stationary_dist(H.get_Stationary());
		cout << "Stationary Distribution: " << stationary_dist << endl;
		
		connx delta(H.get_aut());
		unsigned int num_states = delta.size();
		unsigned int alphabet_size  = delta[1].size();

		cout << "Transition Matrix: " << endl;
		for (size_t i = 0; i < num_states; i++)
		{
			cout << "\t";
			for (symbol s(0); s < alphabet_size; s++)
			{
				cout << delta[i][s] << " ";
			}
			cout << endl;
		}

		pitilde kai(H.get_Xpit());
		cout << "Transition Probability Matrix: " << endl;
		for (size_t i = 0; i < alphabet_size; i++)
		{
			cout << "\t" << kai[i] << endl;
		}	
		cout << endl;
	}
	
  	return 0;
}
