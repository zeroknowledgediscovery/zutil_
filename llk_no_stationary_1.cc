#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include <time.h>
#include <omp.h>
#include <fstream>
#include "semantic.h"

#define DEBUG_ 0

using namespace boost::program_options;
//------------------------------------
//
bool is_positive(PFSA& G)
{
	connx aut = G.get_aut();
	pitilde pit = G.get_pit();
  	const unsigned int numstates = aut.size();
	const unsigned int alphabet = aut[0].size();

	for (size_t i = 0; i < numstates; i++)
	{
		map_sym_state aut_row = aut[i];
		vector<double> pit_row = pit[i];
		for (symbol symb(0); symb < alphabet; symb++)
		{
			state s = aut_row[symb];
			double pr = pit_row[symb];
			if ((s == -1) || (pr < 1e-10))
			{
				return false;
			}
		}
	}
	return true;
}

/*
 * Calculate the llk using the matrix form with a initial distribution
 */
double log_likelihood_matrix(
	const symbol_list_ & s,
	map<symbol, vector<vector<double>>> & Gamma, 
	map<symbol, vector<double>> & pitcol,
	const vector<double> & init_state_distr)
{
	if(s.empty())
	{
    	return -1;
	}

	double llk = 0;

	vector<double> curr_state(init_state_distr);
	size_t numstates = curr_state.size();

	for(unsigned int i = 0; i < s.size(); i++)
	{
		double pr = 0;
	  	for(unsigned int st = 0; st < numstates; st++)
		{
			pr += pitcol[s[i]][st] * curr_state[st];
		}
	  	llk += log(pr);

	  	vector <double> state_vec_tmp(curr_state);
	  	for (unsigned int k = 0; k < numstates; k++)
		{
	  		double V = 1e-7;
	  		for (unsigned int j = 0; j < numstates; j++)
			{
	    		V += Gamma[s[i]][j][k] * state_vec_tmp[j];
			}
	  		curr_state[k] = V;
		}
	  
	  	double S = 0.0;
	  	for (unsigned int i = 0; i < numstates; i++)
		{
			S += curr_state[i];
		}

	  	if(S > 0.0)
		{
			S = 1./S;
	  	}
	  	for (unsigned int i = 0; i < numstates; i++)
		{
			curr_state[i] *= S;
		}	
	}
	return -llk / s.size();
}


double log_likelihood_state(
	const symbol_list_ & s,
	connx & aut,
	pitilde & pit,
	int init_state)
{
	double llk = 0;
	int curr_state = init_state;
	for(unsigned int j = 0; j < s.size(); j++)
	{
		symbol symb(s[j]);
		double pr = pit[curr_state][symb];
		curr_state = aut[curr_state][symb];
		if ((curr_state == -1) || (pr == 0))
		{
			return -1;
		}
		else
		{
	  		llk += log(pr);
		}
	}	
	return -llk / s.size();
}


/*
 * There are four choices for method:
 * 	- "stationary": use the stationary distribution;
 * 	- "uniform": use uniform distribution;
 * 	- "random": use a randomly generated list of states as entering state;
 * 	- "fixed": use [0,...,num_states - 1] as the list of entering state;  
 */
vector<double> log_likelihood(
	PFSA& G, 
	const vector<symbol_list_> & seqs, 
	string method="stationary",
	size_t nProcessors=1)
{
	connx aut = G.get_aut();
	pitilde pit = G.get_pit();
	size_t numstates = aut.size();
	map<symbol, vector<vector<double>>> Gamma = G.get_Gamma(); 
	map<symbol, vector<double>> pitcol;
  	for(size_t st = 0; st < numstates; st++)
    	for(symbol symb(0); symb < aut[0].size(); symb++)
    	  	pitcol[symb].push_back(pit[st][symb]);
	vector<double> stationary(G.get_Stationary());
	vector<double> uniform(numstates, 1. / (double)numstates); 
	
	vector<double> llk(seqs.size(), 0);

	omp_set_num_threads(nProcessors);

	if ((method.compare("stationary") == 0) || (method.compare("uniform") == 0))
	{
		vector<double> init_state_distr = (method.compare("stationary") == 0)? stationary : uniform;
    	#pragma omp parallel
		{
			#pragma omp for
  			for(size_t i = 0; i < seqs.size(); i++)
			{
    			llk[i] = log_likelihood_matrix(seqs[i], Gamma, pitcol, init_state_distr);
			}
		}
	}
	else if ((method.compare("fixed") == 0) || (method.compare("random") == 0))
	{
		vector<state> state_list(numstates, 0);
		for (size_t i = 0; i < numstates; i++)
			state_list[i] = i;

		if (method.compare("random") == 0)
			random_shuffle(state_list.begin(), state_list.end());

    	#pragma omp parallel
		{
			#pragma omp for
  			for(size_t i = 0; i < seqs.size(); i++)
			{	
				bool is_successful = false;
				for (size_t j = 0; j < numstates; j++)
				{
					llk[i] = log_likelihood_state(seqs[i], aut, pit, state_list[j]);
					if (llk[i] > 0)
					{
						is_successful = true;
						break;
					}
				}
				if (!is_successful)
					llk[i] = log_likelihood_matrix(seqs[i], Gamma, pitcol, stationary);
			}
		}
	}
	else
	{
		cout << "method error! Choose from \"stationary\", \"uniform\", \"fixed\", and \"random\"" << endl;
	}
	return llk;
}

//------------------------------------
int main(int argc, char *argv[])
{
	const string version = "Log-Likelihood v0.9 2018 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE = "Exiting. Type -h or --help for usage";

  	string seqfile;
	string pfsafile;
  	symbol_list_ seq;
  	string DATA_DIR="across";
  	unsigned int len=10000000;
	string output="llk";
	size_t nProcessors;
	string method;
  
  	options_description desc( "### Loglikelihood zed.uchicago.edu 2018 ###\n\
		--------------------------\n\
		Note: Multiple input sequences can be given,\n\
		one in each new line of the file (if data in rows) named with option -s\n\
		If data is in columns, then each column is read as a new sequence\n\
		datadir is specified with option -D (deafult row)\n\
		--------------------------\n\
		Example (in testsuite directory): ../bin/llk -f S2.cfg -s seq.dat -x 100\n Usage");
  	desc.add_options()
        ("help,h", "print help message.")
        ("version,V", "print version number")
        ("seq,s", value<string>(&seqfile), "input sequence file")
        ("datadir,D", value<string>(&DATA_DIR),"data direction: row or column [row]")
        ("datalen,x", value<unsigned int>(&len),"length max for input sequence [10000000]")
        ("pfsafile,f", value<string>(&pfsafile), "pfsa file")
        ("method,m", value<string>(&method)->default_value("stationary"), "method of llk [stationary, uniform, fixed, random]")
        ("num_processors,n", value<size_t>(&nProcessors), "number of processors")
        ("output,o", value<string>(&output), "output file. If missing print to stdout.");
  
	positional_options_description p;
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
	if (vm.count("num_processors") == 0)
	{
		nProcessors = omp_get_max_threads();
	}


  	if (DATA_DIR=="row")
	{
    	        DATA_DIR="across";
	}
  	if (DATA_DIR=="column")
	{
    	        DATA_DIR="up";
	}

	srand(time(NULL));
  
  	PFSA G = SCC_UTIL__::read_mc(pfsafile, "PFSA");

  	data_reader *R;
  	R = new data_reader(seqfile, DATA_DIR, len);
  	vector<Symbolic_string_> Svec(R->getsymbolic_string_vector());
	vector<symbol_list_> seqs;
	for (size_t i = 0; i < Svec.size(); i++)
		seqs.push_back(Svec[i].get_symbol_list());

  	vector<double> llk = log_likelihood(G, seqs, method, nProcessors);

        if (vm.count("output"))
        {
	        ofstream outFile(output);
	        for (double l : llk)
	        {
		        outFile << l << " ";
	        }
	        outFile.close();
        }
        else
        {
                for (size_t i = 0; i < llk.size(); i++)
                {
                        cout << llk[i];
                        (i < llk.size() - 1)? cout << " ": cout << endl;
                }
        }

  	return 0;
}
