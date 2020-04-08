#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

#include <dirent.h>
#include <time.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

const string currentDateTime() 
{
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%m%d-%X", &tstruct);
    
	return buf;
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

void save(PFSA& G, string filename, string type)
{
	connx aut = G.get_aut();
	pitilde pit;
	if (type.compare("PFSA") == 0)
	{
		pit = G.get_pit();	
	}
	else if (type.compare("XPFSA") == 0)
	{
		pit = G.get_Xpit();
	}
	else
	{
		cout << "Incorrect save type!" << endl;
		exit(1);
	}

	ofstream output(filename);
	output << "#CONNX" << endl;
	for (size_t i = 0; i < aut.size(); i++)
	{
		for (symbol s(0); s < aut[i].size(); s++)
		{
			output << aut[i][s];
			(s < aut[i].size() - 1) ? output << " " : output << endl;
		}
	}
	
	output << "\n#PITILDE" << endl;
	for (size_t i = 0; i < pit.size(); i++)
	{
		for (size_t j = 0; j < pit[i].size(); j++)
		{
			output << pit[i][j];
			(j < pit[i].size() - 1) ? output << " " : output << endl;
		}
	}
	output.close();
}

void get_llk_matrix(
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


void get_dist_matrix(
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


void get_HMM_llk_matrix(
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
  	const string version="TransducerStudy-comprehensive v0 2020 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
  
	string PFSA_filename;
	string XPFSA_filename;
	size_t length;
	string random_pfsa_parameters;
	string random_xpfsa_parameters;
	string result_folder;

	double min;
	double max;
	size_t num_samples;

	// Parameters for weighted finite-depth distance
	size_t depth;
	size_t num_runs;
 
	// Inferred the parameters
	PFSA G;
	PFSA H; 
	size_t num_states;
	size_t alphabet;
	size_t x_num_states;
	size_t x_alphabet;
	connx aut;
	pitilde pit;
	connx x_aut;
	pitilde Xpit;
 	
	options_description desc(
"\n\
-----------------------### TransducerStudy-comprehensive zed.uchicago.edu 2020 ###---------------------------------\n\
---------------------------------------------------------------------------------------------------------\n\
This is code is used to calculate the llk matrix so as to study the eigen structure of XPFSA transducer.\n\n\
The input includes:\n\
\t1) A PFSA G, given either by a config filename, or a specification of number of states and alphabet size;\n\
\t2) An XPFSA transducer H, given either by a config filename, or a specification of number of states and alphabet size;\n\
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
./bin_practice/test_zbase_comprehensive -p cfgfiles/PFSA_M2_0.cfg -x cfgfiles/XPFSA_2_3_2.cfg -l 1000 -n 100\n\
---------------------------------------------------------------------------------------------------------\n"
);
  	desc.add_options()
    ("help,h", "Print help message.")
    ("version,V", "Print version number.")
    ("pfsa,p", value<string>(&PFSA_filename), "PFSA filename.")
    ("random_pfsa_parameters,q", value<string>(&random_pfsa_parameters), "Parameters for randomly \
generated PFSA. Input in the form \"input_num_states input_alphabet\"")
    ("xpfsa,x", value<string>(&XPFSA_filename), "XPFSA filename.")
    ("random_xpfsa_parameters,y", value<string>(&random_xpfsa_parameters), "Parameters for randomly \
generated XPFSA. Input in the form \"output_num_states output_alphabet\"")
    ("length,l", value<size_t>(&length)->default_value(1000), "Length of the sequence.")
    ("min,a", value<double>(&min)->default_value(-2), "Minimum value of the range.")
    ("max,b", value<double>(&max)->default_value(2), "Maximum value of the range. Must be strictly bigger than min.")
    ("num_samples,n", value<size_t>(&num_samples)->default_value(2), "Number of sample points.\
 The first sample point is min while the last sample point is max.")
	("result_folder,f", value<string>(&result_folder), "If no result folder specified\
a folder ./TransducerStudy_result_mmdd-HH:MM:SS will be used.")
	("depth,d", value<size_t>(&depth)->default_value(5), "The depth to evaluate the weighted finite-depth distance (WFDD).")
	("num_run,r", value<size_t>(&num_runs)->default_value(200), "The number of runs to evaluate the mean and standard deviation of WFDD.");
  	positional_options_description p;
  	variables_map vm;
  	if (argc == 1)
	{
    	cout <<"empty arg, type -h or --help for manual" << endl;
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

  	if ((vm.count("pfsa") == 0) && (vm.count("random_pfsa_parameters") == 0))
	{
    	cout << "Either specify a PFSA file or specify parameters for randomly generated PFSA"<< endl;
		return 1;
	}
	else if ((vm.count("pfsa") == 1))
	{
		if (vm.count("random_pfsa_parameters") == 1)
		{
			cout << "parameter --pfsa(-p) will override parameter --random_pfsa_parameters(-q)";
		}
		G = SCC_UTIL__::read_mc(PFSA_filename, "PFSA");
	
		aut = G.get_aut();
		num_states = aut.size();
		alphabet = aut[0].size();
	}
	else
	{
		stringstream ss;
		ss << random_pfsa_parameters;
		string a;
		ss >> a;
		num_states = (size_t)stoi(a);
		ss >> a;
		alphabet = (size_t)stoi(a);
		
		aut = SCC_UTIL__::get_random_aut(num_states, alphabet);
		pit = SCC_UTIL__::get_random_pit(aut);
		G = PFSA(pit, aut);
		// PFSA_filename = "PFSA.cfg";
	}
  	
	if ((vm.count("xpfsa") == 0) && (vm.count("random_xpfsa_parameters") == 0))
	{
    	cout << "Either specify a XPFSA file or specify parameters for randomly generated PFSA" << endl;
		return 0;
	}
	else if ((vm.count("xpfsa") == 1))
	{
		if (vm.count("random_xpfsa_parameters") == 1)
		{
			cout << "parameter xpfsa(x) will override parameter random_xpfsa_parameters(y)";
		}
		H = SCC_UTIL__::read_mc(XPFSA_filename, "XPFSA");
		
		x_aut = H.get_aut();
		x_num_states = x_aut.size();
		x_alphabet = x_aut[0].size();
	}
	else
	{
		stringstream ss;
		ss << random_xpfsa_parameters;
		string a;
		ss >> a;
		x_num_states = (size_t)stoi(a);
		ss >> a;
		x_alphabet = (size_t)stoi(a);
		
		x_aut = SCC_UTIL__::get_random_aut(x_num_states, alphabet);
		pitilde x_pit = SCC_UTIL__::get_uniform(x_num_states, alphabet);
		Xpit = SCC_UTIL__::get_random_pit(x_num_states, x_alphabet);
		
		H = PFSA(x_pit, x_aut); 
		H.set_Xpit(Xpit);
		// XPFSA_filename = "XPFSA.cfg";
	}
	if (vm.count("result_folder") == 0)
	{
		result_folder = "./TransducerStudy_result_" + currentDateTime();
		boost::filesystem::create_directory(result_folder);
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

	G.describe("PFSA G");
	H.describe("Transducer");
	save(G, result_folder + "/PFSA.cfg", "PFSA");
	save(H, result_folder + "/XPFSA.cfg", "XPFSA");

	
	double range = max - min;
	double step_size = range / (num_samples - 1);
 	vector<double> scale(num_samples, 0.);

	for (size_t i = 0; i < num_samples; i++)
	{
	 	scale[i] = min + i * step_size;
	}
	 
	clock_t start, end; 

	cout << "\nLog-likelihood" << endl;

	start = clock();
	get_llk_matrix(G, H, scale, length, result_folder + "/llk");
	end = clock();
	cout << "\ttime = " << double(end - start) / double(CLOCKS_PER_SEC) << endl; 
	
	cout << "\nWeighted Finite-Depth Distance" << endl;
	start = clock();
	get_dist_matrix(G, H, scale, length, depth, num_runs, result_folder + "/dist");
	end = clock();
	cout << "\ttime = " << double(end - start) / double(CLOCKS_PER_SEC) << endl; 

	cout << "\nHMM Log-likelihood" << endl;
	start = clock();
	get_HMM_llk_matrix(G, H, scale, length, result_folder + "/HMM_llk");
	end = clock();
	cout << "\ttime = " << double(end - start) / double(CLOCKS_PER_SEC) << endl; 

  	return 0;
}
