#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

#include <dirent.h>
#include <time.h>

//#include <gsl/gsl_statistics.h>
#include "semantic.h"
//#include <boost/filesystem.hpp>

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

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
size_t rand_int(size_t min, size_t max)
{
	srand(time(NULL));
	return min + rand() % (max + 1 - min);
}

void save_sequences(const vector<symbol_list_>& sls, string fname)
{
	ofstream output(fname);
	for (size_t i = 0; i < sls.size(); i++)
	{
		symbol_list_ sl = sls[i];
		for (size_t j = 0; j < sl.size(); j++)
		{
			symbol symb = sl[j];
			output << symb;
			if (j < sl.size() - 1)
				output << " ";  
		}
		output << endl;
	}
	output.close();
}

int main(int argc, char** argv)
{
  	const string version="Base set choice v0 2020 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
  
	string result_folder;
	
	size_t alphabet = 2; // Keep simple and start with binary.
	size_t num_states_min;
	size_t num_states_max;
	// bool use_transducer;
	size_t x_num_states_min;
	size_t x_num_states_max;
	size_t length;
	size_t num_samples_each_class;
	//size_t num_datasets;

	string PFSA_fname = "PFSA.cfg";
	string XPFSA_fname = "XPFSA.cfg";
	string seq_fname = "seq.dat";

 	
	options_description desc(
"\n\
-----------------------### TransducerStudy-comprehensive zed.uchicago.edu 2020 ###---------------------------------\n\
---------------------------------------------------------------------------------------------------------\n\
This is code is used to test the choice for the base PFSA in unsupervised llk calculation.\n\n\
The input includes:\n\
\t1) minimum number of states;\n\
\t2) maximum number of states;\n\
\t3) minimum number of transducer states, if not given do not use transducer;\n\
\t4) maximum number of transducer states, if not given do not use transducer;\n\
\t3) sequence length;\n\
\t4) sample size of each class;\n\
\t5) number of random datasets;\n\
Example (in zutil_ folder):\n\
./bin_practice/test_zbase_comprehensive -p cfgfiles/PFSA_M2_0.cfg -x cfgfiles/XPFSA_2_3_2.cfg -l 1000 -n 100\n\
---------------------------------------------------------------------------------------------------------\n"
);
  	desc.add_options()
    ("help,h", "Print help message.")
    ("version,V", "Print version number.")
    ("min_states,a", value<size_t>(&num_states_min)->default_value(2), "minimum number of PFSA states.")
    ("max_states,b", value<size_t>(&num_states_max)->default_value(5), "maximum number of PFSA states.")
    // ("use_transducer,T", value<bool>(&use_transducer)->default_value(false), "whether to use transducer.")
    ("min_states_x,c", value<size_t>(&x_num_states_min)->default_value(2), "minimum number of transducer states.")
    ("max_states_x,d", value<size_t>(&x_num_states_max)->default_value(5), "maximum number of transducer states.")
    ("length,l", value<size_t>(&length)->default_value(2000), "sequence length.")
    ("num_samples,n", value<size_t>(&num_samples_each_class)->default_value(25), "Number of samples in each class.")
	//("num_datasets,r", value<size_t>(&num_datasets), "number of random datasets generated.")
	("PFSA_fname,p", value<string>(&PFSA_fname)->default_value("PFSA.cfg"), "PFSA filename.")
	("XPFSA_fname,x", value<string>(&XPFSA_fname)->default_value("XPFSA.cfg"), "XPFSA filename.")
	("seq_fname,s", value<string>(&seq_fname)->default_value("seq.dat"), "Sequence filename.")
	("result_folder,f", value<string>(&result_folder), "folder name.");
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


	srand (time(NULL));
	vector<symbol_list_> sls;

	size_t num_states = rand_int(num_states_min, num_states_max);

	connx aut = SCC_UTIL__::get_random_aut(num_states, alphabet);
	pitilde pit = SCC_UTIL__::get_random_pit(aut);
	PFSA G = PFSA(pit, aut);
	G.describe("PFSA");
	save(G, result_folder + "/" + PFSA_fname, "PFSA");

	size_t use_transducer = rand() % 2;	
	if (use_transducer == 1)
	{
		size_t x_num_states = rand_int(x_num_states_min, x_num_states_max);
		connx x_aut = SCC_UTIL__::get_random_aut(x_num_states, alphabet);
		pitilde x_pit = SCC_UTIL__::get_uniform(x_num_states, alphabet);
		pitilde Xpit = SCC_UTIL__::get_random_pit(x_num_states, alphabet);
		PFSA H = PFSA(x_pit, x_aut); 
		H.set_Xpit(Xpit);
		H.describe("XPFSA");
		save(H, result_folder + "/" + XPFSA_fname, "XPFSA");
		
		for (size_t i = 0; i < num_samples_each_class; i++)
		{
			symbol_list_ sl = G.gen_transduced_data(H, length);
			sls.push_back(sl);
		}
	}
	else
	{
		for (size_t i = 0; i < num_samples_each_class; i++)
		{
			symbol_list_ sl = G.gen_data(length).get_symbol_list();
			sls.push_back(sl);
		}
	}

	string fname = result_folder + "/" + seq_fname;
	save_sequences(sls, fname);

	return 0;
}
