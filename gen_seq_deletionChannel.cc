#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>
#include <random>

#include <dirent.h>
#include <time.h>

#include "semantic.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;


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
  	const string version = "Generate sequences for the deletion channel project";
  	const string EMPTY_ARG_MESSAGE = "Exiting. Type -h or --help for usage";

	vector<string> PFSA_fnames;
        size_t length;
        size_t num_samples;
        double deletion_rate;
	string seq_fname;
 	
	options_description desc(
"\n\
-----------------------### Simple sequence generation zed.uchicago.edu 2020 ###---------------------------------\n\
---------------------------------------------------------------------------------------------------------\n\
Example (in zutil_ folder):\n\
./bin_practice/gen_seq_simple -p cfgfiles/PFSA_M2_0.cfg -l 1000 -n 100\n\
---------------------------------------------------------------------------------------------------------\n"
);
  	desc.add_options()
        ("help,h", "Print help message.")
        ("version,V", "Print version number.")
	("PFSA_fnames,p", value<vector<string>>(&PFSA_fnames)->multitoken(), "PFSA filenames.")
        ("length,l", value<size_t>(&length)->default_value(2000), "sequence length.")
        ("num_samples,n", value<size_t>(&num_samples)->default_value(25), "Number of samples in each class.")
        ("deletion_rate,d", value<double>(&deletion_rate), "Deletion rate.")
	("seq_fname,s", value<string>(&seq_fname), "Sequence filename.");
  	positional_options_description p;
  	variables_map vm;
  	if (argc == 1)
	{
    	        cout << EMPTY_ARG_MESSAGE << endl;
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


        default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0, 1.0);

        vector<symbol_list_> sls;
        for (string PFSA_fname : PFSA_fnames)
        {
	        PFSA G = SCC_UTIL__::read_mc(PFSA_fname, "PFSA");
		for (size_t i = 0; i < num_samples; i++)
		{
			symbol_list_ sl = G.gen_data(length).get_symbol_list();
                        symbol_list_ sl_del;
                        for (symbol s: sl)
                        {
                                double key = distribution(generator);
                                // cout << key << " ";
                                if (key > deletion_rate)
                                        sl_del.push_back(s);
                        }
			sls.push_back(sl_del);
                        // cout << endl;
		}
	}
	save_sequences(sls, seq_fname);

	return 0;
}
