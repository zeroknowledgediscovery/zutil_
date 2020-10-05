#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <math.h>

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
  	const string version = "Generate sequences v0 2020 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE = "Exiting. Type -h or --help for usage";

	string PFSA_fname;
	string XPFSA_fname;
        size_t length;
        size_t num_samples;
	string seq_fname;
 	
	options_description desc(
"\n\
-----------------------### Simple sequence generation zed.uchicago.edu 2020 ###---------------------------------\n\
---------------------------------------------------------------------------------------------------------\n\
The input includes:\n\
\t1) PFSA filename;\n\
\t2) XPFSA filename (Optional);\n\
\t3) sequence length;\n\
\t4) number of sequences;\n\
\t5) sequence filename.\n\
Example (in zutil_ folder):\n\
./bin_practice/gen_seq_simple -p cfgfiles/PFSA_M2_0.cfg -x cfgfiles/XPFSA_2_3_2.cfg -l 1000 -n 100\n\
---------------------------------------------------------------------------------------------------------\n"
);
  	desc.add_options()
        ("help,h", "Print help message.")
        ("version,V", "Print version number.")
	("PFSA_fname,p", value<string>(&PFSA_fname), "PFSA filename.")
	("XPFSA_fname,x", value<string>(&XPFSA_fname), "XPFSA filename.")
        ("length,l", value<size_t>(&length)->default_value(2000), "sequence length.")
        ("num_samples,n", value<size_t>(&num_samples)->default_value(25), "Number of samples in each class.")
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
	PFSA G = SCC_UTIL__::read_mc(PFSA_fname, "PFSA");
        G.describe("PFSA");
        PFSA H;
        if (vm.count("XPFSA"))
        {
	        H = SCC_UTIL__::read_mc(XPFSA_fname, "XPFSA");
                H.describe("XPFSA");
        }
        
        srand(time(NULL));
	
        vector<symbol_list_> sls;
        if (vm.count("XPFSA"))
        {
	        for (size_t i = 0; i < num_samples; i++)
	        {
		        symbol_list_ sl = G.gen_transduced_data(H, length);
		        sls.push_back(sl);
	        }
        }
	else
	{
		for (size_t i = 0; i < num_samples; i++)
		{
			symbol_list_ sl = G.gen_data(length).get_symbol_list();
			sls.push_back(sl);
		}
	}

	save_sequences(sls, seq_fname);

	return 0;
}
