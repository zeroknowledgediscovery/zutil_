#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include "semantic.h"

#define DEBUG_ 0

using namespace boost::program_options;
//------------------------------------

double log_likelihood(PFSA& G,
		      const symbol_list_& s)
{
  if(s.empty())
    return 0.0;
  
  double llk=0.0;
  const unsigned int numstates=G.get_aut().size();

  map<symbol,vector<double> > map_pitcol;
  for(unsigned int st=0;st<numstates;++st)
    for(symbol i(0); i < G.get_aut()[0].size();++i)
      map_pitcol[i].push_back(G.get_pit()[st][i]);

  vector<double> curr_state(G.get_Stationary());

  for(unsigned int i=0;i<s.size();++i)
    {
      double pr=0.0;
      for(unsigned int st=0;st<numstates;++st)
	pr+=map_pitcol[s[i]][st]*curr_state[st];

      llk+=log2(pr);

      vector <double> state_vec_tmp(curr_state);
      for (unsigned int k=0; k < numstates; k++)
	{
	  double V=0.0000001;
	  for (unsigned int j=0; j < numstates; j++)
	    V+=G.get_Gamma()[s[i]][j][k]*state_vec_tmp[j];
	  curr_state[k]=V;
	}
      
      double S=0.0;
      for (unsigned int i=0;i<numstates;i++)
	S+=curr_state[i];

      if(S>0.0)
	S=1/S;
      
      for (unsigned int i=0;i<numstates;i++)
	curr_state[i] *= S;	
    }

  return -llk/s.size();
}

//------------------------------------
int main(int argc, char *argv[])
{
  const string version="Log-Likelihood v0.9 2018 zed.uchicago.edu";
  const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

  string seqfile="seqfile.dat",pfsafile="pfsa.mod";
  symbol_list_ seq;
  string DATA_DIR="across";
  unsigned int len=10000000;
  
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
    ("seq,s",value<string>(&seqfile), "input sequence file")
    ("datadir,D",value<string>(&DATA_DIR),"data direction: row or column [row]")
    ("datalen,x",value<unsigned int>(&len),"length max for input sequence [10000000]")
    ("pfsafile,f",value<string>(&pfsafile), "pfsa file");
  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    cout <<"empty arg, type -h or --help" << endl;
  try
    {
      store(command_line_parser(argc, argv)
	    .options(desc)
	    .run(), vm);
      notify(vm);
    } 
  catch (std::exception &e)
    {
      cout << endl << e.what() 
	   << endl << desc << endl;
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

  if (DATA_DIR=="row")
    DATA_DIR="across";

  if (DATA_DIR=="column")
    DATA_DIR="up";

  
  PFSA G = SCC_UTIL__::read_mc(pfsafile, "PFSA");
  vector<double> stationary_dist(G.get_Stationary());
  map < symbol, vector < vector < double > > > Gamma(G.get_Gamma());
  vector <vector <double> >  PI_(G.get_PI());

  data_reader *R;
  R = new data_reader(seqfile,DATA_DIR,len);
  vector < Symbolic_string_ > Svec(R->getsymbolic_string_vector());

  vector<double> llk;
  for(unsigned int i=0;i<Svec.size();++i)
    llk.push_back(log_likelihood(G,Svec[i].get_symbol_list()));

  cout << llk << endl;

  if(DEBUG_)
    {
      cout << "Stationary Distribution: " << stationary_dist << endl;
      
      cout << "Transition Probability Matrix: " << endl;
      for (unsigned int i=0; i < G.get_aut().size();++i)
	cout << PI_[i] << endl;

      cout << "Gamma Matrices: " << endl;
      for (symbol i(0); i < G.get_aut()[0].size();++i)
	{
	  cout << "symbol --> " << i << endl;
	  for (unsigned int j=0; j < G.get_aut().size(); j++)
	    cout << Gamma[i][j] << endl;
	}
      cout<< seq<< endl;
    }
 
  
  return 0;
}

