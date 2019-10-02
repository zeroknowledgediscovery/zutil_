#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include <boost/timer/timer.hpp>
#include <random>

#include "semantic.h"
#include "config.h"

using namespace boost::program_options;


vector<option> ignore_numbers(vector<string>& args)
{
  vector<option> result;
  int pos = 0;
  while(!args.empty())
    {
      const auto& arg = args[0];
      double num;
      if(boost::conversion::try_lexical_convert(arg, num))
	{
	  result.push_back(option());
	  option& opt = result.back();

	  opt.position_key = pos++;
	  opt.value.push_back(arg);
	  opt.original_tokens.push_back(arg);

	  args.erase(args.begin());
	}
      else
	break;
    }

  return result;
}
//------------------------------------------------------

//------------------------------------
int main(int argc, char *argv[])
{

  const string version="Run PFSA Model  v1.0 2019 zed.uchicago.edu";
  const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
  string MC_PRINT_FILE="";
  const gsl_rng_type * T;
  gsl_rng * r;  
  gsl_rng_env_setup(); 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set (r, random_seed());

  unsigned int len=1000;
  unsigned int repeat=1;
  string graphpref="pfsa", graphtype="png";
  bool DRAW_GRAPH=false, MC_PRINT=false;
  bool SHOW_STATIONARY=false, SHOW_PI=false, SHOW_GAMMA=false;
  string configfile="config.cfg";
  map <state, map <symbol, string> > REPLACE_SYM;

  connx aut;
  map < state, map < symbol, double > > pitilde_map;
  symbol_list_ omega;


  options_description desc( "### PFSA run zed.uchicago.edu 2019 ###\n\
--------------------------\n\
Example Usage:\n\
./prun -f pfsa.cfg\n\
./prun -f pfsa.cfg -l 1000\n\
./prun -f pfsa.cfg -l 1000 -o out.txt\n\
./prun -f pfsa.cfg -l 1000 -o out.txt -n 10\n\
 Usage");
  desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("mod,f",value<string>(&configfile), "model file")
    ("datalen,l",value<unsigned int>(&len),"output sequence length [1000]")
    ("numrepeat,n",value< unsigned int >(&repeat), "number of runs [1] ")
    ("dfile,o",value< string >(&MC_PRINT_FILE), "output file []")
    ("gamma,G",value< bool >(&SHOW_GAMMA), "print gamma [off]")
    ("pi,P",value< bool >(&SHOW_PI), "print pi [off]")
    ("stationary,S",value< bool >(&SHOW_STATIONARY), "print stationary [off]")
    ("machine,M",value< bool >(&MC_PRINT), "print machine [off]")
    ("graphpref,R",value< string >(&graphpref), "graph name [pfsa]")
    ("graphtype,T",value< string >(&graphtype), "graph type [png]")
    ("graph,g",value< bool >(&DRAW_GRAPH), "draw graph [off]");
  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    {
      cout <<"empty arg, type -h or --help" << endl;
      exit(0);
    }
  try
    {
      store(command_line_parser(argc, argv)
	    .extra_style_parser(&ignore_numbers)
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
      exit(0);
      return 1;
    }
  if (vm.count("version"))
    {
      cout << version << endl;
      exit(0);
      return 1;
    }


  CONFIG cfg(configfile);
  cfg.set_map<state,symbol,state>(aut,"CONNX");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");
  cfg.set_vector<symbol>(omega,"OMEGA");
  cfg.set_map<state,symbol,string>(REPLACE_SYM,"REPLACE_SYM");

  /** pitilde matrix for PFSAs */ 
  pitilde pit;  
  state c_state=0;
  for ( map < state, map < symbol, double > >::iterator iti=pitilde_map.begin();
	iti != pitilde_map.end(); iti++)
    {
      vector <double> vec_tmp;
      for (map <symbol, double>::iterator itj=iti->second.begin();
	   itj != iti->second.end();
	   itj++)
	vec_tmp.push_back(itj->second);
      pit[c_state++] = vec_tmp;
    }
  string graphname=graphpref;
  PFSA G(pit, aut);

  if (DRAW_GRAPH)
    drawGraph(aut,pit,graphname,graphtype);
  if (MC_PRINT) 
    G.mc_print();
  if (MC_PRINT_FILE == "")
    for(unsigned int i=0;i<repeat;++i)
      {
	symbol_list_ s=G.gen_data(len).get_symbol_list();
	cout << s << endl;
      }
  else
    {
      ofstream MOUT(MC_PRINT_FILE.c_str());
      for(unsigned int i=0;i<repeat;++i)
	{
	  symbol_list_ s=G.gen_data(len).get_symbol_list();
	  MOUT <<  s << endl;
	}
      MOUT.close();
    }

  if (SHOW_STATIONARY)
    cout << "Stationary distribution: "<<  G.get_Stationary() << endl;

  if (SHOW_PI)
    {
      cout << "Transition Probability Matrix: " << endl;
      vector < vector <double> > V(G.get_PI());
      for (unsigned int i=0; i < aut.size();++i)
	cout << V[i] << endl;
    }

  if (SHOW_GAMMA)
    {
      cout << "Gamma Matrices: " << endl;
      map < symbol, vector < vector < double > > > Gamma(G.get_Gamma());
      for (symbol i(0); i < aut[0].size();++i)
	{
	  cout << "symbol --> " << i << endl;
	  for (unsigned int j=0; j < aut.size(); j++)
	    cout << Gamma[i][j] << endl;
	}
    }

  if (!REPLACE_SYM.empty())
    {
      Symbolic_string_  sG(G.gen_data(len));
      symbol_list_ sg = sG.get_symbol_list();
      for (unsigned int i=0;i<sg.size();++i)
	cout << REPLACE_SYM[0][sg[i]];
      cout << endl;
    }

  //------------------------------------ 
  //------------------------------------ 
  return 0;
}
