#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <boost/timer/timer.hpp>

#include "semantic.h"

#define DEBUG 0

using namespace boost::program_options;

const string VERSION="drawpfsa v1.0 \n Copyright Ishanu Chattopadhyay 2018";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";


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
//------------------------------------
int main(int argc, char *argv[])
{
  unsigned int len=500;
  string graphpref="pfsa", dotcfg="dotII.cfg";
  bool  MC_PRINT=false, 
    MC_STRING=false, INVERSE=false;
  bool SHOW_STATIONARY=false, SHOW_PI=false, SHOW_GAMMA=false;
  string configfile="config.cfg";
  map <state, map <symbol, string> > REPLACE_SYM;
  string title="";
  int DRAW_GRAPH=1;


 options_description infor( "Program information");
  infor.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number");


  options_description usg( "Usage");
  usg.add_options()
    ("configfile,c",value<string>(&configfile), "config file [default: config.cfg]")
    ("datalen,x",value< unsigned int >(&len), "data length max for input sequence")
    ("graphpref,N",value< string>(&graphpref), "preface of outfile [pfsa]")
    ("dotname,d",value< string>(&dotcfg), "dot file cfg [dotII.cfg]")
    ("graphtype,D",value<int>(&DRAW_GRAPH), "model type, 1 for pfsa, 2 for xpfsa [1]")
    ("mcprint,P",value<bool>(&MC_PRINT), "print machine [off]")
    ("mcinverse,I",value<bool>(&INVERSE), "use inverse machine [off]")
    ("stationary,S",value<bool>(&SHOW_STATIONARY), "show stationary [off]")
    ("pi,p",value<bool>(&SHOW_PI), "show transition matrix [off]")
    ("gamma,G",value<bool>(&SHOW_GAMMA), "show Gamma matrix [off]");

  options_description desc( "Draw PFSA/XPFSA");
  desc.add(infor).add(usg);

  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    {
      cout << EMPTY_ARG_MESSAGE << endl;
      return 1;
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
      return 1;
    }
  if (vm.count("help"))
    {
      cout << desc << endl;
      return 1;
    }
  if (vm.count("version"))
    {
      cout << VERSION << endl; 
      return 1;
    }

  connx aut;
  map < state, map < symbol, double > > pitilde_map;
  symbol_list_ omega;

  CONFIG cfg(configfile);
  cfg.set_map<state,symbol,state>(aut,"CONNX");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");


  /** pitilde matrix for PFSAs */ 
  pitilde pit;  
  //if (!RAND_PIT)
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
  PFSA G(pit, aut,0,len);

  //cout << dotcfg << endl;
  if(DEBUG)
    cout << pit.size() << endl;

  if (DRAW_GRAPH==1)
    G.drawGraph(graphname,dotcfg,title);
  if (DRAW_GRAPH==2)
    G.drawGraphX(graphname,dotcfg,title);


  
  return 0;
}

 
