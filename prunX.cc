#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include "semantic.h"
#include "config.h"

using namespace boost::program_options;




PFSA read_mc(string filename, string TYPE__)
{

  connx aut;
  pitilde pit, Xpit;
  map < state, map < symbol, double > > pitilde_map;

  CONFIG modfile(filename);
  modfile.set_map<state,symbol,state>(aut,"CONNX");
  modfile.set_map<state,symbol,double>(pitilde_map,"PITILDE");


  PFSA *G;
  state c_state=0;
  for ( map < state, map < symbol, 
	  double > >::iterator iti=pitilde_map.begin();
	iti != pitilde_map.end();
	iti++)
    {
      vector <double> vec_tmp;
      for (map <symbol, double>::iterator itj=iti->second.begin();
	   itj != iti->second.end();
	   itj++)
	vec_tmp.push_back(itj->second);
      pit[c_state++] = vec_tmp;
    }

  if(TYPE__=="XPFSA")
    {
      pitilde pit2;
      for (state  j=0; j < (int)aut.size();++j)
	pit2[j]= vector <double> (aut[0].size(),1.0/(aut[0].size()+0.0));
      G = new PFSA(pit2,aut);
    }
  else
    G = new PFSA (pit,aut);
  G->set_Xpit(pit);

  return *G;
};




//-----------------------------------
//-----------------------------------
int main(int argc,char* argv[])
{ 
  const string version="Causality Networks v1.0 copyright Ishanu Chattopadhyay 2015";
  const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

  string configfile="config.cfg",pfsafile="pfsa.mod", xpfsafile="xpfsa.mod";
  vector <double> partition;
  unsigned int len=1000;
  string DATA_DIR="across";
  string DATA_TYPE="continuous";

  options_description desc( "Usage: ");
  desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("datalen,x",value< unsigned int >(), "data length maximum for input sequence")
    ("pfsafile,v",value<string>(), "pfsa file")
    ("xpfsafile,w",value<string>(), "xpfsa file")
    ("configfile,c",value<string>(), "config file");
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

  if (vm.count("datalen"))
    len=vm["datalen"].as<unsigned int>();
  if (vm.count("configfile"))
    configfile=vm["configfile"].as<string>();
  if (vm.count("pfsafile"))
    pfsafile=vm["pfsafile"].as<string>();
  if (vm.count("xpfsafile"))
    xpfsafile=vm["xpfsafile"].as<string>();

  PFSA G = read_mc(pfsafile, "PFSA");
  PFSA X = read_mc(xpfsafile, "XPFSA");


  Symbolic_string_ data = G.gen_data(len);

  cout << data << endl;

  Symbolic_string_ dataX = data.operator*(X);


  cout << dataX << endl;


}
