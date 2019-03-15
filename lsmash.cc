#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include <boost/timer/timer.hpp>

#include "semantic.h"

#define DEBUG_ 0

using namespace boost::program_options;
//------------------------------------

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

connx autM2={{0,{{symbol(0),0},{symbol(1),1}}},
	     {1,{{symbol(0),0},{symbol(1),1}}}};
pitilde pitM2={{0,{0.3,0.7}},{1,{0.7,0.3}}};

connx autS2={{0,{{symbol(0),0},{symbol(1),1}}},{1,{{symbol(0),1},{symbol(1),0}}}};
pitilde pitS2={{0,{0.3,0.7}},{1,{0.7,0.3}}};

connx autT3={{0,{{symbol(0),1},{symbol(1),2}}},{1,{{symbol(0),2},{symbol(1),0}}},{2,{{symbol(0),0},{symbol(1),1}}}};
pitilde pitT3={{0,{0.3,0.7}},{1,{0.7,0.3}},{2,{0.6,0.4}}};

connx autM4={{0,{{symbol(0),0},{symbol(1),1}}},{1,{{symbol(0),2},{symbol(1),3}}},{2,{{symbol(0),0},{symbol(1),1}}},{3,{{symbol(0),2},{symbol(1),3}}}};
pitilde pitM4={{0,{0.3,0.7}},{1,{0.7,0.3}},{2,{0.8,0.2}},{3,{0.2,0.8}}};

//------------------------------------
int main(int argc, char *argv[])
{
  const string version="Log-Likelihood Smash v0.9 2019 zed.uchicago.edu";
  const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

  string seqfile="seqfile.dat",ofile="L.dst";
  vector<string> pfsafile;
  symbol_list_ seq;
  string DATA_DIR="across";
  unsigned int len=10000000;
  vector <double> partition;
  string DATA_TYPE="continuous";
  bool DERIVATIVE=false;
  bool TIMER=true;

  options_description desc( "### Loglikelihood zed.uchicago.edu 2018 ###\n\
--------------------------\n\
Example Usage:\n\
../bin/lsmash -s seq.dat\n\
../bin/lsmash -s seq.dat -x 100 (restrict length of data read)\n\
../bin/lsmash -s seq.dat -x 100 -o L.dst (specify output file)\n\
../bin/lsmash -f S2.cfg M2.cfg T3.cfg -s seq.dat -x 100 (specify PFSA projectors)\n\
 Usage");
  desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("seq,s",value<string>(&seqfile), "input sequence file")
    ("datadir,D",value<string>(&DATA_DIR),"data direction: row or column [row]")
    ("datatype,T",value< string>(&DATA_TYPE), "data type: continous or symbolic")
    ("datalen,x",value<unsigned int>(&len),"length max for input sequence [10000000]")
    ("partition,P",value< vector<double> >(&partition)->multitoken(), "partition")
    ("use_derivative,u",value<bool>(&DERIVATIVE), "use derivative [false]")
    ("pfsafile,f",value< vector<string> >(&pfsafile)->multitoken(), "pfsa files")
    ("timer,t",value< bool >(&TIMER), "display timer [1 (true)] ")
    ("dfile,o",value< string >(&ofile), "output file [L.dst]");
  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    cout <<"empty arg, type -h or --help" << endl;
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

  if (DATA_DIR=="row")
    DATA_DIR="across";

  if (DATA_DIR=="column")
    DATA_DIR="up";

  vector<PFSA> G;
  if(pfsafile.empty())
    {
      G.push_back(PFSA (pitM2,autM2));
      G.push_back(PFSA (pitM4,autM4));
      G.push_back(PFSA (pitS2,autS2));
      G.push_back(PFSA (pitT3,autT3));
    }
  else
    for(unsigned int i=0;i<pfsafile.size();++i)
      G.push_back(SCC_UTIL__::read_mc(pfsafile[i], "PFSA"));
  
  data_reader *R;
  if (DATA_TYPE=="continuous")
    R = new data_reader(seqfile,DATA_DIR,partition,len,false,DERIVATIVE);
  else
    R = new data_reader(seqfile,DATA_DIR,len);


  //cout << G.size() << endl;
  matrix_dbl D;

  vector <symbol_list_> S= R->getlist_vector();

  for(unsigned int i=0;i<S.size();++i)
    {
      Symbolic_string_ s(S[i]);
      S[i]=(~s).get_symbol_list();
    }
  
  if (TIMER)
    {
      timer::auto_cpu_timer t;
      D= SCC_UTIL__::llk_distance(S,G);
    }
  else
    D= SCC_UTIL__::llk_distance(R->getlist_vector(),G);
    
  // cout << D << endl;

  ofstream out(ofile.c_str());
  out << D;
  out.close();
  
  return 0;
}

