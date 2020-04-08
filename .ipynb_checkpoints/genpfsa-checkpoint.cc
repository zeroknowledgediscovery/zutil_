#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include <list>
#include <boost/timer/timer.hpp>
#include <random>

#include "semantic.h"

#define DEBUG_ 0

using namespace boost::program_options;
//------------------------------------
template<typename Numeric, typename Generator = std::mt19937>
Numeric gen_random(Numeric from, Numeric to)
{
    thread_local static Generator gen(std::random_device{}());

    using dist_type = typename std::conditional
    <
        std::is_integral<Numeric>::value
        , std::uniform_int_distribution<Numeric>
        , std::uniform_real_distribution<Numeric>
    >::type;

    thread_local static dist_type dist;

    return dist(gen, typename dist_type::param_type{from, to});
}
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



//------------------------------------
int main(int argc, char *argv[])
{
  const string version="Log-Likelihood Smash v0.9 2019 zed.uchicago.edu";
  const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

  string seqfile="",ofile="L.dst";
  string MC_TYPE="M";
  bool PRINT_MC=true;
  unsigned int alphabet=2;
  unsigned int depth=2;
  unsigned int numstates=0;
  
  options_description desc( "### M/C Gen zed.uchicago.edu 2018 ###\n\
--------------------------\n\
Example Usage:\n\\n\
 Usage");
  desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("depth,D",value<unsigned int>(&depth),"depth [2]")
    ("numstates,N",value<unsigned int>(&numstates),"numstates [2], overrides depth")
    ("alphabet,A",value<unsigned int>(&alphabet),"alphabet [2]")
    ("mctype,T",value< string>(&MC_TYPE), "mc type: M | S | T")
    ("ofile,o",value< string >(&ofile), "output file [L.dst]")
    ("print,p",value< bool >(&PRINT_MC), "print PFSAs [on]");
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

  if(MC_TYPE=="M")
    {
      if(numstates >0)
	MESSAGE("WARNING: redefining numstates to accomodate MC_TYPE");
      numstates=pow(alphabet,depth);
    }

  if(numstates==0)
    numstates=pow(alphabet,depth);

  //depth d-markov m/c
  PFSA G=SCC_UTIL__::generate_mc(alphabet,
				 numstates,
				 MC_TYPE);


  if(PRINT_MC)
    G.mc_print();
  if(ofile!="")
    {
      ofstream out(ofile.c_str());
      G.mc_print(out);
      out.close();
    }
  return 0;
}

