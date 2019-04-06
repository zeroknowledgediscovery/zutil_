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

  string seqfile="",ofile="L.dst";
  vector<string> pfsafile;
  symbol_list_ seq;
  string DATA_DIR="across";
  unsigned int len=10000000;
  vector <double> partition;
  string DATA_TYPE="symbolic";
  bool DERIVATIVE=false;
  bool TIMER=true, PRINT_MC=false;
  unsigned int RANDOM_MC=0;

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
    ("dfile,o",value< string >(&ofile), "output file [L.dst]")
    ("randomproj,R",value<unsigned int >(&RANDOM_MC), "no. of random machines to use [0]")
    ("machines,m",value< bool >(&PRINT_MC), "print PFSAs used [off]");
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

  if(seqfile=="")
      MESSAGE("ERROR: empty seq file");
  if(partition.empty())
    DATA_TYPE="symbolic";
  
  if (DATA_DIR=="row")
    DATA_DIR="across";

  if (DATA_DIR=="column")
    DATA_DIR="up";

    
  data_reader *R;
  if (DATA_TYPE=="continuous")
    R = new data_reader(seqfile,DATA_DIR,partition,len,false,DERIVATIVE);
  else
    R = new data_reader(seqfile,DATA_DIR,len);


  //cout << G.size() << endl;
  matrix_dbl D;

  vector <symbol_list_> S= R->getlist_vector();
  
  unsigned int alphabet=0;
  for(unsigned int i=0;i<S.size();++i)
    {
      Symbolic_string_ s(S[i]);
      symbol alph(s.get_alphabet());
      if(s.get_alphabet()>alphabet)
	alphabet=s.get_alphabet();
      //S[i]=(~s).get_symbol_list();
    }
  
  
  vector<PFSA> G;
  if(pfsafile.empty())
    {
      G.push_back(SCC_UTIL__::generate_mc(alphabet,
					  alphabet,
					  "M"));
      G.push_back(SCC_UTIL__::generate_mc(alphabet,
					  alphabet*alphabet,
					  "M"));
      G.push_back(SCC_UTIL__::generate_mc(alphabet,
					  alphabet,
					  "S"));
      G.push_back(SCC_UTIL__::generate_mc(alphabet,
					  alphabet*2,
					  "T"));
      /*      G.push_back(PFSA (pitM2,autM2));
      G.push_back(PFSA (pitM4,autM4));
      G.push_back(PFSA (pitS2,autS2));
      G.push_back(PFSA (pitT3,autT3));
      */
    }
  else
    for(unsigned int i=0;i<pfsafile.size();++i)
      G.push_back(SCC_UTIL__::read_mc(pfsafile[i], "PFSA"));

  if(RANDOM_MC>0)
    {
      size_t numG=G.size();
      for(unsigned int i=0;i<RANDOM_MC;++i)
	{
	  PFSA G_=SCC_UTIL__::FWN(G[0].get_aut()[0].size());
	  for(unsigned int j=0;j<numG;++j)
	      G_=G_ + (G[j] * gen_random<int>(0, 1)
		       * gen_random<double>(-2.0, 2.0));
	  G.push_back(G_);
	}
      
    }
  
  if(PRINT_MC)
    for(unsigned int i=0;i<G.size();++i)
      G[i].mc_print();


  if (TIMER)
    {
      timer::auto_cpu_timer t;
      D= SCC_UTIL__::llk_distance(S,G);
    }
  else
    D= SCC_UTIL__::llk_distance(S,G);
    
  // cout << D << endl;

  ofstream out(ofile.c_str());
  out << D;
  out.close();
  
  return 0;
}

