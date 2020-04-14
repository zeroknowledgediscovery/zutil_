#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <set>
#include <boost/timer/timer.hpp>
#include <random>
#include <cmath>

#include "semantic.h"

#define DEBUG_ 0

using namespace boost::program_options;
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
ostream& operator << (ostream &out, const vector<vector<double>>& matrix)
{
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			out << matrix[i][j];
			(j < matrix[i].size() - 1) ? out << " " : out << endl;
		}
	}
   	return out;
}
//------------------------------------------------------------------------------
connx autM2={{0,{{symbol(0),0},{symbol(1),1}}},
	     {1,{{symbol(0),0},{symbol(1),1}}}};
pitilde pitM2={{0,{0.3,0.7}},{1,{0.7,0.3}}};

connx autS2={{0,{{symbol(0),0},{symbol(1),1}}},{1,{{symbol(0),1},{symbol(1),0}}}};
pitilde pitS2={{0,{0.3,0.7}},{1,{0.7,0.3}}};

connx autT3={{0,{{symbol(0),1},{symbol(1),2}}},{1,{{symbol(0),2},{symbol(1),0}}},{2,{{symbol(0),0},{symbol(1),1}}}};
pitilde pitT3={{0,{0.3,0.7}},{1,{0.7,0.3}},{2,{0.6,0.4}}};

connx autM4={{0,{{symbol(0),0},{symbol(1),1}}},{1,{{symbol(0),2},{symbol(1),3}}},{2,{{symbol(0),0},{symbol(1),1}}},{3,{{symbol(0),2},{symbol(1),3}}}};
pitilde pitM4={{0,{0.3,0.7}},{1,{0.7,0.3}},{2,{0.8,0.2}},{3,{0.2,0.8}}};

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  	const string version="Log-Likelihood Smash Classification v0 2020 zed.uchicago.edu";
  	const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

  	string seq_train_file = "";
  	string seq_test_file = "";
	string label_train_file = "";
	string label_test_file = "";
	string coor_train_file = "";
	string coor_test_file = "";
  	vector<string> pfsafile;

  	unsigned int len;
  	vector<double> partition;
	string DATA_DIR;
  	string DATA_TYPE;
  	bool DERIVATIVE;

  	vector<PFSA> PFSA_vec;

  	options_description desc( "### Loglikelihood zed.uchicago.edu 2018 ###\n\
--------------------------\n\
Example Usage:\n\
../bin/lsmash -f seq.dat\n\
../bin/lsmash -f seq.dat -x 100 (restrict length of data read)\n\
../bin/lsmash -f seq.dat -x 100 -o L.dst (specify output file)\n\
../bin/lsmash -F S2.cfg M2.cfg T3.cfg -f seq.dat -x 100 (specify PFSA projectors)\n\
 Usage");
  	desc.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number")
    ("seq_train,f",value<string>(&seq_train_file), "train sequence file")
    ("seq_test,g",value<string>(&seq_test_file), "test sequence file")
    ("label_train,l",value<string>(&label_train_file), "train label file")
    ("label_test,k",value<string>(&label_test_file), "test label file")
    ("coor_train,c",value<string>(&coor_train_file), "train coordinate file")
    ("coor_test,d",value<string>(&coor_test_file), "test coordinate file")
    ("datadir,D",value<string>(&DATA_DIR)->default_value("row"),"data direction: row or column")
    ("datatype,T",value<string>(&DATA_TYPE)->default_value("symbolic"), "data type: continuous or symbolic")
    ("datalen,x",value<unsigned int>(&len)->default_value(10000000),"length max for input sequence")
    ("partition,P",value< vector<double>>(&partition)->multitoken(), "partition")
    ("use_derivative,u",value<bool>(&DERIVATIVE)->default_value(false), "use derivative")
    ("pfsafile,F",value<vector<string>>(&pfsafile)->multitoken(), "pfsa files");
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
      	cout << endl << e.what() << endl << desc << endl;
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
	if (vm.count("pfsafile"))
	{
	  	for (string PFSA_filename : pfsafile)
    	{
			PFSA G = SCC_UTIL__::read_mc(PFSA_filename, "PFSA");
			PFSA_vec.push_back(G);
		}
	}
	else
	{
		PFSA M2(pitM2, autM2);
		PFSA_vec.push_back(M2);
		PFSA M4(pitM4, autM4);
		PFSA_vec.push_back(M4);
		PFSA T3(pitT3, autT3);
		PFSA_vec.push_back(T3);
		PFSA S2(pitS2, autS2);
		PFSA_vec.push_back(S2);
	}
  	if (seq_train_file == "")
      	MESSAGE("ERROR: empty train sequence file");
  	if (seq_test_file == "")
      	MESSAGE("ERROR: empty test sequence file");
  	if (label_train_file == "")
      	MESSAGE("ERROR: empty train label file");
  	if (label_test_file == "")
      	MESSAGE("ERROR: empty train label file");

  	if (partition.empty())
    	DATA_TYPE="symbolic";
  	if (DATA_DIR=="row")
    	DATA_DIR="across";
  	if (DATA_DIR=="column")
    	DATA_DIR="up";

    
	timer::auto_cpu_timer t;

	// Get train, test, and label
	data_reader *R_train, *R_test;
  	if (DATA_TYPE == "continuous")
  	{
    	R_train = new data_reader(seq_train_file,DATA_DIR,partition,len,false,DERIVATIVE);
		R_test = new data_reader(seq_test_file,DATA_DIR,partition,len,false,DERIVATIVE);
  	}
  	else
  	{
    	R_train = new data_reader(seq_train_file,DATA_DIR,len);
    	R_test = new data_reader(seq_test_file,DATA_DIR,len);
  	}
	
	int l;
	ifstream input;
	vector<int> labels_train;
	input.open(label_train_file);
	while (input >> l)
		labels_train.push_back(l);	
	input.close();
	vector<int> labels_test;
	input.open(label_test_file);
	while (input >> l)
		labels_test.push_back(l);	
	input.close();

	set<int> class_labels;
	for (int l : labels_train)
	{
		class_labels.insert(l);
	}
	

	// ---------------------Get llk featurization------------------------
	// ------------------------------------------------------------------
  	vector <symbol_list_> S_train = R_train->getlist_vector();
  	vector <symbol_list_> S_test = R_test->getlist_vector();
    
	//llk coordinates of train and test;
	vector<vector<double>> C_train(S_train.size(), vector<double>(PFSA_vec.size(), 0));	
	vector<vector<double>> C_test(S_test.size(), vector<double>(PFSA_vec.size(), 0));
	for(size_t i = 0; i < S_train.size(); i++)
		for(size_t j = 0; j < PFSA_vec.size(); j++)
			C_train[i][j] = PFSA_vec[j].log_likelihood(S_train[i]);
	
	for(size_t i = 0; i < S_test.size(); i++)
		for(size_t j = 0; j < PFSA_vec.size(); j++)
			C_test[i][j] = PFSA_vec[j].log_likelihood(S_test[i]);
	// ------------------------------------------------------------------
	// ------------------------------------------------------------------



	// --------------------Get distance and class------------------------
	// ------------------------------------------------------------------
	// Distance matrix is of dimension |test| \times |number_of_class|
	// vector<vector<double>> Dist(C_test.size(), vector<double>(C_train.size(), 0)); 
	size_t correct = 0;
	map<int, double> class_distance;
	for (size_t i = 0; i < C_test.size(); i++)
	{
		int test_label = labels_test[i];
		for(size_t j = 0; j < C_train.size(); j++)
		{
			int l = labels_train[j];
			double d = 0;
			for(size_t k = 0; k < PFSA_vec.size(); k++)
			{	
				d += fabs(C_test[i][k] - C_train[j][k]);
			} 
			class_distance[l] += d;
		}
		map<int, double>::iterator itr;
		double min_dist = 1e20; 
		int predict_label;
		for (itr = class_distance.begin(); itr != class_distance.end(); ++itr)
		{
			int l = itr->first;
			double d = itr->second;
			if (d < min_dist)
			{
				min_dist = d;
				predict_label = l;
			}
		}
		class_distance.clear();
		if (predict_label == test_label)
			correct += 1;
	}
	double correct_rate = (double)correct / (double)C_test.size();
	cout << "number of classes = " << class_labels.size() << endl;
	cout << "correct_rate = " << correct_rate << endl;
	
	// ------------------------------------------------------------------
	// ------------------------------------------------------------------
	
	ofstream coor_train(coor_train_file.c_str());
	coor_train << C_train;
	coor_train.close();
	
	ofstream coor_test(coor_test_file.c_str());
	coor_test << C_test;
	coor_test.close();
  
  	return 0;
}

