#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <tuple>

#include <dirent.h>

#include <gsl/gsl_statistics.h>
#include "semantic.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

const string VERSION="\nQuantizer with log-likelihood v1.3 \nSimply distance matrix!";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

const bool DEBUG=false;

vector<PFSA> read_PFSA(string PFSA_dir)
{
	DIR *pDir;
	pDir = opendir(PFSA_dir.c_str());

	if (pDir == NULL)
	{
		cout << "Cannot open " << PFSA_dir << "!";
		exit(1);
	}

	vector<PFSA> PFSA_vec;
	struct dirent *pDirent;

	while ((pDirent = readdir(pDir)) != NULL)
	{
		string filename(pDirent->d_name);
		if (filename.compare(".") * filename.compare("..") != 0)
		{
			string PFSA_file = PFSA_dir + "/" + filename;
			PFSA G = SCC_UTIL__::read_mc(PFSA_file, "PFSA");
			PFSA_vec.push_back(G);
			
			if (DEBUG)
			{	
				cout << PFSA_file << endl;
  				vector<double> stationary_dist(G.get_Stationary());
      			cout << "Stationary Distribution: " << stationary_dist << endl;
      			
				connx delta(G.get_aut());
				unsigned int num_states = delta.size();
				unsigned int alphabet_size  = delta[1].size();

				cout << "Transition Matrix: " << endl;
      			for (size_t i = 0; i < num_states; i++)
				{
					cout << "\t";
					for (symbol s(0); s < alphabet_size; s++)
					{
						cout << delta[i][s] << " ";
					}
					cout << endl;
				}

  				vector<vector<double>>  pitilde(G.get_PI());
				cout << "Transition Probability Matrix: " << endl;
      			for (size_t i = 0; i < num_states; i++)
				{
					cout << "\t" << pitilde[i] << endl;
				}	
				cout << endl;
			}
		}
	}

	closedir(pDir);

	return PFSA_vec;
}

double get_log_likelihood(PFSA& G, 
						const vector<double>& stationary, 
						const symbol_list_& seq)
{
  	if(seq.empty())
	{
    	return 0.0;
	}
  
  	double llk = 0.0;
  	size_t num_states = G.get_aut().size();
 
 	vector<double> curr_distr(stationary);

  	for(symbol s : seq)
    {
      	double pr = 0.0;
      	for(unsigned int st = 0; st < num_states; st++)
		{
			pr += curr_distr[st] * G.get_pit()[st][s];
		}

      	llk += log2(pr);

      	vector <double> nxt_distr(num_states, 0.0);
	  	double sum = 0.0;

      	for (size_t st = 0; st < num_states; st++)
		{
			int nxt_state = G.get_aut()[st][s];
			
			if (nxt_state != -1)
			{
				double temp = curr_distr[st] * G.get_pit()[st][s];
				nxt_distr[nxt_state] += temp; 
				sum += temp;
			}
		}
      
      	for (size_t st = 0; st < num_states; st++)
		{
			curr_distr[st] = nxt_distr[st] / sum;	
		}
	}
	return -llk / seq.size();
}


double total_variation_distance(const vector<double>& vec1, const vector<double>& vec2)
{
	double distance = 0;
	for (size_t i = 0; i < vec1.size(); i++)
	{
		distance += abs(vec1[i] - vec2[i]);
	}
	return distance;
}


vector<vector<double>> get_llk_distance_matrix(vector<PFSA>& PFSA_vec, const vector<symbol_list_>& seqs)
{
	vector<vector<double>> stationary_distrs;
	for (PFSA G: PFSA_vec)
	{
		stationary_distrs.push_back(G.get_Stationary());
	}

	vector<vector<double>> llk_matrix;
	for (size_t i = 0; i < seqs.size(); i++)
	{
		vector<double> llk_vec;
		for (size_t j = 0; j < PFSA_vec.size(); j++)
		{
			llk_vec.push_back(get_log_likelihood(PFSA_vec[j], stationary_distrs[j], seqs[i]));
		}
		llk_matrix.push_back(llk_vec);
	}
	
	vector<vector<double>> distance_matrix;
	for (size_t i = 0; i < seqs.size(); i++)
	{
		vector<double> distance_vec;
		for(size_t j = 0; j < seqs.size(); j++)
		{
			distance_vec.push_back(total_variation_distance(llk_matrix[i], llk_matrix[j]));
		}
		distance_matrix.push_back(distance_vec);
	}
	return distance_matrix;
}


unsigned int get_rand_seed()
{
  chrono::time_point<chrono::system_clock> now = chrono::system_clock::now();
  unsigned long long seed = (unsigned long int)chrono::duration_cast<chrono::nanoseconds>(now.time_since_epoch()).count();		
  return (unsigned int)seed % INT_MAX;
}

double random_key()
{
  srand(get_rand_seed());
  double key = (double)rand() / RAND_MAX;
  return key;
}

unsigned int random_int(unsigned int N)
{
  srand(get_rand_seed());
  return rand() % N;
}

vector<unsigned int> choose(unsigned int N, unsigned int k)
{
  if (k > N)
    {
      cerr << "ERROR: there is no way to choose "
	   << k << " number from a set of size " << N << ".\n";
      exit(1);
    }
  set<unsigned int> indices;
  for (size_t i = 0; i < N; i++)
    {
      indices.insert(i);
    }
  vector<unsigned int> chosen;
  for (size_t i = 0; i < k; i++)
    {
      set<unsigned int>::iterator itr = indices.begin();
      unsigned int offset = random_int(indices.size());
      advance(itr, offset);
      unsigned int index = *itr;
      indices.erase(itr);
      chosen.push_back(index);
    }
  sort(chosen.begin(), chosen.end());
  return chosen;
}

void detrend_inplace(vector<double> & vec)
{
  for (size_t i = 0; i < vec.size() - 1; i++)
    {
      vec[i] = vec[i + 1] - vec[i];
    }
  vec.pop_back();
}

void detrend_inplace(vector<double> & vec, unsigned int d)
{
  for (size_t i = 0; i < d; i++)
    {
      detrend_inplace(vec);
    }
}

void detrend_inplace(vector<vector<double>> & vvec, unsigned int d)
{
  for (size_t s = 0; s < vvec.size(); s++)
    {
      detrend_inplace(vvec[s], d);
    }
}

vector<double> detrend(const vector<double> & vec, unsigned int d)
{
  vector<double> new_vec(vec);
  detrend_inplace(new_vec, d);
  return new_vec;
}

vector<vector<double>> detrend(const vector<vector<double>> & vvec, unsigned int d)
{
  vector<vector<double>> new_vvec(vvec);
  detrend_inplace(new_vvec, d);
  return new_vvec;
}

void normalize_inplace(vector<double> & vec)
{
  double mean = gsl_stats_mean(&vec[0], 1, vec.size());
  double sd = gsl_stats_sd(&vec[0], 1, vec.size());
  for (size_t s = 0; s < vec.size(); s++)
    {
      vec[s] = (vec[s] - mean) / sd;
    }
}

void normalize_inplace(vector<vector<double>> & vvec)
{
  for (size_t s = 0; s < vvec.size(); s++)
    {
      normalize_inplace(vvec[s]); 
    }
}

vector<double> normalize(const vector<double> vec)
{
  vector<double> new_vec(vec);
  normalize_inplace(new_vec);
  return new_vec;
}

vector<vector<double>> normalize(const vector<vector<double>> & vvec)
{
  vector<vector<double>> new_vvec(vvec);
  normalize(new_vvec);
  return new_vvec;
}

pair<bool, pair<double, double>> get_freq_val_range(const vector<double> & vec, double _A, double _B)
{
  unsigned int _C = 10;
  unsigned int num_bins = (unsigned int)(_C / _A) + 1;
	
  unsigned int size = vec.size();
  double min = vec.front();
  double max = vec.front();

  for (size_t i = 1; i < size; i++)
    {
      double cur = vec[i];
      min = (cur < min) ? cur : min;
      max = (cur > max) ? cur : max;
    }
	
  double bin_width = (max - min) / num_bins;

  typedef pair<unsigned int, unsigned int> bin;
  vector<bin> bins;
  for (size_t i = 0; i < num_bins; i++)
    {
      bins.push_back(bin(i, 0));
    }
  for (size_t i = 0; i < size; i++)
    {
      unsigned int idx = int(abs(vec[i] - min) / bin_width);
      if (idx == num_bins)
	{
	  idx--;
	}
      bins[idx].second += 1;
    }

  struct
  {
    bool operator()(bin a, bin b) const
    {
      return a.second > b.second ;
    }
  }bin_size_comparator;

  sort(bins.begin(), bins.end(), bin_size_comparator);

  unsigned int min_bin = bins[0].first;
  unsigned int max_bin = min_bin;
  unsigned int bin_sum = bins[0].second;
  bool to_prune = false;
  for (size_t i = 1; i < _C; i++)
    {
      if ((double)bin_sum / size >= _B)
	{
	  to_prune = true;
	  break;
	}
      else
	{
	  unsigned int cur_bin = bins[i].first;
	  if ((cur_bin < min_bin + _C) && (max_bin < cur_bin + _C))
	    {
	      min_bin = (cur_bin < min_bin) ? cur_bin : min_bin;
	      max_bin = (cur_bin > max_bin) ? cur_bin : max_bin;
	      bin_sum += bins[i].second;
	    }
	  else
	    {
	      break;
	    }
	}
    }
  if (to_prune)
    {
      double lower = min + min_bin * bin_width;
      double upper = min + (max_bin  + 1) * bin_width;
      cout << "lower = " << lower <<endl;
      cout << "upper = " << upper <<endl;
      cout << "concentration = " << (double)bin_sum / size <<endl;
      pair<double, double> p(lower, upper);
      return pair<bool, pair<double, double>>(true, p);
    }	
  else
    {
      pair<double, double> p(0, 0);
      return pair<bool, pair<double, double>>(false, p);
    }
}


pair<bool, pair<double, double>> get_freq_val_range(const vector<vector<double>> & vvec, double _A, double _B)
{
  vector<double> vec;
  for (size_t i = 0; i < vvec.size(); i++)
    {
      vec.insert(vec.end(), vvec[i].begin(), vvec[i].end());
    }
  return get_freq_val_range(vec, _A, _B);
}

vector<double> remove_freq_val(const vector<double> & vec, pair<double, double> range, double connectiveness)
{
  unsigned int size = vec.size();
  typedef	pair<unsigned int, unsigned int> segment;
  vector<segment> segments;
  bool is_recording = false;
  for (size_t i = 0; i < size; i++)
    {
      if (vec[i] >= range.first && vec[i] < range.second)
	{
	  if (!is_recording)
	    {
	      segment seg(i, i);
	      segments.push_back(seg);
	      is_recording = true;
	    }
	}
      else
	{
	  if(is_recording)
	    {
	      segments[segments.size() - 1].second = i;
	      is_recording = false;
	    }
	}
    }
  if(is_recording)
    {
      segments[segments.size() - 1].second = size;
    }

  struct
  {
    bool operator()(segment s, segment t) const
    {
      return (s.second - s.first) > (t.second - t.first);
    }
  }size_comparator;
  sort(segments.begin(), segments.end(), size_comparator);
	
  for (size_t i = 0; i < segments.size(); i++)
    {
      if ((segments[i].second - segments[i].first) < vec.size() * connectiveness)
	{
	  segments.resize(i);
	  break;
	}
    }

  struct
  {
    bool operator()(segment s, segment t) const
    {
      return s.first < t.first;
    }
  }start_comparator;
  sort(segments.begin(), segments.end(), start_comparator);
	
  vector<double> pruned;
  int start = 0;
  for (size_t i = 0; i < segments.size(); i++)
    {
      for (size_t j = start; j < segments[i].first; j++)
	{
	  pruned.push_back(vec[j]);
	}
      start = segments[i].second;
    }
  for (size_t j = start; j < size; j++)
    {
      pruned.push_back(vec[j]);
    }
  return pruned;
}

vector<vector<double>> remove_freq_val(const vector<vector<double>> vvec, pair<double, double> range, double connectiveness)
{
  vector<vector<double>> vpruned;
  for (size_t i = 0; i < vvec.size(); i++)
    {
      vector<double> pruned = remove_freq_val(vvec[i], range, connectiveness);
      vpruned.push_back(pruned);
    }
  return vpruned;
}

pair<bool, pair<double, double>> prune_inplace(vector<double> & vec, double _A, double _B, double connectiveness)
{
  pair<bool, pair<double, double>> p = get_freq_val_range(vec, _A, _B);
  if (p.first)
    {
      vec = remove_freq_val(vec, p.second, connectiveness);
    }
  return p;
}

pair<bool, pair<double, double>> prune_inplace(vector<vector<double>> & vvec, double _A, double _B, double connectiveness)
{
  pair<bool, pair<double, double>> p = get_freq_val_range(vvec, _A, _B);
  if (p.first)
    {
      vvec = remove_freq_val(vvec, p.second, connectiveness);
    }
  return p;
}

vector<double> remove_freq_val(const vector<double> & vec, pair<double, double> range)
{
  vector<double> pruned;
  for (double v : vec)
    {
      if (v < range.first || v > range.second)
	{
	  pruned.push_back(v);
	}
    }
  return pruned;
}

vector<vector<double>> remove_freq_val(const vector<vector<double>> & vvec, pair<double, double> range)
{
  vector<vector<double>> vpruned;
  for (size_t i = 0; i < vvec.size(); i++)
    {
      vector<double> pruned = remove_freq_val(vvec[i], range);
      vpruned.push_back(pruned);
    }
  return vpruned;
}

class epsilonComp
{
  double epsilon;
public:
  epsilonComp() :epsilon(.1) {};
  epsilonComp(double e) :epsilon(e) {};
  bool operator()(double l, double r) const
  {
    double d = abs(l - r);
    return (l < r) && (d > epsilon);
  }
};

vector<unsigned int> get_counts(const vector<double> & data, double epsilon)
{
  map<double, unsigned int, epsilonComp> countMap(epsilonComp((const double)epsilon));
  for (size_t i = 0; i < data.size(); i++)
    {
      countMap[data[i]]++;
    }
	
  map<double, unsigned int, epsilonComp>::iterator itr;
  vector<unsigned int> counts;
  for (itr = countMap.begin(); itr != countMap.end(); ++itr)
    {
      unsigned int count = itr->second;
      counts.push_back(count);
    }
	
  return counts;
}

// Function to calculate entropy from quantile vector
double entropy_quantile(const vector<double> & quantile)
{
  unsigned int dim = quantile.size() + 1;
  double lower = 0;
  double ent = 0;
  for (size_t i = 0; i < dim; i++)
    {
      double upper = (i == dim - 1) ? 1. : quantile[i];
      double p = upper - lower;
      ent += -log(p) * p;
      lower = quantile[i];
    }
  return ent;
}

// Function to calculate entropy from quantile vector and a cut vector of it
// For example the quantile vector is (.1, .4, .6, .8)
// and the cut is (0, 2), 
// then we calculate the entropy of the quantile vector (.1, .6)
double entropy_quantile(const vector<double> quantile, const vector<unsigned int> cut)
{
  double lower = 0;
  double ent = 0;
  for (size_t i = 0; i < cut.size(); i++)
    {
      double p = quantile[cut[i]] - lower;
      ent += -log(p) * p;
      lower = quantile[cut[i]];
    }
  ent += -log(1. - lower) * (1. - lower);
  return ent;
}

pair<bool, vector<double>> get_partition(
					 const vector<double> & data, 
					 const vector<double> & quantile, 
					 double epsilon)
{
  vector<double> partition;
	
  unsigned int alphabetSize = quantile.size() + 1;
  // Get the epsilon-counts from data
  vector<unsigned int> counts = get_counts(data, epsilon);
  if (counts.size() < alphabetSize)
    {
      return pair<bool, vector<double>>(false, partition);
    }

  // Calculate the cumulative fractions
  vector<unsigned int> partial_sums;
  unsigned int ps = 0;
  for (size_t i = 0; i < counts.size(); i++)
    {	
      partial_sums.push_back(ps + counts[i]);
      ps = partial_sums[i];
    }
  vector<double> fractions;
  for (size_t i = 0; i < counts.size() - 1; i++)
    {
      double f = partial_sums[i] / (double)ps;
      fractions.push_back(f);
    }
	
  // Find the cuts to the cumulative fractions
  // to get ones that are the most similar to 
  // the given quantile.
  vector<pair<int, int>> bounds;
  int j = 0;
  for (size_t i = 0; i < quantile.size(); i++)
    {
      while (((unsigned int)j < fractions.size()) & (fractions[j] < quantile[i]))
	{
	  j++;	
	}
      pair<int, int> bound_pair(j - 1, j);
      bounds.push_back(bound_pair);
    }
	
  bool is_successful = true;	
  std::queue<vector<unsigned int>> cuts;
  for (size_t i = 0; i < quantile.size(); i++)
    {
      int lower = bounds[i].first;
      int upper = bounds[i].second;
      unsigned int cur_size = cuts.size();
      if (cur_size > 0)
	{
	  for (size_t j = 0; j < cur_size; j++)
	    {
	      vector<unsigned int> cut = cuts.front();
	      cuts.pop();
	      if ((int)cut.back() < lower){
		vector<unsigned int> new_cut(cut);
		new_cut.push_back(lower);
		cuts.push(new_cut);
	      }
	      if (((int)cut.back() < upper) && ((unsigned int)upper < fractions.size()))
		{
		  vector<unsigned int> new_cut(cut);
		  new_cut.push_back(upper);
		  cuts.push(new_cut);
		}
	    }
	}
      else
	{
	  if (lower >= 0)
	    {
	      vector<unsigned int> new_cut(1, lower);
	      cuts.push(new_cut);
	    }
	  if ((unsigned int)upper < fractions.size())
	    {
	      vector<unsigned int> new_cut(1, upper);
	      cuts.push(new_cut);
	    }
	}
      if (cuts.size() == 0)
	{
	  is_successful = false;
	  break;
	}
    }
	
  if (is_successful)
    {
      vector<unsigned int> bestCut;
      double diff = log(alphabetSize);
      double target_ent = entropy_quantile(quantile);
      while (!cuts.empty())
	{
	  vector<unsigned int> cut = cuts.front();
	  cuts.pop();
	  double ent = entropy_quantile(fractions, cut);
	  if (abs(ent - target_ent) < diff)
	    {
	      bestCut = cut;
	      diff = abs(ent - target_ent);
	    }
	}
      for (size_t c : bestCut)
	{
	  partition.push_back(data[partial_sums[c]]);	
	}
    }
  return pair<bool, vector<double>>(is_successful, partition);
}

void print_vec_with_bracket(vector<double> vec)
{
  cout << "[";
  for (size_t i = 0; i < vec.size(); i++)
    {
      printf("%.4f", vec[i]);
      (i < vec.size() - 1)? cout << ", " : cout << "]\n";
    }
}

vector<vector<double>> get_data(string data_filename) 
{
	string line;
	vector<vector<double>> data;
	ifstream data_file(data_filename);	

	while (getline(data_file, line))
	{
		stringstream ss;
		ss << line;
		vector<double> datum;
		double data_point;
	
		while (ss >> data_point)
		{
			datum.push_back(data_point);
		}
	
		data.push_back(datum);
	}

	data_file.close();

	return data;
}


void save_data(string filename, const vector<vector<double>> & data)
{
	ofstream file(filename);
    for (size_t i = 0; i < data.size(); i++)
    {
		for (size_t j = 0; j < data[i].size(); j++)
		{
			file << data[i][j];
			(j < data[i].size() - 1) ? file << " " : file << "\n";
		}
    }
	file.close();
}


double highlight_one_entropy(unsigned int alphabetSize, double p)
{
  return -p * log(p) - (1 - p) * log((1 - p) / (alphabetSize - 1.));
}

double bin_search(unsigned int alphabetSize, double entropy_level)
{
  double max_ent = log(alphabetSize);
  double target_ent = entropy_level * max_ent;
	
  double lower = 1. / (double)alphabetSize;
  double upper = 1.;
  double middle;
  while (true)
    {
      middle = (lower + upper) / 2.;
      double ent = highlight_one_entropy(alphabetSize, middle);
      if (abs(ent - target_ent) < 1e-4)
	{
	  break;
	}
      else
	{
	  if (ent > target_ent)
	    {
	      lower = middle;
	    }
	  else
	    {
	      upper = middle;
	    }
	}
    }	
  return middle;
}

vector<vector<double>> get_quantiles(
				     unsigned int ASLow, 
				     unsigned int ASHigh, 
				     double entLow, 
				     double entHigh, 
				     unsigned int num)
{
  vector<vector<double>> quantiles;
  double step_size = (entHigh - entLow) / (num - 1.);
	
  for (size_t as = ASLow; as <= ASHigh; as++)
    {
      for (size_t n = 0; n < num; n++)
	{
	  double entropy_level = entLow + step_size * (double)n;
	  if (abs(entropy_level - 1.) < 1e-3)
	    {
	      vector<double> evenly;
	      double quant = 0;
	      for (size_t i = 0; i < as - 1; i++)
		{
		  quant += 1. / as;
		  evenly.push_back(quant);
		}
	      quantiles.push_back(evenly);
	      break;
	    }
	  double p = bin_search(as, entropy_level);
			
	  double q = (1. - p) / (as - 1.);
	  for (size_t i = 0; i < as; i++)
	    {
	      vector<double> quantile;
	      double quant = 0;
	      for (size_t j = 0; j < as - 1; j++)
		{
		  quant += (j == i) ? p : q;
		  quantile.push_back(quant);
		}
	      quantiles.push_back(quantile);	
	    }
	}	
    }
  return quantiles;
}

const double EPSILON = 1e-2;

class Quantizer
{
protected:
  	// TRAIN + TEST
  	vector<vector<double>> data;

	// Model PFSA
	vector<PFSA> PFSA_vec;

  	// The parameter for epsilonComp class.
  	// Roughly speaking, a pair of number that don't differ by epsilon will be consider the same. 
  	double epsilon;
  	void _partition(const vector<vector<double>> & quantiles, bool verbose);
  	pair<bool, vector<double>> _partition(const vector<double> & quantile, const vector<double> & pooledData, bool verbose);	
	
  	// Convert a continuous data point to a symbol according given partition vector.
  	symbol _get_symbol(double dataPoint, const vector<double> & partition);
  	// Convert continuous data streams to symbolic ones using partition got.
  	pair<bool, vector<symbol_list_>> _fit(const vector<double> & partition);
  	
	// Calculate the distance matrix, sae, and discrimination of the quantization.
  	void _calc_statistics();
  	pair<bool, vector<vector<double>>> _calc_statistics(const vector<double> & partition);

public:
  	map<unsigned int, pair<bool, vector<double>>> partitions;
  	map<unsigned int, pair<bool, vector<vector<double>>>> dist_matrices;
  	// Constructor:
  	Quantizer(
		vector<vector<double>> data,
  		string PFSA_directory,
		double epsilon,
		vector<vector<double>> quantiles,
		bool verbose);
};


void Quantizer::_partition(const vector<vector<double>> & quantiles, bool verbose)
{
  	vector<double> pooledData;
  	for (size_t n = 0; n < data.size(); n++)
    {
      	pooledData.insert(pooledData.end(), data[n].begin(), data[n].end());
    }
  	sort(pooledData.begin(), pooledData.end());

  	map<vector<double>, bool> partition_existence_map;	
  	map<vector<double>, bool>::iterator itr;
  	vector<double> empty;	
  	for (size_t i = 0; i < quantiles.size(); i++)
    {
      	pair<bool, vector<double>> partition_pair = _partition(quantiles[i], pooledData, verbose);
      	if (partition_pair.first)
		{
	  		itr = partition_existence_map.find(partition_pair.second);
	  		if (itr != partition_existence_map.end())
	    	{
	      		partitions[i] = pair<bool, vector<double>>(false, empty);
	    	}
	  		else
	    	{
	      		partition_existence_map[partition_pair.second] = true;
	      		partitions[i] = partition_pair;
	    	}
		}
    }	
}

pair<bool, vector<double>> Quantizer::_partition(const vector<double> & quantile, const vector<double> & pooledData, bool verbose)
{
  	bool is_successful = true;
  	vector<double> partition;

  	if (epsilon < 0)
    {
      	for (size_t i = 0; i < quantile.size(); i++)
		{
	  		double p = gsl_stats_quantile_from_sorted_data(&pooledData[0], 1, pooledData.size(), quantile[i]);
	  		if (partition.size() != 0 && p == partition.back())
	    	{
	      		is_successful = false;
	      		break;
	    	}
	  		partition.push_back(p);
		}
      	if (!is_successful)
		{
	  		if (verbose)
	    	{
	      		cout << "Data has to grouped to find partition!" << endl;
	    	}
	  		pair<bool, vector<double>> partition_pair = get_partition(pooledData, quantile, EPSILON);
			
	  		if (partition_pair.first)
	    	{
	      		partition = partition_pair.second;
	      		is_successful = true;
	    	}
		}	
    }
  	else
    {
      	pair<bool, vector<double>> partition_pair = get_partition(pooledData, quantile, epsilon);
      	if (partition_pair.first)
		{
	  		partition = partition_pair.second;
		}
      	else
		{
	  		is_successful = false;
		}	
    }
  	return pair<bool, vector<double>>(is_successful, partition);
}


symbol Quantizer::_get_symbol(double dataPoint, const vector<double> & partition)
{
  	symbol sym(0);
  	while (sym < partition.size() && dataPoint >= partition[sym])
    {
      	sym++;
    }
  	return sym;
}

pair<bool, vector<symbol_list_>> Quantizer::_fit(const vector<double> & partition)
{
  vector<symbol_list_> quantizedData;
  unsigned int alphabetSize = partition.size() + 1;
  set<symbol> full;
  for (symbol sym(0); sym < alphabetSize; sym++)
    {
      full.insert(sym);
    }
  bool is_successful = true;
  for (size_t n = 0; n < data.size(); n++)
    {
      set<symbol> scratch(full);
      symbol_list_ quantizedDatum;
      for (size_t i = 0; i < data[n].size(); i++)
	{
	  symbol sym = _get_symbol(data[n][i], partition);
	  quantizedDatum.push_back(sym);
	  scratch.erase(sym);
	}
      if (scratch.size() != 0)
	{
	  is_successful = false;
	  quantizedData.clear();
	  break;
	}
      else
	{
	  quantizedData.push_back(quantizedDatum);
	}
    }
  return pair<bool, vector<symbol_list_>>(is_successful, quantizedData);
}


void Quantizer::_calc_statistics()
{
  	for (size_t i = 0; i < partitions.size(); i++)
  	{
      	if (partitions[i].first)
		{
	  		dist_matrices[i] = _calc_statistics(partitions[i].second);
		}
      	else
		{
			vector<vector<double>> empty;
	  		dist_matrices[i] = make_pair(false, empty);
		}
    }
}


pair<bool, vector<vector<double>>> Quantizer::_calc_statistics(const vector<double> & partition)
{
	pair<bool, vector<symbol_list_>> quantized_pair = _fit(partition);
	
	if (quantized_pair.first)
    {
      	vector<vector<double>> dist_matrix = get_llk_distance_matrix(PFSA_vec, quantized_pair.second);
		return make_pair(true, dist_matrix);
   	}	
  	else
	{
		vector<vector<double>> dist_matrix;
		return make_pair(false, dist_matrix);
	}
}


Quantizer::Quantizer(
	vector<vector<double>> data,
	string PFSA_directory,
	double epsilon,
	vector<vector<double>> quantiles,
	bool verbose)
{
  	this->data = data;

	this->PFSA_vec = read_PFSA(PFSA_directory);

  	this->epsilon = epsilon;

  	_partition(quantiles, verbose);
  	_calc_statistics();
}

int main(int argc, char** argv)
{ 
	try
    {	
		string data_filename;
		string PFSA_directory;	
		bool pruned;
		unsigned int detrend;
		bool normalized;
		double epsilon;
		unsigned int verbose;

		options_description infor( "Program information");
		infor.add_options()
		("help,h", "print help message.")
		("version,V", "print version number");
		
		options_description usg("Allowed options");
		usg.add_options()
    	("data-filename,D", value<string>(&data_filename), "filename of the dataset")
    	("pfsa-directory,P", value<string>(&PFSA_directory), "directory of the model PFSA")
    	("pruned,r", value<bool>(&pruned), "Whether to pruned out frequent values. Try both options when not specified")
    	("detrended,d", value<unsigned int>(&detrend), "The number of times the raw data is differentiated. Choose between 0 and 1 if not specified")
    	("normalized,n", value<bool>(&normalized), "Whether to normalized the data. Try both options if not specified")
    	("epsilon,e", value<double>(&epsilon)->default_value(-1), "-1 for not grouping the data for partition, a positive value gives the resolution for grouping.")
    	("verbose,v", value<unsigned int>(&verbose)->default_value(0), "0 for displaying nothing, 1 for displaying all information.");

      	options_description desc( "\n### Log-likelihood Quantizer ### (yi huang 2019 June)");
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
      		store(command_line_parser(argc, argv).options(desc).run(), vm);
			notify(vm);
    	} 
      	catch (std::exception &e)
    	{
      		cout << endl << e.what() << endl << desc << endl;
      		return 1;
    	}

      	notify(vm);
    	
      	if (vm.count("version"))
    	{
      		cout << VERSION << endl; 
      		return 0;
    	}
      	if (vm.count("help"))
    	{
      		cout << "Usage: options description\n";
      		cout << desc;
      		return 0;
    	}
    	
      	vector<bool> _pruned;
      	if (vm.count("pruned"))
    	{
      		_pruned.push_back(vm["pruned"].as<bool>());
    	}
      	else
    	{
      		_pruned.push_back(false);
      		_pruned.push_back(true);
    	}
    	
      	vector<bool> _detrended;
      	if (vm.count("detrended"))
    	{
      		_detrended.push_back(vm["detrended"].as<unsigned int>());
    	}
      	else
    	{
      		_detrended.push_back(0);
      		_detrended.push_back(1);
    	}
    	
      	vector<bool> _normalized;
      	if (vm.count("normalized"))
    	{
      		_normalized.push_back(vm["normalized"].as<bool>());
    	}
      	else
    	{
      		_normalized.push_back(false);
      		_normalized.push_back(true);
    	}
    	
    	
      	vector<vector<double>> data = get_data(data_filename);

		// with 90%, 95%, and 100% of the maximum entropy
      	vector<vector<double>> quantiles = get_quantiles(2, 2, .90, 1., 3);

      	// Figure out whether pruning is needed. 
      	pair<bool, pair<double, double>> prune_pair;
      	if ((_pruned.size() == 2) || ((_pruned.size() == 1) & _pruned[0]))
    	{
      		prune_pair = get_freq_val_range(data, .05, .6);	
      		if (prune_pair.first == false)
        	{
          		_pruned = vector<bool>(1, false);
          		if (verbose)
    			{
    	  			cout << "Pruning is requested but not necessary!" << endl;
    			}
        	}
    	}

      	map<unsigned int, string> parameter_map;
      	map<unsigned int, vector<vector<double>>> dist_map;
      
		unsigned int count = 0;
      	for (bool r : _pruned)
    	{
      		for (size_t d : _detrended)
        	{
          		for (bool n : _normalized)
    			{
    	  			if (verbose)
    	    		{
    	      			cout << "pruned = " << r
    		   				<< "\tdetrended = " << d
    		   				<< "\tnormalized = " << n << endl;
    	    		}
    	  			
					string prefix;	
					vector<vector<double>> processed_data;
    	  			if (r)
    	    		{
    	      			processed_data = remove_freq_val(data, prune_pair.second);
    	      			double lower = prune_pair.second.first;
    	      			double upper = prune_pair.second.second;
    	      			prefix = "R[" + to_string(lower) + " " + to_string(upper) + "]";
    	    		}
    	  			else
    	    		{
    	      			processed_data = data;
    	    		}
    	  			detrend_inplace(processed_data, d);
    	  			if (n)
    	    		{
    	      			normalize_inplace(processed_data);
    	    		}
    	  			prefix += "D" + to_string(d)+ "N" + to_string(n); 	
    				
    	  			Quantizer q(processed_data, PFSA_directory, epsilon, quantiles, verbose);
    				
    	  			for (size_t i = 0; i < q.dist_matrices.size(); i++)
    	    		{
    	      			if (q.dist_matrices[i].first)
    					{
    		  				string parameter = prefix + "[";
    		  				vector<double> partition = q.partitions[i].second;
    		  				for (size_t j = 0; j < partition.size(); j++)
    		    			{
    		      				parameter += to_string(partition[j]);
    		      				if (j < partition.size() - 1)
    							{
    			  					parameter += " ";
    							}
    		      				else
    							{
    			  					parameter += "]";
    							}
    		    			}
    						
    		  				parameter_map[count] = parameter;
    		  				dist_map[count] = q.dist_matrices[i].second;
							
							if(verbose)
							{
								cout << count << ": " << parameter << endl;
							}					

    		  				count++;
    					}
    	      			else
    					{
    		  				if (verbose)
    		    			{
    		      				cout << "\t Missing symbol" << endl;
    		    			}
    					}
    	    		}			
    			}
        	}
    	}
      	
		if (verbose)
    	{
      		cout << "Found " << count << " valid parameters." << endl;
    	}

      	ofstream parameter_file("valid_parameter");
      	for (size_t i = 0; i < parameter_map.size(); i++)
    	{
      		parameter_file << parameter_map[i] << endl;
          	save_data(parameter_map[i] + "_dist", dist_map[i]);
    	}
      	parameter_file.close();	
    }
  	catch (std::exception & err)
    {
      	cerr << err.what() << endl;
      	return 1;
    }

  	return 0;
}
