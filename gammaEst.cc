/* 
 * Estimate gamma, the coefficient causal dependence, without models
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

// #include "semantic.h"
// #include "config.h"

using namespace std;

typedef map<unsigned int, vector<unsigned int>> DataPhi;

vector<unsigned int> compose(
	const vector<unsigned int>& input, 
	const vector<unsigned int>& output,
	size_t ab_out)
{
 	size_t length = (input.size() < output.size())? input.size() : output.size();
	
	vector<unsigned int> composed;
	for (size_t i = 0; i < length; i++)
	{	
		composed.push_back(input[i] * ab_out + output[i]);
	}
	return composed;
}


double entropy(const vector<double>& prob)
{
  	double e = 0.;
  	for (double p : prob)
    	if (p > 0.)
      		e -= p * log2(p);
  	return e;
}

double entropy(const vector<unsigned int>& count)
{
	unsigned int sum = 0;
	double entropy = 0.;
	for (unsigned int c : count)
	{
		sum += c;
		if (c > 0)
			entropy -= c * log2(c);
	}
	return entropy / sum + log2(sum);
}

DataPhi get_count(
	const vector<unsigned int>& seq, 
	size_t alphabet,
	size_t seq_len)
{
	size_t length = seq.size();

	vector<unsigned int> init;
	for (unsigned int i = 0; i < length; i++)
	{
		init.push_back(i);
	}

	DataPhi odp; 
	odp[0] = init;
	
	DataPhi::iterator itr;

	for (unsigned int sl = 0; sl < seq_len; sl++)
	{
		DataPhi ndp;
		
		for (itr = odp.begin(); itr != odp.end(); ++itr)
		{
			unsigned int wordSig = itr->first;
			for (unsigned int index : itr->second)
			{
				unsigned int i = index + sl;
				if (i < length)
				{
					ndp[wordSig * alphabet + seq[i]].push_back(index);
				}
			}
		}
		odp = ndp;
	}
	
	DataPhi counts;
	for (itr = odp.begin(); itr != odp.end(); ++itr)
	{	
		vector<unsigned int> count(alphabet, 0);

		for (unsigned int index : itr->second)
		{
			unsigned int i = index + seq_len;
			if (i < length)
			{
				count[seq[i]] += 1;
			}
		}
		counts[itr->first] = count;
	}
	return counts;
}


double get_entropy( 
	const vector<unsigned int>& seq, 
	size_t alphabet,
	size_t seq_len)
{
	DataPhi counts = get_count(seq, alphabet, seq_len);

	double entropy_avg = 0.;
	
	DataPhi::iterator itr;
	for (itr = counts.begin(); itr != counts.end(); ++itr)
	{
		unsigned int sum = 0;
		double entropy = 0.;
		for (unsigned int count: itr->second)
		{
			sum += count;
			if (count > 0)
				entropy -= count * log2(count);
		}
		entropy = entropy / sum + log2(sum);
		entropy_avg += entropy * sum;
	}
	return entropy_avg / (seq.size() - seq_len);
}

pair<double, double> get_proj_entropy( 
	const vector<unsigned int>& seq_first, 
	const vector<unsigned int>& seq_second, 
	size_t alphabet_first,
	size_t alphabet_second,
	size_t seq_len)
{
	vector<unsigned int> seq = compose(seq_first, seq_second, alphabet_second);
	size_t alphabet = alphabet_first * alphabet_second;
	DataPhi counts = get_count(seq, alphabet, seq_len);
	

	double entropy_avg_first = 0.;
	double entropy_avg_second = 0.;

	DataPhi::iterator itr;
	for (itr = counts.begin(); itr != counts.end(); ++itr)
	{
		vector<unsigned int> vec = itr->second;
		vector<unsigned int> vec_first(alphabet_first, 0);
		vector<unsigned int> vec_second(alphabet_second, 0);

		unsigned int sum = 0;

		for (size_t i = 0; i < vec.size(); i++)
		{			
			sum += vec[i];
			vec_first[i / 2] += vec[i];
			vec_second[i % 2] += vec[i];
		}
		entropy_avg_first += entropy(vec_first) * sum;
		entropy_avg_second += entropy(vec_second) * sum;
		
	}
	size_t total = seq.size() - seq_len;
	return make_pair(entropy_avg_first /total, entropy_avg_second / total);
}


int main(int argc, char* argv[])
{
	vector<unsigned int> x = {0, 1, 2, 0, 0, 1, 0, 2, 1, 0, 0};
	vector<unsigned int> y = {1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1};
	vector<unsigned int> z = compose(x, y, 2);	


	unsigned int seq_len = stoi(argv[1]);
	
	double composed_entropy = get_entropy(z, 6, seq_len);
	cout << "composed entropy = " << composed_entropy << endl;


 	pair<double, double> proj_entropy = get_proj_entropy(x, y, 3, 2, seq_len);
	cout << "first projected entropy = " << proj_entropy.first << endl;
	cout << "second projected entropy = " << proj_entropy.second << endl;

	double entropy_first = get_entropy(x, 3, seq_len);
	double gamma_first = (entropy_first - proj_entropy.first) / entropy_first;
	cout << "first gamma = " << gamma_first << endl;
	
	double entropy_second = get_entropy(y, 2, seq_len);
	double gamma_second = (entropy_second - proj_entropy.second) / entropy_second;
	cout << "second gamma = " << gamma_second << endl;

	return 0;
}
