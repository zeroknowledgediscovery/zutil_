#include <cstdio>
#include <math.h>
#include <vector> 
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>


#define INF 1e20       //Pseudo Infitinte number for this code
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) sqrt((x-y)*(x-y))

using namespace std;

double dtw(const vector<double> &A, const vector<double> &B, int r, double bsf = INF)
{
    int m = A.size();
    vector<double> cost (2*r+1,INF);

    vector<double> cost_prev (2*r+1,INF);

    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
			{
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost) {   
                min_cost = cost[k];
            }
        }
        
        if (min_cost>bsf) {
            return INF;
        }

        cost.swap(cost_prev);
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    return final_dtw;
}


class reader
{
public:
	vector<double> get_data()
	{
		return data;
	}
	void readNextRow(istream& str)
	{
		data.clear();
		string line;
		getline(str, line);
		stringstream lineStream(line);
        string cell;
       	while(getline(lineStream, cell, ' '))
        {
			double val = stod(cell);
            data.push_back(val);
        }
	}
private:
    vector<double> data;
};

istream& operator>>(istream& str, reader& row)
{
    row.readNextRow(str);
    return str;
}   


int main(int argc, char *argv[])
{
    vector<vector<double>> train;  
    
    ifstream file_TRAIN(argv[1]); 
    reader row;

    while(file_TRAIN >> row)
	{
        train.push_back(row.get_data());
	}
	// cout << "number of sequence = " <<  train.size() << endl;
	// cout << "sequence length = " << train[0].size() << endl;
	// for (size_t i = 0; i < train.size(); i++)
	// {
	// 	for(size_t j = 0; j < train[i].size(); j++)
	// 	{
	// 		cout << train[i][j];
	// 		(j < train[i].size() - 1)? cout << " " : cout << endl;
	// 	}
	// }
	

    // #if defined(_OPENMP)    
    //     omp_set_num_threads(28);
    //     #pragma omp parallel
    //     #pragma omp for reduction(+:correctly_predicted)
    // #endif
    
	size_t dim = train.size();
	int window_size = stoi(argv[2]);

	vector<vector<double>> matrix(dim, vector<double>(dim, 0));
    for (size_t i = 0; i < dim; i++) 
	{
		double min_dist = INF;
        for (size_t j = i + 1; j < dim; j++) 
		{
            double dist = dtw(train[i], train[j], window_size, min_dist);
			matrix[i][j] = dist;
			matrix[j][i] = dist;
			// cout << i << ", " << j << " " << dist << endl;
        }
    }
	
	string output_file = argv[3];
	ofstream output(output_file.c_str());
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for(size_t j = 0; j < matrix[i].size(); j++)
		{
			output << matrix[i][j];
			(j < matrix[i].size() - 1)? output << " " : output << endl;
		}
	}
	output.close();

    return 0;
}
