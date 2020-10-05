#include <cstdio>
#include <math.h>
#include <iostream>
#include <vector> 
#include <omp.h>
#include <fstream>
#include <sstream>
#include <boost/timer/timer.hpp>

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


class UCRRow
{
    public:
        vector<double>&& get_vec()
        {
            return move(m_data);
        }
    
        void normalize() 
        {   
            double e_sq = this->sq_sum/this->m_data.size();
            double e_sum = this->sum/this->m_data.size();
            double stdev = sqrt(e_sq-e_sum*e_sum);
            for ( int i = 0; i < m_data.size(); i++ ) { 
                double val = m_data[i];           
                m_data[i] = (val-e_sum)/stdev;
            } 
        }
    
        double get_label() const
        {
            return label;
        }
    
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;
        
            vector<double>().swap(m_data);
            bool first = true;
            this->sq_sum = 0;
            this->sum = 0;
            while(std::getline(lineStream, cell, '\t'))
            {
                double val =  std::stod(cell);
                if (first) {
                    first = false;
                    label = val;
                } else {
                    this->sq_sum = this->sq_sum + val*val;
                    this->sum = this->sum + val;
                    this->m_data.push_back(val);
                }
            }
        }
    
    private:
        std::vector<double>    m_data;
        double label;
        double sq_sum;
        double sum;
};

std::istream& operator>>(std::istream& str, UCRRow& data)
{
    data.readNextRow(str);
    return str;
}   


int main(int argc, char *argv[])
{
    vector<vector<double>> train;  
    vector<double> lable_train;
    vector<vector<double>> test;  
    vector<double> lable_test;   
    
    std::ifstream       file_TRAIN(argv[1]); 
    UCRRow row;

    while(file_TRAIN >> row)
    {
        //row.normalize();
        train.push_back(row.get_vec());
        lable_train.push_back(row.get_label());
    }

    std::ifstream       file_TEST(argv[2]);
    
    while(file_TEST >> row)
    {   
        //row.normalize();
        test.push_back(row.get_vec());
        lable_test.push_back(row.get_label());
    }
    vector<double> test_predicted (test.size(),0);
    double correctly_predicted = 0.0;
    
	boost::timer::auto_cpu_timer t;
    for ( int tst = 0; tst < test.size(); tst = tst + 1 ) {
        double min_dist = INF;
        for ( int tes = 0; tes < train.size(); tes = tes + 1 ) {
            double dist =  dtw(test[tst], train[tes], train[tes].size(), min_dist);
            if (dist < min_dist) {
                min_dist = dist;
                test_predicted[tst] = lable_train[tes];
            }
        } 
            if (lable_test[tst] == test_predicted[tst]) {
                correctly_predicted = correctly_predicted + 1;
        }
    } 

    std::cout << "accuracy = " << correctly_predicted/lable_test.size() << endl;
    return 0;
}
