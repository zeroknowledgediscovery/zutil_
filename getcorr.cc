#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>

using namespace std;

//--------------------------------
void getcorr(unsigned int row,
	     unsigned int col,
	     size_t n,
	     vector<double>& v1, 
	     vector<double>& v2, 
	     map<unsigned int, 
	     map<unsigned int,
	     double> >&M)
{
  //vector <double> tmp(2*n,0.0);
  double val=gsl_stats_correlation(&v1[0],1, &v2[0],1,n);

  val = 1-((val+1)/(2.0));

  M[row][col]=val;
  M[col][row]=  M[row][col];
  return ;
};
//--------------------------------
ostream& operator << (ostream &out, map<unsigned int, 
		      map<unsigned int,
		      double> >&M)
{
  for( map<unsigned int, 
	 map<unsigned int,
	 double> >::iterator itr=M.begin();
       itr!=M.end();
       ++itr)
    {
      for(map<unsigned int,
	    double>::iterator itr_=itr->second.begin();
	  itr_!=itr->second.end();
	  ++itr_)
	out << itr_->second << " " ;
      out << endl;
    }
  return out;
};
//---------------------------------
int main(int argc, char* argv[])
{
  unsigned int NUM=0;
  string datfile="data.dat",ofile="outres.txt";
  if(argc>1)
    datfile=argv[1];
  if(argc>2)
    ofile=(argv[2]);
  if(argc>3)
    NUM=atoi(argv[3]);

  unsigned int count=0;
  map<unsigned int, vector<double> >H;
  string line;
  double val;

  ifstream IN(datfile.c_str());
  while(getline(IN,line))
    {
      stringstream ss(line);
      while(ss>>val)
	H[count].push_back(val);
      count++;
    }
  IN.close();

  if(NUM==0)
    NUM=H.size();

  size_t n = H[0].size();

  map<unsigned int, map<unsigned int,double> > M;
  for(unsigned int i=0;i<NUM;++i)
    M[i][i]=0.0;


  for(unsigned int i=0; i<NUM;++i)
    for(unsigned int j=0; j<i;++j)
      getcorr(i,j,n,H[i],H[j],M);

  ofstream OUT(ofile.c_str());
  OUT << M ;
  OUT.close();

  return 0;
}
