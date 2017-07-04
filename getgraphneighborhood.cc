// Implementation: Testprogram for 2-dimensional Range Trees
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.
#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <utility>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>
#include <CGAL/squared_distance_2.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <stdlib.h>

using namespace std;


typedef CGAL::Cartesian<double> Representation;
typedef CGAL::Point_3< Representation> Point_3; 
typedef CGAL::Vector_3< Representation> Vector_3; 
typedef CGAL::Direction_3< Representation> Direction_3; 

//typedef CGAL::Range_tree_map_traits_3 <Representation, Point_3  > Traits;

typedef CGAL::Range_tree_map_traits_3 <Representation, unsigned int  > Traits;

typedef CGAL::Range_tree_3 <Traits> Range_tree_3_type;

//typedef CGAL::Range_tree_3 <TraitsII> Range_tree_II_type;

typedef Traits::Key Key;
typedef Traits::Interval Interval;


//typedef TraitsII::Key KeyII;
//typedef TraitsII::Interval IntervalII;


/*
  nodedef>name VARCHAR,label VARCHAR,class VARCHAR, visible BOOLEAN,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,x DOUBLE,y DOUBLE,color VARCHAR
  s1,'Hello "world" !',type1,true,true,10.0,10.0,-52.11296,-25.921143,'114,116,177'
  s2,'Well, this is',type1,true,true,10.986123,10.986123,-20.114172,25.740356,'219,116,251'
  s3,'A correct 'GDF' file',type1,true,true,10.986123,10.986123,8.598924,-26.867584,'192,208,223'
  edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR
  s1,s2,true,'114,116,177'
  s2,s3,true,'219,116,251'
  s3,s2,true,'192,208,223'
  s3,s1,true,'192,208,223'
*/

void writeGDF(string gdfFile,vector <vector <Key> >& G, 
	      map<unsigned int,vector<string> >& DESC )
{
  string COMMA = ",";
  string NODEDEF = "nodedef>name VARCHAR,label VARCHAR, color VARCHAR, weight DOUBLE, family VARCHAR";
  string EDGEDEF = "edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR, weight DOUBLE";
  string COLOR = "'100,100,100'";
  string COLOR_V = "'240,100,100'";

  ofstream OUT(gdfFile.c_str());

  OUT << NODEDEF << endl;

  if(DESC.empty())
    for(unsigned int i=0; i<G.size();++i)
	OUT << i << COMMA << i << COMMA 
	    <<  COLOR_V << COMMA << 0 << COMMA << "x" << endl;
  else
    for(unsigned int i=0; i<G.size();++i)
	OUT << i << COMMA 
	    << DESC[i][0] << COMMA 
	    <<  COLOR_V << COMMA 
	    << DESC[i][3] << COMMA << DESC[i][2] <<endl;

  
  OUT << EDGEDEF << endl;

  for(unsigned int i=0; i<G.size();++i)
    for(unsigned int j=0; j<G[i].size();++j)
      OUT << i << COMMA << G[i][j].second << COMMA << "false" << COMMA << COLOR << COMMA<< 1.0 << endl;;

  OUT.close();

}

//------------------------------------------------
/*!
  \b \f$ \epsilon\f$-resolution comparison class 
*/
class eps_compare_Pts { 

public:
  eps_compare_Pts(){epsilon=0.1;};		
  eps_compare_Pts(double s_){epsilon=s_;};		
  bool operator()(const Point_3 &lhs,const Point_3 &rhs) const 
  { 
    pair<double, double> lhs_p = make_pair(lhs.x(),lhs.y());
    pair<double, double> rhs_p = make_pair(rhs.x(),rhs.y());

    double S = sqrt((lhs.x()-rhs.x())*(lhs.x()-rhs.x()) 
		    + (lhs.y()-rhs.y())*(lhs.y()-rhs.y()));
    return (S>epsilon) && (lhs_p < rhs_p); 
  }

private:
  double epsilon;
};
//------------------------------------------------
ostream& operator << (ostream &out, Point_3  P)
{
  out << P.x() << " " << P.y() << " " << P.z();
  return out;
};
ostream& operator << (ostream &out, Vector_3  P)
{
  out << P.x() << " " << P.y()<< " " << P.z() ;
  return out;
};
ostream& operator << (ostream &out, Direction_3  P)
{
  out << P.dx() << " " << P.dy() << " " << P.dz() ;
  return out;
};
ostream& operator << (ostream &out, vector <Point_3> P)
{
  for(vector <Point_3>::iterator itr=P.begin();
      itr != P.end();
      ++itr)
    out << itr->x() << " " << itr->y()<< " " << itr->z() << endl ;
  return out;
};
ostream& operator << (ostream &out, 
		      map <Point_3,double,eps_compare_Pts> P)
{
  for(map <Point_3,double>::iterator itr=P.begin();
      itr != P.end();
      ++itr)
    out << (itr->first).x() << " " 
	<< (itr->first).y() << " " 
	<< (itr->first).z() << " " 
	<< itr->second  <<  endl ;
  return out;
};
//------------------------------------------------



int main(int argc, char* argv[])
{

  vector<Key> InputList;

  string data_file, outfile="out.dat";
  double EPS=0.00001;
  string gdfFile="net.gdf";
  string statFile="";

  if (argc > 1)
    data_file=argv[1];
  if (argc > 2)
    EPS=atof(argv[2]);
  if (argc > 3)
    gdfFile=argv[3];
  if (argc > 4)
    statFile=argv[4];

  map<unsigned int,vector<string> > DESC;

  if(statFile!="")
    {
      ifstream ST(statFile.c_str());
      string ln;
      unsigned int count=0;
      while(getline(ST,ln))
	{
	  string token;
	  vector <string> tokens;
	  stringstream ss(ln);
	  while(getline(ss,token,','))
	    tokens.push_back(token);
	  DESC[count++]=tokens;
	}
      ST.close();
    }

  //vector <Point_3> points_;
  Point_3 ZERO(0.0,0.0,0.0);
  Vector_3 eps_vectorP(EPS/2,EPS/2,EPS/2);
  Vector_3 eps_vectorM(-EPS/2,-EPS/2,-EPS/2);

  ifstream IN(data_file.c_str());
  string line;
  unsigned int count=0;
  double val;
  vector< vector<double> > Val(3);
  while(getline(IN,line))
    {
      stringstream ss(line);
      while(ss>>val)
	Val[count].push_back(val);
      count++;
    }
  for(unsigned int i=0; i< Val[0].size();++i)
    InputList.push_back(Key(Point_3 (Val[0][i],Val[1][i],Val[2][i]), i) );

   
  Range_tree_3_type Range_tree_3(InputList.begin(),InputList.end());  
  vector <vector <Key> > OUTPUTS;

  for(vector<Key>::iterator thisP=InputList.begin();
      thisP!=InputList.end();
      ++thisP)
    {

      Point_3 curr(thisP->first.x(),thisP->first.y(),thisP->first.z());
      Interval win=Interval(Key(curr+eps_vectorM,0),Key(curr+eps_vectorP,0));
      vector<Key> OutputList;

      Range_tree_3.window_query(win, std::back_inserter(OutputList));
      OUTPUTS.push_back(OutputList);

    }


  for(unsigned int i=0; i<OUTPUTS.size();++i)
    {
      cout << i << " :  ";
      for(unsigned int j=0; j<OUTPUTS[i].size();++j)
	cout << OUTPUTS[i][j].second << " ";
      cout << endl;

    }

  writeGDF(gdfFile, OUTPUTS,DESC);

  return 0;

}
