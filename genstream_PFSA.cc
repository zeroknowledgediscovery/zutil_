#include "semantic.h"
#include "config.h"

//------------------------------------
//------------------------------------
/*
ostream& operator << (ostream &out, connx &M)
{
  for(connx::iterator itr=M.begin();
      itr!=M.end();
      ++itr)
    {
      for(map_sym_state::iterator itr1=itr->second.begin();
	  itr1!=itr->second.end();
	  ++itr1)
	out  << itr1->second << " "; 
      out << endl;
    }
  return out;
}
//------------------------------------
//------------------------------------
ostream& operator << (ostream &out, symbol_list_ &s)
{
  for (unsigned int i=0; i<s.size(); i++)
    out <<  (s[i]) << "|";
  return out;
}
//------------------------------------
//------------------------------------

ostream& operator << (ostream &out, vector<double> v)
{
  for (unsigned int i=0; i<v.size(); i++)
    out <<  (v[i]) << " ";
  return out;
}

//------------------------------------
//------------------------------------
ostream& operator << (ostream &out,  
		      map < symbol, unsigned int > &M)
{
  for(map < symbol, 
	unsigned int >::iterator itr=M.begin(); 
      itr!=M.end();
      ++itr)
    out  << itr->second << " "; 

  return out;
}
//------------------------------------
//------------------------------------
ostream& operator << (ostream &out, phi_data_type_ &Phi)
{
  for(phi_data_type_::iterator itr=Phi.begin(); 
      itr!=Phi.end();
      ++itr)
    {
      symbol_list_ key(itr->first);
      out << key  << ": " << itr->second << endl; 
    }

  return out;
}
//------------------------------------
//------------------------------------
ostream& operator << (ostream &out, stoch_phi_data_type_ &sPhi)
{
  for(stoch_phi_data_type_::iterator itr=sPhi.begin(); 
      itr!=sPhi.end();
      ++itr)
    {
      symbol_list_ key(itr->first);
      out << key  << "-->" ; 
      for(vector <double>::iterator itr1=itr->second.begin();
	  itr1!=itr->second.end();
	  ++itr1)
	out << " " << *itr1 ;
      out << endl;
    }
  return out;
}
//------------------------------------
*/
//------------------------------------
int main(int argc, char *argv[])
{
  const gsl_rng_type * T;
  gsl_rng * r;  
  gsl_rng_env_setup(); 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set (r, random_seed());

  int len=500;
  string graphpref="pfsa", graphtype="png";
  bool DRAW_GRAPH=false, MC_PRINT=false, 
    MC_STRING=false, RAND_PIT=true, INVERSE=false;
  string configfile="config.cfg";

  if (argc > 1)
    configfile = argv[1];
  connx aut;
  map < state, map < symbol, double > > pitilde_map;
  symbol_list_ omega, omega1;

  CONFIG cfg(configfile);
  cfg.set_map<state,symbol,state>(aut,"CONNX");
  cfg.set(len,"DATA_LENGTH");
  cfg.set(DRAW_GRAPH,"DRAW_GRAPH");
  cfg.set(graphpref,"GRAPHNAME_PREF");
  cfg.set(graphtype,"GRAPH_TYPE");
  cfg.set(MC_PRINT,"MC_PRINT");
  cfg.set(INVERSE,"USE_INVERSE_MACHINE");
  cfg.set(MC_STRING,"MC_STRING");
  cfg.set(RAND_PIT,"RANDOM_PITILDE");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");
  cfg.set_vector<symbol>(omega,"OMEGA");

  /*
    cfg.set_vector<symbol>(omega,"OMEGA");
    cfg.set_vector<symbol>(omega1,"OMEGAx");
  */

  /** pitilde matrix for PFSAs */
  pitilde pit;  
  if (!RAND_PIT)
    {
      state c_state=0;
      for ( map < state, map < symbol, double > >::iterator iti=pitilde_map.begin();
	    iti != pitilde_map.end(); iti++)
	{
	  vector <double> vec_tmp;
	  for (map <symbol, double>::iterator itj=iti->second.begin();
	       itj != iti->second.end();
	       itj++)
	    vec_tmp.push_back(itj->second);
	  pit[c_state++] = vec_tmp;
	}
    }
  string graphname=graphpref;
  PFSA G(pit, aut);
  
  if (DRAW_GRAPH==1)
    drawGraph(aut,pit,graphname,graphtype);
  if (MC_PRINT) 
    G.mc_print();
  //------------------------------------ 
  //------------------------------------ 
 
  Symbolic_string_ Sa, Sb;
  Sa = G.gen_data(len);   
  Sb = G.gen_data(len); 
  symbol_list_ data(Sa.get_symbol_list());
  
  genESeSS X(data, Sa.get_alphabet());


  //  X.get_mc().mc_print();

  //  PFSA GG(X.get_mc());


  cout << "ann_err: " << X.annihilation_error << " eps: " << X.eps << endl;

  //PFSA Ginfer(X.get_mc());
  //cout << Ginfer << endl;

  X.get_mc().drawGraph("Inf", "dotI.cfg");


  cout << "phi_lambda:  " <<  X.phi(omega) << endl;
  cout << "sync string: " << X.sync_string() << endl;

  //  X.get_mc().mc_print();
  cout << X << endl;

  /*
  symbol_list_ data1((~Sb).get_symbol_list());

  vector <symbol_list_> Vdata;
  Vdata.push_back(data1); 
  // Vdata.push_back(data1);
  genESeSS_multistream Y(Vdata, Sa.get_alphabet());
  Y.get_mc().drawGraph("multistream", "dotI.cfg");
  cout << Y << endl;
  */
  // X.get_mc().mc_print();

  /*

  
  vector <symbol_list_> Vdata;
  Vdata.push_back(data); 
  Vdata.push_back(data);

  genESeSS_multistream Y(Vdata, Sa.get_alphabet());

  Y.get_mc().mc_print();
  cout << Y.DEV_ANN_ << endl;

  PFSA Ginfer1(Y.get_mc());

  Y.get_mc().drawGraph("multistream", "dotI.cfg");
  cout << Y.phi(omega) << endl;

  cout << Y.sync_string() << endl;
  */




  return 0;
}

 


//------------------- print aut --------------
/*cout << epsilon << "-----------------" << endl;
for (unsigned int i=0; 
     i < autR.size(); 
     i++)
  {
    for (unsigned int j=0; 
	 j < autR[0].size(); 
	 j++)
      cout << autR[i][j] << " " ;
    cout << endl;
  }
*/
//------------------- print aut --------------


/*

  phi_data_type_ Phi; 
  pin_data_type_ Datapin;
*/
//getPhi(omega1,data,Phi,Datapin);
//getPhi(omega,data,Phi,Datapin);

//cout << dec<< Phi << endl;

/*stoch_phi_data_type_ sPhi;
  mk_stochastic(Phi,sPhi);
  cout << sPhi << endl;
*/
/*
  vector < symbol_list_  > Vdata;
  vector <pin_data_type_>  VDatapin;
  Vdata.push_back(Sa.get_symbol_list());
*/
// Vdata.push_back((!Sb).get_symbol_list());
  
/*  stoch_phi_data_type_ sPhi1;
    phi_data_type_ Phi1;
    
    getPhi(omega,Vdata,Phi1,VDatapin);
    getPhi(omega1,Vdata,Phi1,VDatapin);

    mk_stochastic(Phi1,sPhi1);
    cout << sPhi1 << endl;
*/
/*
  alphabet=Sa.get_alphabet();
  symbol_list_ wsync(getSync(Vdata,Phi,VDatapin)); 
  cout <<  "vector sync:  " << wsync << endl;
*/
/*
  symbol_list_ omega_syn(getSync(data,Phi,Datapin)); 
  cout << "sync:  " <<  omega_syn << endl;
*/
/*
  stoch_phi_data_type_ sPhi2;
  mk_stochastic(Phi,sPhi2);
  cout << sPhi2 << endl;
*/
/*
  cout << wsync1 << endl;
  cout << inc(wsync1) << endl;
  cout   << inc(inc(wsync1)) << endl;
  cout << inc(wsync1) << endl;
*/
  



//------------------------------------
/*
  double entropy(vector <double> &v)
  {
  double S=0.0;
  for(unsigned int i=0; i<v.size();++i)
  if (v[i] > 0)
  S-=v[i]*log(v[i]);
  return S;
  }
  //------------------------------------
  double entropy(map <symbol, unsigned int > &H)
  {
  vector <double> v(mk_stochastic(H));
  return entropy(v);
  }
  //------------------------------------
  */
