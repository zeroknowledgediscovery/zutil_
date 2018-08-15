#include "semantic.h"
#include "config.h"

//------------------------------------
int main(int argc, char *argv[])
{
  string MC_PRINT_FILE="";
  const gsl_rng_type * T;
  gsl_rng * r;  
  gsl_rng_env_setup(); 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set (r, random_seed());

  int len=500;
  string graphpref="pfsa", graphtype="png";
  bool DRAW_GRAPH=false, MC_PRINT=false, 
    MC_STRING=true, INVERSE=false;
  bool SHOW_STATIONARY=false, SHOW_PI=false, SHOW_GAMMA=false;
  string configfile="config.cfg";
  map <state, map <symbol, string> > REPLACE_SYM;

  if (argc > 1)
    configfile = argv[1];
  connx aut;
  map < state, map < symbol, double > > pitilde_map;
  symbol_list_ omega;

  CONFIG cfg(configfile);
  cfg.set_map<state,symbol,state>(aut,"CONNX");
  cfg.set(len,"DATA_LENGTH");
  cfg.set(DRAW_GRAPH,"DRAW_GRAPH");
  cfg.set(graphpref,"GRAPHNAME_PREF");
  cfg.set(graphtype,"GRAPH_TYPE");
  cfg.set(MC_PRINT,"MC_PRINT");
  cfg.set(INVERSE,"USE_INVERSE_MACHINE");
  cfg.set(MC_STRING,"MC_STRING");
  cfg.set(SHOW_STATIONARY,"SHOW_STATIONARY");
  cfg.set(SHOW_PI,"SHOW_PI");
  cfg.set(SHOW_GAMMA,"SHOW_GAMMA");
  cfg.set(MC_PRINT_FILE,"MC_PRINT_FILE");



  //############ test #############

  //  cfg.set(RAND_PIT,"RANDOM_PITILDE");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");
  cfg.set_vector<symbol>(omega,"OMEGA");
  cfg.set_map<state,symbol,string>(REPLACE_SYM,"REPLACE_SYM");

  if (argc > 2)
    len = atoi(argv[2]);
  if (argc > 3)
    MC_PRINT_FILE = argv[3];

  /** pitilde matrix for PFSAs */ 
  pitilde pit;  
  //if (!RAND_PIT)
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
  string graphname=graphpref;
  PFSA G(pit, aut,0,len);


  if (DRAW_GRAPH==1)
    drawGraph(aut,pit,graphname,graphtype);
  if (MC_PRINT) 
    G.mc_print();
  if (MC_STRING)
    {
      if (MC_PRINT_FILE == "")
	cout << G << endl;
      else
	{
	  ofstream MOUT(MC_PRINT_FILE.c_str());
	  MOUT <<  G;
	  MOUT.close();
	}
    }

  if (SHOW_STATIONARY)
    cout << "Stationary distribution: "<<  G.get_Stationary() << endl;

  if (SHOW_PI)
    {
      cout << "Transition Probability Matrix: " << endl;
      vector < vector <double> > V(G.get_PI());
      for (unsigned int i=0; i < aut.size();++i)
	cout << V[i] << endl;
    }

  if (SHOW_GAMMA)
    {
      cout << "Gamma Matrices: " << endl;
      map < symbol, vector < vector < double > > > Gamma(G.get_Gamma());
      for (symbol i(0); i < aut[0].size();++i)
	{
	  cout << "symbol --> " << i << endl;
	  for (unsigned int j=0; j < aut.size(); j++)
	    cout << Gamma[i][j] << endl;
	}
    }

  /* for (unsigned int i=0;i<REPLACE_SYM[0].size();++i)
    cout << REPLACE_SYM[0][symbol (i)];
  */

  if (!REPLACE_SYM.empty())
    {
      Symbolic_string_  sG(G.gen_data(len));
      symbol_list_ sg = sG.get_symbol_list();
      for (unsigned int i=0;i<sg.size();++i)
	cout << REPLACE_SYM[0][sg[i]];
      cout << endl;
    }


  //------------------------------------ 
  //------------------------------------ 
 
  //  Symbolic_string_ Sa, Sb;
  
  return 0;
}

 
