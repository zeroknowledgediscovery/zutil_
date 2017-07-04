#include "semantic.h"
#include "config.h"

//------------------------------------
int main(int argc, char *argv[])
{
  int len=500;
  string graphpref="pfsa", dotcfg="dotII.cfg";
  bool  MC_PRINT=false, 
    MC_STRING=false, INVERSE=false;
  bool SHOW_STATIONARY=false, SHOW_PI=false, SHOW_GAMMA=false;
  string configfile="config.cfg";
  map <state, map <symbol, string> > REPLACE_SYM;
  string title="";
  int DRAW_GRAPH=1;

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
  cfg.set(MC_PRINT,"MC_PRINT");
  cfg.set(INVERSE,"USE_INVERSE_MACHINE");
  cfg.set(MC_STRING,"MC_STRING");
  cfg.set(SHOW_STATIONARY,"SHOW_STATIONARY");
  cfg.set(SHOW_PI,"SHOW_PI");
  cfg.set(SHOW_GAMMA,"SHOW_GAMMA");



  //  cfg.set(RAND_PIT,"RANDOM_PITILDE");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");
  cfg.set_vector<symbol>(omega,"OMEGA");
  cfg.set_map<state,symbol,string>(REPLACE_SYM,"REPLACE_SYM");

  if (argc > 2)
    DRAW_GRAPH = atoi(argv[2]);
  if (argc > 3)
    graphpref = argv[3];
  if (argc > 4)
    dotcfg = argv[4];
  if (argc > 5)
    title = argv[5];

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

  //cout << dotcfg << endl;
  

  if (DRAW_GRAPH==1)
    G.drawGraph(graphname,dotcfg,title);
  if (DRAW_GRAPH==2)
    G.drawGraphX(graphname,dotcfg,title);


  
  return 0;
}

 
