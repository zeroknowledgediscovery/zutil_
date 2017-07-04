#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <utility>
#include <sys/time.h>
#include <unistd.h>

#include "semantic.h"
#include "config.h"

int main(int argc, char *argv[])
{
  const gsl_rng_type * T;
  gsl_rng * r;   /* create  generator chosen by  env  var GSL_RNG_TYPE */   
  gsl_rng_env_setup(); 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set (r, random_seed());

  /** Length of data stream from each PFSA*/
  int len=500;
  /** Number of PFSAs considered */
  int NUM_MC=10, NUM_EACH=1;
  string graphpref="pfsa", graphtype="png";
  bool DRAW_GRAPH=false, MC_PRINT=false, MC_STRING=false, RAND_PIT=false, INVERSE=false;
  string configfile_master="config.cfg", configfile_prefix="cfg_MC",
    outfileM="mapH.txt", 
    outfileE="outE.txt", 
    outfileD="dim.txt";

  if (argc > 1)
    configfile_master = argv[1];
  /** Connection matrix for each automaton*/
  /** Read Config file*/
  CONFIG cfg0(configfile_master);
  cfg0.set(len,"DATA_LENGTH");
  cfg0.set(NUM_EACH,"NUM_EACH");
  cfg0.set(DRAW_GRAPH,"DRAW_GRAPH");
  cfg0.set(graphpref,"GRAPHNAME_PREF");
  cfg0.set(graphtype,"GRAPH_TYPE");
  cfg0.set(MC_PRINT,"MC_PRINT");
  cfg0.set(INVERSE,"USE_INVERSE_MACHINE");
  cfg0.set(MC_STRING,"MC_STRING");
  cfg0.set(outfileM,"OUTFILE_MAP");
  cfg0.set(outfileD,"OUTFILE_DIM_ERROR");
  cfg0.set(outfileE,"OUTFILE_EMBEDDING_COORDINATES");
  cfg0.set(configfile_prefix,"PREF_STRING");
  cfg0.set(NUM_MC,"NUM_MC");

  vector <PFSA> G;
  /** pitilde matrix for PFSAs */
  for (int num_pfsa=0; num_pfsa < NUM_MC; num_pfsa++)
    {
      char buff[1024];
      sprintf(buff,"%d",num_pfsa);
      string graphname=graphpref+buff;
      string configfile_counter=configfile_prefix + buff;

      connx aut;
      map < state, map < symbol, double > > pitilde_map;

      CONFIG cfg(configfile_counter);
      cfg.set_map<state,symbol,state>(aut,"CONNX");
      //      cfg.set(len,"DATA_LENGTH");
      //      cfg.set(NUM_EACH,"NUM_EACH");
      cfg.set(DRAW_GRAPH,"DRAW_GRAPH");
      //      cfg.set(graphpref,"GRAPHNAME_PREF");
      //      cfg.set(graphtype,"GRAPH_TYPE");
      //    cfg.set(MC_PRINT,"MC_PRINT");
      cfg.set(INVERSE,"USE_INVERSE_MACHINE");
      //      cfg.set(MC_STRING,"MC_STRING");
      cfg.set(RAND_PIT,"RANDOM_PITILDE");
      cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");

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
      unsigned int alphabet =aut[0].size();
      unsigned int NUM_STATE=aut.size();

      PFSA G_tmp(pit, aut);
      G.push_back(G_tmp);
      if (DRAW_GRAPH==1)
	drawGraph(aut,pit,graphname,graphtype);
      if (MC_PRINT)
	G_tmp.mc_print();
      if (MC_STRING)
	{
	  Symbolic_string_ Stmp;
	  if (!INVERSE)
	    Stmp = G_tmp.gen_data(len);      
	  else
	    Stmp = (!G_tmp).gen_data(len);      
	  cout << Stmp << endl;
	}
    }

  vector <Symbolic_string_> SS;

  for (int num_pfsa_i=0; num_pfsa_i < NUM_MC; num_pfsa_i++)
    {
      Symbolic_string_ Si = G[num_pfsa_i].gen_data(len);      
      SS.push_back(Si);
    }
  Set_symbolic_string_ DSS(SS,NUM_EACH);
  matrix_dbl M1= DSS.distance_matrix() ;
  vector <double> E(DSS.embedding_error());
  matrix_dbl M2= DSS.embedding_coordinates(3000);


  ofstream outM(outfileM.c_str());
  outM << M1 << endl;
  ofstream outD(outfileD.c_str());
  outD << E << endl;
  ofstream outE(outfileE.c_str());
  outE  << M2 << endl;
  outM.close();
  outD.close();
  outE.close();

}


