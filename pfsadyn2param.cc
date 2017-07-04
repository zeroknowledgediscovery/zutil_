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
  bool DRAW_GRAPH=false, MC_PRINT=false, MC_STRING=false, RAND_PIT=true, INVERSE=false;
  string configfile="config.cfg",
    outfileM="mapH.txt", 
    outfileE="outE.txt", 
    outfileD="dim.txt";

  if (argc > 1)
    NUM_MC = atoi(argv[1]);
  if (argc > 2)
    configfile = argv[2];
  /** Connection matrix for each automaton*/
  connx aut;
  map < state, map < symbol, double > > pitilde_map;

  /** Read Config file*/
  CONFIG cfg(configfile);
  cfg.set_map<state,symbol,state>(aut,"CONNX");
  cfg.set(len,"DATA_LENGTH");
  cfg.set(NUM_EACH,"NUM_EACH");
  cfg.set(DRAW_GRAPH,"DRAW_GRAPH");
  cfg.set(graphpref,"GRAPHNAME_PREF");
  cfg.set(graphtype,"GRAPH_TYPE");
  cfg.set(MC_PRINT,"MC_PRINT");
  cfg.set(INVERSE,"USE_INVERSE_MACHINE");
  cfg.set(MC_STRING,"MC_STRING");
  cfg.set(RAND_PIT,"RANDOM_PITILDE");
  cfg.set(outfileM,"OUTFILE_MAP");
  cfg.set(outfileD,"OUTFILE_DIM_ERROR");
  cfg.set(outfileE,"OUTFILE_EMBEDDING_COORDINATES");
  cfg.set_map<state,symbol,double>(pitilde_map,"PITILDE");

  vector <PFSA> G;
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

  
  
  unsigned int alphabet =aut[0].size();

  unsigned int NUM_STATE=aut.size();
  for (int num_pfsa=0; num_pfsa < NUM_MC; num_pfsa++)
    {
	if (RAND_PIT)
	{
		for (unsigned int num_state=0;num_state<NUM_STATE; num_state++)
		{
			vector <double> pitrow;
			double remaining_prob=1.0;
			for (symbol s(0); s < alphabet-1; s++)
			{
				double probval= gsl_ran_flat (r, 0.0,remaining_prob);
				remaining_prob -= probval;
				pitrow.push_back(probval);
			}
		pitrow.push_back(remaining_prob);
		pit[num_state] = pitrow;
		}
	  }
      char buff[1024];
      sprintf(buff,"%d",num_pfsa);
      string graphname=graphpref+buff;
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
      // cout << "LENGTH " << dec<<  len << endl;
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

/*
  Symbolic_string_ T, TT(!Si), TTT(~Si), TTTT(!!!Si);	      
  T=!!Si;
  Si.get_norm();
  cout << "Si " << Si.norm << endl;
  TTT.get_norm();
  cout << "~Si " <<  TTT.norm << endl;
  TT.get_norm();
  cout << "!Si " << TT.norm << endl;
  T.get_norm();
  cout << "!!Si " << T.norm << endl;
  TTTT.get_norm();
  cout << "!!!Si " <<  TTTT.norm << endl;
	      
*/
