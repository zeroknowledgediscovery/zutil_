#include <stdlib.h>
#include "semantic.h"
#include "config.h"

int main(int argc, char *argv[])
{
  int NUM_EACH=100;
  unsigned int len=1000,BEG=1,END=150,NAMEWIDTH=3;
  unsigned int COLUMNNUM=1;
  unsigned int EVALUATION_DEPTH=2;
  vector <double> partition;
  bool VERBOSE_=false, ONLY_SAE=false, ONLY_PAST=false;
  int HIST=1;

  //matrix_dbl partition_m;
  //vector<unsigned int> COLUMNNUM_M;
  
  string configfile="config.cfg", 
    datafile="data.txt",
    outfileM="mapH.txt", 
    outfileE="outE.txt", 
    outfileS="outS.txt", 
    outfileD="dim.txt";
  string DATA_DIR="up";
  string DATA_TYPE="continuous";

  if (argc > 1)
    configfile = argv[1];
  /*
  if (argc > 2)
    datafile = argv[2];
  */

  /** Read Config file*/
  CONFIG cfg(configfile);
  cfg.set(BEG,"BEG");
  cfg.set(END,"END");
  cfg.set(NAMEWIDTH,"NAMEWIDTH");
  cfg.set(COLUMNNUM,"COLUMNNUM");
  cfg.set(outfileM,"OUTFILE_MAP");
  cfg.set(outfileD,"OUTFILE_DIM_ERROR");
  cfg.set(outfileE,"OUTFILE_EMBEDDING_COORDINATES");
  cfg.set(outfileS,"OUTFILE_SAE");
  cfg.set(len,"DATA_LENGTH");
  cfg.set(NUM_EACH,"NUM_EACH");
  cfg.set(VERBOSE_,"VERBOSE");
  cfg.set(ONLY_SAE,"ONLY_SAE");
  cfg.set(ONLY_PAST,"ONLY_PAST");
  cfg.set(HIST,"ONLY_PAST_HISTORY");
  cfg.set(EVALUATION_DEPTH,"EVALUATION_DEPTH");
  cfg.set_vector<double>(partition,"PARTITION");
  cfg.set(DATA_TYPE,"DATA_TYPE");
  cfg.set(DATA_DIR,"DATA_DIR");
  /*
    cfg.set_vector<unsigned int>(COLUMNNUM_M,"COLUMNNUM_M");
    cfg.set_map<unsigned int, unsigned int,double>(partition_m,"PARTITION_M");
  */

  if (argc > 2)
    END = atoi(argv[2]);

  if (argc > 3)
    datafile = argv[3];


  
  vector < Symbolic_string_ > Svec;

  if (BEG > 0)
    get_continuous_DataMatrix(BEG,END,NAMEWIDTH,
			      COLUMNNUM,partition,
			      len,Svec); //SUEL3
  else
    {
      data_reader *R;
      if (DATA_TYPE=="continuous")
	R = new data_reader(datafile,DATA_DIR,partition,len);
      else
	R = new data_reader(datafile,DATA_DIR,len);
      Svec = R->getsymbolic_string_vector();
    }

  /*
    get_continuous_DataMatrix(BEG,END,NAMEWIDTH,
    COLUMNNUM_M,partition_m,
    len,Svec); //SUEL3
  */
  Set_symbolic_string_ *DSS;
  if (!ONLY_PAST)
    {
      DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH, ONLY_SAE);
  
      if (!ONLY_SAE)
	{
	  matrix_dbl M1= DSS->distance_matrix() ;
	  vector <double> E(DSS->embedding_error());
	  matrix_dbl M2= DSS->embedding_coordinates(3000);

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
      else
	{
	  ofstream outS(outfileS.c_str());
	  vector <double> sae(DSS->get_sae());
	  outS << sae << endl;
	  outS.close();
	}
    }
  else
    {
      DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH,ONLY_PAST, HIST);
      ofstream outS(outfileS.c_str());
      vector <double> hdiff(DSS->get_hdiff());
      outS << hdiff << endl;
      outS.close();
    }
}

