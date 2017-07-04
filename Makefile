CC = g++ 
CFLAGS = -std=c++11 -static -O3 -Wall -Wextra -Wunused -fopenmp -Wl,--as-needed
CFLAGS_DYN = -std=c++11  -O3 -Wall -Wextra -Wunused -fopenmp -Wl,--as-needed
ZBASE=./zbase

LIBSO=  -static-libstdc++ -static-libgcc -lgomp  -lm 
LIBS=   -lboost_system -lboost_thread -lboost_program_options -lboost_timer  -lboost_chrono -Bdynamic -lgsl -lgslcblas
LIBPATH= $(ZBASE)/lib

DEPS = 

INCLUDES = -I$(ZBASE)


OBJ =   prun prunX drawpfsa pfsadyn2param  computepfsadistance

all:	prun prunX drawpfsa pfsadyn2param  computepfsadistance  clear mvbin


# compile libraries --------------------------------
semantic.o: semantic.cc $(DEPS)
	$(CC) -c -o $@ $< -fopenmp $(CFLAGS) $(INCLUDES); ar rcs libsemcrct.a semantic.o; rm semantic.o; mv libsemcrct.a ./lib

config.o: config.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS); ar rcs libconfigfile.a config.o; rm config.o; mv libconfigfile.a ./lib

# compile binaries --------------------------------


prun.o: prun.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)

prun: prun.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lconfigfile -lsemcrct $(LIBSO) $(LIBS) 


drawpfsa.o: drawpfsa.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)

drawpfsa: drawpfsa.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lconfigfile -lsemcrct $(LIBSO) $(LIBS) 


prunX.o: prunX.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)

prunX: prunX.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lconfigfile -lsemcrct $(LIBSO) $(LIBS) 


pfsadyn2param.o: pfsadyn2param.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)


computepfsadistance.o: computepfsadistance.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)

pfsadyn2param:	pfsadyn2param.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lconfigfile -lsemcrct $(LIBSO) $(LIBS) 


computepfsadistance:	computepfsadistance.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lconfigfile -lsemcrct $(LIBSO) $(LIBS) 



# utility commands --------------------------------

clear:	
	rm -rf *.o
clean:	
	rm -rf *.o *~ $(OBJ) *.tgz

mvbin:
	rm -rf ./bin >&/dev/null; mkdir bin; mv $(OBJ) ./bin
src:
	tar -czvf src.tgz *.cc *.h Makefile Doxyfile
ref:
	doxygen Doxyfile; cd doc/latex; make; make

 # EOF --------------------------------
