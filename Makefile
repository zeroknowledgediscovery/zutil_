CC = g++ -std=c++11

CFLAGS = -O3 -Wall -Wextra -Wunused -fopenmp -Wl,--as-needed
CFLAGS =  -static -O3 -Wall -Wextra -Wunused -fopenmp -Wl,--as-needed

ZBASE=./zbase

LIBSO=  -lgomp  -lm 
LIBS=   -lboost_system -lboost_thread -lboost_program_options -lboost_timer  -lboost_chrono -Bdynamic -lgsl -lgslcblas

LIBSO=  -static-libstdc++ -static-libgcc -lgomp  -lm 
LIBS=   -lboost_system -lboost_thread -lboost_program_options -lboost_timer  -lboost_chrono -Bdynamic -lgsl -lgslcblas
LIBPATH= $(ZBASE)/lib

DEPS = 
INCLUDES = -I$(ZBASE)


OBJ =   llk prun prunX drawpfsa pfsadyn2param  computepfsadistance

all:	$(OBJ)  clear mvbin


# compile binaries --------------------------------


%.o :	%.cc
	$(CC) $(INCLUDES)  -c -o $@ $< $(CFLAGS)

llk: llk.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lsemcrct -lconfigfile  $(LIBSO) $(LIBS) 

prun: prun.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH)  -lsemcrct -lconfigfile $(LIBSO) $(LIBS) 


drawpfsa: drawpfsa.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lsemcrct -lconfigfile $(LIBSO) $(LIBS) 


prunX: prunX.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lsemcrct -lconfigfile  $(LIBSO) $(LIBS) 



pfsadyn2param:	pfsadyn2param.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lsemcrct -lconfigfile $(LIBSO) $(LIBS) 


computepfsadistance:	computepfsadistance.o
	$(CC)  $(CFLAGS) -o $@ $^  -L$(LIBPATH) -lsemcrct -lconfigfile $(LIBSO) $(LIBS) 



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
