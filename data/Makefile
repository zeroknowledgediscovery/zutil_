run: compile
	./dtw_ucr_compiled Trace_TRAIN.tsv Trace_TEST.tsv	
compile: dtw_ucr.cpp
	g++ -O3 -fopenmp dtw_ucr.cpp -o dtw_ucr_compiled 
valgrind: dtw_ucr.cpp
	g++ -g -O3 -fopenmp dtw_ucr.cpp -o dtw_ucr_compiled 
	valgrind  --leak-check=full --track-origins=yes ./dtw_ucr_compiled Adiac_TRAIN.tsv Adiac_TEST.tsv	
