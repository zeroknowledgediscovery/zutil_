#!/bin/bash

for i in {1..100}
do
	echo $i
	cmd="./llk_baseChoice -f seqDatasets/ -p PFSA_$i.cfg -x XPFSA_$i.cfg -s seq_$i.dat"
	eval $cmd
	sleep 1
done
