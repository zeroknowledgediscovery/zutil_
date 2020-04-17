#!/bin/bash

cat ../seqDatasets/seq_${1}.dat ../seqDatasets/seq_${2}.dat > seq.dat
../lsmash -f seq.dat -T symbolic -o lsmash.dst -S 0
../dtw seq.dat $3 dtw.dst
