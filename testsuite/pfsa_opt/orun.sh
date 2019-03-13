#!/bin/bash

PFSA="M2.cfg M4.cfg S2.cg T3.cfg"

coeff0=$1
coeff1=$2
coeff2=$3
coeff3=$4
data=$5

PATHA="/home/ishanu/ZED/Research/structural_algebra_/bin/psolve -o tmp.cfg -x "
EQ=" ($coeff0 *M2.cfg) + ($coeff1 *M4.cfg) + ($coeff2 *S2.cfg) + ($coeff3 *T3.cfg)  "

$PATHA " $EQ "

PATHB="/home/ishanu/Dropbox/ZED/Research/zutil_/bin/llk -f tmp.cfg -T symbolic  -s $data "

$PATHB
 
