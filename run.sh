# XPFSA="./cfgfiles/XPFSA_2_3_2.cfg"
for PFSA in ./cfgfiles/PFSA_M*
do
	# PFSA=./cfgfiles/PFSA_M2_9.cfg
	for XPFSA in ./cfgfiles/XPFSA_2*
	do
		cmd="./bin_practice/test_zbase_llk -p $PFSA -x $XPFSA -l 1000 -n 100"
		echo $cmd
		eval $cmd	
	done
done
