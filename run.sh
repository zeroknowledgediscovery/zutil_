XPFSA="./cfgfiles/XPFSA_2_3_2.cfg"
for PFSA in ./cfgfiles/PFSA_M*
do
	# PFSA=./cfgfiles/PFSA_M2_9.cfg
	# for XPFSA in ./cfgfiles/XPFSA_2*
	# do
		cmd="./bin_practice/practice $PFSA $XPFSA 5000" # llk_matrices/llk_matrix"
		echo $cmd
		eval $cmd	
	# done
done
