#!/bin/bash

get_name () {
	file=$1
	index=`echo $file | awk -F"." '{print length($0)-length($NF)}'`
	if [[ $index == 0 ]]
	then
		label_file=${file}_$2
	else
		prefix=${file:0:(($index-1))}
		suffix=${file:$index:((${#file} - $index))}
		label_file=${prefix}_$2.${suffix}
	fi
	echo $label_file
}
eval ./dtw_ucr_compiled $1 $2
train_label_file=`get_name $1 label`
test_label_file=`get_name $2 label`

train_seq_file=`get_name $1 seq`
test_seq_file=`get_name $2 seq`

train_coord_file=`get_name $1 coord`
test_coord_file=`get_name $2 coord`

awk -F"\t" '{print $1}' $1 > $train_label_file
awk -F"\t" '{print $1}' $2 > $test_label_file
awk -F"\t" '{for(i=2;i<=NF;i++){ printf("%s",( (i>2) ? " " : "" ) $i) } ; print ;}' $1 > $train_seq_file
awk -F"\t" '{for(i=2;i<=NF;i++){ printf("%s",( (i>2) ? " " : "" ) $i) } ; print ;}' $1 > $test_seq_file

# cmd="./lsmash_classification -f $train_seq_file -g $test_seq_file -l $train_label_file -k $test_label_file -c $train_coord_file -d $test_coord_file -u 1 -P 0 -T continuous -F baseChoice_cfgfiles/PFSA_0.cfg baseChoice_cfgfiles/PFSA_1.cfg baseChoice_cfgfiles/PFSA_2.cfg baseChoice_cfgfiles/PFSA_3.cfg baseChoice_cfgfiles/PFSA_4.cfg"

cmd="./lsmash_classification -f $train_seq_file -g $test_seq_file -l $train_label_file -k $test_label_file -c $train_coord_file -d $test_coord_file -u 1 -P 0 -T continuous"
num_arg=$#
if [[ $num_arg == 3 ]]
then
	if [[ -f $3/PFSAfiles ]]
	then
		rm $3/PFSAfiles
	fi
	
	for file in `ls $3`
	do
		echo $3/${file} >> $3/PFSAfiles
	done
	sed -i ':a;N;$!ba;s/\n/ /g' $3/PFSAfiles
	PFSAfiles=`cat $3/PFSAfiles`
	cmd="${cmd} -F $PFSAfiles"
fi
echo $cmd
eval $cmd

rm $train_seq_file $test_seq_file $train_label_file $test_label_file 
if [[ $num_arg == 3 ]]
then
	rm $3/PFSAfiles
fi
