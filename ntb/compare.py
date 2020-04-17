import numpy as np
import time
import os
import sys
import subprocess
import pandas as pd

def separation(D, labels):
	dim = len(D)
	A, B = 0, 0
	ca, cb = 0, 0
	for i in range(dim):
		for j in range(dim):
			k, l = labels[i], labels[j]
			if k == l:
				A += D[i][j]
				ca += 1
			else:
				B += D[i][j]
				cb += 1
	r = (A / ca) / (B / cb)
	return r


if __name__ == '__main__':
	
	lsmash = '/home/yhuang10/D3M/zutil_/lsmash'
	smash = "/home/yhuang10/D3M/zutil_/smash"
	dtw = '/home/yhuang10/D3M/zutil_/dtw'
	datasets = '/home/yhuang10/D3M/zutil_/seqDatasets/'

	window_sizes = [5, 10, 20, 30, 40, 50, 100]
	sequence_length = int(sys.argv[1])


	labels = [0] * 25 + [1] * 25
	num_datasets = 200

	r_lsmash, r_smash = [], [] # the ratio (mean within-class distance) / (mean between-class distance)
	t_lsmash, t_smash = [], [] # running time 
	r_dtw = {'r_dtw_{}'.format(ws): [] for ws in window_sizes}
	t_dtw = {'t_dtw_{}'.format(ws): [] for ws in window_sizes}
	I, J = [], []

	n = 0
	while n < num_datasets:
		i = np.random.randint(100) + 1
		j = np.random.randint(100) + 1
		if i == j:
			# Do not calculate separtion if two classes are identical
			# yes, it is rare, but it does happen.
			continue
		
		n += 1
		I.append(i)
		J.append(j)
		
		# Form random bi-classification datasets  
		os.system('cat {}/seq_{}.dat {}/seq_{}.dat | cut -f 1-{} -d \" \" > seq.dat'.format(datasets, i, datasets, j, sequence_length))

		# lsmash
		t0 = time.time()
		status_lsmash = os.system('{} -f seq.dat -T symbolic -o lsmash.dst -S 0'.format(lsmash))
		t1 = time.time()
		t_lsmash.append(t1 - t0)
		dist_lsmash = np.genfromtxt('lsmash.dst')
		r = separation(dist_lsmash, labels)
		r_lsmash.append(r)

		# smash
		t0 = time.time()
		status_smash = os.system('{} -f seq.dat -n 2 -T symbolic -D row -o smash.dst'.format(smash))
		t1 = time.time()
		t_smash.append(t1 - t0)
		dist_smash = np.genfromtxt('smash.dst')
		r = separation(dist_smash, labels)
		r_smash.append(r)
		
		# dtw
		for ws in window_sizes:
			t0 = time.time()
			status_dtw = os.system('{} seq.dat {} dtw.dst'.format(dtw, ws))
			t1 = time.time()
			t_dtw['t_dtw_{}'.format(ws)].append(t1 - t0)
			dist_dtw = np.genfromtxt('dtw.dst')
			r = separation(dist_dtw, labels)
			r_dtw['r_dtw_{}'.format(ws)].append(r)
		
		print('lsmash time = {}, smash time = {}, '.format(t_lsmash[-1], t_smash[-1])),
		for ws in window_sizes:
			print('dtw time {} = {}'.format(ws, t_dtw['t_dtw_{}'.format(ws)][-1])),
		print
		
		print('round {}, ({}, {}), lsmash = {}, smash = {}, '.format(n, i, j, r_lsmash[-1], r_smash[-1])),
		for ws in window_sizes:
			print('dtw {} = {}'.format(ws, r_dtw['r_dtw_{}'.format(ws)][-1])),
		print

	data = {
		'class_0': I, 
		'class_1': J, 
		'r_lsmash': r_lsmash, 
		'r_smash': r_smash,  
		't_lsmash': t_lsmash,
		't_smash': t_smash,
	}

	data.update(r_dtw)
	data.update(t_dtw)	
	df = pd.DataFrame(data=data)
	print(df)
	df.to_csv('performance_compare_sl{}.csv'.format(sequence_length))
