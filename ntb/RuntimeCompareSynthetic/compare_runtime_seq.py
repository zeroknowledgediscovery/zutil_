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
	
	lsmash = 'bin/lsmash'
	smash = "bin/smash"
	dtw = 'bin/dtw'
	datasets = 'seqDatasets/'

	sequence_length = 200 * np.arange(1, 11)

	window_size = int(sys.argv[1])

	labels = [0] * 25 + [1] * 25
	num_datasets = 3

	m_lsmash, m_smash, m_dtw = [], [], []
	s_lsmash, s_smash, s_dtw = [], [], []
	for sl in sequence_length:
		
		print('sequence length = {}'.format(sl))
		n = 0
		t_lsmash, t_smash, t_dtw = [], [], [] # running time 
		while n < num_datasets:
			i = np.random.randint(100) + 1
			j = np.random.randint(100) + 1
			if i == j:
				# Do not calculate separtion if two classes are identical
				# yes, it is rare, but it does happen.
				continue
		
			n += 1
		
			# Form random bi-classification datasets  
			os.system('cat {}/seq_{}.dat {}/seq_{}.dat | cut -f 1-{} -d \" \" > seq.dat'.format(datasets, i, datasets, j, sl))

			# lsmash
			t0 = time.time()
			status_lsmash = os.system('{} -f seq.dat -T symbolic -o lsmash.dst -S 0'.format(lsmash))
			t1 = time.time()
			t_lsmash.append(t1 - t0)

			# smash
			t0 = time.time()
			status_smash = os.system('{} -f seq.dat -n 2 -T symbolic -D row -o smash.dst'.format(smash))
			t1 = time.time()
			t_smash.append(t1 - t0)
		
			# dtw
			t0 = time.time()
			status_dtw = os.system('{} seq.dat {} dtw.dst'.format(dtw, window_size))
			t1 = time.time()
			t_dtw.append(t1 - t0)
		
			print('\nlsmash time = {}, smash time = {}, dtw time = {}'.format(t_lsmash[-1], t_smash[-1], t_dtw[-1]))
		
		m_lsmash.append(np.mean(t_lsmash))
		m_smash.append(np.mean(t_smash))
		m_dtw.append(np.mean(t_dtw))
		s_lsmash.append(np.std(t_lsmash))
		s_smash.append(np.std(t_smash))
		s_dtw.append(np.std(t_dtw))

	data = {
		'sl': sequence_length,
		'm_lsmash': m_lsmash, 
		'm_smash': m_smash,  
		'm_dtw': m_dtw,
		's_lsmash': s_lsmash,
		's_smash': s_smash,
		's_dtw': s_dtw
	}

	df = pd.DataFrame(data=data)
	print(df)
	df.to_csv('runtime_seqlen_{}.csv'.format(window_size))
