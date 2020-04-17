import numpy as np
import time
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import scipy.stats as stats
import glob
import os

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
    dtw = '/home/yhuang10/D3M/zutil_/dtw'
    datasets = '/home/yhuang10/D3M/zutil_/seqDatasets/'
    
    labels = [0] * 25 + [1] * 25
    num_datasets = 20
    window_size = 10

    r_lsmash, r_dtw = [], [] # the ratio (mean within-class distance) / (mean between-class distance)

    n = 0
    while n < num_datasets:
        i = np.random.randint(100) + 1
        j = np.random.randint(100) + 1
        if i == j:
            # Do not calculate separtion if two classes are identical
            # yes, it is rare, but it does happen.
            continue
        
        n += 1
        os.system('cat {}/seq_{}.dat {}/seq_{}.dat > seq.dat'.format(datasets, i, datasets, j))

        status_lsmash = os.system('{} -f seq.dat -T symbolic -o lsmash.dst -S 0'.format(lsmash))
        status_dtw = os.system('{} seq.dat {} dtw.dst'.format(dtw, window_size))
        
        print('lsmash status = {}, dtw status = {}'.format(status_lsmash, status_dtw))
        
        dist_lsmash = np.genfromtxt('lsmash.dst')
        r = separation(dist_lsmash, labels)
        r_lsmash.append(r)
        
        dist_dtw = np.genfromtxt('dtw.dst')
        r = separation(dist_dtw, labels)
        r_dtw.append(r)

    #     # randomly generated base PFSA
    #     cmd = '../lsmash -f seq.dat -T symbolic -o random.dst -S 0 -F {}'.format(' '.join(glob.glob('../baseChoice_cfgfiles/*')))
    #     ! {cmd}
    #     dist_random = np.genfromtxt('random.dst')
    #     r = separation(dist_random, labels)
    #     r_random.append(r)

        print('round {}, ({}, {}), lsmash = {}, dtw={}'.format(n, i, j, r_lsmash[-1], r_dtw[-1]))
