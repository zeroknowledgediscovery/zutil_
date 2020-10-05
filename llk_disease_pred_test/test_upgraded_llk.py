import pandas as pd
import numpy as np
from glob import glob
from time import time
from subprocess import check_output

from scipy.spatial import distance
import numpy.linalg as LA


def run(method):

	print('\n{}'.format(method))
	time0 = time()

	for pfsa in sorted(glob('*pfsa')):
		prefix = pfsa.rstrip('.pfsa')
		print('\t{}'.format(prefix))
		
		out = prefix + '.out'
		check_output('./llk_v1 -s testfile.dat -f {} -o {} -m {}'.format(pfsa, out, method), shell=True)
	
	time1 = time()
	time_elapsed = (time1 - time0) / 60.
	print('time elapsed = {} min.s'.format(time_elapsed))
	
	
	result = []
	index = []
	for pfsa in sorted(glob('*pfsa')):
		prefix = pfsa.rstrip('.pfsa')
		out = prefix + '.out'
		res = np.genfromtxt(out)
		check_output('rm {}'.format(out), shell=True)
		index.append(prefix)
		result.append(res)
		
	df = pd.DataFrame(index=index, data=result)

	df.to_csv('llk_result_{}.csv'.format(method))

	return time_elapsed
	

if __name__ == '__main__':

	results = {}
	for method in ["stationary", "uniform", "fixed", "random"]:
		run(method)
		results[method] = pd.read_csv('llk_result_{}.csv'.format(method), index_col=0)

	indices = results["stationary"].index
	standard = results["stationary"].values

	for method in ["uniform", "fixed", "random"]:
		print('\n{}'.format(method))
		variation = results[method].values
		diff  = standard - variation

		for index, s, v, d in zip(indices, standard, variation, diff):
			ns = LA.norm(s)
			nv = LA.norm(v)
			nd = LA.norm(d)
			cosine = distance.cosine(s, v)
			print('\t'.format(index))
			print('\tlength of standard = {:.6f}'.format(ns))
			print('\tlength of variation = {:.6f}'.format(nv))
			print('\tdistance = {:.6f}'.format(nd))
			print('\tcosin distance = {:.6f}\n'.format(cosine))
