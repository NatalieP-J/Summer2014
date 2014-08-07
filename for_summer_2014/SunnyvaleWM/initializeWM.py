import multiprocessing
import subprocess
from numpy import *

def LoadDataTab(fname):
    f=open(fname,'r')
    data=[]
    for line in f.readlines():
        data.append(line.replace('\n','').split('\t'))
    f.close()
    return data

WM = array(LoadDataTab('WM04.dat'))[:,][:-10]
names = WM[:,0]
dists = array([float(i) for i in WM[:,2]])
rbs = 10**array([float(i) for i in WM[:,3]])
mubs = array([float(i) for i in WM[:,4]])
alphas = array([float(i) for i in WM[:,5]])
betas = array([float(i) for i in WM[:,6]]) + 1
gammas = array([float(i) for i in WM[:,7]]) + 1
M2Ls = array([float(i) for i in WM[:,8]])
MBH1s = 10**array([float(i) for i in WM[:,10]])
MBH2s = 10**array([float(i) for i in WM[:,12]])

def worker(i):
	subprocess.Popen(['python','WMrateget.py',names[i],alphas[i],betas[i],gammas[i],M2Ls[i],MBH1s[i],rbs[i],mubs[i],MBH2s[i]])

if __name__ == '__main__':
	for i in range(8):
		p = multiprocessing.Process(target=worker(i))
		p.start()
