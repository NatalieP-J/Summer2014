from numpy import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
from construction import fromfileplot
from rhomodels import NukerModelGenRho
import sys


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

GENERATE = True

name = names[i]
alpha = alphas[i]
beta = betas[i]
gamma = gammas[i]
M2L = M2Ls[i]
MBH_Msun = MBH1s[i]
rb = rbs[i]
mub = mubs[i]
rho0 = findrho0(rb,M2L,mub)
model1 = NukerModelGenRho('{0}_1'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
model1.getrho()
Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood = getrate(model1)
MBH_Msun = MBH2s[i]
model2 = NukerModelGenRho('{0}_2'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
model2.getrho()
Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood = getrate(model2)



