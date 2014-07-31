from numpy import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
from construction import fromfileplot

GENERATE = True

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

from rhomodels import NukerModelGenRho
for galaxy in range(15,16):
    print 'galaxy ',galaxy+1, ' of ',len(WM)
    name = names[galaxy]
    alpha = alphas[galaxy]
    beta = betas[galaxy]
    gamma = gammas[galaxy]
    M2L = M2Ls[galaxy]
    MBH_Msun = MBH1s[galaxy]
    rb = rbs[galaxy]
    mub = mubs[galaxy]
    rho0 = findrho0(rb,M2L,mub)
    model1 = NukerModelGenRho('{0}'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
    model1.getrho()
    Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrate(model1)
    utest = arange(-7,1,0.05)
    utest = 10**utest
    d = dgdlnrp(utest,model1)
    pklwrite('{0}/rrate.pkl'.format(model1.directory),utest)
    pklwrite('{0}/rate.pkl'.format(model1.directory),d)
    #fromfileplot(model1.name,model1.directory)
    #fromfileplot(model2.name,model2.directory)





