from numpy import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
from construction import fromfileplot

GENERATE = False
Menc,psi,Jc2,g,G,f,rate = 0,1,2,3,4,5,6
partial = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF",rate:"ON"}

Gtemp = 6.67e-11
ctemp = 3e8
Msuntemp = 2e30

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
Rs1s = []
Rs2s = []
from rhomodels import NukerModelGenRho
for galaxy in range(6):#len(WM)):
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
    model1 = NukerModelGenRho('{0}_1'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
    Rs1 = 2*Gtemp*MBH_Msun*Msuntemp/(ctemp**2)
    Rs1s.append(Rs1/model1.rT)
    #model1.getrho()
    #Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood1 = getrate(model1,partial)

    MBH_Msun = MBH2s[galaxy]
    model2 = NukerModelGenRho('{0}_2'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
    Rs2 = 2*Gtemp*MBH_Msun*Msuntemp/(ctemp**2)
    Rs2s.append(Rs2/model2.rT)
    #model2.getrho()
    #Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood2 = getrate(model2,partial)
'''
    utest = arange(-2,0)
    utest = 10**utest
    try:
        plt.figure()
        plt.loglog(utest**2, 10**rategood1(log10(utest)))
        plt.loglog(utest**2, 10**rategood1(log10(utest)))
    except (ValueError,TypeError):
        pass
'''



