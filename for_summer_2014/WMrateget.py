from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from manage import *

try:
    call(['rm -f', 'rhoratefcns.pyc'])
    call(['rm -f', 'construction.pyc'])
    call(['rm -f', 'models.pyc'])
except OSError:
    pass

GENERATE = True


MsunV = 4.83
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600

rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

def getrate(model):
    
    Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
    seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"ON"}
    verbosity = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}

    if model.generate == True:
        call(["mkdir","{0}".format(model.directory)])
        seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}
    
    rarray,rchange,rstart = rgrid([model],4,-6,0.03)
    Mencgood = compute([model],["Menc",Menc],funcMenc,rtest,[4,-6,0.03],rgrid,[3-model.g,0],[['r','M'],False],[seton[Menc],verbosity[Menc]])

    if Mencgood == 0:
        print 'Menc failed'
        return 0,0,0,0,0,0

    if Mencgood != 0:

        psigood = compute([model,'Model'],["psi",psi],funcpsi,rtest,[4.3,-6,0.03],rgrid,[-1,-1],[['r','$\psi$'],False],[seton[psi],verbosity[psi]])
        
        if psigood == 0:
            print 'psi failed'
            return Mencgood,0,0,0,0,0

        if psigood != 0:

            Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
            
            Jc2good = compute(Jprereqs,["Jc2",Jc2],funcJc2,rtest,[3,-4,0.01],Egrid,[-1,-1],[['E','Jc2'],False],[seton[Jc2],verbosity[Jc2]])

            if Jc2good == 0:
                print 'Jc2 failed'
                return Mencgood,psigood,0,0,0,0

            if Jc2good != 0:

                lgprereqs = [model,'Model',psigood,"psi"]

                ggood = compute(lgprereqs,["g",g],funclg,rtest,[3,-3,0.1],Egrid,[model.b-0.5,model.g-0.5],[['E','g'],False],[seton[g],verbosity[g]])

                if ggood == 0:
                    print 'g failed'
                    return Mencgood,psigood,Jc2good,0,0,0

                if ggood != 0:

                    bGprereqs = [model,'Model',psigood, "psi",ggood,"g"]

                    Gtest = arange(-5,5,0.01)
                    Gtest = append(Gtest,40)
                    Gtest = insert(Gtest,0,-40)
                    Gtest = 10**Gtest

                    Ggood = compute(bGprereqs,["G",G],funcbG,Gtest,[3,-3,0.1],Egrid,[model.b-4,model.g-4],[['E','G'],False],[seton[G],verbosity[G]])
                    
                    psibG_memo = {}
                    part2bG_memo = {}
                    part3bG_memo = {}

                    if Ggood == 0:
                        print 'mathcalG failed'
                        return Mencgood,psigood,Jc2good,ggood,0,0

                    if Ggood != 0:

                        fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
                        
                        ftest = arange(-3,5,0.01)
                        ftest = append(ftest,40)
                        ftest = insert(ftest,0,-40)
                        ftest = 10**ftest

                        fgood = compute(fprereqs,["f",f],funcf,ftest,[5,-3,0.03],Egrid,[model.b-1.5,model.g-1.5],[['E','f'],False],[seton[f],verbosity[f]])

                        print('\a')

                        if fgood == 0:
                            print 'f failed'
                            return Mencgood,psigood,Jc2good,ggood,Ggood,0

                        if fgood != 0:

                            return Mencgood,psigood,Jc2good,ggood,Ggood,fgood

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

from models import NukerModelRho
for galaxy in range(4,len(WM)):
    print galaxy+1, ' of ',len(WM)
    name = names[galaxy]
    alpha = alphas[galaxy]
    beta = betas[galaxy]
    gamma = gammas[galaxy]
    M2L = M2Ls[galaxy]
    MBH_Msun = MBH1s[galaxy]
    rb = rbs[galaxy]
    mub = mubs[galaxy]
    rho0 = findrho0(rb,M2L,mub)
    model = NukerModelRho('{0}1'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE)
    Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrate(model)
    MBH_Msun = MBH2s[galaxy]
    model = NukerModelRho('{0}2'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE)
    Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrate(model)





