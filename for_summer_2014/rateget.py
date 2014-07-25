from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages

DISPLAY = True
GENERATE = False

alpha = 7.52
beta = 3.13
gamma = 1.98
r0pc = 1.
rb = 10**2.38
r0pc = rb
mub = 19.98
M2L = 6.27
MsunV = 4.83
rho0 = findrho0(rb,M2L,mub)
MBH_Msun = 10**6.04
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
    seton = {Menc:"ON",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
    verbosity = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}

    if model.generate == True:
        call(["mkdir","{0}".format(model.directory)])
        seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}

    #pp = PdfPages('{0}/Master.pdf'.format(model.directory))
    
    rarray,rchange,rstart = rgrid([model],4,-6,0.03)
    Mencgood = compute([model],["Menc",Menc],funcMenc,rtest,[4,-6,0.03],rgrid,[3-model.g,0],[['r','M'],False],[seton[Menc],verbosity[Menc]])

    #pp.savefig()

    psigood = compute([model,'Model'],["psi",psi],funcpsi,rtest,[4.3,-6,0.03],rgrid,[-1,-1],[['r','$\psi$'],False],[seton[psi],verbosity[psi]])

    #pp.savefig()
    
    Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
    
    Jc2good = compute(Jprereqs,["Jc2",Jc2],funcJc2,rtest,[3,-4,0.01],Egrid,[-1,-1],[['E','Jc2'],False],[seton[Jc2],verbosity[Jc2]])

    #pp.savefig()

    lgprereqs = [model,'Model',psigood,"psi"]

    ggood = compute(lgprereqs,["g",g],funclg,rtest,[3,-3,0.1],Egrid,[model.b-0.5,model.g-0.5],[['E','g'],False],[seton[g],verbosity[g]])
    
    #pp.savefig()

    bGprereqs = [model,'Model',psigood, "psi",ggood,"g"]

    Gtest = arange(-5,5,0.01)
    Gtest = append(Gtest,40)
    Gtest = insert(Gtest,0,-40)
    Gtest = 10**Gtest

    Ggood = compute(bGprereqs,["G",G],funcbG,Gtest,[3,-3,0.1],Egrid,[model.b-4,model.g-4],[['E','G'],False],[seton[G],verbosity[G]])
    
    psibG_memo = {}
    part2bG_memo = {}
    part3bG_memo = {}

    #pp.savefig()

    fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
    
    ftest = arange(-3,5,0.01)
    ftest = append(ftest,40)
    ftest = insert(ftest,0,-40)
    ftest = 10**ftest

    fgood = compute(fprereqs,["f",f],funcf,ftest,[5,-3,0.03],Egrid,[model.b-1.5,model.g-1.5],[['E','f'],False],[seton[f],verbosity[f]])

    #pp.savefig()

    print('\a')
    
    #pp.close()

    if DISPLAY == True:
        plt.show()

    return Mencgood,psigood,Jc2good,ggood,#Ggood,fgood


#choose model from rho or genrho options
#for example:
from models import NukerModelRho
model = NukerModelRho('NGC4467',alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE)
Mencgood,psigood,Jc2good,ggood = getrate(model)



