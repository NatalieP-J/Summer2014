from numpy import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
from construction import fromfileplot
from rhomodels import NukerModelGenRho
import sys

GENERATE = True


name = sys.argv[1]
alpha = sys.argv[2]
beta = sys.argv[3]
gamma = sys.argv[4]
M2L = sys.argv[5]
MBH_Msun = sys.argv[6]
rb = sys.argv[7]
mub = sys.argv[8]
rho0 = findrho0(rb,M2L,mub)
model1 = NukerModelGenRho('{0}_1'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
model1.getrho()
Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood = getrate(model1)
MBH_Msun = sys.argv[9]
model2 = NukerModelGenRho('{0}_2'.format(name),alpha,beta,gamma,rb,rho0,MBH_Msun,GENERATE,memo = False)
model2.getrho()
Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood = getrate(model2)




