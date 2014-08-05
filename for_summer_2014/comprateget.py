from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
import matplotlib.pyplot as plt
plt.ion()

#create piecewise fcn for genrho

def rarray(bottom,top,step,edgy1,edgy2):
	vals = arange(bottom,top,step)
	vals = append(vals,edgy1)
	vals = insert(rtest,0,edgy2)
	return 10**vals

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

alpha = 1.0
beta = 4.0
gamma = 1.5
r0pc = 1.0
rho0 = 1e5
MBH_Msun = 1e3
name = 'finaltestform'

from rhomodels import NukerModelRho
model = NukerModelRho('{0}1'.format(name),alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE,memo = True)
name1 = str(model).split(' ')[0][11:]
Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrate(model)
from rhomodels import NukerModelGenRho
model1 = NukerModelGenRho('{0}2'.format(name),alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE,memo = True)
name2 = str(model1).split(' ')[0][11:]
model1.getrho()
Mencgood1,psigood1,Jc2good1,ggood1,Ggood1,fgood1 = getrate(model1)

rtest2 = arange(-20,20,0.01)
rtest2 = 10**rtest2

plt.figure()
plt.loglog(rtest2,model.rho(rtest2))
plt.loglog(rtest2,model1.rho(rtest2))
plt.ylabel(r'$\rho$')
plt.xlabel(r'$r$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_rho.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,abs(model.drhodr(rtest2)))
plt.loglog(rtest2,abs(model1.drhodr(rtest2)))
plt.ylabel(r'$\rho$')
plt.xlabel(r'$r$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_drhodr.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,model.d2rhodr2(rtest2))
plt.loglog(rtest2,model1.d2rhodr2(rtest2))
plt.ylabel(r'$\rho$')
plt.xlabel(r'$r$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_d2rhodr2.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**Mencgood(log10(rtest2)))
plt.loglog(rtest2,10**Mencgood1(log10(rtest2)))
plt.ylabel(r'$M_{enc}$')
plt.xlabel(r'$r$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_Menc.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**psigood(log10(rtest2)))
plt.loglog(rtest2,10**psigood1(log10(rtest2)))
plt.ylabel(r'$\psi$')
plt.xlabel(r'$r$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_psi.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**Jc2good(log10(rtest2)))
plt.loglog(rtest2,10**Jc2good1(log10(rtest2)))
plt.ylabel(r'$J_c^2$')
plt.xlabel(r'$E$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_Jc2.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**ggood(log10(rtest2)))
plt.loglog(rtest2,10**ggood1(log10(rtest2)))
plt.ylabel(r'$g$')
plt.xlabel(r'$E$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_lg.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**Ggood(log10(rtest2)))
plt.loglog(rtest2,10**Ggood1(log10(rtest2)))
plt.ylabel(r'$G$')
plt.xlabel(r'$E$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_bg.png'.format(name1,name2))

plt.figure()
plt.loglog(rtest2,10**fgood(log10(rtest2)))
plt.loglog(rtest2,10**fgood1(log10(rtest2)))
plt.ylabel(r'$f$')
plt.xlabel(r'$E$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_f.png'.format(name1,name2))

plt.show()

