from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from rateget import *
import matplotlib.pyplot as plt
plt.ion()

u = array([  0.00000000e+00,   1.00000000e-04,   4.00000000e-04,
         9.00000000e-04,   1.60000000e-03,   2.50000000e-03,
         3.60000000e-03,   4.90000000e-03,   6.40000000e-03,
         8.10000000e-03,   1.00000000e-02,   1.21000000e-02,
         1.44000000e-02,   1.69000000e-02,   1.96000000e-02,
         2.25000000e-02,   2.56000000e-02,   2.89000000e-02,
         3.24000000e-02,   3.61000000e-02,   4.00000000e-02,
         4.41000000e-02,   4.84000000e-02,   5.29000000e-02,
         5.76000000e-02,   6.25000000e-02,   6.76000000e-02,
         7.29000000e-02,   7.84000000e-02,   8.41000000e-02,
         9.00000000e-02,   9.61000000e-02,   1.02400000e-01,
         1.08900000e-01,   1.15600000e-01,   1.22500000e-01,
         1.29600000e-01,   1.36900000e-01,   1.44400000e-01,
         1.52100000e-01,   1.60000000e-01,   1.68100000e-01,
         1.76400000e-01,   1.84900000e-01,   1.93600000e-01,
         2.02500000e-01,   2.11600000e-01,   2.20900000e-01,
         2.30400000e-01,   2.40100000e-01,   2.50000000e-01,
         2.60100000e-01,   2.70400000e-01,   2.80900000e-01,
         2.91600000e-01,   3.02500000e-01,   3.13600000e-01,
         3.24900000e-01,   3.36400000e-01,   3.48100000e-01,
         3.60000000e-01,   3.72100000e-01,   3.84400000e-01,
         3.96900000e-01,   4.09600000e-01,   4.22500000e-01,
         4.35600000e-01,   4.48900000e-01,   4.62400000e-01,
         4.76100000e-01,   4.90000000e-01,   5.04100000e-01,
         5.18400000e-01,   5.32900000e-01,   5.47600000e-01,
         5.62500000e-01,   5.77600000e-01,   5.92900000e-01,
         6.08400000e-01,   6.24100000e-01,   6.40000000e-01,
         6.56100000e-01,   6.72400000e-01,   6.88900000e-01,
         7.05600000e-01,   7.22500000e-01,   7.39600000e-01,
         7.56900000e-01,   7.74400000e-01,   7.92100000e-01,
         8.10000000e-01,   8.28100000e-01,   8.46400000e-01,
         8.64900000e-01,   8.83600000e-01,   9.02500000e-01,
         9.21600000e-01,   9.40900000e-01,   9.60400000e-01,
         9.80100000e-01])

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


d = ([model,Jc2good,Ggood,fgood],u)
d1 = ([model1,Jc2good1,Ggood1,fgood1],u)

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

plt.figure()
plt.loglog(u**2,d)
plt.loglog(u**2,d1)
plt.ylabel(r'$\frac{d\gamma}{dlnr_p}$')
plt.xlabel(r'$u^2$')
plt.title('{0} vs {1}'.format(name1,name2))
plt.savefig('NukerRhoGals/compimages/newpsi_{0}_{1}_rate.png'.format(name1,name2))

plt.show()

