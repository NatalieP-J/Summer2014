import numpy as np
nsum = np.sum
from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import time
import datetime
from subprocess import call
import pickle
from besselgen import BesselGen
import matplotlib
plt.ion()
MsunV = 4.83
alpha =7.52
beta = 2.13
gamma = 0.98
r0pc = 1#10**2.38
mub = 19.98
M2L = 6.27
rho0 = 1e5#(1./r0pc)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 
MBH_Msun = 1e3#10**6.04
galname = 'testing'#'NGC4467'
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
generate = True
rtest = arange(-2,2,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

def integrator(vals,fcn,downlim,uplim,tol=1.49e-7,args = [],verbose = True):
    try:
        problems = []
        results = []
        for i in range(len(vals)):
            if args != []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],args = args[i],epsabs = tol,full_output = 1)
            elif args == []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],epsabs = tol,full_output = 1)
            try:
                if temp[3] != '':
                    problems.append(i)
                    if verbose == True:
                        print '{0},\t i = {1},\t message = {2}'.format(fcn[1],i,temp[3])
            except IndexError:
                pass
            results.append(temp[0])
        return array(results),problems
    except TypeError:
        problems = []
        if args != []:
            temp = intg.quad(fcn[0],downlim,uplim,args = args,epsabs = tol,full_output = 1)
        elif args == []:
            temp = intg.quad(fcn[0],downlim,uplim,epsabs = tol,full_output = 1)
        try:
            if temp[3] != '':
                problems.append(i)
                if verbose == True:
                    print '{0},\t i = {1},\t message = {2}'.format(fcn[1],i,temp[3])
        except IndexError:
            pass
        return temp[0],problems

def dblintegrator(vals,fcn,lims,tol=1.49e-7,args = [],verbose = True):
    a,b,c,d = lims
    try:
        results = []
        for i in range(len(vals)):
            if verbose == True:
                print i+1, ' of ',len(vals)
            if args != []:
                temp = intg.dblquad(fcn,a[i],b[i],c,d,args = args[i],epsabs = tol)
            elif args == []:
                temp = intg.dblquad(fcn,a[i],b[i],c,d,epsabs = tol)
            results.append(temp[0])
        return array(results)
    except TypeError:
        if args != []:
            temp = intg.dblquad(fcn,a,b,c,d,args = args,epsabs = tol)
        elif args == []:
            temp = intg.dblquad(fcn,a,b,c,d,epsabs = tol)
        return temp[0]


class NukerModel:
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=False):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.b = beta
        self.g = gamma
        #starting radius
        self.r0 = r0pc
        #starting density
        self.rho0 = rho0
        #black holes mass in units of Msun
        self.MBH = MBH_Msun
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale (currently unused)
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #prefactor in I
        self.factor = 2**((self.b-self.g)/self.a)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'Nuker/{0}_a{1}_b{2}_g{3}_r{4}_rho{5}_MBH{6}'.format(self.name,self.a,self.b,self.g,self.r0,self.rho0,self.MBH)
    
    def I(self,r):
        return self.factor*(r**-self.g)*(1+(r**self.a))**(-(self.b-self.g)/self.a)

    def dIdR_2(self,r):
        return self.I(r)*((-self.g/r)-(self.b-self.g)*(r**(self.a-1))*(1+(r**self.a))**-1)

    def d2IdR2_2(self,r):
        part1 = (-self.g/r)*(self.dIdR_2(r) - (self.I(r)/r))
        part2a = (self.g-self.b)*(r**(self.a-1))*(1+(r**self.a))**-1
        part2b = (self.dIdR_2(r) + (1./r)*((self.a-1) - self.a*(r**self.a)*(1+(r**self.a))**1)*self.I(r))
        part2 = part2a*part2b
        return part1 + part2
    
    def d3IdR3_2(self,r):
        A = self.a
        R = r**self.a
        Rm = r**(self.a-1)
        Rexp = (1+(r**self.a))**-1
        part1 = (-self.g/r)*(self.d2IdR2_2(r) - (2./r)*self.dIdR_2(r) + (2./r**2)*self.I(r))
        part2fac = (self.g-self.b)*(Rm*Rexp)
        part2a = ((A-1)/r - Rm*Rexp*A)*(self.dIdR_2(r) + (1./r)*((A-1)-A*R*Rexp)*self.I(r))
        part2b = self.d2IdR2_2(r)
        part2c = (-1./r**2)*((A-1)-A*R*Rexp)*self.I(r)
        part2d = (A**2/r)*Rm*Rexp*(-1 + R*Rexp)*self.I(r)
        part2e = (1./r)*((A-1)-A*R*Rexp)*self.dIdR_2(r)
        part2 = part2fac*(part2a+part2b+part2c+part2d+part2e)
        return part1 + part2

    #compute surface density first derivative
    def dIdR(self,r):
        part1 = r**(-self.g-1)
        part2 = (1+(r**self.a))**(-((self.b-self.g)/self.a)-1)
        part3 = (self.g - (self.b*(r**self.a)))
        return -self.factor*part1*part2*part3
        #return -self.factor*(r**(-self.g-1))*((1+(r**self.a))**(((-self.b+self.g)/self.a)-1))*(self.g - self.b*(r**self.a))
    
    #and its second
    def d2IdR2(self,r):
        part1 = r**(-2-self.g)
        part2 = (1+(r**self.a))**((-2*self.a - self.b + self.g)/self.a)
        part3 = (r**(2*self.a))*self.b*(1+self.b) + self.g + (self.g**2) + (r**self.a)*(self.b - self.a*self.b + self.g + self.a*self.g + 2*self.b*self.g)
        return self.factor*part1*part2*part3

    #and third derivatives
    def d3IdR3(self,r):
        part1 = (r**(-3-self.g))*((1+(r**self.a))**((-3*self.a - self.b + self.g)/self.a))
        part2a = (-r**(3*self.a))*self.b*(1+self.b)*(2+self.b) - self.g*(1+self.g)*(2+self.g)
        part2b = (r**(2*self.a))*((-1+self.a)*self.b*(4+self.a + 3*self.b) - (2+(self.a**2) + 3*self.a*(1+self.b) + 3*self.b*(2+self.b))*self.g)
        part2c = (r**self.a)*((1+self.a)*(-4+self.a-3*self.g)*(self.g-self.b)*(2+(self.a**2) - 3*self.a*(1+self.g) + 3*self.g*(2+self.g)))
        part2 = part2a + part2b + part2c
        return self.factor*part1*part2
    
    def rhointerior(self,theta,r):
        return self.dIdR_2(r/cos(theta))/cos(theta)

    def funcrho(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2.   
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,verbose = verbose)[0]
        
    def drhodrinterior(self,theta,r):
        return self.dI2dR2_2(r/cos(theta))/cos(theta)**2

    def funcdrhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -(1./pi)*integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,verbose = verbose)[0]
    
    def d2rhodr2interior(self,theta,r):
        return self.dI3dR3_2(r/cos(theta))/cos(theta)**3

    def funcd2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return -integrator(r,[self.rhointerior,'rho'],lows,highs,args = rargs,verbose = verbose)[0]
    
    def oldrho(self,r):
        G = self.g + 1
        B = self.b + 1
        return (r**-G)*(1+r**self.a)**((G-B)/self.a)

    def olddrhodr(self,r):
        G = self.g + 1
        B = self.b + 1
        return (-r**(-1-G))*((1+r**self.a)**((G-self.a-B)/self.a))*(G+B*r**self.a)
    #and second derivatives
    def oldd2rhodr2(self,r):
        G = self.g + 1
        B = self.b + 1
        part1 = r**(-2-G)
        part2 = (1+r**self.a)**((G-B-2*self.a)/self.a)
        part3a = B*(1+B)*r**(2*self.a)
        part3b = G + G**2
        part3c = (B - (self.a*B) + G + (self.a*G) + (2*B*G))*r**self.a
        part3 = part3a + part3b + part3c
        return part1*part2*part3
    
    def getrho(self):
        rtest = arange(-7,7,0.01)
        rtest = append(rtest,40)
        rtest = insert(rtest,0,-40)
        rtest = 10**rtest
        if self.generate == True:
            call(['mkdir','{0}'.format(self.directory)])
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rrho'),"wb")
            pickle.dump(rtest,pklrfile)
            pklrfile.close()
            tab1 = self.funcrho(rtest)
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rho'),"wb")
            pickle.dump(tab1,pklrfile)
            pklrfile.close()
            tab2 = self.funcdrhodr(rtest)
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"wb")
            pickle.dump(tab2,pklrfile)
            pklrfile.close()
            tab3 = self.funcd2rhodr2(rtest)
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"wb")
            pickle.dump(tab3,pklrfile)
            pklrfile.close()
            self.inter1 = interp1d(log10(rtest),log10(tab1))
            self.inter2 = interp1d(log10(rtest),log10(tab2))
            self.inter3 = interp1d(log10(rtest),log10(tab3))
            print 'Saved densities and interpolated'
            
        elif self.generate != True:
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rrho'),"rb")
            rtest = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'rho'),"rb")
            tab1 = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"rb")
            tab2 = pickle.load(pklrfile)
            pklrfile.close()
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"rb")
            tab3 = pickle.load(pklrfile)
            pklrfile.close()
            self.inter1 = interp1d(log10(rtest),log10(tab1))
            self.inter2 = interp1d(log10(rtest),log10(tab2))
            self.inter3 = interp1d(log10(rtest),log10(tab3))
            print 'Loaded densities and interpolated'
        
    def rho(self,r):
        return 10**self.inter1(log10(r))
    
    def drhodr(self,r):
        return 10**self.inter2(log10(r))
    
    def d2rhodr2(self,r):
        return 10**self.inter3(log10(r))


########******************* CONSTRUCT MODEL *******************########

model = NukerModel(galname,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate)

model.getrho()

rhovals = model.funcrho(rtest)
rhovals2 = model.oldrho(rtest)
plt.figure()
plt.title(r'$\alpha = {0},\beta = {1},\gamma = {2}$'.format(model.a,model.b,model.g))
plt.xlabel('r')
plt.ylabel(r'$\rho$')
plt.loglog(rtest[1:-1],rhovals[1:-1],label = r'New $\rho$')
plt.loglog(rtest[1:-1],rhovals2[1:-1],label = r'Old $\rho$')
plt.legend(loc = 'best')
plt.show()


rhovals = abs(model.drhodr(rtest))
rhovals2 = abs(model.olddrhodr(rtest))
plt.figure()
plt.title(r'$\alpha = {0},\beta = {1},\gamma = {2}$'.format(model.a,model.b,model.g))
plt.xlabel('r')
plt.ylabel(r'$\frac{d\rho}{dr}$')
plt.loglog(rtest[1:-1],rhovals[1:-1],label = r'New $\frac{d\rho}{dr}$')
plt.loglog(rtest[1:-1],rhovals2[1:-1],label = r'Old $\frac{d\rho}{dr}$')
plt.legend(loc = 'best')
plt.show()

rhovals = model.d2rhodr2(rtest)
rhovals2 = model.oldd2rhodr2(rtest)
plt.figure()
plt.title(r'$\alpha = {0},\beta = {1},\gamma = {2}$'.format(model.a,model.b,model.g))
plt.xlabel('r')
plt.ylabel(r'$\frac{d^2\rho}{dr^2}$')
plt.loglog(rtest[1:-1],rhovals[1:-1],label = r'New $\frac{d^2\rho}{dr^2}$')
plt.loglog(rtest[1:-1],rhovals2[1:-1],label = r'Old $\frac{d^2\rho}{dr^2}$')
plt.legend(loc = 'best')
plt.show()

