from numpy import *
import math
from subprocess import call
import pickle
from scipy.interpolate import interp1d

Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
MsunV = 4.83

from construction import integrator
 
class NukerModeldIdR:
    """
    Call <modelname>.getrho() after model has been initialized as <modelname>
    """
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=False):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.B = beta
        self.b = beta+1
        self.G = gamma
        self.g = gamma+1
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
        self.directory = 'NukerIGals/{0}_dIdR_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
    
    #compute luminosity density
    def I(self,r):
        return self.factor*(r**-self.g)*(1+(r**self.a))**(-(self.b-self.g)/self.a)
    #and its first
    def dIdR(self,r):
        return -(self.I(r)/(r*(1+r**self.a)))*(self.g + self.b*r**self.a)

    #second
    def d2IdR2(self,r):
        return (self.I(r)/(r*(1+r**self.a))**2)*((r**(2*self.a))*self.b*(1+self.b) + self.g + (self.g**2) + (r**self.a)*(self.b-self.a*self.b + self.g + 2*self.b*self.g))
    
    #and third derivatives
    def d3IdR3(self,r):
        return (self.I(r)/(r*(1+r**self.a))**3)*((-r**(3*self.a))*self.b*(1+self.b)*(2+self.b) - self.g*(1 + self.g)*(2+self.g) + (r**(2*self.a))*((-1+self.a)*self.b*(4+self.a+3*self.b) - (2+(self.a**2) + 3*self.a*(1+self.b) + 3*self.b*(2+self.b))*self.g) + (r**self.a)*((1+self.a)*(-4+self.a-3*self.g)*self.g - self.b*(2+(self.a**2) - 3*self.a*(1+self.g) + 3*self.g*(2+self.g))))

class SersicModeldIdR:
    #initialize variables that constitute the model
    def __init__(self,model_name,n,rho0,Re,M2L,I0,MBH_Msun,r0,generate):
        #model name
        self.name = model_name
        #Sersic index
        self.n = n
        #Effective radius (half-light radius)
        self.Re = Re
        #Scale radius
        self.r0 = r0
        #Mass to light ratio
        self.M2L = M2L
        #Intensity at r0
        self.I0 = I0
        #Mass of the black hole in units of Msun
        self.MBH = MBH_Msun
        #determine other parameters based on n
        if self.n <10. and self.n > 0.6:
            self.b = 2*self.n - (1./3.) + 0.009876/self.n
            self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
        #density at r0
        self.rho0 = M2L*I0*(self.b**(self.n*(1-self.p)))*(math.gamma(2*self.n)/(2*Re*math.gamma(self.n*(3-self.p))))
        #Coulomb logarithm
        self.Lam = self.MBH*0.4
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = 'SersicRhoGals/{0}_GenRho_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)

    def I(self,r):
        return exp(-self.b*(r**(1./self.n)))

    def dIdR(self,r):
        return -self.I(r)*(self.b/self.n)*r**(-1.+(1./n))

    def dI2dR2(self,r):
        return self.I(r)*(self.b/self.n**2)*(-1+self.n+self.b*r**(1./n))*r**(-2.+(1./n))
    
    def dI3dR3(self,r):
        return self.I(r)*(self.b/self.n**3)*(-1+3*self.n-2*self.n**2-3*self.b*(-1+n)*r**(1./n)-k**2*r**(2./n))*r**(-3+(1./n))
