from numpy import *
import math
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
MsunV = 4.83

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
        self.directory = '{0}_Nuker_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
    
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
    
    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def rho(self,r,verbose = False):
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
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

    def drhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,verbose = verbose)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

    def d2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,verbose = verbose)[0]


class NukerModelGenRho:
    """
    Call <modelname>.getrho() after model has been initialized as <modelname>
    """
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
        self.directory = '{0}_Nuker_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)
    
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
    
    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

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
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

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
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,verbose = verbose)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

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
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,verbose = verbose)[0]
    
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
            tab2 = abs(self.funcdrhodr(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"wb")
            pickle.dump(tab2,pklrfile)
            pklrfile.close()
            tab3 = abs(self.funcd2rhodr2(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"wb")
            pickle.dump(tab3,pklrfile)
            pklrfile.close()
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
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
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
            print 'Loaded densities and interpolated'
        
    def rho(self,r):
        return 10**self.inter0(log10(r))
    
    def drhodr(self,r):
        return -10**self.inter1(log10(r))
    
    def d2rhodr2(self,r):
        return 10**self.inter2(log10(r))


class NukerModelRho:
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
        #dynamical timescale
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #start a new directory?
        self.generate = generate
        #directory name
        self.directory = '{0}_Nuker_a{1}_b{2}_g{3}_MBH{4}'.format(self.name,self.a,self.b,self.g,self.MBH)

    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        return (-r**(-1-self.g))*((1+r**self.a)**((self.g-self.a-self.b)/self.a))*(self.g+self.b*r**self.a)
    #and second derivatives
    def d2rhodr2(self,r):
        part1 = r**(-2-self.g)
        part2 = (1+r**self.a)**((self.g-self.b-2*self.a)/self.a)
        part3a = self.b*(1+self.b)*r**(2*self.a)
        part3b = self.g + self.g**2
        part3c = (self.b - (self.a*self.b) + self.g + (self.a*self.g) + (2*self.b*self.g))*r**self.a
        part3 = part3a + part3b + part3c
        return part1*part2*part3

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
        self.directory = '{0}_Sersic_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)

    def I(self,r):
        return exp(-self.b*(r**(1./self.n)))

    def dIdR(self,r):
        return -self.I(r)*(self.b/self.n)*r**(-1.+(1./n))

    def dI2dR2(self,r):
        return self.I(r)*(self.b/self.n**2)*(-1+self.n+self.b*r**(1./n))*r**(-2.+(1./n))
    
    def dI3dR3(self,r):
        return self.I(r)*(self.b/self.n**3)*(-1+3*self.n-2*self.n**2-3*self.b*(-1+n)*r**(1./n)-k**2*r**(2./n))*r**(-3+(1./n))

    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

    def rho(self,r,verbose = False):
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
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

    def drhodr(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,verbose = verbose)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

    def d2rhodr2(self,r,verbose = False):
        try:
            r.shape
            rargs = [tuple((i,)) for i in r]
            lows = zeros(len(r))
            highs = zeros(len(r)) + pi/2. 
        except (AttributeError,TypeError):
            rargs = tuple((r,))
            lows = 0
            highs = pi/2.
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,verbose = verbose)[0]

class SersicModelGenRho:
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
        self.directory = '{0}_Sersic_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)

    def I(self,r):
        return exp(-self.b*(r**(1./self.n)))

    def dIdR(self,r):
        return -self.I(r)*(self.b/self.n)*r**(-1.+(1./n))

    def dI2dR2(self,r):
        return self.I(r)*(self.b/self.n**2)*(-1+self.n+self.b*r**(1./n))*r**(-2.+(1./n))
    
    def dI3dR3(self,r):
        return self.I(r)*(self.b/self.n**3)*(-1+3*self.n-2*self.n**2-3*self.b*(-1+n)*r**(1./n)-k**2*r**(2./n))*r**(-3+(1./n))

    #use luminosity density to compute rho
    def rhointerior(self,theta,r):
        return self.dIdR(r/cos(theta))/cos(theta)

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
    
    #and its first
    def drhodrinterior(self,theta,r):
        return self.d2IdR2_2(r/cos(theta))/cos(theta)**2

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
        return (1./pi)*integrator(r,[self.drhodrinterior,'drhodr'],lows,highs,args = rargs,verbose = verbose)[0]
    
    #and second derivative
    def d2rhodr2interior(self,theta,r):
        return self.d3IdR3_2(r/cos(theta))/cos(theta)**3

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
        return (1./pi)*integrator(r,[self.d2rhodr2interior,'d2rhodr2'],lows,highs,args = rargs,verbose = verbose)[0]
    
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
            tab2 = abs(self.funcdrhodr(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'drhodr'),"wb")
            pickle.dump(tab2,pklrfile)
            pklrfile.close()
            tab3 = abs(self.funcd2rhodr2(rtest))
            pklrfile = open('{0}/{1}.pkl'.format(self.directory,'d2rhodr2'),"wb")
            pickle.dump(tab3,pklrfile)
            pklrfile.close()
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
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
            self.inter0 = interp1d(log10(rtest),log10(tab1))
            self.inter1 = interp1d(log10(rtest),log10(tab2))
            self.inter2 = interp1d(log10(rtest),log10(tab3))
            print 'Loaded densities and interpolated'
        
    def rho(self,r):
        return 10**self.inter0(log10(r))
    
    def drhodr(self,r):
        return -10**self.inter1(log10(r))
    
    def d2rhodr2(self,r):
        return 10**self.inter2(log10(r))

class SersicModelRho:
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
        self.directory = '{0}_Sersic_n{1}_MBH{2}'.format(self.name,self.n.self.MBH)
    
    #Compute density
    def rho(self,r):
        return (r**-self.p)*exp(-self.b*(r**(1./self.n)))
    #and its first
    def drhodr(self,r):
        pre = (r**-self.p)*exp(-self.b*(r**(1./self.n)))*((self.n*r)**-1)
        post = (self.n*self.p)+self.b*(r**(1./n))
        return pre*post
    #and second derivatives
    def d2rhodr2(self,r):
        pre = (r**-p)*exp(-b*(r**(1./self.n)))*((self.n*r)**-2)
        post = (self.p*(1+self.p)*self.n**2) + self.b*(-1 + self.n + 2*self.n*self.p)*(r**(1./n)) + (b**2)*(r**-p)
        return pre*post
