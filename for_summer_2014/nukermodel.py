from numpy import *
#from astropy import constants as const
import scipy.integrate as intg
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar
from scipy.optimize import golden
from scipy.optimize import brenth
from scipy.optimize import newton
from scipy.optimize import broyden1
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
G = 6.67e-11#const.G
realMsun = 2e30#const.M_sun
Rsun = 7e8#const.R_sun
pc = 3e16#const.pc
#pi = np.pi
#exp = np.exp
#log = np.log
class MakeModel:
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun):
        #name of model
        self.name = model_name
        #Nuker fit parameters
        self.a = alpha
        self.b = beta
        self.g = gamma
        #starting radius ************************************************
        self.r0 = r0pc
        #starting density
        self.rho0 = rho0
        #black holes mass in units of Msun
        self.MBH = MBH_Msun
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale (currently unused)
        self.tdyn0 = ((G*self.rho0*realMsun)/pc**3)**(-1./2)
    
    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        return (self.g+r**(self.a*self.b))*(r**(-1-self.g))*(1+r**self.a)**((self.g-self.a-self.b)/self.a)
    #and second derivatives
    def drho2dr2(self,r):
        return(r**(-2-self.g))*((1+r**self.a)**((-2*self.a-self.b+self.g)/self.a))*(self.b*(1+self.b)*r**(2*self.a)+self.g+self.g**2+(self.b-self.a*self.b+self.g+self.a*self.g+2*self.b*self.g)*r**self.a)
    
    #a function describing interior mass when integrated
    def Minterior(self,r):
        return self.rho(r)*r**2
    
    #computes enclosed mass - right now this has to cycle through r values, can this be more efficient?
    def Menc(self,r,verbose=False):
        try:
            t = r.shape
            Mencs = []
            for i in range(len(r)):
                temp = intg.quad(self.Minterior,0,r[i],full_output=1)
                try:
                    if temp[3]!='' and verbose==True:
                        print 'Menc, r = ',r[i],'message = ',temp[3],'\n'
                except IndexError:
                    pass
                Mencs.append(4*pi*temp[0])
            return array(Mencs)
        except AttributeError:
            return 4*pi*intg.quad(self.Minterior,0,r)[0]
    
    #function to go in the solver for finding rH (its implicit definition)
    def rHimplicit(self,r):
        return self.Mnorm-self.Menc(r)
    #solve for rH
    def rH(self):
        rresult=root(self.rHimplicit,1e-4)
        return rresult.x
    
    #def psi2interior(self,r):
    #    return self.rho(exp(r))*exp(2*r)
    
    #LINDA'S USES LOGS, IS THIS EQUIVALENT? (no warnings generated this way, unlike with logs)
    
    #compute part 2 of psi - right now this has to cycle through r values, can this be more efficient?
    def psi2(self,r,verbose=False):
        try:
            t = r.shape
            psi2s = []
            for i in range(len(r)):
                temp=intg.quad(self.Minterior,r[i],inf,full_output=1)
                try:
                    if temp[3]!='' and verbose==True:
                        print 'psi2, index =',i,'r = ',r[i],'message = ',temp[3],'\n'
                except IndexError:
                    pass
                psi2s.append(4*pi*temp[0])
            return array(psi2s)
        except AttributeError:    
            return 4*pi*intg.quad(self.Minterior,r,inf)[0]
    
    #compute psi (potential)
    def psi(self,r):
        return (self.Mnorm/r) + (self.Menc(r)/r) + self.psi2(r)
    
    #generate rgrid
    def rgrid(self,upstep=5,downstep=-5,step=0.03):
        rmin = min([self.rH(),[1.]])
        rmax = max([self.rH(),[1.]])
        rimin = log10(rmin)+downstep
        rimax = log10(rmax)+upstep
        dri = step
        rarray = arange(rimin,rimax,dri)
        rarray = 10**rarray
        #rarray = rarray[:389] #Come up with better handling for this!
        rchange = rarray[len(rarray)-1]
        rstart = rarray[0]
        return rarray,rchange,rstart

    def Egrid(self,upstep=5,downstep=-3,step=0.1):
        rmin = min([self.rH(),[1.]])
        rmax = max([self.rH(),[1.]])
        rmin = rmin[0] #mystery problem here ****************************
        rmax = rmax[0] #why do they need to be manually unpacked?
        eimin = log10(self.Menc(rmax)/rmax) + downstep
        eimax = log10(self.Mnorm/rmin) + upstep
        dei = step
        Earray = arange(eimin,eimax,dei)
        Earray = 10**Earray
        Echange = Earray[len(Earray)-1]
        Estart = Earray[0]
        return Earray,Echange,Estart
    
    #generate a piecewise expression via interpolation and power law approximations at extremes
    def piecewise2(self,r,inter,start,end,lim1,lim2,smallrexp,largerexp,conds=False):
        if conds==False:
        #identify three domains
            set1 = r[(r<lim1)]
            set2 = r[(r>=lim1)&(r<=lim2)]
            set3 = r[(r>lim2)]
        #describe the function on each domain
            piece1 = start*(set1/lim1)**smallrexp
            piece2 = 10**(inter(log10(set2)))
            piece3 = end*(set3/lim2)**largerexp
        #return function across the whole array
            return concatenate((piece1,piece2,piece3))
        if conds!=False:
            print 'your piecewise function has further conditions'
            set1 = r[(r<lim1)]
            set2 = r[(r>=lim1)&(r<=lim2)]
            set3 = r[(r>lim2)]
            piece1 = start*(set1/lim1)**smallrexp
            piece2 = 10**(inter(log10(set2)))
            tip = conds[0]
            if self.b>tip:
                piece3 = end*(set3/lim2)**conds[1]
            elif self.b<tip:
                piece3 = end*(set3/lim2)**conds[2]
            elif self.b == tip:
                piece3 = end + conds[3]*log(set3/lim2)
            return concatenate((piece1,piece2,piece3))
        
    def psigood(self,r,smallrexp=-1,largerexp = -1,plotting=False):
        #generate rgrid and its outer boundaries
        rarray,rchange,rstart = self.rgrid(5,-5,0.03)
        #generate psi values based on rgrid
        psitab = self.psi(rarray)
        #right now this is a useless construction, but may need it, depends on interpolator
        construct = column_stack((log10(rarray),log10(psitab)))
        #interpolate over given range
        psiinter = interp1d(log10(rarray),log10(psitab)) #check this does what you think it does!!!! (when documentation comes back online)
        #identify psi boundaries for power law approximations
        start = psitab[0]
        end = psitab[len(rarray)-1]
        m = self.piecewise2(r,psiinter,start,end,rstart,rchange,smallrexp,largerexp)

        if plotting==True:
            plt.clf()
            plt.loglog(r,m,'.')
            plt.ylabel(r'$\psi$')
            plt.xlabel('r')
            plt.xlim(min(r),max(r))
            plt.ylim(min(m),max(m))
            #interr=53231.5270085 #identifies r value at which interpolation errors began
            #plt.axvline(interr,color='g',label='Start of quad errors')
            plt.axvline(rstart, color='r',label='Limits of interpolation')
            plt.axvline(rchange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m
    
    def Mencgood(self,r,smallerexp=0,largerexp = 0,plotting=False):
        rarray,rchange,rstart = self.rgrid(5,-5,0.03)
        Menctab = self.Menc(rarray)
        Mencinter = interp1d(log10(rarray),log10(Menctab))
        start = Menctab[0]
        end = Menctab[len(rarray)-1]
        smallrexp=3-self.g
        conds = [2,0,3-self.b,4*pi*self.rho(rchange)*(rchange**3)]
        m = self.piecewise2(r,Mencinter,start,end,rstart,rchange,smallrexp,largerexp,conds=conds)
        if plotting==True:
            plt.clf()
            plt.loglog(r,m,'.')
            plt.ylabel(r'M$_{enc}$')
            plt.xlabel('r')
            plt.xlim(min(r),max(r))
            plt.ylim(min(m),max(m))
            plt.axvline(rstart, color='r',label='Limits of interpolation')
            plt.axvline(rchange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m

    def rapoimplicit(self,r,E):
        return self.psi(abs(r))-E #needs abs, b/c occasionally guesses negative

    def rapo(self,E):
        rresult = root(self.rapoimplicit,0.01*E**-1,args=E)
        return abs(rresult.x) #see rapoimplicit

    def Jc2implicit(self,r,E):
        return self.psigood(r)-E-((self.Mencgood(r)+self.Mnorm)/(2*r))
    
    def Jc2(self,E):
        rresult = root(self.Jc2implicit,0.01*E**-1,args=E)
        return (self.Mencgood(rresult.x)+self.Mnorm)*rresult.x
        
    def ginterior(self,r,E):
        return (self.drhodr(1./r)/(r**2))*(E-self.psi(1./r))**-0.5
      
    def funcg(self,E,verbose = False):
        print 'starting g evaluation'
        try:
            t = E.shape
            gans = []
            for i in range(len(E)):
                rapoval = self.rapo(E[i]) #THIS STEP IS SLOW
                print 'rapoval=',rapoval
                print i+1, 'of', len(E), '\n'
                temp = intg.quad(self.ginterior,0,1./rapoval,args = E[i],full_output=1)
                t = temp[0]
                try:
                    if temp[3]!='' and verbose==True:
                        print 'g, E = ',E[i],'message = ',temp[3],'\n'
                        t = 1
                except IndexError:
                    pass
                #print t
                gans.append(-pi*t)
            return array(gans)
        except AttributeError:
            rapoval = self.rapo(E)
            print rapoval
            return -pi*intg.quad(self.ginterior,0,1./rapoval,args = E)[0]

    def ggood(self,E,smallerexp = 0,largerexp = 0,plotting=False):
        smallrexp = self.b-0.5
        largerexp = self.g-0.5
        Earray,Echange,Estart = self.Egrid(5,-3,0.1)
        gtab = self.funcg(Earray,verbose=True)
        gint = interp1d(log10(Earray),log10(gtab))
        start = gtab[0]
        end = gtab[len(Earray)-1]
        #print gtab
        m = self.piecewise2(E,gint,start,end,Estart,Echange,smallrexp,largerexp)
        #print m
        if plotting==True:
            plt.clf()
            plt.loglog(E,m,'.')
            plt.ylabel(r'g')
            plt.xlabel('E')
            plt.xlim(min(E),max(E))
            plt.ylim(min(m),max(m))
            plt.axvline(Estart, color='r',label='Limits of interpolation')
            plt.axvline(Echange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m

    #******************************* 
    #******************************* 
    #******************************* 
    #******************************* def mathcalG(self,E,psigood,ggood):
    #******************************* def distribution(self,E,Mencgood,psigood):
                                

model = MakeModel('testing',1.,4.,1.5,1.,1.e5,1000)
#test1 = 0  #relied on inital dictionary definition in gen_params, now defunct
#test2 = model.rho(1.) #calculate rho
#test3 = model.drho2dr2(1.)
#test4 = model.Menc(0.1)
#test5 = model.psi2(1)
rtest = arange(-12,12,0.01)
rtest = 10**rtest
#test6 = model.psigood(rtest,plotting=True)
#test7 = model.Mencgood(rtest,plotting=True)
test8 = model.ggood(rtest,plotting=True) #still broken, debugging
#test9 = model.funcg(10) #still broken
#test10 = model.Jc2(1.)
