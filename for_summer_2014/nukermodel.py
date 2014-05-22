from numpy import *
#from astropy import constants as const
import scipy.integrate as intg
from scipy.optimize import root
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
        self.name = model_name
        self.a = alpha
        self.b = beta
        self.g = gamma
        self.r0 = r0pc
        self.rho0 = rho0
        self.MBH = MBH_Msun
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        self.rT = Rsun*(self.MBH)**(1./3)
        self.r0_rT=(self.r0*pc)/self.rT
        self.tdyn0 = ((G*self.rho0*realMsun)/pc**3)**(-1./2)
    #generate some additional parameters
    #def gen_params(self):
        
        #param_dict = {'MBHnorm':MBHnorm,
        #              'rT':rT,
        #              'r0_rT':r0_rT,
        #              'tdyn0':tdyn0}
        #params = [MBHnorm, rT, r0_rT,tdyn0]
        #return params,param_dict
    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        return (self.g+r**(self.a*self.b))(r**(-1-self.g))*(1+r**self.a)**((self.g-self.a-self.b)/self.a)
    #and second derivatives
    def drho2dr2(self,r):
        return(r**(-2-self.g))*((1+r**self.a)**((-2*self.a-self.b+self.g)/self.a))*(self.b*(1+self.b)*r**(2*self.a)+self.g+self.g**2+(self.b-self.a*self.b+self.g+self.a*self.g+2*self.b*self.g)*r**self.a)
    #a function describing interior mass when integrated
    def Minterior(self,r):
        return self.rho(r)*r**2
    #computes enclosed mass
    def Menc(self,r):
        try:
            t = r.shape
            Mencs = []
            for i in range(len(r)):
                Mencs.append(4*pi*intg.quad(self.Minterior,0,r[i])[0])
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
    #compute part 2 of psi
    def psi2(self,r):
        try:
            t = r.shape
            psi2s = []
            for i in range(len(r)):
                psi2s.append(4*pi*intg.quad(self.Minterior,r[i],inf)[0])
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
        #print 'rmin,rmax=',rmin,rmax
        rimin = log10(rmin)+downstep
        rimax = log10(rmax)+upstep
        dri = step
        rarray = arange(rimin,rimax,dri)
        rarray = 10**rarray
        rchange = rarray[len(rarray)-1]
        rstart = rarray[0]
        return rarray, rchange,rstart
    #generate psigood via interpolation and power law approximations at extremes
    def piecewise2(self,r,inter,start,end,lim1,lim2,smallrexp,largerexp,conds=False):
        set1 = r[(r<lim1)]
        set2 = r[(r>=lim1)&(r<=lim2)]
        set3 = r[(r>lim2)]
        piece1 = start*(set1/lim1)*smallrexp
        piece2 = 10**(inter(log10(set2)))
        piece3 = end*(set3/lim2)**largerexp
        return concatenate((piece1,piece2,piece3))
        
    def psigood(self,r,smallrexp=-1,largerexp = -1):
        rarray,rchange,rstart = self.rgrid()
        #print 'start=',rstart,'end=',rchange
        print r
        print r[0]
        psitab = self.psi(rarray)
        construct = column_stack((log10(rarray),log10(psitab)))
        psiinter = interp1d(log10(rarray),log10(psitab))
        start = psitab[0]
        end = psitab[len(rarray)-1]
        m = self.piecewise2(r,psiinter,start,end,rstart,rchange,smallrexp,largerexp)
        print m
        return m
        

    #******************************* def rapo(self,E,psigood):
    #******************************* def Jc2(self,E,Mencgood,psigood):
    #******************************* def funcg(self,E,psigood):
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
test6 = model.psigood(rtest)
plt.plot(rtest,test6,'.')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\psi$')
plt.xlabel('r')
