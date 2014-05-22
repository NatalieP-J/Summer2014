from numpy import *
from astropy import constants as const
import scipy.integrate as intg
from scipy.optimize import root
G = const.G
realMsun = const.M_sun
Rsun = const.R_sun
pc = const.pc
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
        return 4*pi*intg.quad(self.Minterior,0,r)[0]
    #function to go in the solver for finding rH (its implicit definition)
    def rHimplicit(self,r):
        return self.Mnorm-Menc(r)
    #solve for rH
    def rH(self):
        rresult=root(self.rHimplicit,1e-4)
        return rresult.x
    
    #def psi2interior(self,r):
    #    return self.rho(exp(r))*exp(2*r)
    
    #LINDA'S USES LOGS, IS THIS EQUIVALENT? (no warnings generated this way, unlike with logs)
    #compute part 2 of psi
    def psi2(self,r):
        return 4*pi*intg.quad(self.Minterior,r,np.inf)[0]
    #compute psi (potential)
    def psi(self,r):
        return (self.Mnorm/r) + (self.Menc(r)/r) + self.psi2(r)
    #generate rgrid
    def rgrid(self,upstep=5,downstep=-5,step=0.03):
        rmin = min([self.rH,1])
        rmax = max([self.rH,1])
        rimin = log10(rmin)+downstep
        rimax = log10(rmax)+upstep
        dri = step
        rarray = np.array(rimin,rimax,dri)
        rarray = 10**rarray
        rchange = rarray(len(rarray))
        return rarray, rchange
    #generate psigood via interpolation and power law approximations at extremes
    def psigood(self,largerexp = -1):
        rarray,rchange = self.rgrid()
        psitab = self.psi(rarray)
        construct = np.column_stack((log10(rarray),log10(psitab)))
        
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
test6 = model.psi(1)
