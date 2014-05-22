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
    
    #computes enclosed mass - right now this has to cycle through r values, can this be more efficient?
    def Menc(self,r):
        try:
            t = r.shape
            Mencs = []
            for i in range(len(r)):
                temp = intg.quad(self.Minterior,0,r[i],full_output=1)
                try:
                    if temp[3]!='':
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
    def psi2(self,r):
        try:
            t = r.shape
            psi2s = []
            for i in range(len(r)):
                temp=intg.quad(self.Minterior,r[i],inf,full_output=1)
                try:
                    if temp[3]!='':
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
        rarray = rarray[:389] #Come up with better handling for this!
        rchange = rarray[len(rarray)-1]
        rstart = rarray[0]
        return rarray, rchange,rstart
    
    #generate a piecewise expression via interpolation and power law approximations at extremes
    def piecewise2(self,r,inter,start,end,lim1,lim2,smallrexp,largerexp,conds=False):
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
        
    def psigood(self,r,smallrexp=-1,largerexp = -1):
        #generate rgrid and its outer boundaries
        rarray,rchange,rstart = self.rgrid()
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
        return m,rstart,rchange
        

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
test6,lim1,lim2 = model.psigood(rtest)
plt.clf()
plt.loglog(rtest,test6,'.')
plt.ylabel(r'$\psi$')
plt.xlabel('r')
plt.xlim(min(rtest),max(rtest))
plt.ylim(min(test6),max(test6))
#interr=53231.5270085 #identifies r value at which interpolation errors began
plt.axvline(lim1, color='r',label='Limits of interpolation')
plt.axvline(lim2, color='r')
#plt.axvline(interr,color='g',label='Start of quad errors')
plt.legend(loc='best')
plt.show()
