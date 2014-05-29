import numpy as np
nround = np.round
nabs = np.abs
from numpy import *
#from astropy import constants as const
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import time
G = 6.67e-11#const.G
realMsun = 2e30#const.M_sun
Rsun = 7e8#const.R_sun
pc = 3e16#const.pc
#pi = np.pi
#exp = np.exp
#log = np.log
class MakeModel:
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=True):
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
        #generate data?
        self.generate = generate
    
    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        return (-r**(-1-self.g))*((1+r**self.a)**((self.g-self.a-self.b)/self.a))*(self.g+self.b*r**self.a)
    #and second derivatives
    def d2rhodr2(self,r):
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
    
    def psi2interior(self,r):
        return self.rho(r)*r
    
    def psi2(self,r,verbose=False):
        try:
            t = r.shape
            psi2s = []
            for i in range(len(r)):
                temp=intg.quad(self.psi2interior,r[i],inf,full_output=1)
                try:
                    if temp[3]!='' and verbose==True:
                        print 'psi2, index =',i,'r = ',r[i],'message = ',temp[3],'\n'
                except IndexError:
                    pass
                psi2s.append(4*pi*temp[0])
            return array(psi2s)
        except AttributeError:    
            return 4*pi*intg.quad(self.psi2interior,r,inf)[0]
    
    #compute psi (potential)
    def psi(self,r):
        return (self.Mnorm/r) + (self.Menc(r)/r) + self.psi2(r,verbose=False)#************************************True)
    
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
        #print psitab
        m = self.piecewise2(r,psiinter,start,end,rstart,rchange,smallrexp,largerexp)

        if plotting==True:
            plt.figure()
            plt.loglog(r,m,'.')
            plt.ylabel(r'$\psi$')
            plt.xlabel('r')
            plt.xlim(min(r),max(r))
            plt.ylim(min(m),max(m))
            interr=53231.5270085 #identifies r value at which interpolation errors began
            plt.axvline(interr,color='g',label='Start of quad errors')
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
        return abs(self.psi(abs(r))-E) #needs abs, b/c occasionally guesses negative

    def rapo(self,E):
        if E**-1 > 0.2:
            rguess = 10*E**-1
        elif E**-1 < 0.2:
            rguess = 0.01*E**-1
        rresult = root(self.rapoimplicit,rguess,args=E)
        if rresult.success == True:
            return abs(rresult.x) #see rapoimplicit
        elif rresult.success == False:
            print 'Failed to evaluate rapo'
            print rresult.message
            return abs(rresult.x)

    def Jc2implicit(self,r,E):
        return self.psigood(r)-E-((self.Mencgood(r)+self.Mnorm)/(2*r))
    
    def Jc2(self,E):
        rresult = root(self.Jc2implicit,0.01*E**-1,args=E)
        return (self.Mencgood(rresult.x)+self.Mnorm)*rresult.x
        
    def ginterior(self,r,E):
        return (self.drhodr(1./r))*(r**-2)*((sqrt(E-self.psi(1./r)))**-1)
    def funcg(self,E,verbose = False):
        try:
            t = E.shape
            gans = []
            for i in range(len(E)):
                #print i+1, 'of', len(E)
                rapoval = self.rapo(E[i]) #THIS STEP IS SLOW
                temp = intg.quad(self.ginterior,0,1./rapoval,args = E[i],full_output=1)
                t = temp[0]
                try:
                    if temp[3]!='' and verbose==True:
                        print 'g, E = ',E[i],'message = ',temp[3],'\n'
                        t = 1
                except IndexError:
                    pass
                gans.append(-pi*t)
            return array(gans)
        except AttributeError:
            rapoval = self.rapo(E)
            return -pi*intg.quad(self.ginterior,0,1./rapoval,args = E)[0]

    def ggood(self,E,smallrexp = 0,largerexp = 0,plotting=False):
        smallrexp = self.b-0.5
        largerexp = self.g-0.5
        Earray,Echange,Estart = self.Egrid(5,-3,0.1)
        gtab = self.funcg(Earray,verbose=True)
        #for i in range(len(gtab)):
        #    print 'E = ',Earray[i],' g = ',gtab[i]
        gint = interp1d(log10(Earray),log10(gtab))
        start = gtab[0]
        end = gtab[len(Earray)-1]
        #print gtab
        m = self.piecewise2(E,gint,start,end,Estart,Echange,smallrexp,largerexp)
        #print m
        if plotting==True:
            plt.clf()
            plt.loglog(E,m,'.')
            #plt.loglog(Earray,gtab,'ro',label = 'Evaluation points')
            plt.ylabel(r'g')
            plt.xlabel('E')
            plt.xlim(min(E),max(E))
            plt.ylim(min(m),max(m))
            plt.axvline(Estart, color='r',label='Limits of interpolation')
            plt.axvline(Echange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m

    def mathcalGinterior(self,theta,r,E):                              
        return (r**2/sqrt(self.psi(r)-E))*(sqrt(theta)**-1 - np.sqrt(theta))*self.funcg(self.psi(r)*(1-theta)+E*theta)
    
    def mathcalG(self,E,verbose = False):
        print 'starting G evaluation'
        try:
            t = E.shape
            Gans = []
            for i in range(len(E)):
                rapoval = self.rapo(E[i])
                temp = intg.dblquad(self.mathcalGinterior,0,rapoval,lambda r:0,lambda r:1,args=(E[i],))
                t = temp[0]
                try:
                    if temp[3]!='' and verbose==True:
                        print 'G, E = ',E[i],'message = ',temp[3]
                        t = 1
                except IndexError:
                    pass
                Gans.append(t)
            return array(Gans)
        except AttributeError:
            rapoval = self.rapo(E)
            print rapoval
            return intg.dblquad(self.mathcalGinterior,0,rapoval,lambda r:0,lambda r:1,args = (E,))[0]

    def mathcalGgood(self,E,smallrexp = 0,largerexp = 0,plotting = False):
        smallrexp = self.b-4
        largerexp = self.g-4
        Earray,Echange,Estart = self.Egrid(3,-4,0.1)
        Gtab = self.mathcalG(Earray,verbose=True)
        Gint = interp1d(log10(Earray),log10(Gtab))
        start = Gtab[0]
        end = Gtab[len(Earray)-1]
        m = self.piecewise2(E,Gint,start,end,Estart,Echange,smallrexp,largerexp)
        if plotting == True:
            plt.clf()
            plt.loglog(E,m,'.')
            plt.ylabel(r'G')
            plt.xlabel('E')
            plt.xlim(min(E),max(E))
            plt.ylim(min(m),max(m))
            plt.axvline(Estart, color='r',label='Limits of interpolation')
            plt.axvline(Echange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m
    
    def finterior1(self,r,E,rapoval):
        var = rapoval/r
        result = (var**2)*(1./sqrt(E-self.psi(var)))*(var/self.d2rhodr2(var))
        return result
    
    def finterior2(self,r,E,rapoval):
        var = rapoval/r
        return (var**2)*(1./sqrt(E-self.psi(var)))*self.drhodr(var)

    def finterior3(self,r,E,rapoval):
        var = rapoval/r
        return (var**2)*(1./sqrt(E-self.psi(var)))*(1./(2*rapoval))*self.drhodr(var)*((self.Mnorm*(r-1)+r*self.Menc(var) - self.Menc(rapoval))/(E-self.psi(var)))
                   
    def funcf(self,E,verbose=False):
        print 'starting f evaluation'
        epsilon = 1e-7
        try:
            t = E.shape
            fans = []
            for i in range(len(E)):
                print i+1, ' of ', len(E)
                rapoval = self.rapo(E[i])
                prefactor = (1./(sqrt(8)*pi**2*(self.Mnorm + self.Menc(rapoval))))
                temp1 = intg.quad(self.finterior1,epsilon,1-epsilon,args=(E[i],rapoval),full_output = 1)
                temp2 = intg.quad(self.finterior2,epsilon,1-epsilon,args=(E[i],rapoval),full_output = 1)
                temp3 = intg.quad(self.finterior3,epsilon,1-epsilon,args=(E[i],rapoval),full_output = 1)
                t = temp1[0] + temp2[0] + temp3[0]
                try:
                    if verbose==True:
                        if temp1[3] != '':
                            print 'f, E = ',E[i],'message = ',temp1[3]
                        elif temp2[3] != '':
                            print 'f, E = ',E[i],'message = ',temp2[3]
                        elif temp3[3] != '':
                            print 'f, E = ',E[i],'message = ',temp3[3]
                        t = -1
                except IndexError:
                    pass
                fans.append(prefactor*t)
            return array(fans)
        except AttributeError:
            rapoval = self.rapo(E)
            prefactor = (1./(sqrt(8)*pi**2*(self.Mnorm + self.Menc(rapoval))))
            temp1 = intg.quad(self.finterior1,0,1,args=(E,rapoval))
            temp2 = intg.quad(self.finterior2,0,1,args=(E,rapoval))
            temp3 = intg.quad(self.finterior3,0,1,args=(E,rapoval))
            return prefactor*(temp1[0] + temp2[0] + temp3[0])
                
    def fgood(self,E,smallrexp = 0, largerexp = 0,plotting = False):
        smallrexp = self.b-1.5
        largerexp = self.g-1.5
        Earray,Echange,Estart = self.Egrid(5,-3,0.03)
        ftab = self.funcf(Earray,verbose=True)
        fint = interp1d(log10(Earray),log10(ftab))
        start = ftab[0]
        end = ftab[len(Earray)-1]
        m = self.piecewise2(E,fint,start,end,Estart,Echange,smallrexp,largerexp)
        if plotting==True:
            plt.clf()
            plt.loglog(E,m,'.')
            plt.loglog(Earray,gtab,'ro',label = 'Evaluation points')
            plt.ylabel(r'g')
            plt.xlabel('E')
            plt.xlim(min(E),max(E))
            plt.ylim(min(m),max(m))
            plt.axvline(Estart, color='r',label='Limits of interpolation')
            plt.axvline(Echange, color='r')
            plt.legend(loc='best')
            plt.show()
        return m


model = MakeModel('testing',1.,4.,1.5,1.,1.e5,1000,generate = False)
rtest = arange(-12,12,0.01)
rtest = 10**rtest
tic = time.clock()
test = model.ggood(rtest,plotting=True)
toc = time.clock()
print 'g runs in ', toc-tic, 's'
'''

if model.generate == True:
    psi = model.psi(rtest)
    savetxt('{0}_psi.dat'.format(model.name),psi)
    print 'made psi'
    g = model.ggood(rtest)
    savetxt('{0}_g.dat'.format(model.name),g)
    print 'made g'

elif model.generate == False:
    psi = loadtxt('{0}_psi.dat'.format(model.name))
    g = loadtxt('{0}_g.dat'.format(model.name))

else:
    print 'Invalid choice for "generate"'

def locate(somelist,checker,someval):
    return somelist[abs(checker)==min(abs(checker))]

'''
def locate(somelist,checker,someval):
    idx = (nabs(checker-someval)).argmin()
    #if someval - checker[idx] < 0:
    #    print 'someval = ', someval, 'checker = ', checker[idx], 'diff = ', someval-checker[idx]
    return somelist[idx]
'''
def mathcalGinterior1(theta,r,E):
    #print 'called integrand',r,theta
    rcheck = rtest-r
    #rcheck = rtest
    psir = locate(psi,rcheck,r)
    #if psir-E <= 0:
    #    print 'psir = ',psir, 'E = ', E, 'diff = ', psir-E
    Echeck = rtest-(psir*(1-theta)+E*theta)
    #Echeck = rtest
    gE = locate(g,Echeck,(psir*(1-theta)+E*theta))
    #print E,psir,r,theta
    part1 = (r**2/sqrt(abs(psir-E)))
    part2 = (sqrt(theta)**-1)
    return part1*part2*gE

def mathcalGinterior2(theta,r,E):
    #print 'called integrand',r,theta
    rcheck = rtest-r
    psir = locate(psi,rcheck,r)
    Echeck = rtest-(psir*(1-theta)+E*theta)
    gE = locate(g,Echeck,(psir*(1-theta)+E*theta))
    #print E,psir,r,theta
    part1 = (r**2/sqrt(psir-E))
    part2 = (-sqrt(theta))
    return part1*part2*gE
    
def mathcalG(E,verbose = False):
    print 'starting G evaluation'
    eps = 1e-7
    try:
        t = E.shape
        Gans = []
        for i in range(len(E)):
            print i+1, ' of ', len(E)
            rapoval = model.rapo(E[i])
            temp1 = intg.dblquad(mathcalGinterior1,0,rapoval,lambda r:0,lambda r:1,args=(E[i],))
            temp2 = intg.dblquad(mathcalGinterior2,0,rapoval,lambda r:0,lambda r:1,args=(E[i],))
            t = temp1[0] + temp2[0]
            try:
                if temp[3]!='' and verbose==True:
                    print 'G, E = ',E[i],'message = ',temp[3]
                    t = 1
            except IndexError:
                pass
            Gans.append(t)
        return array(Gans)
    except AttributeError:
        rapoval = model.rapo(E)
        #print rapoval
        temp1 = intg.dblquad(mathcalGinterior1,0,rapoval,lambda r:0,lambda r:1,args=(E,))
        print 'done 1'
        temp2 = intg.dblquad(mathcalGinterior2,0,rapoval,lambda r:0,lambda r:1,args=(E,))
        print 'done 2'
        t = temp1[0] + temp2[0]
        return t

print mathcalG(1e3,verbose=True)
'''
