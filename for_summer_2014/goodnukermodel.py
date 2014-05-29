from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import time
import math
Gconst = 6.67e-11
realMsun = 1.99e30
Rsun = 6.9e8
pc = 3.1e16
Menc = 0
psi = 1
Jc2 = 2
g = 3
G = 4
f = 5
run = {Menc:"ON",psi:"ON",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}

########******************* CONSTRUCTION FUNCTIONS *******************########

def piecewise2(r,inter,start,end,lim1,lim2,smallrexp,largerexp,conds=False):
    set1 = r[(r<lim1)]
    set2 = r[(r>=lim1)&(r<=lim2)]
    set3 = r[(r>lim2)]
    piece1 = start*(set1/lim1)**smallrexp
    piece2 = 10**(inter(log10(set2)))
    if conds==False:
        piece3 = end*(set3/lim2)**largerexp
        return concatenate((piece1,piece2,piece3))
    if conds!=False:
        tip = conds[0]
        if model.b>tip:
            piece3 = end*(set3/lim2)**conds[1]
        elif model.b<tip:
            piece3 = end*(set3/lim2)**conds[2]
        elif model.b == tip:
            piece3 = end + conds[3]*log(set3/lim2)
        return concatenate((piece1,piece2,piece3))

def plotter(r,m,rstart,rchange,labels):
    plt.figure()
    plt.loglog(r,m,'.')
    plt.ylabel(r'{0}'.format(labels[1]))
    plt.xlabel('{0}'.format(labels[0]))
    plt.xlim(min(r),max(r))
    plt.ylim(min(m),max(m))
    plt.axvline(rstart, color='r',label='Limits of interpolation')
    plt.axvline(rchange, color='r')
    plt.legend(loc='best')
    plt.show()

def makegood(func,r,size,grid,smallrexp,largerexp,verbose = False,conds = False,plotting=False):
    rarray,rchange,rstart = grid(size[0],size[1],size[2])
    tab,problems = func(rarray,verbose)
    tab = [i for j, i in enumerate(tab) if j not in problems]
    rarray = [i for j, i in enumerate(rarray) if j not in problems]
    inter = interp1d(log10(rarray),log10(tab))
    start = tab[0]
    end = tab[len(rarray)-1]
    m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp,conds)
    if plotting != False:
        plotter(r,m,rstart,rchange,plotting)
    return m,inter


########******************* MODEL FRAMEWORK *******************########
class MakeModel:
    #initialize variables that constitute our model
    def __init__(self,model_name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate=True):
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
        #black hole mass normalized to galaxy density and radius
        self.Mnorm = self.MBH/(self.rho0*(self.r0)**3)
        #tidal disruption radius
        self.rT = Rsun*(self.MBH)**(1./3)
        #number of tidal radii to span galaxy
        self.r0_rT=(self.r0*pc)/self.rT
        #dynamical timescale (currently unused)
        self.tdyn0 = ((Gconst*self.rho0*realMsun)/pc**3)**(-1./2)
        #generate data?
        self.generate = generate
    
    #compute density
    def rho(self,r):
        return (r**-self.g)*(1+r**self.a)**((self.g-self.b)/self.a)
    #and its first
    def drhodr(self,r):
        return (-r**(-1-self.g))*((1+r**self.a)**((self.g-self.a-self.b)/self.a))*(self.g+self.b*r**self.a)
    #and second derivatives
    '''
    def d2rhodr2(self,r):
        return(r**(-2-self.g))*((1+r**self.a)**((-2*self.a-self.b+self.g)/self.a))*(self.b*(1+self.b)*r**(2*self.a)+self.g+self.g**2+(self.b-self.a*self.b+self.g+self.a*self.g+2*self.b*self.g)*r**self.a)
    '''
    def d2rhodr2(self,r):
        part1 = r**(-2-self.g)
        part2 = (1+r**self.a)**((self.g-self.b-2*self.a)/self.a)
        part3a = self.b*(1+self.b)*r**(2*self.a)
        part3b = self.g + self.g**2
        part3c = (self.b - (self.a*self.b) + self.g + (self.a*self.g) + (2*self.b*self.g))*r**self.a
        part3 = part3a + part3b + part3c
        return part1*part2*part3
                                
########******************* CONSTRUCT MODEL *******************########

model = MakeModel('testing',1.,4.,1.5,1.,1.e5,1.e3,generate = False)
rtest = arange(-12,12,0.01)
rtest = 10**rtest

########******************* ENCLOSED MASS *******************########

def Minterior(r):
    return model.rho(r)*r**2

def funcMenc(r,verbose=False):
    try:
        problems = []
        t = r.shape
        Mencs = []
        for i in range(len(r)):
            temp = intg.quad(Minterior,0,r[i],full_output=1)
            try:
                if temp[3]!='' and verbose==True:
                    print 'Menc, r = ',r[i],'message = ',temp[3],'\n'
                    problems.append(i)
            except (IndexError, TypeError):
                pass
            if temp[0] >= 0:
                Mencs.append(4*pi*temp[0])
            elif temp[0] < 0:
                Mencs.append(Mencs[i-1])
        return array(Mencs),array(problems)
    except AttributeError:
        problem = []
        temp = 4*pi*intg.quad(Minterior,0,r)[0]
        try:
            if temp[3]!='' and verbose==True:
                print 'Menc, r = ',r,'message = ',temp[3],'\n'
                problem = [r]
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t,problem

########******************* RADIUS OF INFLUENCE *******************######## 

def rHimplicit(r):
    return abs(model.Mnorm-Menc(r)[0])

def rH():
    rresult=root(rHimplicit,1e-4)
    if rresult.success == True:
        return rresult.x
    elif rresult.success == False:
        print 'Failed to evaluate rH'
        print rrseult.message
        return rresult.x

########******************* CONSTRUCTION FUNCTIONS *******************########

def rgrid(upstep=5,downstep=-5,step=0.03):
    rmin = min([rH(),[1.]])
    rmax = max([rH(),[1.]])
    rimin = log10(rmin) + downstep
    rimax = log10(rmax) + upstep
    dri = step
    rarray = arange(rimin,rimax,dri)
    rarray = 10**rarray
    rchange = rarray[len(rarray)-1]
    rstart = rarray[0]
    return rarray,rchange,rstart

rarray,rchange,rstart = rgrid(10,-10,0.03)

def Egrid(upstep=5,downstep=-3,step=0.1):
    rmin = min([rH(),[1.]])[0]
    rmax = max([rH(),[1.]])[0]
    eimin = log10(Menc(rmax)[0]/rmax) + downstep
    eimax = log10(model.Mnorm/rmin) + upstep
    dei = step
    Earray = arange(eimin,eimax,dei)
    Earray = 10**Earray
    Echange = Earray[len(Earray)-1]
    Estart = Earray[0]
    return Earray,Echange,Estart

########******************* COMPUTE MENC *******************######## 

if run[Menc] == "ON":
    tic = time.clock()
    Mencvals,Mencgood = makegood(Menc,rtest,[10,-10,0.03],rgrid,3-model.g,0,conds = [2,0,3-model.b,4*pi*model.rho(rchange)*(rchange**3)])#,plotting = ['r','M'])
    toc = time.clock()
    print 'Menc ran in ',toc-tic, 's'

########******************* APOHELION RADIUS (???) *******************######## 

def rapoimplicit(r,E):
    return abs(psi(abs(r),verbose=False)[0]-E)

def rapo(E):
    if E**-1 > 0.2:
        rguess = 10*E**-1
    elif E**-1 < 0.2:
        rguess = 0.01*E**-1
    rresult = root(rapoimplicit,rguess,args=E)
    if rresult.success == True:
        return abs(rresult.x)
    elif rresult.success == False:
        print 'Failed to evaluate rapo'
        print rresult.message
        return abs(rresult.x)

########******************* POTENTIAL *******************######## 
        
def psi2interior(r):
    return model.rho(r)*r

def psi2(r,verbose=False):
    try:
        problems = []
        t = r.shape
        psi2s = []
        for i in range(len(r)):
            temp=intg.quad(psi2interior,r[i],inf,full_output=1)
            try:
                if temp[3]!='' and verbose==True:
                    print 'psi2, index =',i,'r = ',r[i],'message = ',temp[3],'\n'
                    problems.append(r)
            except (IndexError,TypeError):
                pass
            psi2s.append(4*pi*temp[0])
        return array(psi2s),problems
    except AttributeError: 
        problem = []
        temp = 4*pi*intg.quad(psi2interior,r,inf)[0]
        try:
            if temp[3]!='' and verbose==True:
                print 'psi2, r = ',r,'message = ',temp[3],'\n'
                problem = [r]
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t, problem

def funcpsi(r,verbose):
    part1 = (model.Mnorm/r)
    part2 = Menc(r,verbose)
    part3 =  psi2(r,verbose)
    problems = array([])
    if part2[1] != []:
        problems = concatenate((problems,array(part2[1])))
    if part3[1] != []:
        print part3[1]
        problems = concatenate((problems,array(part3[1])))
    return part1 + (part2[0]/r) + part3[0],problems

########******************* COMPUTE PSI *******************######## 

if run[psi] == "ON":
    if run[Menc] == "ON":
        tic = time.clock()
        psivals,psigood = makegood(psi,rtest,[10,-10,0.03],rgrid,-1,-1)#, plotting = ['r','$\psi$'])
        toc = time.clock()
        print 'psi ran in ',toc-tic, 's'
    elif run[Menc] != "ON":
        print 'To compute psi please turn Menc ON'

########******************* CIRCULAR ANGULAR MOMENTUM *******************######## 

def Jc2implicit(r,E,verbose):
    return abs(psigood(abs(r))-E-((Mencgood(abs(r))+model.Mnorm)/(2*r)))

def funcJc2(E,verbose):
    try:
        t = E.shape
        Jcs = []
        problems = []
        for i in range(len(E)):
            rresult = root(Jc2implicit,0.01*E[i]**-1,args=(E[i],verbose))
            if rresult.success == True:
                Jcs.append(((Menc(rresult.x,verbose)[0]+model.Mnorm)*rresult.x)[0])
            elif rresult.success == False:
                print 'Failed to evaluate Jc2'
                print rresult.message
                Jcs.append(((Menc(rresult.x,verbose)[0]+model.Mnorm)*rresult.x)[0])
                problems.append(i)
        return array(Jcs),problems
    except AttributeError:
        rresult = root(Jc2implicit,0.01*E[i]**-1,args=(E[i],verbose))
        problem = []
        if rresult.success == True:
            Jc = ((Menc(rresult.x,verbose)[0]+model.Mnorm)*rresult.x)[0]
        elif rresult.success == False:
            print 'Failed to evaluate Jc2'
            print rresult.message
            Jc = ((Menc(rresult.x,verbose)[0]+model.Mnorm)*rresult.x)[0]
            problem = [E]
        return Jc, problem

########******************* COMPUTE Jc2 *******************######## 

if run[Jc2] == "ON":
    if run[psi] == "ON" and run[Menc] == "ON":
        tic = time.clock()
        Jc2good = makegood(Jc2,rtest,[3,-3,0.01],Egrid,-1,-1,verbose = True,plotting = ['E','Jc2'])
        toc = time.clock()
        print 'Jc2good ran in ',toc-tic, 's'
    elif run[psi] != "ON":
        print 'To compute Jc2 please turn psi ON'
    elif run[Menc] != "ON":
        print 'To compute Jc2 please turn Menc ON'

########******************* g *******************######## 

def ginterior(r,E):
    return (model.drhodr(1./r))*(r**-2)*((sqrt(E-psigood((1./r))))**-1)
    
def funcg(E,verbose=False):
    try:
        t = E.shape
        gans = []
        problems = []
        for i in range(len(E)):
            rapoval = rapo(E[i])
            temp = intg.quad(ginterior,0,1./rapoval,args = E[i],full_output = 1)
            t = temp[0]
            try:
                if temp[3] != '' and verbose == True:
                    print 'g, E = ',E[i], 'message = ', temp[3]
                    problems.append(i)
            except (IndexError,TypeError):
                pass
            gans.append(-pi*t)
        return array(gans),problems
    except (AttributeError,TypeError):
        problem = []
        rapoval = rapo(E)
        temp = intg.quad(ginterior,0,1./rapoval,args = E,full_output = 1)
        t = temp[0]
        try:
            if temp[3] != '' and verbose == True:
                print 'g, E = ',E, 'message = ', temp[3]
                problem = [E]
        except (IndexError,TypeError):
            pass
        return -pi*t, problem

########******************* COMPUTE g *******************######## 

if run[g] == "ON":
    if run[psi] == "ON":
        tic = time.clock()
        gvals,ggood = makegood(funcg,rtest,[5,-3,0.1],Egrid,model.b-0.5,model.g-0.5)#,plotting = ['E','g'])
        toc = time.clock()
        print 'g ran in ',toc-tic, 's'
    elif run[psi] != "ON":
        print 'To compute g, please turn psi ON'

########******************* mathcalG *******************######## 

def Ginterior(theta,r,E,verbose):
    part1 = (r**2)/sqrt(psigood(r)-E)
    part2 = ggood(psigood(r)*(1-theta) + E*theta)
    part3 = (1./sqrt(theta))-sqrt(theta)
    return part1*part2*part3

def funcG(E,verbose = False):
    print E
    print E.shape
    try:
        t = E.shape
        Gans = []
        problems = []
        for i in range(len(E)):
            print i+1, 'of', len(E)
            rapoval = rapo(E[i])
            temp = intg.dblquad(Ginterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E[i],verbose))
            try:
                if temp[3] != '' and verbose == True:
                    print 'G, E = ', E[i], 'message = ', temp[3]
                    problems.append(i)
            except IndexError:
                pass
            Gans.append(temp)
        return array(Gans),problems
    except AttributeError:
        rapoval = rapo(E)
        temp = intg.dblquad(Ginterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E,verbose))
        problem = []
        try:
            if temp[3] != '' and verbose == True:
                print 'G, E = ', E, 'message = ', temp[3]
                problem = [E]
        except IndexError:
            pass
        return temp,problem

########******************* COMPUTE G *******************######## 
if run[G] == "ON":
    if run[psi] == "ON" and run[g] == "ON":
        tic = time.clock()
        Ggood = makegood(funcG,rtest,[3,-4,0.1],Egrid,model.b-4,model.g-4,verbose=False,plotting = ['E','G'])
        toc = time.clock()
        print 'Ggood ran in ', toc-tic, 's'
    elif run[psi] != "ON":
        print 'To compute G, please turn psi ON'
    elif run[g] != "ON":
        print 'To compute G, please turn g ON'

########******************* DISTRIBUTION FUNCTION *******************######## 

def finterior1(r,E,rapoval,verbose):
    var = rapoval/r
    psi = (10**psigood(log10(var)))[0]
    result = (var**3)*(1./sqrt(abs(E-psi)))*model.d2rhodr2(var)
    return result

def finterior2(r,E,rapoval,verbose):
    var = rapoval/r
    psi = (10**psigood(log10(var)))[0]
    return (var**2)*(1./sqrt(abs(E-psi)))*model.drhodr(var)

def finterior3(r,E,rapoval,verbose):
    var = rapoval/r
    psi = (10**psigood(log10(var)))[0]
    Mencvar = (10**Mencgood(log10(var)))[0]
    Mencrap = (10**Mencgood(log10(rapoval)))[0]
    return -(var**2)*(1./sqrt(abs(E-psi)))*(1./(2*rapoval))*model.drhodr(var)*((model.Mnorm*(r-1) + r*Mencvar - Mencrap)/abs(E-psi))
                   
def funcf(E,verbose=False):
    print 'starting f evaluation'
    epsilon = 0
    try:
        t = E.shape
        fans = []
        problems = []
        for i in range(len(E)):
            rapoval = rapo(E[i])
            prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + Menc(rapoval,verbose)[0])))
            temp1 = intg.quad(finterior1,epsilon,1-epsilon,args=(E[i],rapoval,verbose),full_output = 1)
            temp2 = intg.quad(finterior2,epsilon,1-epsilon,args=(E[i],rapoval,verbose),full_output = 1)
            temp3 = intg.quad(finterior3,epsilon,1-epsilon,args=(E[i],rapoval,verbose),full_output = 1)
            t = temp1[0] + temp2[0] + temp3[0]
            try:
                if verbose==True:
                    if temp1[3] != '':
                        print 'f, E = ',E[i],'message = ',temp1[3]
                        problems.append(i)
                    elif temp2[3] != '':
                        print 'f, E = ',E[i],'message = ',temp2[3]
                        problems.append(i)
                    elif temp3[3] != '':
                        print 'f, E = ',E[i],'message = ',temp3[3]
                        problems.append(i)                    
            except IndexError:
                pass
            fans.append((prefactor*t)[0])
        return array(fans),problems
    except AttributeError:
        rapoval = rapo(E)
        prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + Menc(rapoval,verbose)[0])))
        temp1 = intg.quad(finterior1,0,1,args=(E,rapoval,verbose),full_output = 1)
        temp2 = intg.quad(finterior2,0,1,args=(E,rapoval,verbose),full_output = 1)
        temp3 = intg.quad(finterior3,0,1,args=(E,rapoval,verbose),full_output = 1)
        t = temp1[0] + temp2[0] + temp3[0]
        problem = []
        try:
            if verbose==True:
                if temp1[3] != '':
                    print 'f, E = ',E,'message = ',temp1[3]
                    problem = [E]
                elif temp2[3] != '':
                    print 'f, E = ',E,'message = ',temp2[3]
                    problem = [E]
                elif temp3[3] != '':
                    print 'f, E = ',E,'message = ',temp3[3]
                    problem = [E]
                t = -1
        except IndexError:
            pass
        return prefactor*t, problem

########******************* COMPUTE f *******************######## 
if run[f] == "ON":
    if run[Menc] == "ON" and run[psi] == "ON":
        tic = time.clock()
        fgood = makegood(funcf,rtest,[5,-3,0.03],Egrid,model.b-1.5,model.g-1.5,verbose=False,plotting = ['E','f'])
        toc = time.clock()
        print 'fgood ran in ', toc-tic, 's'
    elif run[Menc] != "ON":
        print 'To compute f, please turn Menc ON'
    elif run[psi] != "ON":
        print 'To compute f, please turn psi ON'



