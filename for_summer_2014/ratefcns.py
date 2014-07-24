import numpy as np
nsum = np.sum
from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time
import datetime
import pickle
from besselgen import BesselGen

alpha = 7.52
beta = 3.13
gamma = 1.98
r0pc = 1.
rb = 10**2.38
r0pc = rb
mub = 19.98
M2L = 6.27
MsunV = 4.83
rho0 = (1./rb)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 
MBH_Msun = 10**6.04
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
generate = False
seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
verbosity = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
plot = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}

rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

#choose model from rho or genrho options
#for example:
from models import NukerModel
model = NukerModel('NGC4467',alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate)

if model.generate == True:
    call(["mkdir","{0}".format(directory)])
    seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}

########******************* ENCLOSED MASS *******************########

def Minterior(r):
    """
    interior of the Menc integral
    """
    return model.rho(r)*r**2

def funcMenc(r,verbose=False):
    """
    functional form of Menc
    relies on Minterior
    returns Menc(r)
    """
    try:
        problems = []
        t = r.shape
        Mencs = []
        for i in range(len(r)):
            temp = intg.quad(Minterior,0,r[i],full_output=1,epsabs = tol,epsrel = tol)
            try:
                if temp[3]!='':
                    problems.append(i)
                    if verbose==True:
                        print 'Menc, r = ',r[i],'message = ',temp[3],'\n'
            except (IndexError, TypeError):
                pass
            Mencs.append(4*pi*temp[0])
        return array(Mencs),array(problems)
    except (AttributeError,TypeError):
        problem = []
        temp = 4*pi*intg.quad(Minterior,0,r)[0]
        try:
            if temp[3]!='':
                problems = [r]
                if verbose==True:
                    print 'Menc, r = ',r,'message = ',temp[3],'\n'
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t,problem

########******************* RADIUS OF INFLUENCE *******************######## 

def rHimplicit(r):
    """
    equation that has its minimum when r = rH
    """
    return abs(model.Mnorm-funcMenc(abs(r))[0])

def rH(verbose=True):
    """
    finds root of rHimplicit
    """
    rresult=root(rHimplicit,1e-4)
    if rresult.success == True:
        return abs(rresult.x)
    elif rresult.success == False:
        if verbose == True:
            print 'Failed to evaluate rH'
            print rresult.message
        return abs(rresult.x)

rH = rH()

########******************* POTENTIAL *******************######## 
        
def psi2interior(r):
    """
    interior of psi part 2 integral
    """
    return model.rho(r)*r

tol = 1e-3
def psi2(r,verbose=False):
    """
    functional form of psi part 2
    relies on psi2interior
    returns psi2(r)
    """
    try:
        problems = []
        t = r.shape
        psi2s = []
        for i in range(len(r)):
            temp=intg.quad(psi2interior,r[i],inf,full_output=1,epsabs = tol,epsrel = tol)
            try:
                if temp[3]!='':
                    if verbose==True:
                        print 'psi2, index =',i,'r = ',r[i],'message = ',temp[3],'\n'
                    problems.append(i)
            except (IndexError,TypeError):
                pass
            psi2s.append(4*pi*temp[0])
        return array(psi2s),problems
    except AttributeError: 
        problem = []
        temp = 4*pi*intg.quad(psi2interior,r,inf)[0]
        try:
            if temp[3]!='':
                if verbose==True:
                    print 'psi2, r = ',r,'message = ',temp[3],'\n'
                problem = [r]
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t, problem

def funcpsi(r,verbose=False):
    """
    returns potential as a function of r
    """
    part1 = (model.Mnorm/r)
    part2 = funcMenc(r,verbose)
    part3 =  psi2(r,verbose)
    problems = array([])
    if part2[1] != []:
        problems = concatenate((problems,array(part2[1])))
    if part3[1] != []:
        problems = concatenate((problems,array(part3[1])))
    return part1 + (part2[0]/r) + part3[0],problems

########******************* APOCENTER RADIUS *******************######## 

def rapoimplicit(r,E):
    """
    function with a minimum at r=rapo
    """
    return abs(10**psigood(log10(abs(r)))-E)

def rapo(E):
    """
    finds root of rapoimplicit
    """
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

########******************* CIRCULAR ANGULAR MOMENTUM *******************######## 

def Jc2implicit(r,E,verbose):
    """
    """
    return abs(10**psigood(log10(abs(r)))-E-((10**Mencgood(log10(abs(r)))+model.Mnorm)/(2*r)))

def funcJc2(E,verbose):
    """
    see Jc2implicit
    """
    try:
        t = E.shape
        Jcs = []
        problems = []
        for i in range(len(E)):
            #print i+1, ' of ', len(E)
            if E[i]**-1 > 0.4:
                rguess = 10*E[i]**-1
            elif E[i]**-1 < 0.4:
                rguess = 0.01*E[i]**-1
            rresult = root(Jc2implicit,rguess,args=(E[i],verbose))
            rresult.x = abs(rresult.x)
            if rresult.success == True:
                Jcs.append(((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0])
            elif rresult.success == False:
                if verbose==True:
                    print 'Failed to evaluate Jc2'
                    print rresult.message
                Jcs.append(((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0])
                problems.append(i)
        return array(Jcs),problems
    except AttributeError:
        rresult = root(Jc2implicit,0.01*E[i]**-1,args=(E[i],verbose))
        problem = []
        rresult.x = abs(rresult.x)
        if rresult.success == True:
            Jc = ((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0]
        elif rresult.success == False:
            if verbose==True:
                print 'Failed to evaluate Jc2'
                print rresult.message
            Jc = ((10**Mencgood(log10(rresult.x))+model.Mnorm)*rresult.x)[0]
            problem = [E]
        return Jc, problem

########******************* g *******************######## 

def lginterior(r,E):
    """
    interior of g integral
    """
    return (model.drhodr(1./r))*(r**-2)*((sqrt(abs(E-10**psigood(log10(1./r)))))**-1)
    
def funclg(E,verbose=False): #removed negative from final answer while incorporating alternate Nuker model
    """
    functional form of g
    relies on ginterior
    returns g(E)
    """
    try:
        t = E.shape
        gans = []
        problems = []
        for i in range(len(E)):
            rapoval = rapo(E[i])
            temp = intg.quad(lginterior,0,1./rapoval,args = E[i],full_output = 1)
            t = temp[0]
            try:
                if temp[3] != '':
                    if verbose == True:
                        print 'g, E = ',E[i], 'message = ', temp[3]
                    problems.append(i)
            except (IndexError,TypeError):
                pass
            gans.append(-pi*t)
        return array(gans),problems
    except (AttributeError,TypeError) as e:
        problem = []
        rapoval = rapo(E)
        temp = intg.quad(lginterior,0,1./rapoval,args = E,full_output = 1)
        t = temp[0]
        try:
            if temp[3] != '':
                if verbose == True:
                    print 'g, E = ',E, 'message = ', temp[3]
                problem = [E]
        except (IndexError,TypeError):
            pass
        return -pi*t, problem

########******************* mathcalG *******************########

psibG_memo = {}
part2bG_memo = {}
part3bG_memo = {}

def bGinterior(theta,r,E):
    """
    interior of G integral
    """
    if not r in psibG_memo:
        psibG_memo[r] = 10**psigood(log10(r))
    psir = psibG_memo[r]
    part1 = (r**2)/sqrt(psir-E)
    if not log10(psir*(1-theta) + E*theta) in part2bG_memo:
        part2bG_memo[log10(psir*(1-theta) + E*theta)] = 10**ggood(log10(psir*(1-theta) + E*theta))
    part2 = part2bG_memo[log10(psir*(1-theta) + E*theta)]
    if not theta in part3bG_memo:
        part3bG_memo[theta]= (1./sqrt(theta))-sqrt(theta)
    part3 = part3bG_memo[theta]
    return part1*part2*part3

def funcbG(E,verbose = False):
    """
    functional form of mathcalG
    relies on Ginterior
    returns mathcalG(E)
    """
    tolerance = 1.49e-8
    try:
        t = E.shape
        Gans = []
        problems = []
        for i in range(len(E)):
            print i+1, 'of', len(E)
            rapoval = rapo(E[i])
            try:
                temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 1e-4, lambda r: 1,args = (E[i],),epsabs = tolerance,epsrel = tolerance)
            except UserWarning as e:
                if verbose == True:
                    print 'G, E = ', E[i], 'message = ', e
                problems.append(i)
            Gans.append(temp[0])
        return array(Gans),problems
    except AttributeError:
        rapoval = rapo(E)
        problem = []
        try:
            temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E,verbose))
        except UserWarning as e:
            if verbose == True:
                print 'G, E = ', E, 'message = ', temp[3]
            problem = [E]
        return temp[0],problem

########******************* DISTRIBUTION FUNCTION *******************######## 

def finterior(r,E,rapoval):
    var = rapoval/r
    psi = (10**psigood(log10(var)))[0]
    Mencvar = (10**Mencgood(log10(var)))[0]
    Mencrap = (10**Mencgood(log10(rapoval)))[0]
    result1 = (var**3)*(1./sqrt(abs(E-psi)))*model.d2rhodr2(var) 
    result2 = (var**2)*(1./sqrt(abs(E-psi)))*model.drhodr(var) 
    result3 = -(var**2)*(1./sqrt(abs(E-psi)))*(1./(2*rapoval))*model.drhodr(var)*((model.Mnorm*(r-1) + r*Mencvar - Mencrap)/abs(E-psi))  
    return result1+result2+result3

def funcf(E,verbose=False):
    """
    functional form of f
    relies on finterior1,finterior2,finterior3
    returns f(E)
    """
    epsilon = 0
    try:
        t = E.shape
        fans = []
        problems = []
        for i in range(len(E)):
            print i+1, ' of ', len(E)
            rapoval = rapo(E[i])
            prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + 10**Mencgood(log10(rapoval)))))
            temp = intg.quad(finterior,epsilon,1-epsilon,args=(E[i],rapoval),full_output = 1)
            t = temp[0]
            try:
                if temp[3] != '':
                    if verbose == True:
                        print 'f, E = ',E[i],'message = ',temp[3]
                    problems.append(i)                  
            except IndexError:
                pass
            fans.append((prefactor*t)[0])
        return array(fans),problems
    except AttributeError:
        rapoval = rapo(E)
        prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + 10**Mencgood(log10(rapoval)))))
        temp = intg.quad(finterior,epsilon,1-epsilon,args=(E,rapoval),full_output = 1)
        t = temp[0]
        problem = []
        try:
            if temp1[3] != '':
                if verbose == True:
                    print 'f, E = ',E,'message = ',temp1[3]
                problem = [E]
        except IndexError:
            pass
        return prefactor*t, problem

########******************* ADDITIONAL FUNCTIONS *******************######## 

def funcq(r):
    return (4./pi)*log(model.Lam)*(model.r0_rT/model.MBH)*10**Ggood(log10(r))

def Rlc(r):
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return -log(interior)
