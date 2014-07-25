import numpy as np
nsum = np.sum
from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
import matplotlib.pyplot as plt
from besselgen import BesselGen
from construction import *

def findrho0(rb,M2L,mub):
    MsunV = 4.83
    return (1./rb)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 

########******************* ENCLOSED MASS *******************########

def Minterior(r,prereqs):
    """
    interior of the Menc integral
    """
    model = prereqs
    return model.rho(r)*r**2

def funcMenc(r,verbose,prereqs):
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
            temp = intg.quad(Minterior,0,r[i],args = prereqs,full_output=1,epsabs = tol,epsrel = tol)
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
                    print 'Menc, r = ',r,args = prereqs,'message = ',temp[3],'\n'
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t,problem

########******************* RADIUS OF INFLUENCE *******************######## 

def rHimplicit(r,prereqs):
    """
    equation that has its minimum when r = rH
    """
    model = prereqs
    return abs(model.Mnorm-funcMenc(abs(r))[0])

def rH(prereqs,verbose=True):
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

########******************* GRID CREATION  *******************######## 

def rgrid(prereqs,upstep=5,downstep=-5,step=0.03):
    """
    constructs a grid in radius and adds one point at each extreme (up and down)
    returns 10**grid
    """
    model = prereqs
    rmin = min([rH(model),[1.]])
    rmax = max([rH(model),[1.]])
    rimin = log10(rmin) + downstep
    rimax = log10(rmax) + upstep
    dri = step
    rarray = arange(rimin,rimax,dri)
    rarray = 10**rarray
    rchange = rarray[len(rarray)-1]
    rstart = rarray[0]
    return rarray,rchange,rstart

rarray,rchange,rstart = rgrid(5,-5,0.03)

def Egrid(prereqs,upstep=5,downstep=-3,step=0.1):
    """
    constructs a grid in energy and adds one point at each extreme (up and down)
    returns 10**grid
    """
    model = prereqs
    rmin = min([rH(model),[1.]])[0]
    rmax = max([rH(model),[1.]])[0]
    eimin = log10(funcMenc(rmax)[0]/rmax) + downstep
    eimax = log10(model.Mnorm/rmin) + upstep
    dei = step
    Earray = arange(eimin,eimax,dei)
    Earray = 10**Earray
    Echange = Earray[len(Earray)-1]
    Estart = Earray[0]
    return Earray,Echange,Estart

########******************* POTENTIAL *******************######## 
        
def psi2interior(r,prereqs):
    """
    interior of psi part 2 integral
    """
    model = prereqs
    return model.rho(r)*r

tol = 1e-3
def psi2(r,verbose,prereqs):
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
            temp=intg.quad(psi2interior,r[i],inf,args = prereqs,full_output=1,epsabs = tol,epsrel = tol)
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
        temp = 4*pi*intg.quad(psi2interior,r,inf,args = prereqs)[0]
        try:
            if temp[3]!='':
                if verbose==True:
                    print 'psi2, r = ',r,'message = ',temp[3],'\n'
                problem = [r]
                t = temp[0]
        except (IndexError,TypeError):
            t = temp
        return t, problem

def funcpsi(r,verbose,prereqs):
    """
    returns potential as a function of r
    """
    model = prereqs
    part1 = (model.Mnorm/r)
    part2 = funcMenc(r,verbose,prereqs)
    part3 =  psi2(r,verbose,prereqs)
    problems = array([])
    if part2[1] != []:
        problems = concatenate((problems,array(part2[1])))
    if part3[1] != []:
        problems = concatenate((problems,array(part3[1])))
    return part1 + (part2[0]/r) + part3[0],problems

########******************* APOCENTER RADIUS *******************######## 

def rapoimplicit(r,E,prereqs):
    """
    function with a minimum at r=rapo
    """
    psigood = prereqs
    return abs(10**psigood(log10(abs(r)))-E)

def rapo(E,prereqs):
    """
    finds root of rapoimplicit
    """
    if E**-1 > 0.2:
        rguess = 10*E**-1
    elif E**-1 < 0.2:
        rguess = 0.01*E**-1
    rresult = root(rapoimplicit,rguess,args=(E,prereqs))
    if rresult.success == True:
        return abs(rresult.x)
    elif rresult.success == False:
        print 'Failed to evaluate rapo'
        print rresult.message
        return abs(rresult.x)

########******************* CIRCULAR ANGULAR MOMENTUM *******************######## 

def Jc2implicit(r,E,verbose,prereqs):
    """
    """
    model,Mencgood,psigood = prereqs
    return abs(10**psigood(log10(abs(r)))-E-((10**Mencgood(log10(abs(r)))+model.Mnorm)/(2*r)))

def funcJc2(E,verbose,prereqs):
    """
    see Jc2implicit
    """
    model,Mencgood,psigood = prereqs
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
            rresult = root(Jc2implicit,rguess,args=(E[i],verbose,prereqs))
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
        rresult = root(Jc2implicit,0.01*E**-1,args=(E,verbose,prereqs))
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

def lginterior(r,E,prereqs):
    """
    interior of g integral
    """
    model,psigood = prereqs
    return (model.drhodr(1./r))*(r**-2)*((sqrt(abs(E-10**psigood(log10(1./r)))))**-1)
    
def funclg(E,verbose,prereqs): #removed negative from final answer while incorporating alternate Nuker model
    """
    functional form of g
    relies on ginterior
    returns g(E)
    """
    model,psigood = prereqs
    try:
        t = E.shape
        gans = []
        problems = []
        for i in range(len(E)):
            rapoval = rapo(E[i],psigood)
            temp = intg.quad(lginterior,0,1./rapoval,args = (E[i],prereqs),full_output = 1)
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
        rapoval = rapo(E,psigood)
        temp = intg.quad(lginterior,0,1./rapoval,args = (E,prereqs),full_output = 1)
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

def bGinterior(theta,r,E,prereqs):
    """
    interior of G integral
    """
    model,psigood,ggood = prereqs
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

def funcbG(E,verbose,prereqs):
    """
    functional form of mathcalG
    relies on Ginterior
    returns mathcalG(E)
    """
    model,psigood,ggood = prereqs
    tolerance = 1.49e-8
    try:
        t = E.shape
        Gans = []
        problems = []
        for i in range(len(E)):
            print i+1, 'of', len(E)
            rapoval = rapo(E[i],psigood)
            try:
                temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 1e-4, lambda r: 1,args = (E[i],prereqs),epsabs = tolerance,epsrel = tolerance)
            except UserWarning as e:
                if verbose == True:
                    print 'G, E = ', E[i], 'message = ', e
                problems.append(i)
            Gans.append(temp[0])
        return array(Gans),problems
    except AttributeError:
        rapoval = rapo(E,psigood)
        problem = []
        try:
            temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E,prereqs))
        except UserWarning as e:
            if verbose == True:
                print 'G, E = ', E, 'message = ', temp[3]
            problem = [E]
        return temp[0],problem

########******************* DISTRIBUTION FUNCTION *******************######## 

def finterior(r,E,rapoval,prereqs):
    model,Mencgood,psigood = prereqs
    var = rapoval/r
    psi = (10**psigood(log10(var)))[0]
    Mencvar = (10**Mencgood(log10(var)))[0]
    Mencrap = (10**Mencgood(log10(rapoval)))[0]
    result1 = (var**3)*(1./sqrt(abs(E-psi)))*model.d2rhodr2(var) 
    result2 = (var**2)*(1./sqrt(abs(E-psi)))*model.drhodr(var) 
    result3 = -(var**2)*(1./sqrt(abs(E-psi)))*(1./(2*rapoval))*model.drhodr(var)*((model.Mnorm*(r-1) + r*Mencvar - Mencrap)/abs(E-psi))  
    return result1+result2+result3

def funcf(E,verbose,prereqs):
    """
    functional form of f
    relies on finterior1,finterior2,finterior3
    returns f(E)
    """
    model,Mencgood,psigood = prereqs
    epsilon = 0
    try:
        t = E.shape
        fans = []
        problems = []
        for i in range(len(E)):
            print i+1, ' of ', len(E)
            rapoval = rapo(E[i],psigood)
            prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + 10**Mencgood(log10(rapoval)))))
            temp = intg.quad(finterior,epsilon,1-epsilon,args=(E[i],rapoval,prereqs),full_output = 1)
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
        rapoval = rapo(E,psigood)
        prefactor = (1./(sqrt(8)*pi**2*(model.Mnorm + 10**Mencgood(log10(rapoval)))))
        temp = intg.quad(finterior,epsilon,1-epsilon,args=(E,rapoval,prereqs),full_output = 1)
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

def funcq(r,Ggood,model):
    return (4./pi)*log(model.Lam)*(model.r0_rT/model.MBH)*10**Ggood(log10(r))

def Rlc(r,Jc2good,model):
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return -log(interior)

########******************* IMPORT DATA TABLES *******************########

bessel = BesselGen(['alpham_table.pkl','xi_table.txt','Bessel_table.txt','mpiece_table.txt'])

########******************* RATE *******************########

def dgdlnrpinterior(E,u,qmin):
    """
    interior of the integral used to calculate rate as function of pericentre
    """
    sumlim = max([200,2*qmin**-0.5])
    ms = arange(1,sumlim,1.)
    alphas = bessel.alpham(ms)
    qval = funcq(E)
    fval = 10**fgood(log10(E))
    xival = bessel.xi(qval)
    bfin = bessel.besselfin(ms,u)
    mpiece = bessel.mpiece(ms)
    part1 = array(fval/(1+(qval**-1)*(xival)*Rlc(E)))
    part2list = exp(-array(matrix(alphas**2).T*matrix(qval/4)))
    part2list = array([(bfin/mpiece)[i]*part2list[i] for i in range(len(alphas))])
    part2 = 1-2*nsum(part2list,axis = 0)
    return part1*part2

def dgdlnrp(u,Emin = 0.01,Emax=100,verbose = False):
    """
    rp - pericentre radius
    Emin, Emax - bounds of the integral
    verbose = True - print error messages from integration
    returns the rate for given rp
    """
    #u = sqrt(rp/model.rT)
    prefactor = (8*pi**2)*model.MBH*(model.r0_rT**-1)*((model.tdyn0/(3600*24*365))**-1)
    #print 'prefactor = ',prefactor
    qmin = funcq(Emax)
    qmax = funcq(Emin)
    try:
    	result_list = []
    	for i in range(len(u)):
    		print i+1, ' of ', len(u)
    		result = intg.romberg(dgdlnrpinterior,Emin,Emax,args = (u[i],qmin),divmax = 20)#,full_output = 1)
    		t = result[0]
    		result_list.append(prefactor*(u[i]**2)*t)
    		try:
        		if result[3] != '':
        			if verbose == True:
        				print 'dgdlnrp, rp = ',rp[i], 'message = ',result[3]
        	except (IndexError,TypeError):
        		pass
        return array(result_list)
    except (AttributeError,TypeError) as e:
    	result = intg.romberg(dgdlnrpinterior,Emin,Emax,args = (u,qmin),divmax = 20)#,full_output = 1)
    	t = result[0]
    	try:
    		if result[3] != '':
    			if verbose == True:
    				print 'dgdlnrp, rp = ',rp, 'message = ',result[3]
    	except (IndexError,TypeError):
    		pass
    	return prefactor*(u**2)*t #units yr^-1

def Ndotinterior(r):
    u = sqrt(r/model.rT)
    return r*dgdlnrp(u)

def Ndot():
    return intg.quad(Ndotinterior,0,model.rT,full_output=1)[0]


def ginterior(E):
    qval = funcq(E)
    return (10**fgood(log10(E))*qval)/((qval/bessel.xi(qval)) + Rlc(E))

def gdirect(Emin = 0.01,Emax = 100,verbose = False):
    prefactor = (8*pi**2)*model.MBH*(model.r0_rT**-1)*(model.tdyn0**-1)
    qmin = funcq(Emax)
    qmax = funcq(Emin)
    result = intg.quad(ginterior,Emin,Emax,full_output = 1)
    t = result[0]
    try:
        if result[3] != '':
            if verbose == True:
                print 'direct rate message = ',result[3]
    except (IndexError,TypeError):
        pass
    return prefactor*t
