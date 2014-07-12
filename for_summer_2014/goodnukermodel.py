import numpy as np
nsum = np.sum
from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time
import datetime
from subprocess import call
import pickle
from besselgen import BesselGen
from scipy.special import hyp2f1
import matplotlib
plt.ion()
try:
    call(['rm -f', 'Imodel_test.pyc'])
except OSError:
    pass
alpha = 7.52 #1.0
beta = 3.13-1#4.0
gamma = 1.98-1#1.5
r0pc = 1.
rb = 10**2.38
r0pc = rb
mub = 19.98
M2L = 6.27
MsunV = 4.83
rho0 = 1e5
rho0 = (1./rb)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 
MBH_Msun = 10**6.04#1e3
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
generate =False
seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
verbosity = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
plot = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}
########******************* MODEL FRAMEWORK *******************########

from Imodel_test import NukerModel
                                
########******************* CONSTRUCT MODEL **3*****************########

model = NukerModel('ImodelNGC4467',alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate)

model.getrho()

rtest = arange(-7,6,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

directory = "{0}_a{1}_b{2}_g{3}_r{4}_rho{5}_MBH{6}".format(model.name,model.a,model.b,model.g,model.r0,model.rho0,model.MBH)
if model.generate == True:
    call(["mkdir","{0}".format(directory)])
    seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
########******************* CONSTRUCTION FUNCTIONS *******************########
def piecewise2(r,inter,start,end,lim1,lim2,smallrexp,largerexp,conds=False):
    """
    r - independent variable
    inter - interpolated object
    start - first computed value of function
    end - last computed value of function
    lim1 - first piecewise break
    lim2 - second piecewise break
    smallrexp - log slope at small r or large E
    largerexp - log slope at large r or small E
    conds - any extra conditions

    returns a value for the function at r
    """
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

def plotter(name,r,inter,rstart,rchange,start,end,smallrexp,largerexp,conds,labels):
    """
    name - name under which figure will be saved
    r - independent variable array
    inter - interpolated object
    rstart - first element of r
    rchange - last element of r
    start - first computed value of function
    end - last computed value of function
    smallrexp - log slope at small r or large E
    largerexp - log slope at large r or small E
    conds - any extra conditions
    labels - plot axis labels, x and then y

    saves a plot of the piecewise function
    """
    m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp,conds)
    plt.figure()
    plt.loglog(r[1:-1],m[1:-1],'c',linewidth = 5)
    plt.ylabel(r'{0}'.format(labels[1]))
    plt.xlabel('{0}'.format(labels[0]))
    plt.xlim(min(r[1:-1]),max(r[1:-1]))
    plt.ylim(min(m[1:-1]),max(m[1:-1]))
    #if rstart>r[0] and rchange<r[len(r)-1]:
        #plt.axvline(rstart, color='r',label='Limits of interpolation')
        #plt.axvline(rchange, color='r')
        #plt.legend(loc='best')
    plt.savefig('{0}/{1}.png'.format(directory,name))
    #plt.show()

def makegood(func,r,size,grid,smallrexp,largerexp,verbose = False,conds = False,plotting=False,problem = True):
    """

    func - function to be evaluated
    r - independent variable array
    size - size of generated independent variable array with format 
    	   [upstep,downstep,max,min,stepsize]
    grid - choice of grid generator function
    smallrexp - log slope at small r or large E
    largerexp - log slope at large r or small E
    verbose = False - suppresses warnings and error messages from 
    	              integration and rootfinders
    conds = False - no special conditions on the piecewise function
    plotting = False - 	do not save plots
    problem = True - eliminate problem points
    
    returns an interpolated object version of the function based 
    computed values
    
    """
    rarray,rchange,rstart = grid(size[0],size[1],size[2])
    tab,problems = func(rarray,verbose)
    print 'fraction reporting a message: {0}'.format(float(len(problems))/float(len(tab)))
    if problem == True:
        tab = [i for j, i in enumerate(tab) if j not in problems]
        rarray = [i for j, i in enumerate(rarray) if j not in problems]
    inter = interp1d(log10(rarray),log10(tab))
    start = tab[0]
    end = tab[len(rarray)-1]
    m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp,conds)
    inter2 = interp1d(log10(r),log10(m))
    pklrfile = open('{0}/r{1}.pkl'.format(directory,str(func)[10:15]),"wb")
    pickle.dump(r,pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,str(func)[10:15]),"wb")
    pickle.dump(m,pklffile)
    pklffile.close()
    if plotting != False:
        plotter(str(func)[10:15],r,inter,rstart,rchange,start,end,smallrexp,largerexp,conds,plotting)
    return inter2

def compute(dependencies,name,function,rtest,size,grid,exps,kwargs):
    """
    dependencies - other functions needed to compute this one, 
                   format [func1, "func1",func2,"func2",...]
    name - name of function in the dictionaries, 
           format ["name",name]
    function - name of the functional form
    rtest - independent variable array
    size - size of generated independent variable array 
    	   with format [upstep,downstep,max,min,stepsize]
    grid - choice of grid generator function
    exps - extreme r or E behaviour, 
           format [smallrexp,largerexp]
    kwargs - additional information used to specify 
            conditions and plotting information, 
            format [conds,plotting,problem]
    
    finds interpolated form based on conditions in the 
    dictionaries and pickles it or unpickles intepolated form
    returns interpolated form
    """
    strname,name = name
    if seton[name] == "ON":
        try:   
            i = 0
            while i<len(dependencies):
                dependencies[i](1)
                i+=2
            smallrexp,largerexp = exps
            conditions,plotdat,prob = kwargs
            if verbosity[name] == "ON" and plot[name] == "ON":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,verbose = True,conds = conditions,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "ON" and plot[name] == "OFF":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,verbose = True,conds = conditions,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "OFF" and plot[name] == "ON":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,conds = conditions,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "OFF" and plot[name] == "OFF":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,conds = conditions,problem = prob)
                toc = time.clock()
                delt = toc-tic
            print '{0}good ran in \t {1}'.format(strname,str(datetime.timedelta(seconds=delt)))
            return good
        except TypeError as e:
            print 'e = ',e
            print 'To compute {0}, please turn {1} ON'.format(strname,dependencies[i+1])
    elif seton[name] != "ON":
        try:
            tic = time.clock()
            pklrfile = open('{0}/r{1}.pkl'.format(directory,str(function)[10:15]),"rb")
            rarray = pickle.load(pklrfile)
            pklrfile.close()
            pklffile = open('{0}/{1}.pkl'.format(directory,str(function)[10:15]),"rb")
            tab = pickle.load(pklffile)
            pklffile.close()
            good =  interp1d(log10(rarray),log10(tab))
            toc = time.clock()
            delt = toc-tic
            print '{0}good loaded in \t {1}'.format(strname,str(datetime.timedelta(seconds=delt)))
            return good
        except IOError:
            print '{0} not yet generated, please turn it ON'.format(strname)
            
########******************* ENCLOSED MASS *******************########

def Minterior(r):
    """
    interior of the Menc integral
    """
    #print r
    return model.rho(r)*r**2

tol = 1.e-3
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

########******************* CONSTRUCTION FUNCTIONS *******************########

def rgrid(upstep=5,downstep=-5,step=0.03):
    """
    constructs a grid in radius and adds one point at each extreme (up and down)
    returns 10**grid
    """
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

rarray,rchange,rstart = rgrid(5,-5,0.03)

def Egrid(upstep=5,downstep=-3,step=0.1):
    """
    constructs a grid in energy and adds one point at each extreme (up and down)
    returns 10**grid
    """
    rmin = min([rH(),[1.]])[0]
    rmax = max([rH(),[1.]])[0]
    eimin = log10(funcMenc(rmax)[0]/rmax) + downstep
    eimax = log10(model.Mnorm/rmin) + upstep
    dei = step
    Earray = arange(eimin,eimax,dei)
    Earray = 10**Earray
    Echange = Earray[len(Earray)-1]
    Estart = Earray[0]
    return Earray,Echange,Estart

########******************* COMPUTE MENC *******************######## 

Mencgood = compute([],["Menc",Menc],funcMenc,rtest,[4,-6,0.03],rgrid,[3-model.g,0],[[2,0,3-model.b,4*pi*model.rho(rchange)*(rchange**3)],['r','M'],False])

def Menc2(r):
    const = -4*pi*((-3+model.g)**-1)
    rfact = r**(3-model.g)
    hyp = hyp2f1(((3-model.g)/model.a),(-(-model.b+model.g)/model.a),(1+((3-model.g)/model.a)),-r**model.a)
    return const*rfact*hyp
'''
#print 'Mdiff = {0}'.format(10**Mencgood(log10(rtest))/Menc2(rtest))
plt.figure()
plt.loglog(rtest,Menc2(rtest),'g',label = 'Analytic')
plt.loglog(rtest,10**Mencgood(log10(rtest)),'r',label = 'Integral')
plt.xlabel('r')
plt.ylabel('Menc')
plt.legend(loc = 'best')
plt.title(r'{0}, $\alpha$ = {1}, $\beta$ = {2}, $\gamma$ = {3}'.format(model.name,alpha,beta,gamma))
plt.show()
plt.savefig('{0}/Mcompare.png'.format(directory))
'''
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

########******************* COMPUTE PSI *******************######## 

psigood = compute([],["psi",psi],funcpsi,rtest,[4.3,-6,0.03],rgrid,[-1,-1],[False,['r','$\psi$'],False])

def funcpsi2(r):
    const = ((-2+model.b)**-1)
    rfact = r**(2-model.b)
    hyp = hyp2f1(((-2+model.b)/model.a),((model.b-model.g)/model.a),((-2+model.a+model.b)/model.a),-r**-model.a)
    part1 = 4*pi*const*rfact*hyp
    part2 = Menc2(r)/r
    part3 = model.Mnorm/r
    return part1 + part2 + part3

'''
#print 'psidiff = {0}'.format(10**psigood(log10(rtest))/funcpsi2(rtest))
plt.figure()
plt.loglog(rtest,funcpsi2(rtest),'g',label = 'Analytic')
plt.loglog(rtest,10**psigood(log10(rtest)),'r',label = 'Integral')
plt.xlabel('r')
plt.ylabel('psi')
plt.legend(loc = 'best')
plt.title(r'{0}, $\alpha$ = {1}, $\beta$ = {2}, $\gamma$ = {3}'.format(model.name,alpha,beta,gamma))
plt.show()
plt.savefig('{0}/psicompare.png'.format(directory))
'''
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
    not actually sure what this solves yet ********************************************************************
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

########******************* COMPUTE Jc2 *******************######## 

prereqs = [Mencgood,"Menc",psigood,"psi"]

Jc2good = compute(prereqs,["Jc2",Jc2],funcJc2,rtest,[3,-4,0.01],Egrid,[-1,-1],[False,['E','Jc2'],False])

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

########******************* COMPUTE g *******************######## 
prereqs = [psigood,"psi"]

ggood = compute(prereqs,["g",g],funclg,rtest,[3,-3,0.1],Egrid,[model.b-0.5,model.g-0.5],[False,['E','g'],False])

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

########******************* COMPUTE G *******************######## 
prereqs = [psigood, "psi",ggood,"g"]

Gtest = arange(-5,5,0.01)
Gtest = append(Gtest,40)
Gtest = insert(Gtest,0,-40)
Gtest = 10**Gtest

Ggood = compute(prereqs,["G",G],funcbG,Gtest,[3,-3,0.1],Egrid,[model.b-4,model.g-4],[False,['E','G'],False])

psibG_memo = {}
part2bG_memo = {}
part3bG_memo = {}

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

########******************* COMPUTE f *******************######## 
prereqs = [Mencgood,"Menc",psigood,"psi"]

ftest = arange(-3,5,0.01)
ftest = append(ftest,40)
ftest = insert(ftest,0,-40)
ftest = 10**ftest

fgood = compute(prereqs,["f",f],funcf,ftest,[5,-3,0.03],Egrid,[model.b-1.5,model.g-1.5],[False,['E','f'],False])

########******************* ADDITIONAL FUNCTIONS *******************######## 

def funcq(r):
    return (4./pi)*log(model.Lam)*(model.r0_rT/model.MBH)*10**Ggood(log10(r))

def Rlc(r):
    interior = 2*(model.Mnorm/model.r0_rT)*(1./10**Jc2good(log10(r)))
    return -log(interior)

etest = 10**arange(-2,2,0.01)
#plt.loglog(etest,funcq(etest))
#plt.show()

########******************* IMPORT DATA TABLES *******************########

tic = time.clock()
bessel = BesselGen(['alpham_table.pkl','xi_table2.txt','Bessel_table.txt','mpiece_table.txt'])
toc = time.clock()
delt = toc-tic
print 'bessel loaded in \t {0}'.format(str(datetime.timedelta(seconds=delt)))

########******************* CALCULATE RATE *******************########

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
 
rps = arange(-5,0,0.1)
#rps = concatenate((rps,arange(-2,0,0.001)))
#rps = arange(-0.01,0,0.001)
rps = 10**rps                                
#rps *= model.rT

rps = array([  0.00000000e+00,   1.00000000e-04,   4.00000000e-04,
         9.00000000e-04,   1.60000000e-03,   2.50000000e-03,
         3.60000000e-03,   4.90000000e-03,   6.40000000e-03,
         8.10000000e-03,   1.00000000e-02,   1.21000000e-02,
         1.44000000e-02,   1.69000000e-02,   1.96000000e-02,
         2.25000000e-02,   2.56000000e-02,   2.89000000e-02,
         3.24000000e-02,   3.61000000e-02,   4.00000000e-02,
         4.41000000e-02,   4.84000000e-02,   5.29000000e-02,
         5.76000000e-02,   6.25000000e-02,   6.76000000e-02,
         7.29000000e-02,   7.84000000e-02,   8.41000000e-02,
         9.00000000e-02,   9.61000000e-02,   1.02400000e-01,
         1.08900000e-01,   1.15600000e-01,   1.22500000e-01,
         1.29600000e-01,   1.36900000e-01,   1.44400000e-01,
         1.52100000e-01,   1.60000000e-01,   1.68100000e-01,
         1.76400000e-01,   1.84900000e-01,   1.93600000e-01,
         2.02500000e-01,   2.11600000e-01,   2.20900000e-01,
         2.30400000e-01,   2.40100000e-01,   2.50000000e-01,
         2.60100000e-01,   2.70400000e-01,   2.80900000e-01,
         2.91600000e-01,   3.02500000e-01,   3.13600000e-01,
         3.24900000e-01,   3.36400000e-01,   3.48100000e-01,
         3.60000000e-01,   3.72100000e-01,   3.84400000e-01,
         3.96900000e-01,   4.09600000e-01,   4.22500000e-01,
         4.35600000e-01,   4.48900000e-01,   4.62400000e-01,
         4.76100000e-01,   4.90000000e-01,   5.04100000e-01,
         5.18400000e-01,   5.32900000e-01,   5.47600000e-01,
         5.62500000e-01,   5.77600000e-01,   5.92900000e-01,
         6.08400000e-01,   6.24100000e-01,   6.40000000e-01,
         6.56100000e-01,   6.72400000e-01,   6.88900000e-01,
         7.05600000e-01,   7.22500000e-01,   7.39600000e-01,
         7.56900000e-01,   7.74400000e-01,   7.92100000e-01,
         8.10000000e-01,   8.28100000e-01,   8.46400000e-01,
         8.64900000e-01,   8.83600000e-01,   9.02500000e-01,
         9.21600000e-01,   9.40900000e-01,   9.60400000e-01,
         9.80100000e-01])

rps = rps[1:]

realrate = array([  0.00000000e+00,   1.67623931e-08,   6.70704114e-08,
         1.50987061e-07,   2.68540076e-07,   4.20107520e-07,
         6.05648287e-07,   8.25473505e-07,   1.07985916e-06,
         1.36912840e-06,   1.69364784e-06,   2.05383223e-06,
         2.45014392e-06,   2.88309310e-06,   3.35324407e-06,
         3.86120693e-06,   4.40765754e-06,   4.99331246e-06,
         5.61896220e-06,   6.28544482e-06,   6.99367337e-06,
         7.74461703e-06,   8.53927080e-06,   9.37882966e-06,
         1.02644506e-05,   1.11975475e-05,   1.21790553e-05,
         1.32108644e-05,   1.42946151e-05,   1.54311459e-05,
         1.66229669e-05,   1.78717320e-05,   1.91801245e-05,
         2.05481363e-05,   2.19811017e-05,   2.34787930e-05,
         2.50435708e-05,   2.66804352e-05,   2.83911928e-05,
         3.01803631e-05,   3.20482782e-05,   3.39998869e-05,
         3.60366259e-05,   3.81665921e-05,   4.03923401e-05,
         4.27188113e-05,   4.51547181e-05,   4.76973813e-05,
         5.03511016e-05,   5.31371028e-05,   5.60480559e-05,
         5.90873746e-05,   6.22843374e-05,   6.56188410e-05,
         6.91278333e-05,   7.27947147e-05,   7.66542472e-05,
         8.06983829e-05,   8.49559547e-05,   8.94176335e-05,
         9.41300019e-05,   9.90694215e-05,   1.04303229e-04,
         1.09792035e-04,   1.15621399e-04,   1.21765590e-04,
         1.28253146e-04,   1.35164856e-04,   1.42468763e-04,
         1.50222354e-04,   1.58510996e-04,   1.67310778e-04,
         1.76697072e-04,   1.86782795e-04,   1.97559741e-04,
         2.09110822e-04,   2.21587562e-04,   2.35068054e-04,
         2.49616163e-04,   2.65405425e-04,   2.82650108e-04,
         3.01498380e-04,   3.22157117e-04,   3.44937632e-04,
         3.70241108e-04,   3.98480030e-04,   4.30110166e-04,
         4.65850896e-04,   5.06669269e-04,   5.53734344e-04,
         6.08614360e-04,   6.73498288e-04,   7.51634395e-04,
         8.47924562e-04,   9.70024338e-04,   1.13091795e-03,
         1.35718390e-03,   1.69284273e-03,   2.28540635e-03,
         3.72563196e-03])

realrate = realrate[1:]

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

d = dgdlnrp(rps)
plt.figure()
plt.plot(rps,d)
plt.xlabel('u')
plt.ylabel(r'$\frac{d\gamma}{d ln r_p}$')
plt.figure()
plt.loglog(rps,d)
plt.xlabel('u')
plt.ylabel(r'$\frac{d\gamma}{d ln r_p}$')

pklrfile = open('{0}/rrate.pkl'.format(directory),"wb")
pickle.dump(rps,pklrfile)
pklrfile.close()
plt.ylabel(r'$\frac{d^2\rho}{dr^2}$')
pklrfile = open('{0}/rate.pkl'.format(directory),"wb")
pickle.dump(d,pklrfile)
pklrfile.close()

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
