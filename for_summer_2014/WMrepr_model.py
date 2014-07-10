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
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600
alpha = 2.94#7.52 #1.0
beta = 2.23#3.13#4.0
gamma = 1.80#1.98#1.5
r0pc = 1.
rb = 10**2.46#10**2.38
r0pc = rb
mub = 18.83#19.98
M2L = 7.25#6.27
MsunV = 4.83
rho0 = 1e5
rho0 = (1./rb)*(1./(10)**2)*(206265**2)*M2L*10**((MsunV-mub)/2.5) 
MBH_Msun = 10**7.11#10**6.04#1e3
masses = [4,6,8,10,12]
galname = 'NGC4551'
rapos = []
fs = []
qs = []
Rlcs = []

plotarrays = [arange(0.9,2.1,0.01),arange(0.9,2.1,0.01),arange(0.9,4,0.01),arange(0.9,4,0.01),arange(0.9,4,0.01)]

generates = [False,False,False,False,False]

for i in range(len(masses)):
    MBH_Msun = 10**masses[i]
    print 'MBH = ',MBH_Msun
    Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
    generate = generates[i]
    seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
    verbosity = {Menc:"ON",psi:"ON",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
    plot = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
    plotarray = 10**plotarrays[i]

########******************* MODEL FRAMEWORK *******************########
    class NukerModel:
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
        #start a new directory?
            self.generate = generate
    
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
                                
########******************* CONSTRUCT MODEL *******************########

    model = NukerModel(galname,alpha,beta,gamma,r0pc,rho0,MBH_Msun,generate)

    rtest = arange(-7,6,0.01)
    rtest = append(rtest,40)
    rtest = insert(rtest,0,-40)
    rtest = 10**rtest

    directory = "{0}_a{1}_b{2}_g{3}_r{4}_rho{5}_MBH{6}".format(model.name,model.a,model.b,model.g,model.r0,model.rho0,model.MBH)
    if model.generate == True:
        call(["mkdir","{0}".format(directory)])
        seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}
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

    def plotter(name,r,inter,rstart,rchange,start,end,smallrexp,largerexp,conds,labels):
        m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp,conds)
        plt.figure()
        plt.loglog(r[1:-1],m[1:-1],'c',linewidth = 5)
        plt.ylabel(r'{0}'.format(labels[1]))
        plt.xlabel('{0}'.format(labels[0]))
        plt.xlim(min(r[1:-1]),max(r[1:-1]))
        plt.ylim(min(m[1:-1]),max(m[1:-1]))
        plt.savefig('{0}/{1}.png'.format(directory,name))

    def makegood(func,r,size,grid,smallrexp,largerexp,verbose = False,conds = False,plotting=False,problem = True):
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
        return model.rho(r)*r**2

    tol = 1.e-3
    def funcMenc(r,verbose=False):
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
        return abs(model.Mnorm-funcMenc(abs(r))[0])

    def rH(verbose=True):
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

########******************* POTENTIAL *******************######## 
        
    def psi2interior(r):
        return model.rho(r)*r

    tol = 1e-3
    def psi2(r,verbose=False):
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

########******************* APOCENTER RADIUS *******************######## 

    def rapoimplicit(r,E):
        return abs(10**psigood(log10(abs(r)))-E)

    def rapo(e):
        try:
            rs = []
            for i in range(len(e)):
                E = e[i]
                if E**-1 > 0.2:
                    rguess = 10*E**-1
                elif E**-1 < 0.2:
                    rguess = 0.01*E**-1
                rresult = root(rapoimplicit,rguess,args=E)
                if rresult.success == True:
                    rs.append(abs(rresult.x))
                elif rresult.success == False:
                    print 'Failed to evaluate rapo'
                    print rresult.message
                    rs.append(abs(rresult.x))
            return array(rs)
        except (TypeError,AttributeError):
            E = e
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
        return abs(10**psigood(log10(abs(r)))-E-((10**Mencgood(log10(abs(r)))+model.Mnorm)/(2*r)))

    def funcJc2(E,verbose):
        try:
            t = E.shape
            Jcs = []
            problems = []
            for i in range(len(E)):
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
        return (model.drhodr(1./r))*(r**-2)*((sqrt(abs(E-10**psigood(log10(1./r)))))**-1)
    
    def funclg(E,verbose=False):
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
        #if not r in psibG_memo:
        #    psibG_memo[r] = 10**psigood(log10(r))
        #psir = psibG_memo[r]
        psir = 10**psigood(log10(r))
        part1 = (r**2)/sqrt(psir-E)
        #if not log10(psir*(1-theta) + E*theta) in part2bG_memo:
        #    part2bG_memo[log10(psir*(1-theta) + E*theta)] = 10**ggood(log10(psir*(1-theta) + E*theta))
        #part2 = part2bG_memo[log10(psir*(1-theta) + E*theta)]
        part2 = 10**ggood(log10(psir*(1-theta) + E*theta))
        #if not theta in part3bG_memo:
        #    part3bG_memo[theta]= (1./sqrt(theta))-sqrt(theta)
        #part3 = part3bG_memo[theta]
        part3 = (1./sqrt(theta))-sqrt(theta)
        return part1*part2*part3

    def funcbG(E,verbose = False):
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
        return interior

    def lnRlc(r):
        return -log(Rlc(r))

    rapos.append(rapo(plotarray))
    fs.append(10**(fgood(log10(plotarray))))
    qs.append(funcq(plotarray))
    Rlcs.append(Rlc(plotarray))

plt.figure()
for i in range(len(rapos)):
    plt.loglog(10**plotarrays[i],rapos[i],label = 'log10(mass) = {0}'.format(masses[i]))
plt.xlabel('E')
plt.ylabel(r'$r_{apo}$')
plt.title(galname)
plt.legend(loc = 'best')
plt.savefig('{0}/rapofig.png'.format(galname))

plt.figure()
for i in range(len(fs)):
    plt.loglog(10**plotarrays[i],fs[i],label = 'log10(mass) = {0}'.format(masses[i]))
plt.xlabel('E')
plt.ylabel('f')
plt.title(galname)
plt.legend(loc = 'best')
plt.savefig('{0}/ffig.png'.format(galname))

plt.figure()
for i in range(len(qs)):
    plt.loglog(10**plotarrays[i],qs[i],label = 'log10(mass) = {0}'.format(masses[i]))
plt.xlabel('E')
plt.ylabel('q')
plt.title(galname)
plt.legend(loc = 'best')
plt.savefig('{0}/qfig.png'.format(galname))

plt.figure()
for i in range(len(Rlcs)):
    plt.loglog(10**plotarrays[i],Rlcs[i],label = 'log10(mass) = {0}'.format(masses[i]))
plt.xlabel('E')
plt.ylabel(r'$R_{lc}$')
plt.title(galname)
plt.legend(loc = 'best')
plt.savefig('{0}/Rlcfig.png'.format(galname))

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
    part1 = array(fval/(1+(qval**-1)*(xival)*lnRlc(E)))
    part2list = exp(-array(matrix(alphas**2).T*matrix(qval/4)))
    part2list = array([(bfin/mpiece)[i]*part2list[i] for i in range(len(alphas))])
    part2 = 1-2*nsum(part2list,axis = 0)
    return part1*part2
 
rps = arange(-5,0,0.1)
#rps = concatenate((rps,arange(-2,0,0.001)))
#rps = arange(-0.01,0,0.001)
rps = 10**rps                                
#rps *= model.rT

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
'''
d = dgdlnrp(rps)
plt.figure()
plt.plot(rps,d)
plt.xlabel('u')
plt.ylabel(r'$\frac{d\gamma}{d ln r_p}$')
plt.figure()
plt.loglog(rps,d)
plt.xlabel('u')
plt.ylabel(r'$\frac{d\gamma}{d ln r_p}$')
'''
def ginterior(E):
    qval = funcq(E)
    return (10**fgood(log10(E))*qval)/((qval/bessel.xi(qval)) + lnRlc(E))

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
