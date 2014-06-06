from numpy import *
import scipy.integrate as intg
from scipy.optimize import root
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time
import datetime
from subprocess import call
import pickle

Lam = exp(1)# ****************************************************************
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e16
km = 10**5
yr = 365*24*3600
Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
generate = True
seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
verbosity = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
plot = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}

rtest = arange(-12,12,0.01)
rtest = 10**rtest

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

########******************* PARAMETER TESTS *******************########

alphas = arange(1.0,10)
betas = arange(1.0,10)
gammas = arange(1.0,10)

for a in range(len(alphas)):
    alpha = alphas[a]
    for b in range(len(betas)):
        beta = betas[b]
        for gam in range(len(gammas)):
            gamma = gammas[gam]
            print 'alpha = ',alpha, ' beta = ',beta, ' gamma = ',gamma
                                
########******************* CONSTRUCT MODEL *******************########

            model = NukerModel('testing',alpha,beta,gamma,1.,1.e5,1.e3,generate)            
            directory = "{0}_a{1}_b{2}_g{3}_r{4}_rho{5}_MBH{6}".format(model.name,model.a,model.b,model.g,model.r0,model.rho0,model.MBH)
            statfile = open('{0}/status.txt'.format(directory),"wb")
            print 'Open statfile'
            statfile.write('{0}\t{1}\t{2}\n'.format(alpha,beta,gamma))
            if model.generate == True:
                call(["mkdir","{0}".format(directory)])
                seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"OFF",f:"OFF"}
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
                nanstrack = []
                for i in range(len(m)):
                    if math.isnan(m[i]) == True:
                        nanstrack.append(i)
                nansfrac = float(len(nanstrack))/float(len(m))
                if nansfrac > 0.3:
                    return nansfrac
                else:
                    plt.figure()
                    if any(t<0 for t in m):
                        plt.loglog(r,abs(m),'.',label = 'Negative values')
                    else:
                        plt.loglog(r,m,'.')
                        plt.ylabel(r'{0}'.format(labels[1]))
                        plt.xlabel('{0}'.format(labels[0]))
                        plt.xlim(min(r),max(r))
                        plt.ylim(min(m),max(m))
                        if rstart>r[0] and rchange<r[len(r)-1]:
                            plt.axvline(rstart, color='r',label='Limits of interpolation')
                            plt.axvline(rchange, color='r')
                        plt.legend(loc='best')
                        plt.savefig('{0}/{1}.png'.format(directory,name))
                    return 0

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
                rarray,rchange,rstart = grid(size[0],size[1],size[2],size[3],size[4])
                tab,problems = func(rarray,verbose)
                print 'fraction reporting a message: {0}'.format(float(len(problems))/float(len(tab)))
                prob = float(len(problems))/float(len(tab))
                if problem == True:
                    tab = [i for j, i in enumerate(tab) if j not in problems]
                    rarray = [i for j, i in enumerate(rarray) if j not in problems]
                #if any(math.isnan(t) == True for t in tab):
                #    return 'Fail'
                #print 'tab = ',tab
                nanstrack = []
                for i in range(len(tab)):
                    if math.isnan(tab[i]) == True:
                        nanstrack.append(i)
                #print 'nanstrack = ',nanstrack
                #print 'nan fraction = ',float(len(nanstrack))/float(len(tab))
                if float(len(nanstrack))/float(len(tab)) > 0.3:
                    return 'Fail'
                else:
                    inter = interp1d(log10(rarray),log10(tab))
                    pklrfile = open('{0}/r{1}.pkl'.format(directory,str(func)[10:15]),"wb")
                    pickle.dump(rarray,pklrfile)
                    pklrfile.close()
                    pklffile = open('{0}/{1}.pkl'.format(directory,str(func)[10:15]),"wb")
                    pickle.dump(tab,pklffile)
                    pklffile.close()
                    start = tab[0]
                    end = tab[len(rarray)-1]
                    if plotting != False:
                       p = plotter(str(func)[10:15],r,inter,rstart,rchange,start,end,smallrexp,largerexp,conds,plotting)
                       if p > 0:
                           return 'Fail'
                       else:
                           return inter
                    if plotting == False:
                        return inter

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
                        temp = intg.quad(Minterior,0,r[i],full_output=1)
                        try:
                            if temp[3]!='':
                                problems.append(i)
                                if verbose==True:
                                    print 'Menc, r = ',r[i],'message = ',temp[3],'\n'
                        except (IndexError, TypeError):
                            pass
                        if r[i] > 10**10:
                            Mencs.append(10.)
                        elif temp[0] >= 0:
                            Mencs.append(4*pi*temp[0])
                        elif temp[0] < 0 or r[i] > 10**10:
                            Mencs.append(10.)
                    return array(Mencs),array(problems)
                except AttributeError:
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
                return abs(model.Mnorm-funcMenc(r)[0])
            
            def rH():
                """
                finds root of rHimplicit
                """
                rresult=root(rHimplicit,1e-4)
                if rresult.success == True:
                    return rresult.x
                elif rresult.success == False:
                    #print 'Failed to evaluate rH'
                    #print rresult.message
                    return rresult.x

########******************* CONSTRUCTION FUNCTIONS *******************########

            def rgrid(upstep=5,downstep=-5,up=12,down=-12,step=0.03):
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
                rarray = append(rarray,up)
                rarray = insert(rarray,0,down)
                rarray = 10**rarray
                rchange = rarray[len(rarray)-1]
                rstart = rarray[0]
                return rarray,rchange,rstart

            rarray,rchange,rstart = rgrid(5,-5,0.03)

            def Egrid(upstep=5,downstep=-3,up=12,down=-12,step=0.1):
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
                Earray = append(Earray,up)
                Earray = insert(Earray,0,down)
                Earray = 10**Earray
                Echange = Earray[len(Earray)-1]
                Estart = Earray[0]
                return Earray,Echange,Estart

########******************* COMPUTE MENC *******************######## 

            Mencgood = compute([],["Menc",Menc],funcMenc,rtest,[3,-3,40,-40,0.03],rgrid,[3-model.g,0],[[2,0,3-model.b,4*pi*model.rho(rchange)*(rchange**3)],['r','M'],False])
            if Mencgood == 'Fail':
                print 'Mencgood = ',Mencgood
                statfile.write('Mencgood = {0}'.format(Mencgood))
                statfile.close()
                continue

########******************* POTENTIAL *******************######## 
        
            def psi2interior(r):
                """
                interior of psi part 2 integral
                """
                return model.rho(r)*r
            
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
                        temp=intg.quad(psi2interior,r[i],inf,full_output=1)
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

            psigood = compute([],["psi",psi],funcpsi,rtest,[3,-3,40,-40,0.03],rgrid,[-1,-1],[False,['r','$\psi$'],False])

            if psigood == 'Fail':
                print 'psigood = ',psigood
                statfile.write('psigood = {0}'.format(psigood))
                statfile.close()
                continue
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
                    #print 'Failed to evaluate rapo'
                    #print rresult.message
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
                        if E[i]**-1 > 0.2:
                            rguess = 10*E[i]**-1
                        elif E[i]**-1 < 0.2:
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
                except IndexError as e:
                    print 'Jc2 err = ', e
                    print E
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

            Jc2good = compute(prereqs,["Jc2",Jc2],funcJc2,rtest,[3,-3,12,-12,0.01],Egrid,[-1,-1],[False,['E','Jc2'],False])

            if Jc2good == 'Fail':
                print 'Jc2good = ',Jc2good
                statfile.write('Jc2good = {0}'.format(Jc2good))
                statfile.close()
                continue

########******************* g *******************######## 

            def lginterior(r,E):
                """
                interior of g integral
                """
                return (model.drhodr(1./r))*(r**-2)*((sqrt(abs(E-10**psigood(log10(1./r)))))**-1)
    
            def funclg(E,verbose=False):
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

            ggood = compute(prereqs,["g",g],funclg,rtest,[3,-3,20,-20,0.1],Egrid,[model.b-0.5,model.g-0.5],[False,['E','g'],False])

            if ggood == 'Fail':
                print 'ggood = ',ggood
                statfile.write('ggood = {0}'.format(ggood))
                statfile.close()
                continue

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
                            temp = intg.dblquad(bGinterior,0,rapoval,lambda r: 0, lambda r: 1,args = (E[i],),epsabs = tolerance,epsrel = tolerance)
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
            
            Ggood = compute(prereqs,["G",G],funcbG,rtest,[2,-2,12,-12,0.1],Egrid,[model.b-4,model.g-4],[False,['E','G'],False])

            psibG_memo = {}
            part2bG_memo = {}
            part3bG_memo = {}

            if Ggood == 'Fail':
                print 'Ggood = ',Ggood
                statfile.write('Ggood = {0}'.format(Ggood))
                statfile.close()
                continue

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
                        except TypeError:
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
                        if temp[3] != '':
                            if verbose == True:
                                print 'f, E = ',E,'message = ',temp1[3]
                            problem = [E]
                    except IndexError:
                        pass
                    return prefactor*t, problem

########******************* COMPUTE f *******************######## 
            prereqs = [Mencgood,"Menc",psigood,"psi"]
            
            fgood = compute(prereqs,["f",f],funcf,rtest,[5,-3,12,-12,0.03],Egrid,[model.b-1.5,model.g-1.5],[False,['E','f'],False])

            if fgood == 'Fail':
                print 'fgood = ',fgood
                statfile.write('fgood = {0}'.format(fgood))
                statfile.close()
                continue
            
            status = 'Pass'
            statfile.write(status)
            statfile.close()
            

########******************* ADDITIONAL FUNCTIONS *******************######## 

            def funcq(r):
                return (4./pi)*log(Lam)*(model.r0_rT/model.MBH)*10**Ggood(log10(r))


'''
def Rlc(r):
    interior = 2*(model.MBHnorm./model.r0_rT)*(1./Jc2good(r))
    return -log(interior)

# dependent on a lot of mystery functions
def dgdlnrp(Emin = 0.01,Emax=100):
    prefactor = (8*pi**2)*model.MBH_Msun*(model.r0_rT**-1)*(model.tdyn0**-1)
    qmin = funcq(Emax)
    qmax = funcq(Emin)
'''
