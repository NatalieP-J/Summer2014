from numpy import *
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
import scipy.integrate as intg
import time
import datetime
import os
from matplotlib.backends.backend_pdf import PdfPages

devnull = open(os.devnull,'w')

def pklread(fname):
    pklffile = open(fname,"rb")
    dat = pickle.load(pklffile)
    pklffile.close()
    return dat

def pklwrite(fname,dat):
    pklffile = open(fname,"wb")
    pickle.dump(dat,pklffile)
    pklffile.close()

def piecewise2(r,inter,start,end,lim1,lim2,smallrexp,largerexp):
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
    piece3 = end*(set3/lim2)**largerexp
    return concatenate((piece1,piece2,piece3))

def makegood(prereqs,func,r,size,grid,smallrexp,largerexp,verbose,plotting,problem = True):
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
    plotting = False - 	do not save plots
    problem = True - eliminate problem points
    
    returns an interpolated object version of the function based 
    computed values
    
    """
    model = prereqs[0]
    rarray,rchange,rstart = grid([model],size[0],size[1],size[2])
    tab,problems = func(rarray,verbose,prereqs)
    frac = float(len(problems))/float(len(tab))
    print 'fraction reporting a message: {0}'.format(frac)
    model.statfile.write('mesg frac = {0}\n'.format(frac))
    if frac != 1.0:
        if problem == True:
            tab = [i for j, i in enumerate(tab) if j not in problems]
            rarray = [i for j, i in enumerate(rarray) if j not in problems]
        lrarray = log10(rarray)
        ltab = log10(tab)
        inter = interp1d(lrarray,ltab)
        start = tab[0]
        end = tab[len(rarray)-1]
        m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp)
        inter2 = interp1d(log10(r),log10(m))
        pklrfile = open('{0}/r{1}.pkl'.format(model.directory,str(func)[10:15]),"wb")
        pickle.dump(r,pklrfile)
        pklrfile.close()
        pklffile = open('{0}/{1}.pkl'.format(model.directory,str(func)[10:15]),"wb")
        pickle.dump(m,pklffile)
        pklffile.close()
        
        if plotting != False:
            xaxis,yaxis = plotting
            plt.figure()
            plt.loglog(r[1:-1],m[1:-1],'c',linewidth = 5)
            plt.loglog(rarray,tab,'.',color = 'DarkOrange')
            plt.ylabel(r'{0}'.format(xaxis))
            plt.xlabel('{0}'.format(yaxis))
            plt.xlim(min(r[1:-1]),max(r[1:-1]))
            plt.ylim(min(m[1:-1]),max(m[1:-1]))
            plt.title(model.name)
            model.pdfdump.savefig()
            plt.close()

        return inter2
    elif frac == 1.0:
        return 0

def compute(dependencies,name,function,rtest,size,grid,exps,kwargs,create):
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
    model = dependencies[0]
    prereqs = dependencies[0::2]
    strname,name = name
    gen,verbosity = create
    if gen == "ON":
        try: 
            if len(dependencies) > 2:
                i = 2
                while i<len(dependencies):
                    dependencies[i](1)
                    i+=2
            smallrexp,largerexp = exps
            plotdat,prob = kwargs
            if verbosity == "ON":
                tic = time.clock()
                good = makegood(prereqs,function,rtest,size,grid,smallrexp,largerexp,verbose = True,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity == "OFF":
                tic = time.clock()
                good = makegood(prereqs,function,rtest,size,grid,smallrexp,largerexp,verbose = False,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            print '{0}good ran in \t {1}'.format(strname,str(datetime.timedelta(seconds=delt)))
            model.statfile.write('{0}good ran in \t {1}\n'.format(strname,str(datetime.timedelta(seconds=delt))))
            return good
        except TypeError as e:
            print 'e = ',e
            print 'To compute {0}, please turn {1} ON'.format(strname,dependencies[i+1])
            model.statfile.write('To compute {0}, please turn {1} ON\n'.format(strname,dependencies[i+1]))
    elif gen != "ON":
        try:
            pklrfile = open('{0}/r{1}.pkl'.format(model.directory,str(function)[10:15]),"rb")
            rarray = pickle.load(pklrfile)
            pklrfile.close()
            pklffile = open('{0}/{1}.pkl'.format(model.directory,str(function)[10:15]),"rb")
            tab = pickle.load(pklffile)
            pklffile.close()
            good =  interp1d(log10(rarray),log10(tab))
            print '{0}good loaded in'.format(strname)
            return good
        except IOError:
            model.statfile.write('{0} not yet generated, please turn it ON\n'.format(strname))
            return 1


def integrator(vals,fcn,downlim,uplim,tol=1.49e-7,args = [],fileobj=devnull,prefactor = 1):
    if isinstance(vals,(list,ndarray))==True:
        problems = []
        results = []
        for i in range(len(vals)):
            if args != []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],args = args[i],epsabs = tol,full_output = 1)
            elif args == []:
                temp = intg.quad(fcn[0],downlim[i],uplim[i],epsabs = tol,full_output = 1)
            try:
                if temp[3] != '':
                    problems.append(i)
                    fileobj.write('{0},\t i = {1},\t message = {2}\n'.format(fcn[1],i,temp[3]))
            except IndexError:
                pass
            if prefactor != 1:
                results.append(prefactor[i]*temp[0])
            elif prefactor == 1:
                results.append(temp[0])
        return array(results),problems
    elif isinstance(vals,(int,float))==True:
        problems = []
        if args != []:
            temp = intg.quad(fcn[0],downlim,uplim,args = args,epsabs = tol,full_output = 1)
        elif args == []:
            temp = intg.quad(fcn[0],downlim,uplim,epsabs = tol,full_output = 1)
        try:
            if temp[3] != '':
                problems.append(vals)
                fileobj.write('\n{0},\t i = {1},\t message = {2}'.format(fcn[1],i,temp[3]))
        except IndexError:
            pass
        return prefactor*temp[0],problems


def dblintegrator(vals,fcn,downlim,uplim,tol=1.49e-7,args = [],fileobj=devnull,prefactor = 1):
    if isinstance(vals,(list,ndarray))==True:
        problems = []
        results = []
        for i in range(len(vals)):
            with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
                if args != []:
                    temp = intg.dblquad(fcn[0],downlim[0][i],uplim[0][i],downlim[1][i],uplim[1][i],args = args[i],epsabs = tol)
                elif args == []:
                    temp = intg.dblquad(fcn[0],downlim[0][i],uplim[0][i],downlim[1][i],uplim[1][i],epsabs = to1)
            if prefactor != 1:
                results.append(prefactor[i]*temp)
            elif prefactor == 1:
                results.append(temp)
        return array(results),problems
    elif isinstance(vals,(int,float))==True:
        problems = []
        with RedirectStdStreams(stdout = fileobj,stderr = fileobj):
            if args != []:
                temp = intg.dblquad(fcn[0],downlim[0],uplim[0],downlim[1],uplim[1],args = args,epsabs = tol)
            elif args == []:
                temp = intg.dblquad(fcn[0],downlim[0],uplim[0],downlim[1],uplim[1],epsabs = tol)
        return prefactor*temp,problems

def fromfileplot(name,directory):
    r = arange(-4,4,0.01)
    r = 10**r
    pp = PdfPages('{0}/{1}_master.pdf'.format(directory,name))
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcM'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcM'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$r$')
    plt.ylabel(r'$M_{enc}$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcp'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcp'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\psi$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcJ'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcJ'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$E$')
    plt.ylabel(r'$J_c^2$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcl'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcl'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$E$')
    plt.ylabel(r'$g$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcb'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcb'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$E$')
    plt.ylabel(r'$G$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pklrfile = open('{0}/r{1}.pkl'.format(directory,'funcf'),"rb")
    rarray = pickle.load(pklrfile)
    pklrfile.close()
    pklffile = open('{0}/{1}.pkl'.format(directory,'funcf'),"rb")
    tab = pickle.load(pklffile)
    pklffile.close()
    Mgood =  interp1d(log10(rarray),log10(tab))
    plt.figure()
    plt.loglog(r,10**Mgood(log10(r)))
    plt.xlabel(r'$E$')
    plt.ylabel(r'$f$')
    plt.title(name)
    pp.savefig()
    plt.close()
    pp.close()
