from numpy import *
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d

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
    if conds==False:
        piece3 = end*(set3/lim2)**largerexp
        return concatenate((piece1,piece2,piece3))

def plotter(name,r,inter,rstart,rchange,start,end,smallrexp,largerexp,labels):
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
    labels - plot axis labels, x and then y

    saves a plot of the piecewise function
    """
    m = piecewise2(r,inter,start,end,rstart,rchange,smallrexp,largerexp)
    plt.figure()
    plt.loglog(r[1:-1],m[1:-1],'c',linewidth = 5)
    plt.ylabel(r'{0}'.format(labels[1]))
    plt.xlabel('{0}'.format(labels[0]))
    plt.xlim(min(r[1:-1]),max(r[1:-1]))
    plt.ylim(min(m[1:-1]),max(m[1:-1]))
    plt.savefig('{0}/{1}.png'.format(directory,name))

def makegood(func,r,size,grid,smallrexp,largerexp,verbose = False,plotting=False,problem = True):
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
            plotdat,prob = kwargs
            if verbosity[name] == "ON" and plot[name] == "ON":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,verbose = True,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "ON" and plot[name] == "OFF":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,verbose = True,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "OFF" and plot[name] == "ON":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,plotting = plotdat,problem = prob)
                toc = time.clock()
                delt = toc-tic
            elif verbosity[name] == "OFF" and plot[name] == "OFF":
                tic = time.clock()
                good = makegood(function,rtest,size,grid,smallrexp,largerexp,problem = prob)
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
