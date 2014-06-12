from numpy import *

#**************THINK ABOUT SORTING NECESSITIES IN CONCATENATE************
#**************ADAPTIVE QUADRATURE AND TOLERANCE?***********************

class Integrator:
    def __init__(self):
        self.memoize = {}

    def quadintegrate(self,function,uplim,downlim,div,addargs):
        pts = arange(downlim,uplim,div)
        ptsold = array([i for i in pts if i in self.memoize])
        ptsnew = array([i for i in pts if i not in self.memoize])
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
        elif addargs != '':
            valsnew = function(ptsnew,addargs)
        memolist = dict(zip(ptsnew,valsnew))
        self.memoize = dict(memolist.items() + self.memoize.items())
        vals = concatenate((valsold,valsnew))
        return sum(vals*div)

    def trapintegrate(self,function,uplim,downlim,div,addargs):
        pts = arange(downlim,uplim+div,div)
        ptsold = array([i for i in pts if i in self.memoize])
        ptsnew = array([i for i in pts if i not in self.memoize])
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
        elif addargs != '':
            valsnew = function(ptsnew,addargs)
        memolist = dict(zip(ptsnew,valsnew))
        self.memoize = dict(memolist.items() + self.memoize.items())
        vals = concatenate((valsold,valsnew))
        vals = zip(vals[::2],vals[1::2]) #makes consecutive pairs of items in the list
        vals = array([sum(i) for i in vals])
        return sum(vals*div)

    def midintegrate(self,function,uplim,downlim,div,addargs):
        pts = arange(downlim,uplim,div)
        pts = zip(pts[::2],pts[1::2])
        pts = array([sum(i) for i in pts])/2.
        ptsold = array([i for i in pts if i in self.memoize])
        ptsnew = array([i for i in pts if i not in self.memoize])
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
        elif addargs != '':
            valsnew = function(ptsnew,addargs)
        memolist = dict(zip(ptsnew,valsnew))
        self.memoize = dict(memolist.items() + self.memoize.items())
        vals = concatenate((valsold,valsnew))
        return sum(2*vals*div)

    def adaptivemid(self,function,uplim,downlim,div,tolerance,addargs):
        pts = arange(downlim,uplim,div)
        mpts = zip(pts[::2],pts[1::2])
        mpts = array([sum(i) for i in mpts])/2.
        ptsold = array([i for i in pts if i in self.memoize])
        ptsnew = array([i for i in pts if i not in self.memoize])
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
            vals = function(pts)
        elif addargs != '':
            valsnew = function(ptsnew,*addargs)
            vals = function(pts,*addargs)
        vals = zip(vals[::2],vals[1::2])
        vals = array([i[1]-i[0] for i in vals])
        err = array((vals*div**2)/3.)
        
    
    def __call__(self,function,uplims,downlims,div,tolerance = 1.49e-7,addargs = '',method = 'quad'):
        intglist = array([])
        if method == 'quad':
            for i in range(len(uplims)):
                intglist = append(intglist,self.quadintegrate(function,uplims[i],downlims[i],div,addargs[i]))
        if method == 'trap':
            for i in range(len(uplims)):
                intglist = append(intglist,self.trapintegrate(function,uplims[i],downlims[i],div,addargs[i]))
        if method == 'mid':
            for i in range(len(uplims)):
                intglist = append(intglist,self.midintegrate(function,uplims[i],downlims[i],div,addargs[i]))
        return intglist
                            
    
def parabola(x,p):
    a,b,c = p
    return a*x**2 + b*x + c


intg = Integrator()
I1 = intg(parabola,[1,2,3,4,5],[0,0,0,0,0],1e-3,addargs = [[1,2,3]]*5,method = 'quad')
I2 = intg(parabola,[1,2,3,4,5],[0,0,0,0,0],1e-3,addargs = [[1,2,3]]*5,method = 'mid')
I3 = intg(parabola,[1,2,3,4,5],[0,0,0,0,0],1e-3,addargs = [[1,2,3]]*5,method = 'trap')
print 'I1 = ',I1,'\nI2 = ',I2,'\nI3 = ',I3
