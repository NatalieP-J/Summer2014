from numpy import *

#**************THINK ABOUT SORTING NECESSITIES IN CONCATENATE************
#**************ADAPTIVE QUADRATURE AND TOLERANCE?***********************

class Integrator:
    def __init__():
        self.memoize = {}

    def quadintegrate(self,function,uplim,downlim,div,addargs):
        pts = arange(downlim,uplim,div)
        ptsold = pts[(not i in self.memoize == False)]
        ptsnew = pts[(not i in self.memoize == True)]
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
        ptsold = pts[(not i in self.memoize == False)]
        ptsnew = pts[(not i in self.memoize == True)]
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
        elif addargs != '':
            valsnew = function(ptsnew,addargs)
        memolist = dict(zip(ptsnew,valsnew))
        self.memoize = dict(memolist.items() + self.memoize.items())
        vals = concatenate((valsold,valsnew))
        vals = zip(vals[::2],vals[1::2]) #makes consecutive pairs of items in the list
        vals = array([sum(i) for i in vals])/2.
        return sum(vals*div)

    def midintegrate(self,function,uplim,downlim,div,addargs):
        pts = arange(downlim,uplim,div)
        pts = zip(pts[::2],pts[1::2])
        pts = array([sum(i) for i in vals])/2.
        ptsold = pts[(not i in self.memoize == False)]
        ptsnew = pts[(not i in self.memoize == True)]
        valsold = array([self.memoize[i] for i in ptsold])
        if addargs == '':
            valsnew = function(ptsnew)
        elif addargs != '':
            valsnew = function(ptsnew,addargs)
        memolist = dict(zip(ptsnew,valsnew))
        self.memoize = dict(memolist.items() + self.memoize.items())
        vals = concatenate((valsold,valsnew))
        return sum(vals*div)
    
    def __call__(self,function,uplims,downlims,div,addargs = '',method):
        intglist = array([])
        if method == 'quad':
            for i in range(len(uplims)):
                intglist = append(intglist,self.quadintegrate(function,uplims[i],downlims[i],div,addargs))
        if method == 'trap':
            for i in range(len(uplims)):
                intglist = append(intglist,self.trapintegrate(function,uplims[i],downlims[i],div,addargs))
        if method == 'mid':
            for i in range(len(uplims)):
                intglist = append(intglist,self.midintegrate(function,uplims[i],downlims[i],div,addargs))
        return intglist
                            
    
