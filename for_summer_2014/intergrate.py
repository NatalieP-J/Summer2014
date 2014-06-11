from numpy import *

#**************THINK ABOUT SORTING NECESSITIES IN CONCATENATE************

class Integrator:
    def __init__():
        self.memoize = {}

    def quadintegrate(self,function,uplim,downlim,tolerance,addargs):
        pts = arange(downlim,uplim,tolerance)
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
        return sum(vals*tolerance)

    def trapintegrate(self,function,uplim,downlim,tolerance,addargs):
        pts = arange(downlim,uplim+tolerance,tolerance)
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
        vals = zip(vals[::2],vals[1::2])
        vals = array([sum(i) for i in vals])/2.
        return sum(vals*tolerance)

    def midintegrate(self,function,uplim,downlim,tolerance,addargs):
        pts = arange(downlim,uplim,tolerance)
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
        return sum(vals*tolerance)
    
    def __call__(self,function,uplims,downlims,tolerance,addargs = '',method):
        intglist = array([])
        if method == 'quad':
            for i in range(len(uplims)):
                intglist = append(intglist,self.quadintegrate(function,uplims[i],downlims[i],tolerance,addargs))
        if method == 'trap':
            for i in range(len(uplims)):
                intglist = append(intglist,self.trapintegrate(function,uplims[i],downlims[i],tolerance,addargs))
        if method == 'mid':
            for i in range(len(uplims)):
                intglist = append(intglist,self.midintegrate(function,uplims[i],downlims[i],tolerance,addargs))
        return intglist
                            
    
