from numpy import *
from scipy.interpolate import interp1d
import time
import pickle
class BesselGen:
    def __init__(self,files):
        datafile = open('{0}'.format(files[0]),"rb")
        data = pickle.load(datafile)
        self.ainter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[1]))
        self.xinter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[2]))
        self.binter = interp1d(data[:,0],data[:,1])
        
    def alpham(self,m):
        try:
            alpha1 = self.ainter(m[(m<10**6)])
            alpha2 = pi*(m[(m>=10**6)] - 0.25)
            return concatenate((alpha1,alpha2))
        except (TypeError,IndexError):
            if m<10**6:
                return self.ainter(m)
            elif m>=10**6:
                return pi*(m-0.25)
        
    def xi(self,q):
        qimin = -10.
        qimax = 1.
        dqi = 0.03
        xichange = 10**(floor((qimax-qimin)/dqi)*dqi + qimin)
        try:
            xi1 = sqrt(4*q[(q<10**qimin)]/pi)
            xi2 = self.xinter(q[(q>=10**qimin)&(q<=xichange)])
            xi3 = q[(q>xichange)]/q[(q>xichange)]
            return concatenate((xi3,xi2,xi1))
        except (TypeError,IndexError):
            if q < 10**qimin:
                return array(sqrt(4*q/pi))
            elif q >= 10**qimin and q < xichange:
                return self.xinter(q)
            elif q >= xichange:
                return array(1.)

    def besselfin(self,m,u):
        zimin = 10**-2
        zimax = 10**2
        dzi = 10**-2
        z = self.alpham(m)*u
        zi = ((z-zimin)/dzi)+1
        try:
            besselfin1 = zi[(zi<zimin)]/zi[(zi<zimin)]
            besselfin2 = self.binter(zi[(zi>=zimin)&(zi<=zimax)])
            besselfin3 = sqrt(2./((m[(zi>zimax)]-0.25)*u*pi**2))*cos(pi*(u*(m[(zi>zimax)]-0.25)-0.25))
            return concatenate((besselfin1,besselfin2,besselfin3))
        except (IndexError,TypeError) as e:
            if zi<zimin:
                print 'in 1'
                return 1
            elif zi>=zimin and zi<=zimax:
                print 'in 2'
                return self.binter(zi)
            elif zi>zimax:
                print 'in 3'
                return sqrt(2./((m-0.25)*u*pi**2))*cos(pi*(u*(m-0.25)-0.25))
