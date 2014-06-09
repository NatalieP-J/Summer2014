from numpy import *
from scipy.interpolate import interp1d

class BesselGen:
    def __init__(self,files):
        data = loadtxt('{0}'.format(files[0]))
        self.ainter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[1]))
        self.xinter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[2]))
        self.binter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[3]))
        self.minter = interp1d(data[:,0],data[:,1])

    def alpham(self,m):
        try:
            alpha1 = self.ainter(m[(m<10**6)])
            alpha2 = pi*(m[(m>=10**6)] - 0.25)
            return concatenate((alpha1,alpha2))
        except TypeError:
            if m<10**6:
                return self.ainter(m)
            elif m>=10**6:
                return pi*(m-0.25)
        
    def xi(self,q):
        qimin = -10.
        qimax = 1.
        dqi = 0.03
        qi = ((log10(q) - qimin)/dqi)+1.
        print qi
        xichange = 10**(floor((qimax-qimin)/dqi)*dqi + qimin)
        print xichange
        print 10**qimin
        try:
            xi1 = sqrt(4*q[(qi<10**qimin)]/pi)
            xi2 = self.xinter(qi[(qi>=10**qimin)&(qi<=xichange)])
            xi3 = qi[(qi>xichange)]/qi[(qi>xichange)]
            return cocatenate((xi1,xi2,xi3))
        except TypeError:
            if qi < 10**qimin:
                return sqrt(4*q/pi)
            if qi >= 10**qimin and qi <= xichange:
                return self.xinter(qi)
            if qi > xichange:
                return 1.


    def besselfin(self,u,m):
        zimin = 10**-2
        zimax = 10**2
        dzi = 10**-2
        z = u*self.alpham(m)
        zi = ((z-zimin)/dzi)+1
        print zi
        try:
            besselfin1 = zi[(z<zimin)]/zi[(z<zimin)]
            besselfin2 = self.binter(zi[(z>=zimin)&(z<=zimax)])
            besselfin3 = sqrt(2./((m-0.25)*u[(z>zimax)]*pi**2))*cos(pi*(u[(z>zimax)]*(m-0.25)-0.25))
            return concatenate((besselfin1,besselfin2,besselfin3))
        except IndexError:
            if z<zimin:
                return 1
            if z>=zimin and z<=zimax:
                return self.binter(zi)
            if z>zimax:
                return sqrt(1./((m-0.25)*u*pi**2))*cos(pi*(u*(m-0.25)-0.25))
    
    def mpiece(self,m):
        mlim = 200.
        try:
            mpiece1 = self.minter(m[(m<=mlim)])
            mpiece2 = (-1)**(m-1)*sqrt(2*m - 0.5)
            return concatenate((mpiece1,mpiece2))
        except TypeError:
            if m<=mlim:
                return self.minter(m)
            if m>mlim:
                return (-1)**(m-1)*sqrt(2*m - 0.5)
                   
        
