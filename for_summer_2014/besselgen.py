from numpy import *
from scipy.interpolate import interp1d

class BesselGen:
    def __init__(self,files):
        data = loadtxt('{0}'.format(files[0]))
        self.ainter = interp1d(data[:,0],data[:,1])
        data = loadtxt('{0}'.format(files[1]))
        self.xinter = interp1d(10**data[:,0],data[:,1])
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
        q = 10**q
        qimin = -10.
        qimax = 1.
        dqi = 0.03
        qi = ((log10(q) - qimin)/dqi)+1.
        xichange = 10**(floor((qimax-qimin)/dqi)*dqi + qimin)
        try:
            xi1 = sqrt(4*q[(q<10**qimin)]/pi)
            xi2 = self.xinter(q[(q>=10**qimin)&(q<=xichange)])
            xi3 = q[(q>xichange)]/qi[(q>xichange)]
            return cocatenate((xi1,xi2,xi3))
        except TypeError:
            if q < 10**qimin:
                return sqrt(4*q/pi)
            elif q >= 10**qimin and q < xichange:
                return self.xinter(q)
            elif q >= xichange:
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
            elif z>=zimin and z<=zimax:
                return self.binter(zi)
            elif z>zimax:
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
            elif m>mlim:
                return (-1)**(m-1)*sqrt(2*m - 0.5)
                   
        
