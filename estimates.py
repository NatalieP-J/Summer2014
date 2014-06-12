from numpy import *

Rsun = 6.95e8
Msun = 1.98e30
MBH_MW = (4e6)*Msun
G = 6.67e-11
c = 3e8
rh_MW = 3e16
Rbulge = 2*(1e3)*(3e16)
Mbulge = (2e10)*Msun
v_typical = 2e5
Nstar_MW = MBH_MW/Msun

tunits = {'s':1.,'min':60.,'hour':3600.,'day':86400,'year':3.1536e7,'decade':3.1536e9,'kyr':3.1536e10,'Myr':3.1536e13,'Gyr':3.1536e16}

def time_units(num,units):
    return num/tunits[units]
def v_inrh(MBH,r):
    return sqrt((G*MBH)/r)
def v_outrh(MBH,rh):
    return sqrt((G*MBH)/rh)
def numdense(rh):
    return 1./((4./3.)*pi*(rh**3))
def Rs(MBH):
    return (2*G*MBH)/c**2
def rT(MBH):
    return Rsun*(MBH/Msun)**(1./3.)
def sigma(MBH):
    return pi*(rT(MBH)**2)
def rmin(v):
    return (2*G*Msun)/v**2
def Coulomb(R,v):
    return log(R/rmin(v))
def trel(v,n,R,rmin):
    return ((v**3)/(((G*m)**2)*n*Coulomb(R,v)))
def basicrate(v,rh,MBH,Nstar,units = 'Gyr'):
    return 1./time_units(1./(Nstar*numdense(rh)*v*sigma(MBH)),units)
def gravrate(v,rh,MBH,Nstar,units = 'kyr'):
    return 1./time_units(1./(Nstar*(numdense(rh)*pi*G*MBH*rT(MBH))/v),units)
def Rlc(r,rh,MBH):
    if r>=rh:
        return (rT(MBH)*rh)/r**2
    if r<rh:
        return (rT(MBH))/r

print 'Basic rate = ',basicrate(v_typical,rh_MW,MBH_MW,Nstar_MW),'per Gyr'
print 'Focus rate = ',gravrate(v_typical,rh_MW,MBH_MW,Nstar_MW),'per kyr'
