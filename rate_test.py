from numpy import *
import matplotlib.pyplot as plt

rT = 1
rH = 1e5
G = 6.67e-11
m = 1.98e30
lnlam = 20
n = 1

def n(r):
    return r**0

def diffusive(r):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = (n(set1)**2)*(set1**4.5)*((G*m)**0.5)*lnlam
    part2 = (n(set2)**2)*(set2**3)*((G*m)**0.5)*lnlam
    return concatenate((part1,part2))

def flc(r):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = n(set1)*(set1**0.5)*rT*((G*m)**0.5)
    part2 = n(set2)*set2*rT*((rH*G*m)**0.5)
    return concatenate((part1,part2))

def funcrp(r,rp):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = n(set1)*((set1*G*m)**0.5)*rp
    part2 = n(set2)*rp*((rH*G*m)**0.5)
    return concatenate((part1,part2))

r = arange(0,10,0.01)
r = 10**r

plt.figure()
plt.title('Rate of Tidal Disruption as a Function of Starting Radius')
plt.loglog(r,funcrp(r,0.90),'m',label = r'$r_p = 0.90$')
plt.loglog(r,funcrp(r,0.80),'DarkOrange',label = r'$r_p = 0.80$')
plt.loglog(r,funcrp(r,0.70),'r',label = r'$r_p = 0.70$')
plt.loglog(r,funcrp(r,0.60),'SpringGreen',label = r'$r_p = 0.60$')
plt.axvline(rH,color = 'b',label = 'Radius of influence')
plt.xlabel(r'$r_0 [AU]$')
plt.ylabel(r'$d\gamma/dlnr_p$')
plt.legend(loc = 'best')
plt.show()
