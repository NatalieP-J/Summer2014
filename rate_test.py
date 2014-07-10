from numpy import *
import matplotlib.pyplot as plt
from estimates import *

rT = 1
rH = 1e5
G = 6.67e-11
m = 1.98e30
MBH = 1e6*m
lnlam = 20
n = 1

def n(r,expo):
    return r**expo

def thetalc(r,expo):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = sqrt(rT/set1)
    part2 = sqrt(rT*rH/set2)
    return concatenate((part1,part2))

def thetadyn(r,expo):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    pre1 = ((m/MBH)**2)*n(set1,expo)*lnlam   
    pre2 = ((m/MBH)**2)*n(set2,expo)*lnlam
    part1 = sqrt(pre1*set1**3)
    part2 = sqrt(pre2*set2*rH**2)
    return concatenate((part1,part2))

def diffusive(r,expo):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = ((n(set1,expo)**2)*(set1**(9./2))*((G*m)**0.5)*lnlam)/((G*MBH)**1.5)
    part2 = ((n(set2,expo)**2)*(rH**1.5)*(set2**3)*((G*m)**0.5)*lnlam)/((G*MBH)**1.5)
    return concatenate((part1,part2))

def flc(r,expo):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = n(set1,expo)*(set1**0.5)*rT*((G*MBH)**0.5)
    part2 = n(set2,expo)*set2*rT*((rH*G*MBH)**0.5)
    return concatenate((part1,part2))

def comborate(r,expo):
    set1 = r[(thetadyn(r,expo) <= thetalc(r,expo))]
    set2 = r[(thetadyn(r,expo) > thetalc(r,expo))]
    return concatenate((diffusive(set1,expo),flc(set2,expo)))


def funcrp(r,rp,expo):
    set1 = r[(r<rH)]
    set2 = r[(r>=rH)]
    part1 = n(set1,expo)*((set1*G*m)**0.5)*rp
    part2 = n(set2,expo)*rp*((rH*G*m)**0.5)
    return concatenate((part1,part2))

r = arange(-10,10,0.05)
r = 10**r
'''
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
'''
expos = arange(-2.0,0.5,0.5)
plt.figure()
plt.title('Rate of Tidal Disruption as a Function of Starting Radius')
for i in range(len(expos)):
    c = comborate(r,expos[i])
    for j in range(len(c)):
        if c[j] == max(c):
            rmax = r[j]
    plt.loglog(r,c,linewidth = 3,label = r'exp = {0}'.format(expos[i]))
    plt.loglog(rmax,max(c),'b*',markersize = 20, color = 'yellow')
plt.xlabel(r'$r_0$ [AU]')
plt.ylabel(r'$d\gamma/dlnr_0$')
plt.axvline(rH,color = 'r',label = 'Radius of Influence')
plt.legend(loc = 'best')
plt.show()
