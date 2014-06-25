from numpy import *
import math
import matplotlib.pyplot as plt
import matplotlib

class SersicModel:
	def __init__(self,model_name,n,rho0,Re,M2L,I0):
		self.name = model_name
		self.n = n
		self.Re = Re
		self.M2L = M2L
		self.I0 = I0
		if self.n <10. and self.n > 0.6:
			self.b = 2*self.n - (1./3.) + 0.009876/self.n
			self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
		self.rho0 = M2L*I0*(self.b**(n*(1-self.p)))*(math.gamma(2*n)/(2*Re*math.gamma(n*(3-self.p))))
		#self.rho0 = rho0

	def rho(self,r):
		return self.rho0*((r/self.Re)**-self.p)*exp(-self.b*((r/self.Re)**(1./self.n)))

	def drhodr(self,r):
		pre = self.rho0*((r/self.Re)**-p)*exp(-b*((r/self.Re)**(1./self.n)))*((self.n*r)**-1)
		post = (self.n*self.p)+self.b*((r/self.Re)**(1./n))
		return pre*post

	def d2rhodr2(self,r):
		pre = self.rho0*((r/self.Re)**-p)*exp(-b*((r/self.Re)**(1./self.n)))*((self.n*r)**-2)
		post = (self.p*(1+self.p)*self.n**2) + self.b*(-1 + self.n + 2*self.n*self.p)*((r/self.Re)**(1./n)) + (b**2)*((r/self.Re)**-p)
		return pre*post


model1 = SersicModel('plot',1,1e5,1.,10.,1e32)
model2 = SersicModel('plot',2,1e5,1.,10.,1e32)
model3 = SersicModel('plot',4,1e5,1.,10.,1e32)
rtest = arange(-3,3,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 25}
matplotlib.rc('font', **font)
plt.figure()
plt.loglog(rtest[1:-1],model1.rho(rtest[1:-1]),'SpringGreen',linewidth = 5.,label = 'n = {0}'.format(1))
plt.loglog(rtest[1:-1],model2.rho(rtest[1:-1]),'Gold',linewidth = 5.,label = 'n = {0}'.format(2))
plt.loglog(rtest[1:-1],model3.rho(rtest[1:-1]),'DeepSkyBlue',linewidth = 5.,label = 'n  = {0}'.format(4))
#plt.axvline(1)
plt.ylim(1e-50,1e40)
plt.xlabel('radius (normalized to effective radius)')
plt.ylabel('density')
plt.legend(loc = 'best')
plt.title('SERSIC MODELS')
