from numpy import *
import math

class SersicModel:
	def __init__(self,model_name,n,Re,M2L,I0):
		self.name = model_name
		self.n = n
		self.Re = Re
		self.M2L = M2L
		self.I0 = I0
		if self.n <10. and self.n > 0.6:
			self.b = 2*self.n - (1./3.) + 0.009876/self.n
			self.p = 1. - 0.6097/self.n + 0.05563/self.n**2
		self.rho0 = M2L*I0*(self.b**(n*(1-self.p)))*(math.gamma(2*n)/(2*Re*math.gamma(n*(3-self.p))))

	def rho(self,r):
		return self.rho0*((r/self.Re)**-p)*exp(-b*((r/self.Re)**(1./self.n)))

	def drhodr(self,r):
		pre = self.rho0*((r/self.Re)**-p)*exp(-b*((r/self.Re)**(1./self.n)))*((self.n*r)**-1)
		post = (self.n*self.p)+self.b*((r/self.Re)**(1./n))
		return pre*post

	def d2rhodr2(self,r):
		pre = self.rho0*((r/self.Re)**-p)*exp(-b*((r/self.Re)**(1./self.n)))*((self.n*r)**-2)
		post = (self.p*(1+self.p)*self.n**2) + self.b*(-1 + self.n + 2*self.n*self.p)*((r/self.Re)**(1./n)) + (b**2)*((r/self.Re)**-p)
		return pre*post
