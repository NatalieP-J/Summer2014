from numpy import *

def HeavisideStep(x):
	try:
		l = len(x)
		part1 = zeros(len(x[(x>0)]))+1.
		part2 = zeros(len(x[(x<=0)]))
		return concatenate((part1,part2))
	except TypeError:
		if x>0:
			return 1.
		elif x<=0:
			return 0.

class TerzicGrahamModel:
	def __init__(self,model_name,n,Re,rb,rhob,gamma,alpha,generate):
		self.name = model_name
		self.n = n
		self.Re = Re
		self.rb = rb
		self.g = gamma
		self.generate = generate
		self.rhob = rhob
		if self.n <10. and self.n > 0.6:
			self.b = 2*self.n - (1./3.) + 0.009876/self.n
			self.p = 1. - 0.6097/self.n + 0.05563/self.n**2

	def rho(self,r):
		rhobar = self.rhob*(2**((self.p-self.g)/self.a))*((self.rb/self.Re)**p)*exp(self.b*((2**(1./self.a))*((self.rb/self.Re)))**(1./self.n))
		rho1 = rhobar*((1+(self.rb/r)**self.a)**(self.g/self.a))
		rho2 = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-self.p/self.a))
		rho3 = exp(-b*((r**self.a + self.rb**self.a)/self.Re**self.a)**(1./(self.n*self.a)))
		return rho1*rho2*rho3

	def drhodr(self,r):
		rhobar = self.rhob*(2**((self.p-self.g)/self.a))*((self.rb/self.Re)**p)*exp(self.b*((2**(1./self.a))*((self.rb/self.Re)))**(1./self.n))
		part1a = -self.g*self.rb*(self.rb/r)**(-1+self.a)*rhobar/r**2
		part1b = exp(-b*((r**self.a + self.rb**self.a)/self.Re**self.a)**(1./(self.n*self.a)))
		part1c = (1+(self.rb/r)**self.a)**(-1+(self.g/self.a))
		part1d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-self.p/self.a))
		part1 = part1a*part1b*part1c*part1d
		part2a = -rhobar*p*(r**(-1+self.a))*self.Re**(-self.a)
		part2b = part1b
		part2c = ((1+(self.rb/r)**self.a)**(self.g/self.a))
		part2d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-1+(-self.p/self.a)))
		part2 = part2a*part2b*part2c*part2d
		part3a = (-1./self.n)*b*rhobar*(r**(-1+self.a))*(self.Re**-self.a)
		part3b = part1b
		part3c = part2c
		part3d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-1+(-self.p/self.a)+(1./(self.a*self.n))))
		part3 = part3a*part3b*part3c*part3d
		return part1+part2+part3

	def d2rhodr2(self,r):
		rhobar = self.rhob*(2**((self.p-self.g)/self.a))*((self.rb/self.Re)**p)*exp(self.b*((2**(1./self.a))*((self.rb/self.Re)))**(1./self.n))
		part1a = rhobar*self.a*self.g*(-1+(self.g/self.a))*(self.rb**2)*((self.rb/r)**(-2+2*self.a))*r**-4
		part1b = exp(-b*((r**self.a + self.rb**self.a)/self.Re**self.a)**(1./(self.n*self.a)))
		part1c = (1+(self.rb/r)**self.a)**(-2+(self.g/self.a))
		part1d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-self.p/self.a))
		part1 = part1a*part1b*part1c*part1d
		part2a = (r**-4)*(-1+self.a)*self.g*(self.rb**2)*((self.rb/r)**(-2*self.a))*rhobar
		part2b = part1b
		part2c = (1+(self.rb/r)**self.a)**(-1+(self.g/self.a))
		part2d = part1d
		part2 = part2a*part2b*part2c*part2d
		part3a = (r**-3)*2*self.g*self.rb*((self.rb/r)**(-1+self.a))*rhobar
		part3b = part1b
		part3c = part2c
		part3d = part1d
		part3 = part3a*part3b*part3c*part3d
		part4a = -self.a*self.p*(-1-(self.p/self.a))*(r**(-2+2*self.a))*(self.Re**(-2*self.a))*rhobar
		part4b = part1b
		part4c = (1+(self.rb/r)**self.a)**(self.g/self.a)
		part4d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-2+(-self.p/self.a)))
		part4 = part4a*part4b*part4c*part4d
		part5a = -(-1+self.a)*self.p*(r**(-2+self.a))*(self.Re**(-self.a))*rhobar
		part5b = part1b
		part5c = part4c
		part5d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-1+(-self.p/self.a)))
		part5 = part5a*part5b*part5c*part5d
		part6a = 2*self.g*self.p*(r**(-3+self.a))*self.rb*((self.rb/r)**(-1+self.a))*(self.Re**(-self.a))*rhobar
		part6b = part1b
		part6c = part2c
		part6d = part5d
		part6 = part6a*part6b*part6c*part6d
		part7a = (1./self.n)*self.b*self.p*(r**(-2+2*self.a))*(self.Re**(-2*self.a))*rhobar
		part7b = part1b
		part7c = part4c
		part7d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-2+(-self.p/self.a)+(1./(self.a*self.n))))
		part7 = part7a*part7b*part7c*part7d
		part8a = -(1./self.n)*self.a*(r**(-2+2*self.a))*(self.Re**(-2*self.a))*rhobar*(-1+(1./(self.a*self.n))-(self.p/self.a))
		part8b = part1b
		part8c = part4c
		part8d = part7d
		part8 = part8a*part8b*part8c*part8d
		part9a = -(1./self.n)*(-1+self.a)*self.b*(r**(-2+self.a))*(self.Re**self.a)*rhobar
		part9b = part1b
		part9c = part4c
		part9d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-1+(-self.p/self.a)+(1./(self.a*self.n))))
		part9 = part9a*part9b*part9c*part9d
		part10a = (1./self.n)*2*self.b*self.g*(r**(-3+self.a))*((self.rb/r)**(-1+self.a))*(self.Re**-self.a)*rhobar
		part10b = part1b
		part10c = part2c
		part10d = part9d
		part10 = part10a*part10b*part10c*part10d
		part11a = (self.n**-2)*(self.b**2)*(r**(-2+2*self.a))*(self.Re**(-2*self.a))*rhobar
		part11b = part1b
		part11c = part4c
		part11d = (((r**self.a + self.rb**self.a)/self.Re**self.a)**(-2+(-self.p/self.a)+(2./(self.a*self.n))))
		part11 = part11a*part11b*part11c*part11d
		return part1+part2+part3+part4+part5+part6+part7+part8+part9+part10+part11














