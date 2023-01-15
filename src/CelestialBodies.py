from src.Orbit import Orbit

class CelestialBodies(Orbit):
	'''Celestial bodies parameters'''
	def __init__(self):
		self.mu = 0
		self.radius = 0
		self.J2 = 0

	def __str__(self):
		return self.mu, self.radius
			

	def earth(self):
		'''Earth parameters'''
		self.mu = 398600.4418  # km^3/s^2
		self.radius = 6378.137 # km
		self.J2 = 1.08262668e-3 # Earth J2
		self.f = 1/298.257223563 # Earth flattening

	def sun(self):
		'''Sun parameters'''
		self.mu = 13271244001.800001 # km^3/s^2
		self.radius = 695700   # km

	def mars(self):
		'''Mars parameters'''
		self.mu = 42828.3	   # km^3/s^2
		self.radius = 3396.19  # km
		self.J2 =  1.96045e-3	   # Mars J2

	def jupiter(self):
		'''Jupiter parameters'''
		self.mu = 126712764.5 # km^3/s^2
		self.radius = 71492   # km
		self.J2 = 14.736e-3  # Jupiter J2

	def saturn(self):
		'''Saturn parameters'''
		self.mu = 37931207.8  # km^3/s^2
		self.radius = 60268   # km
		self.J2 = 16.298e-3	   # Saturn J2

	def uranus(self):
		'''Uranus parameters'''
		self.mu = 5793939.6   # km^3/s^2
		self.radius = 25559   # km
		self.J2 =  3.34343e-3		   # Uranus J2

	def neptune(self):
		'''Neptune parameters'''
		self.mu = 6836527.4   # km^3/s^2
		self.radius = 24764   # km
		self.J2 =   3.411e-3		   # Neptune J2

	def pluto(self):
		'''Pluto parameters'''
		self.mu = 977         # km^3/s^2
		self.radius = 1151    # km
	
	def moon(self):
		'''Moon parameters'''
		self.mu = 4902.801    # km^3/s^2
		self.radius = 1737.4  # km
		self.J2 = 202.7e-6 # Moon J2

	def mercury(self):
		'''Mercury parameters'''
		self.mu = 22032		  # km^3/s^2
		self.radius = 2439.7  # km
		self.J2 = 60.0e-6	  # Mercury J2

	def venus(self):
		'''Venus parameters'''
		self.mu = 324858      # km^3/s^2
		self.radius = 6051.8  # km
		self.J2 = 4.458e-6	  # Venus J2

	def personal(self, mu, radius):
		'''Personal parameters'''
		self.mu = mu
		self.radius = radius


	




		