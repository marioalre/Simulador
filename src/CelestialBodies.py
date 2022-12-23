from src.Orbit import Orbit

class CelestialBodies(Orbit):
	def __init__(self):
		self.mu = 0
		self.radius = 0

	def __str__(self):
		return self.mu, self.radius
			

	def earth(self):
		'''Earth parameters'''
		self.mu = 398600.4418  # km^3/s^2
		self.radius = 6378.137 # km

	def sun(self):
		'''Sun parameters'''
		self.mu = 1.32712440018e20 # km^3/s^2
		self.radius = 695700   # km

	def mars(self):
		'''Mars parameters'''
		self.mu = 42828.3	   # km^3/s^2
		self.radius = 3396.19  # km

	def jupiter(self):
		'''Jupiter parameters'''
		self.mu = 126712764.5 # km^3/s^2
		self.radius = 71492   # km

	def saturn(self):
		'''Saturn parameters'''
		self.mu = 37931207.8  # km^3/s^2
		self.radius = 60268   # km

	def uranus(self):
		'''Uranus parameters'''
		self.mu = 5793939.6   # km^3/s^2
		self.radius = 25559   # km

	def neptune(self):
		'''Neptune parameters'''
		self.mu = 6836527.4   # km^3/s^2
		self.radius = 24764   # km

	def pluto(self):
		'''Pluto parameters'''
		self.mu = 977         # km^3/s^2
		self.radius = 1151    # km
	
	def moon(self):
		'''Moon parameters'''
		self.mu = 4902.801    # km^3/s^2
		self.radius = 1737.4  # km

	def mercury(self):
		'''Mercury parameters'''
		self.mu = 22032		  # km^3/s^2
		self.radius = 2439.7  # km

	def venus(self):
		'''Venus parameters'''
		self.mu = 324858      # km^3/s^2
		self.radius = 6051.8  # km

	def call(self):	
		if self.name == 'earth':
			self.earth()
		elif self.name == 'sun':
			self.sun()
		elif self.name == 'mars':
			self.mars()
		elif self.name == 'jupiter':
			self.jupiter()
		elif self.name == 'saturn':
			self.saturn()
		elif self.name == 'uranus':
			self.uranus()
		elif self.name == 'neptune':
			self.neptune()
		elif self.name == 'pluto':
			self.pluto()
		elif self.name == 'moon':
			self.moon()
		elif self.name == 'mercury':
			self.mercury()
		elif self.name == 'venus':
			self.venus()
		else:
			print('Error: Invalid celestial body')

	




		