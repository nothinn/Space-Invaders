import functions

import math

import numpy as np

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.examples import iss

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.examples import iss

from poliastro.plotting import plot

from poliastro.plotting import OrbitPlotter

import matplotlib.pyplot as plt


def time_range(period, samples):
	i = 0
	while i < period:
		yield i
		i += period/samples

def main():
	earth_radius = 6371 #km

	#Create an orbit object, which a debris is.
	r = [-6045, -3490, 0]
	v = [-3457, 6618, 0]
	rtest, vtest = functions.get_random_ellipse_orbit()
	ss_orig = Orbit.from_vectors(Earth, r * u.km, v * u.m / u.s)
	ss_orig = ss_orig.propagate(0.001*u.s)

	#ss= iss

	
	best_i = 0
	best_j = 0
	best_angle = 0
	best_periapsis = 100000 *u.km

	superiour_periapsis = best_periapsis
	best_ss = ss_orig

	op = OrbitPlotter()
	
	#Plot original orbit
	#op.plot(ss,label = "Original orbit")
	#plt.show()


	period = ss_orig.state.period
	

	for time in time_range(period,30):
		ss_prop = ss_orig.propagate(time)

		#Reset best periapsis:
		best_periapsis = 10000000 *u.km

		#We try different angles but always with the same magnitude
		#The angle is relative to the orbit.
		#This means that 0 degrees is a head on collision
		#90 degrees is a collision so that it will line up towards earth
		#180 degrees is a collision from behind, which is hard to accomplish.
		for angle in range(-180,180):

			orbit_vector = functions.get_vector_orbit(ss_prop)

			new_angle = functions.vec_to_angle_2d(orbit_vector[0].value,orbit_vector[1].value)

			new_angle = new_angle -180 + angle

			#The angle is normalized, meaning it will have the same length no matter what.
			x,y = functions.angle_to_vec_2d(new_angle)


			ss = functions.orbit_impulse(ss_prop,[x*10,y*10,0])

			#The periapsis tells us the closest point to earth.
			periapsis_radius = ss.r_p

			if periapsis_radius < best_periapsis:
				best_periapsis = periapsis_radius
				best_x = x
				best_y = y
				best_angle = angle
			#If this is within the karman line, the debris will be destroyed
			#print(periapsis_radius - earth_radius* u.km)

		print()
		x = best_x
		y = best_y
		print(x)
		print(y)
		print(best_angle)

		ss = functions.orbit_impulse(ss_prop,[x*10,y*10,0])


		#The periapsis tells us the closest point to earth.
		periapsis_radius = ss.r_p

		if periapsis_radius < superiour_periapsis:
			superiour_periapsis = periapsis_radius
			best_ss = ss


		#If this is within the karman line, the debris will be destroyed
		print(periapsis_radius - earth_radius* u.km)

		#Plot the best result
		label = "angle: {}, length: {:4.4f}".format(best_angle, periapsis_radius - earth_radius* u.km)
		op.plot(ss,label = label)


	op2 = OrbitPlotter()

	op2.plot(ss_orig,label = "original")
	op2.plot(best_ss, label = "Best outcome")
	plt.show()

	True