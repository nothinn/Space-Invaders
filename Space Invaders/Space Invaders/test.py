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

import matplotlib.pylab as plt2



import seaborn as sns

from matplotlib import patches

import csv


def time_range(period, samples):
	i = 0
	while i < period:
		yield i
		i += period/samples

def main():
	#Create an orbit object, which a debris is.
	r = [-6045, -3490, 0]
	v = [-3457, 6618, 0]
	rtest, vtest = functions.get_random_ellipse_orbit()
	ss_orig = Orbit.from_vectors(Earth, r * u.km, v * u.m / u.s)

	print(optimal_impact(ss_orig))

	heatmap_orbit(ss_orig,100,30)



def min_distance(ss_orig):
		#We find the angle to the debris

		pos = functions.orbit_to_position(ss_orig)
		first_angle = functions.vec_to_angle_2d(pos[0],pos[1])

		#Reset best periapsis:
		best_periapsis = 10000000 *u.km

		#We try different angles but always with the same magnitude
		#The angle is relative to the orbit.
		#This means that 0 degrees is a head on collision
		#90 degrees is a collision so that it will line up towards earth
		#180 degrees is a collision from behind, which is hard to accomplish.




		#We apply gradient descent to this.
		#for angle in range(-180,180,1):

		#	orbit_vector = functions.get_vector_orbit(ss_orig)
		#	new_angle = functions.vec_to_angle_2d(orbit_vector[0].value,orbit_vector[1].value)
		#	new_angle = new_angle -180 + angle

		#	#The angle is normalized, meaning it will have the same length no matter what.
		#	x,y = functions.angle_to_vec_2d(new_angle)


		#	ss = functions.orbit_impulse(ss_orig,[x*10,y*10,0])

		#	#The periapsis tells us the closest point to earth.
		#	periapsis_radius = ss.r_p

		#	if periapsis_radius < best_periapsis:
		#		best_periapsis = periapsis_radius
		#		best_x = x
		#		best_y = y
		#		best_angle = angle
		#	#If this is within the karman line, the debris will be destroyed
		#	#print(periapsis_radius - earth_radius* u.km)


		current_periapsis = 10000000*u.km
		last_periapsis = 1*u.km
		angle = -179
		change = 40
		swapped_last = False
		direction = 1
		#While the change is larger than 1 meter

		iterations = 0

		while abs(last_periapsis - current_periapsis) > 1*u.m:
			iterations += 1
			angle += change * direction

			last_periapsis = current_periapsis
			orbit_vector = functions.get_vector_orbit(ss_orig)
			new_angle = functions.vec_to_angle_2d(orbit_vector[0].value,orbit_vector[1].value)
			new_angle = new_angle -180 + angle

			#The angle is normalized, meaning it will have the same length no matter what.
			x,y = functions.angle_to_vec_2d(new_angle)
			ss = functions.orbit_impulse(ss_orig,[x*10,y*10,0])

			#The periapsis tells us the closest point to earth.
			current_periapsis = ss.r_p

			#The bad case
			if current_periapsis > last_periapsis:
				#We move in the other direction
				direction *= -1

				#If we did not change direction last time, we half the distance.
				if not swapped_last:
					change /= 2
				last_periapsis = 0*u.km
				swapped_last = True
			else:
				swapped_last = False

		print("Iterations:")
		print(iterations)
		return current_periapsis, angle






def optimal_impact(ss_orig):
	propagated_time = 0*u.h

	best_angle = 0
	best_periapsis = 100000 *u.km

	superiour_periapsis = best_periapsis
	best_ss = ss_orig

	period = ss_orig.state.period
	
	#original_distance = ss_orig.r_p.value - earth_radius

	#We do a gradient descent on the time.
	#First we test when moving backwards
	#Then test when moving forward.
	#We then see what improves the result and go from that.

	ss_prop_front = ss_orig.propagate(1*u.min)
	ss_prop_back = ss_orig.propagate(-1*u.min)


	best_periapsis_front, best_angle1 = min_distance(ss_prop_front)
	best_periapsis_back, best_angle2 = min_distance(ss_prop_back)

	if best_periapsis_back < best_periapsis_front:
		direction = -1
		last_periapsis = best_periapsis_back
	else:
		direction = 1
		last_periapsis = best_periapsis_front

	current_periapsis = 1000000000*u.km

	current_ss = ss_orig

	

	prop_time = current_ss.state.period / 10

	swapped_last = False

	iterations = 0

	while(abs(last_periapsis - current_periapsis)> 100*u.m):
		iterations += 1
		last_periapsis = current_periapsis
		old_ss = current_ss
		current_ss = current_ss.propagate(prop_time*direction)
		propagated_time += prop_time

		current_periapsis,best_angle = min_distance(current_ss)
		#If we do not improve, we swap direction and half the propagation time
		if current_periapsis > last_periapsis:

			propagated_time -= prop_time
			
			direction *= -1
			if not swapped_last:
				prop_time /= 2

			current_ss = old_ss #current_ss.propagate(prop_time*direction)
			#current_periapsis,best_angle = min_distance(current_ss.propagate(prop_time*direction))
			last_periapsis = 0*u.km
			swapped_last = True
			
		else:
			swapped_last = False



	pos = functions.orbit_to_position(current_ss)
	first_angle = functions.vec_to_angle_2d(pos[0],pos[1])
	angle_of_debris = first_angle

	print(current_periapsis)

	print("Iterations last")
	print(iterations)
	return propagated_time,best_angle, angle_of_debris



def heatmap_orbit(ss_orig,xticks,yticks):
	ss_orig = ss_orig.propagate(0.001*u.s)
	earth_radius = 6371 #km

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
	
	myData = []
	
	top_line = []

	original_distance = ss_orig.r_p.value - earth_radius

	top_line.append(ss_orig.r_p)



	angle_myData = []
	
	for time in time_range(period,yticks):
		row = []


		ss_prop = ss_orig.propagate(time)


		#We find the angle to the debris

		pos = functions.orbit_to_position(ss_prop)
		first_angle = functions.vec_to_angle_2d(pos[0],pos[1])

		angle_myData.append(first_angle)


		#Reset best periapsis:
		best_periapsis = 10000000 *u.km

		#We try different angles but always with the same magnitude
		#The angle is relative to the orbit.
		#This means that 0 degrees is a head on collision
		#90 degrees is a collision so that it will line up towards earth
		#180 degrees is a collision from behind, which is hard to accomplish.
		xlabels = []
		for angle in range(-180,180,int(360/yticks)):
			xlabels.append(angle)

			orbit_vector = functions.get_vector_orbit(ss_prop)

			new_angle = functions.vec_to_angle_2d(orbit_vector[0].value,orbit_vector[1].value)

			new_angle = new_angle -180 + angle

			#The angle is normalized, meaning it will have the same length no matter what.
			x,y = functions.angle_to_vec_2d(new_angle)


			ss = functions.orbit_impulse(ss_prop,[x*10,y*10,0])

			#The periapsis tells us the closest point to earth.
			periapsis_radius = ss.r_p

			#We insert the radius
			row.append(periapsis_radius.value - earth_radius)

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

		#Save the row into the data.
		myData.append(row)


	op2 = OrbitPlotter()



	
	pos = functions.orbit_to_position(best_ss)
	first_angle = functions.vec_to_angle_2d(pos[0],pos[1])
	print(first_angle)
	print(best_ss.r_p.value - earth_radius)

	op2.plot(ss_orig,label = "original")
	op2.plot(best_ss, label = "Best outcome")



	plt.show()


	#Sort myData:
	new_data = [x for _,x in sorted(zip(angle_myData,myData))]

	ylabels = sorted(angle_myData)




	ax = sns.heatmap(new_data, linewidth=0,cmap="RdYlGn_r",center = original_distance,cbar_kws={'label': 'Distance to periapsis minus radius of earth [km]'})

	xticks = ax.get_xticks()

	xlabels2 = ax.get_xticklabels()

	for i, value in enumerate(xticks):
		xlabels2[i] = str(round(float(xlabels[int(value)]), 1))
	ax.set_xticklabels(xlabels2)


	yticks = ax.get_yticks()

	ylabels2 = ax.get_yticklabels()

	for i, value in enumerate(yticks):
		ylabels2[i] = str(round(float(ylabels[int(value)]), 1))
	ax.set_yticklabels(ylabels2)


	ax.set_xlabel("Angle of impact [Degrees]")
	ax.set_ylabel("Angle of debris [Degrees]")


	for tick in ax.get_xticklabels():
		tick.set_rotation(45)


	for tick in ax.get_yticklabels():
		tick.set_rotation(0)


	plt2.show()


	myFile = open('saved_values.csv', 'w')

	with myFile:
		writer = csv.writer(myFile)
		writer.writerows(myData)
