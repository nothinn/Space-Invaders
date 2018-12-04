import math
import random

import Classes

import numpy as np

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.maneuver import Maneuver

from poliastro.plotting import plot

from astropy import units as u

import matplotlib.pyplot as plt

def collision(projectile, debris):

	if abs(projectile.model_x - debris.model_x) < 10000 * u.m:
		if abs(projectile.model_y - debris.model_y) < 10000 * u.m:
			return True
	return False

def velocity_change(projectile,debris):

	debris_vel_vector = get_vector_orbit(debris.orbit)


	vel_x = (((projectile.mass*projectile.vel_vector[0])+(debris.mass*debris_vel_vector[0]))/(projectile.mass+debris.mass))-debris_vel_vector[0]
	vel_y = (((projectile.mass*projectile.vel_vector[1])+(debris.mass*debris_vel_vector[1]))/(projectile.mass+debris.mass))-debris_vel_vector[1]

	return [vel_x, vel_y]

def vec_to_angle_2d(x,y):
	if type(x) == type(1*u.km) :
		x = x.value
		y = y.value
	
	radians = math.atan2(y,x)
	degrees = math.degrees(radians)
	
	return degrees

def orbit_to_position(orbit):
	#a = orbit.state.a
	#e = orbit.state.ecc
	#theta = orbit.state.nu / u.rad

	#Equation from https://en.wikipedia.org/wiki/Kepler_orbit#Johannes_Kepler
	#The radius equation
	#upper_part = a*(1-e*e)
	#lower_part = 1+e*math.cos(theta )

	#distance = upper_part/lower_part

	#x,y = angle_to_vec_2d(math.degrees(theta))

	#x *= distance
	#y *= distance

	#print("x:{0:08.2f} y:{1:08.2f}".format(x,y))

	#Multiply by 1000 to get it in meters

	x = orbit.r[0]
	y = orbit.r[1]

	return x, y

def get_vector_orbit(orbit_element):
	return orbit_element.state.v
	
def get_angle_between_two_orbits(orbit1, orbit2):
	x1, y1 = orbit_to_position(orbit1)
	x2, y2 = orbit_to_position(orbit2)

	dist_x = x2 - x1
	dist_y = y2 - y1

	return np.arctan2(dist_y, dist_x).to(u.deg)

	


def orbit_impulse(orbit, vector):
	
	dv = [vector[0].value,vector[1].to(vector[0].unit).value, 0]*vector[0].unit


	#if type(vector) != type([0,0,0]*u.m/u.s):
	#	dv = vector *u.m/u.s

	man = Maneuver.impulse(dv)

	return orbit.apply_maneuver(man)


def angle_to_vec_2d(angle):
	x = math.cos(math.radians(angle))
	y = math.sin(math.radians(angle))
	return x, y


def get_net_angle_immediate(abs_vel_net, satellite, debris):
	pos_x_satalite = satellite.model_x
	pos_y_satalite = satellite.model_y
	vel_x_debris = debris.vel_vector[0]
	vel_y_debris = debris.vel_vector[1]
	pos_x_debris = debris.model_x
	pos_y_debris = debris.model_y

    #First we move the satalite to the center of coordinate system
	pos_x_debris_center = pos_x_debris - pos_x_satalite
	pos_y_debris_center = pos_y_debris - pos_y_satalite
    
    # Now we rotate the system such that the debris is traveling parallel to the y axis
	abs_vel_debris = np.sqrt(vel_x_debris**2 + vel_y_debris**2)
	angle_before_rotate = vec_to_angle_2d(vel_x_debris,vel_y_debris)
	degrees_of_rotation = 0

	if angle_before_rotate >= 0:
	    degrees_of_rotation = 90 - angle_before_rotate
	else:
	    degrees_of_rotation = -1*angle_before_rotate + 90

	degrees_of_rotation *= -1

	pos_x_debris_after_rotation = pos_y_debris_center*math.sin(math.radians(degrees_of_rotation)) + pos_x_debris_center*math.cos(math.radians(degrees_of_rotation))
	pos_y_debris_after_rotation = pos_y_debris_center*math.cos(math.radians(degrees_of_rotation)) - pos_x_debris_center*math.sin(math.radians(degrees_of_rotation))

	#We test that the value is not negative. If assertion fails, it means there is a negative number.
	assert((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2 >= 0)
	assert(((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) >= 0)

    #calculate the two possiple angles to aim for in rotated system
	rotated_aim_1 = np.arctan2(((pos_x_debris_after_rotation**2) * abs_vel_debris + np.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2),
	                          (-abs_vel_debris * pos_y_debris_after_rotation + np.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2))
    
	'''rotated_aim_2 = math.atan2(((pos_x_debris_after_rotation**2) * abs_vel_debris - math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2), 
                               -((abs_vel_debris * pos_y_debris_after_rotation + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2)))
    '''
	collision_time = pos_x_debris_after_rotation / (abs_vel_net * np.cos(rotated_aim_1)) # the time from fireing to the collision in seconds

    # get the angle in the non-roted system
	aim_1 = rotated_aim_1.to(u.deg) + degrees_of_rotation*u.deg
    #aim_2 = math.degrees(rotated_aim_2) + degrees_of_rotation

	return (aim_1.value%360)*u.deg, collision_time.to(u.s) #, aim_2



   #this functions returns the the time to wait before shooting if we want to hit the debris at a 90 degree angle
def get_time_to_shoot(abs_vel_net, satellite, debris):
	pos_x_satalite = satellite.center_x
	pos_y_satalite = satellite.center_y
	vel_x_debris = debris.vel_vector[0]
	vel_y_debris = debris.vel_vector[1]
	pos_x_debris = debris.center_x
	pos_y_debris = debris.center_y

	#First we move the satalite to the center of coordinate system
	pos_x_debris_center = pos_x_debris - pos_x_satalite
	pos_y_debris_center = pos_y_debris - pos_x_satalite
    
    # Now we rotate the system such that the debris is traveling parallel to the y axis
	abs_vel_debris = math.sqrt(vel_x_debris**2 + vel_y_debris**2)
	angle_before_rotate = vec_to_angle_2d(vel_x_debris, vel_y_debris)
	degrees_of_rotation = 0

	if angle_before_rotate >= 0:
		degrees_of_rotation = 90 - angle_before_rotate
	else:
		degrees_of_rotation = -1*angle_before_rotate + 90

	degrees_of_rotation *= -1 # we rotate such that the debris is moving in the posative y direction

	pos_x_debris_after_rotation = pos_y_debris_center*math.sin(math.radians(degrees_of_rotation)) + pos_x_debris_center*math.cos(math.radians(degrees_of_rotation))
	pos_y_debris_after_rotation = pos_y_debris_center*math.cos(math.radians(degrees_of_rotation)) - pos_x_debris_center*math.sin(math.radians(degrees_of_rotation))

	if pos_y_debris_after_rotation > 0:
		return -1 # it is too late to hit the debris
	
	# calculate the time for the net to reach collition point :
	if pos_x_debris_after_rotation>=0:
		rotated_aim_angle = 0 #debris is moving in y direction only on posative site of x, then the 0 degrees is where we should aim to hit with 90 degrees
	else:
		rotated_aim_angle = math.pi #debris is moving in y direction only on posative site of x, then the 180 degrees is where we should aim to hit with 90 degrees

	t_n = pos_x_debris_after_rotation / (abs_vel_net*math.cos(rotated_aim_angle))

	# calculate time for the debris to reach collition point :
	coll_p_y = 0 #collision point in y direction
	
	t_d = (coll_p_y - pos_y_debris_after_rotation) / abs_vel_debris

	#the time that we must wait are the difference between the two times

	t_wait = t_d - t_n

	aim = math.degrees(rotated_aim_angle) + degrees_of_rotation

	if t_wait < 0:
		return -1 # the time to wait is negative which means we do not have time.
	
	return t_wait/100, aim

	

def get_canvas_pos(x, y, canvas_info):
	ref_x = canvas_info[0]
	ref_y = canvas_info[1]
	zoom_mult = canvas_info[2]

	if(type(x) == type(1*u.km)):
		x = x.value * 1000
		y = y.value * 1000

	return (x*zoom_mult - ref_x*zoom_mult + 300), (y*zoom_mult - ref_y*zoom_mult + 300)


def get_random_circular_orbit():
	Er = 6371 # earth radius in km
	r_length = random.uniform(160.0 + Er, 2000.0 + Er) # km
	angle = random.uniform(0.0, math.pi)

	r = [math.cos(angle)*r_length, math.sin(angle)*r_length, 0.0]

	G = 6.67408*10**-11 #Gravitational constant m^3*kg^-1*s^-2
	M_e = 5.9722*10**24 #mass of earth kg 
	cirular_abs_vel = math.sqrt((G*M_e) / (r_length*10**3))

	ran_seed = random.uniform(-1.0, 1.0)

	if(ran_seed < 0):
		vel_angle = angle + (math.pi * -1)/2
	else:
		vel_angle = angle + math.pi/2

	vel = [math.cos(vel_angle) * cirular_abs_vel, math.sin(vel_angle) * cirular_abs_vel, 0]

	return r, vel

def get_random_ellipse_orbit():
	Er = 6371 # earth radius in km
	r_length =  random.uniform(200.0 + Er, 2000.0 + Er) # km
	angle = random.uniform(0.0, math.pi)

	r = [math.cos(angle)*r_length, math.sin(angle)*r_length, 0.0]

	G = 6.67408*10**-11 #Gravitational constant m^3*kg^-1*s^-2
	M_e = 5.9722*10**24 #mass of earth kg 
	cirular_abs_vel = math.sqrt((G*M_e) / (r_length*10**3))


	ran_add_angle = random.uniform(-0.0125, 0.0125)*((r_length-Er)/200) #different amount of ellipticallity is allowed for different radiuses
	ran_seed = random.uniform(-1.0, 1.0)

	if(ran_seed < 0):
		vel_angle = angle + ((math.pi/2) * -1) + ran_add_angle

	else:
		vel_angle = angle + (math.pi/2) + ran_add_angle

	vel = [math.cos(vel_angle) * cirular_abs_vel, math.sin(vel_angle) * cirular_abs_vel, 0]

	return r, vel

def rotate_satellite(satallite, angle_goal, start_time): 
	# This function calculates the parameters for a rotation.
	# The function assumes a initial angular velocity of 0.
	# As there is constant acceleration we can calculate the time to get halfway to the goal, 
	# then deaccelerate for the same time.
	m = satallite.mass #kg
	r = satallite.radius_to_thruster #m
	F = satallite.rotate_thrust_force #N

	degrees_of_rotation = angle_goal%360 - satallite.angle
	if  degrees_of_rotation < 180 and degrees_of_rotation >= 0:
		aa = F/(m*r) #rad/s^2 angular acceleration
	elif degrees_of_rotation > 180:
		degrees_of_rotation = (degrees_of_rotation%180 - 180)
		aa = -(F/(m*r))
	elif degrees_of_rotation < -180:
		degrees_of_rotation = (degrees_of_rotation%180)
		aa = F/(m*r)
	else:
		aa = -(F/(m*r))


	#calculate time to get to half the goal angle
	t = math.sqrt((math.radians(degrees_of_rotation))/aa) #

	return [True, t, degrees_of_rotation, aa, start_time, satallite.angle]

def update_satellite_rotation(satellite, current_time):
	# This function updates the angle of the satillite during rotation
	t = satellite.rotation_info[1]
	degrees_of_rotation = satellite.rotation_info[2]
	aa = satellite.rotation_info[3]
	start_time = satellite.rotation_info[4]
	start_angle = satellite.rotation_info[5]

	if current_time > 2*t+start_time: # If the rotation is over set the angle to the goal and stop the rotation
		return (start_angle + degrees_of_rotation) % 360, [False, 0, 0, 0, 0, 0]

	time_passed = current_time - start_time
	
	if time_passed <= t: # Set angle doing acceleration
		rotated = (aa/2)*(time_passed**2)
		new_angle = math.degrees(rotated) + start_angle
		
	else: # Set angle doing deacceleration
		rotated =  math.radians(degrees_of_rotation)/2 + aa*t*(time_passed-t) + (-aa/2)*((time_passed-t)**2) 
		new_angle = math.degrees(rotated) + start_angle

	new_angle = new_angle % 360

	return new_angle, [True, t, degrees_of_rotation, aa, start_time, start_angle]

def distance_distance_two_objects(ax, ay, bx, by):
	return np.sqrt((bx-ax)**2 + (by-ay)**2)

def orbit_direction(orbit):
	initial_x, initial_y = orbit_to_position(orbit)
	later_x, later_y = orbit_to_position(orbit.propagate(60 * u.s))

	initial_angle = vec_to_angle_2d(initial_x, initial_y)
	later_angle = vec_to_angle_2d(later_x, later_y)

	angle_difference = (later_angle%360) - (initial_angle%360) 

	# we assume that all periods are well above 30min thats means that we can check the boundery condition as below
	if angle_difference >= 0 and angle_difference < 180:
		res = -1
	elif angle_difference < 0 and angle_difference > -180:
		res = 1
	elif angle_difference <= 0:
		res = -1
	else:
		res = 1
	
	return res #-1 counter clockwise, 1: clockwise

def find_crossing_times(satellite, debris_list, seek_time):
	import copy
	#We make a copy to have the same orbit afterwards.
	old_version = satellite.orbit.propagate(0.000001*u.s)

	# start by figuring out of satellite are moving clockwise or counter clockwise
	far_distance = 3000000*u.m
	close_distance = 2000000*u.m
	accept_angle = 15

	satelite_rot_dir = orbit_direction(satellite.orbit)
	result_list = list()
	index_count = -1

	for debris in debris_list: # we are testing for all debris
		#We make a copy of the debris to return to original state.
		old_debris = debris.orbit.propagate(0.000001*u.s)
		if satelite_rot_dir == orbit_direction(debris.orbit): # if the satelite travels in same direction we wont have enough impulse and we wont calculate on this debris.
			result_list.append(False)
			index_count += 1
			continue
		result_list.append(list())
		index_count += 1
		seek_time_unit = seek_time * u.s
		debris_period = debris.orbit.state.period
		satellite_period = satellite.orbit.state.period

		nr_of_crossings = seek_time_unit/debris_period + seek_time_unit/satellite_period
		

		for i in range(1, int(nr_of_crossings.value)*2): #searching at different samples - nyquist style amount of samples
			start_time = (seek_time_unit/(int(nr_of_crossings.value)*2)) * i 
			time_increments = seek_time_unit/(int(nr_of_crossings.value)*2)

			gradiant_time = start_time
			gradiant_diff = time_increments/2

			#First we find slobe of the distance function at sample point
			sat_orbit_copy = satellite.orbit.propagate(gradiant_time)
			deb_orbit_copy = debris.orbit.propagate(gradiant_time)

			sat_x, sat_y = orbit_to_position(sat_orbit_copy)
			deb_x, deb_y = orbit_to_position(deb_orbit_copy)

			dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

			sat_orbit_copy = satellite.orbit.propagate(gradiant_time + 0.1*u.s)
			deb_orbit_copy = debris.orbit.propagate(gradiant_time + 0.1*u.s)

			sat_x, sat_y = orbit_to_position(sat_orbit_copy)
			deb_x, deb_y = orbit_to_position(deb_orbit_copy)

			dist2 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

			slobe = dist2 - dist1 #a negative slope: the we seek a point that are in positive time direction.

			if dist1 < far_distance and slobe < 0*u.m: # see if this is a valid point already
				j = 0
				if dist1 < close_distance: #we are too close to get good angle
					while True:
						j += 1
						sat_orbit_copy = satellite.orbit.propagate(gradiant_time - j*1*u.s)
						deb_orbit_copy = debris.orbit.propagate(gradiant_time - j*1*u.s)

						sat_x, sat_y = orbit_to_position(sat_orbit_copy)
						deb_x, deb_y = orbit_to_position(deb_orbit_copy)

						dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

						if dist1 > close_distance:
							break
				
				satellite_temp = copy.deepcopy(satellite)
				debris_temp = copy.deepcopy(debris)

				satellite_temp.orbit = sat_orbit_copy
				debris_temp.orbit = deb_orbit_copy

				satellite_temp.set_vel_vector()
				debris_temp.set_vel_vector()

				aim_angle, time_to_collision = get_net_angle_immediate(np.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

				debris_temp.orbit = debris_temp.orbit.propagate(time_to_collision)
				debris_temp.set_vel_vector()

				debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

				aim_angle = get_angle_between_two_orbits(satellite_temp.orbit, debris_temp.orbit)

				satellite_temp.angle =  aim_angle.value

				collision_angle = aim_angle.value - debris_angle_flip

				satellite_vec_angle = (np.arctan2(satellite_temp.vel_vector[1], satellite_temp.vel_vector[0]).to(u.deg).value % 360) * u.deg

				aim_dif = aim_angle - satellite_vec_angle
				if aim_dif.value < 0:
					aim_dif = (aim_dif.value%360) * u.deg
				elif aim_dif.value > 180:
					aim_dif = -(aim_dif.value - 360) * u.deg

				if collision_angle < -180:
					collision_angle = -(collision_angle%360)
				
				if collision_angle < accept_angle and collision_angle > -accept_angle and aim_dif < 90 * u.deg and aim_dif > -90 * u.deg: #succes
					#We can now calculate the weight needed to make it hit the Karman line at this position:
					weight = weight_needed(debris_temp,satellite_temp)
					result_list[index_count].append((gradiant_time - j*1*u.s , collision_angle, aim_angle, time_to_collision, weight))#SAVE SUCCES
				else:		#fail
					result_list[index_count].append(False)
				
				continue

			#We test edge case for this samples search interval
			if slobe <= 0*u.m: 
				sat_orbit_copy = satellite.orbit.propagate(gradiant_time + gradiant_diff)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time + gradiant_diff)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

				sat_orbit_copy = satellite.orbit.propagate(gradiant_time + gradiant_diff + 0.1*u.s)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time + gradiant_diff + 0.1*u.s)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist2 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

				slobe = dist2 - dist1

				if dist1 < far_distance and slobe < 0*u.m:
					j = 0
					if dist1 < close_distance: #we are too close to get good angle
						while True:
							j += 1
							sat_orbit_copy = satellite.orbit.propagate(gradiant_time + gradiant_diff - j*1*u.s)
							deb_orbit_copy = debris.orbit.propagate(gradiant_time + gradiant_diff - j*1*u.s)

							sat_x, sat_y = orbit_to_position(sat_orbit_copy)
							deb_x, deb_y = orbit_to_position(deb_orbit_copy)

							dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

							if dist1 > close_distance:
								break
					satellite_temp = copy.deepcopy(satellite)
					debris_temp = copy.deepcopy(debris)

					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(np.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit = debris_temp.orbit.propagate(time_to_collision)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					aim_angle = get_angle_between_two_orbits(satellite_temp.orbit, debris_temp.orbit)

					satellite_temp.angle =  aim_angle.value

					collision_angle = aim_angle.value - debris_angle_flip

					satellite_vec_angle = (np.arctan2(satellite_temp.vel_vector[1], satellite_temp.vel_vector[0]).to(u.deg).value % 360) * u.deg

					aim_dif = aim_angle - satellite_vec_angle
					if aim_dif.value < 0:
						aim_dif = (aim_dif.value%360) * u.deg
					elif aim_dif.value > 180:
						aim_dif = -(aim_dif.value - 360) * u.deg

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < accept_angle and collision_angle > -accept_angle and aim_dif < 90*u.deg and aim_dif > -90 * u.deg:#SAVE SUCCES
						weight = weight_needed(debris_temp,satellite_temp)
						result_list[index_count].append((gradiant_time + gradiant_diff - j*1*u.s, collision_angle, aim_angle, time_to_collision,weight))
					else: # fail
						result_list[index_count].append(False)
					continue

				if slobe <= 0*u.m: #failure
					continue
				else: # sets intial parameters for gradient decent - binary search style
					#Gradiant decent in posative direction
					gradiant_diff *= 0.5
					gradiant_time += gradiant_diff
			else: # test the other edge case
				sat_orbit_copy = satellite.orbit.propagate(gradiant_time - gradiant_diff)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time - gradiant_diff)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

				sat_orbit_copy = satellite.orbit.propagate(gradiant_time - gradiant_diff + 0.1*u.s)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time - gradiant_diff + 0.1*u.s)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist2 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

				slobe = dist2 - dist1

				if dist1 < far_distance and slobe < 0*u.m:
					j = 0
					if dist1 < close_distance: #we are too close to get good angle
						while True:
							j += 1
							sat_orbit_copy = satellite.orbit.propagate(gradiant_time - gradiant_diff - j*1*u.s)
							deb_orbit_copy = debris.orbit.propagate(gradiant_time - gradiant_diff - j*1*u.s)

							sat_x, sat_y = orbit_to_position(sat_orbit_copy)
							deb_x, deb_y = orbit_to_position(deb_orbit_copy)

							dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

							if dist1 > close_distance:
								break
					satellite_temp = copy.deepcopy(satellite)
					debris_temp = copy.deepcopy(debris)

					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(np.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit = debris_temp.orbit.propagate(time_to_collision)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					aim_angle = get_angle_between_two_orbits(satellite_temp.orbit, debris_temp.orbit)

					satellite_temp.angle =  aim_angle.value

					collision_angle = aim_angle.value - debris_angle_flip

					satellite_vec_angle = (np.arctan2(satellite_temp.vel_vector[1], satellite_temp.vel_vector[0]).to(u.deg).value % 360) * u.deg

					aim_dif = aim_angle - satellite_vec_angle
					if aim_dif.value < 0:
						aim_dif = (aim_dif.value%360) * u.deg
					elif aim_dif.value > 180:
						aim_dif = -(aim_dif.value - 360) * u.deg

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < accept_angle and collision_angle > -accept_angle and aim_dif < 90 * u.deg and aim_dif > -90 * u.deg: #save sucess
						weight = weight_needed(debris_temp,satellite_temp)
						result_list[index_count].append((gradiant_time - gradiant_diff - j*1*u.s, collision_angle, aim_angle, time_to_collision,weight))						
					else: # fail
						result_list[index_count].append(False)
					continue

				if slobe > 0*u.m: #failure
					continue
				else:# sets intial parameters for gradient decent - binary search style
					#Gradiant decent in negative direction
					gradiant_diff *= 0.5
					gradiant_time -= gradiant_diff

			
			while True: #Binary search style
				sat_orbit_copy = satellite.orbit.propagate(gradiant_time)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

				sat_orbit_copy = satellite.orbit.propagate(gradiant_time + 0.1*u.s)
				deb_orbit_copy = debris.orbit.propagate(gradiant_time + 0.1*u.s)

				sat_x, sat_y = orbit_to_position(sat_orbit_copy)
				deb_x, deb_y = orbit_to_position(deb_orbit_copy)

				dist2 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)
				
				slobe = dist2 - dist1

				if dist1 < far_distance and slobe < 0*u.m:
					j = 0
					if dist1 < close_distance: #we are too close to get good angle
						while True:
							j += 1
							sat_orbit_copy = satellite.orbit.propagate(gradiant_time - j*1*u.s)
							deb_orbit_copy = debris.orbit.propagate(gradiant_time - j*1*u.s)

							sat_x, sat_y = orbit_to_position(sat_orbit_copy)
							deb_x, deb_y = orbit_to_position(deb_orbit_copy)

							dist1 = distance_distance_two_objects(sat_x, sat_y, deb_x, deb_y)

							if dist1 > close_distance:
								break
					satellite_temp = copy.deepcopy(satellite)
					debris_temp = copy.deepcopy(debris)
					
					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(np.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit = debris_temp.orbit.propagate(time_to_collision)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					aim_angle = get_angle_between_two_orbits(satellite_temp.orbit, debris_temp.orbit)

					satellite_temp.angle =  aim_angle.value

					collision_angle = aim_angle.value - debris_angle_flip

					satellite_vec_angle = (np.arctan2(satellite_temp.vel_vector[1], satellite_temp.vel_vector[0]).to(u.deg).value % 360) * u.deg

					aim_dif = aim_angle - satellite_vec_angle
					if aim_dif.value < 0:
						aim_dif = (aim_dif.value%360) * u.deg
					elif aim_dif.value > 180:
						aim_dif = -(aim_dif.value - 360) * u.deg

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < accept_angle and collision_angle > -accept_angle and aim_dif < 90 * u.deg and aim_dif > -90 * u.deg: #SAVE SUCCES
						weight = weight_needed(debris_temp,satellite_temp)
						result_list[index_count].append((gradiant_time - j*1*u.s, collision_angle, aim_angle, time_to_collision, weight))
					
						break
					else: #failure
						result_list[index_count].append(False)
						break
				
				if slobe <= 0*u.m: # we move in positive time direction
					gradiant_diff *= 0.5
					gradiant_time += gradiant_diff

				else: # we move in negative time direction
					gradiant_diff *= 0.5
					gradiant_time -= gradiant_diff

		debris.orbit = old_debris

	
	# The result is a list of list the outer list has a entry for every piece of debris. It has the value False, if this debris orbits in the same direction as the satellite.
	# If the debris orbits in the oppeset direction this index is a list of sampling point where the code has tried to find an point in the orbit where i can hit the debris.
	# this outer list has an entry for every sampling point. An entry is set to False if no point could be found in search space. If a point is found this entry will have the
	# value (delta time, collision angle, aim angle, time_to_collision) -> where the delta time, is the time from the search started to the satellite should shoot. collision angle is the angle
	# collison where 0 degrees is head on. the aim angle is the angle that the satellite should shoot the debris. Time_to_collision is how long it takes from shooting to hitting.

	satellite.orbit = old_version
	return result_list

def get_first_shoot(search_list):
	for debris_info in search_list:
		if debris_info == False:
			continue
		else:
			for entry in debris_info:
				if entry == False:
					continue
				else:
					return entry
	return False

def get_lightest_shot(search_list):
	lowest_weight = 10000*u.kg
	best = False
	for debris_info in search_list:
		if debris_info != False:
			for entry in debris_info:
				if entry != False:
					#If the fourth entry weighs less than the best, we take that instead.
					if entry[4] < lowest_weight:
						lowest_weight = entry[4]
						best = entry
	return best

def print_best_shots(search_list):
	for debris_info in search_list:
		if debris_info == False:
			print("No hit")
		else:
			for entry in debris_info:
				if entry != False:
					print(entry)

			
	
def weight_needed(debris,satellite):

	#Goal is the karman line and 5 km within for full impact
	goal = 6371 *u.km+ 100*u.km - 5*u.km

	#We start with a 10 gram projectile
	weight = 10*u.g

	#We use the mass of the debris as a rule for how much to change the weight of the projectile.
	weight_change = debris.mass / 100
	
	#The first distance is the original periapsis of the debris
	distance = debris.orbit.r_p

	iterations = 0
	#We continue until we reach within the karman line.
	while(True):
		if iterations > 100:
			weight = 100000000*u.kg
			break

		last_distance = distance
		iterations += 1

		projectile = satellite.get_projectile(weight)

		vel_vec = velocity_change(projectile,debris)
	
		new_orbit = orbit_impulse(debris.orbit,vel_vec)

		distance = new_orbit.r_p

		#We do a sort of binary search
		#When too far from earth, we add weight
		if distance > goal:
			weight += weight_change

		#We stop when we hit within 5 km of the karman line.
		elif distance > goal - 15*u.km:
			break
		#When it is too close to earth, we remove weight instead.
		else:
			weight -= weight_change
			weight_change /= 2

	print(iterations)
	return weight


from numpy.random import rand

def print_debris(debris_list):

	print()
	print("Debris | Eccentricity |   Period   |   Weight  | periapsis |")
	for number, debris in  enumerate(debris_list):
		print("{:6} |    {:#9.2} | {:6.5} | {:6.4} | {:6.5} |".format(
			number, 
			debris.orbit.ecc.value,
			debris.orbit.period.to(u.min),
			debris.mass,
			debris.orbit.r_p - 6371*u.km
			))

def plot_result(crossing_times):
	
	plt.ion()
	fig, ax = plt.subplots()

	
	labels = np.zeros((len(crossing_times),1))
	for count, debris in enumerate( crossing_times):
		#What to write in the label
		labels[count] =str(count)


	for count, debris in enumerate( crossing_times):
		if debris != False:
			x = []
			y = []
			for collision in debris:
				if collision != False:
					x.append(collision[0].value) # Time to shoot
					y.append(collision[4].value) # Weight needed


			if len(x) > 0:
				print(x,y)
				ax.scatter(np.asarray(x),np.asarray(y), label = labels[count])

	ax.legend()
	ax.set_xlabel("Time to shoot [seconds]")
	ax.set_ylabel("Weight needed [g]")

	ax.grid(True)

	plt.pause(0.001)

	
def take_time(shoot_info):
	return shoot_info[0]

def get_ordered_target_list(search_list):
	light_list = list()
	for debris_info in search_list: #This for loop makes a list of lightest shots at every debris.
		if debris_info != False:
			lowest_weight = 1000 * u.kg
			for entry in debris_info:
				if entry != False:
					#If the fourth entry weighs less than the best, we take that instead.
					if entry[4] < lowest_weight:
						lowest_weight = entry[4]
						best = entry
			if lowest_weight != 1000 * u.kg:
				light_list.append(best)

	light_list.sort(key = take_time) #Sorts the list on the time to shoot

	return light_list

def plot_percentage(tot_debris, ordered_list):
	x_axis = list()
	y_axis = list()
	count = 0
	tot_mass = 0 * u.kg

	x_axis.append(0)
	y_axis.append(0)
	
	for entry in ordered_list:
		count += 1
		x_axis.append(entry[0].to(u.hour).value)
		y_axis.append((count/tot_debris) * 100)
		tot_mass += entry[4]

	plt.plot(x_axis, y_axis)
	plt.plot(x_axis, y_axis, 'ro')
	plt.axis([0, x_axis[-1], 0, 100])
	plt.xlabel("Time [Hours]")
	plt.ylabel("Debris Removed [%]")
	plt.grid(True)

	plt.show()

	print("average mass: {}".format(tot_mass/count) )

