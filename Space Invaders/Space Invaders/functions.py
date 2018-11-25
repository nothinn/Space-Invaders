import math
import random

import numpy as np

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.maneuver import Maneuver

from poliastro.plotting import plot

from astropy import units as u

def collision(projectile, debris):

	if abs(projectile.center_x - debris.center_x) < 10:
		if abs(projectile.center_y - debris.center_y) < 10:
			return True
	return False

def velocity_change(projectile,debris):

	vel_x = (((projectile.mass*projectile.vel_vector[0])+(debris.mass*debris.vel_vector[0]))/(projectile.mass*debris.mass))-debris.vel_vector[0]
	vel_y = (((projectile.mass*projectile.vel_vector[1])+(debris.mass*debris.vel_vector[1]))/(projectile.mass*debris.mass))-debris.vel_vector[1]

	return [vel_x, vel_y]

def vec_to_angle_2d(x,y):
	
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

	return x.value * 1000, y.value*1000

def get_vector_orbit(orbit_element):
	return orbit_element.state.v
	

def orbit_impulse(orbit, vector):
	
	dv = vector * u.m / u.s
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
	abs_vel_debris = math.sqrt(vel_x_debris**2 + vel_y_debris**2)
	angle_before_rotate = vec_to_angle_2d(vel_x_debris,vel_y_debris)
	degrees_of_rotation = 0

	if angle_before_rotate >= 0:
	    degrees_of_rotation = 90 - angle_before_rotate
	else:
	    degrees_of_rotation = -1*angle_before_rotate + 90

	degrees_of_rotation *= -1

	pos_x_debris_after_rotation = pos_y_debris_center*math.sin(math.radians(degrees_of_rotation)) + pos_x_debris_center*math.cos(math.radians(degrees_of_rotation))
	pos_y_debris_after_rotation = pos_y_debris_center*math.cos(math.radians(degrees_of_rotation)) - pos_x_debris_center*math.sin(math.radians(degrees_of_rotation))

    #calculate the two possiple angles to aim for in rotated system
	rotated_aim_1 = math.atan2(((pos_x_debris_after_rotation**2) * abs_vel_debris + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2),
	                          (-abs_vel_debris * pos_y_debris_after_rotation + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2))
    
	'''rotated_aim_2 = math.atan2(((pos_x_debris_after_rotation**2) * abs_vel_debris - math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2), 
                               -((abs_vel_debris * pos_y_debris_after_rotation + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2)))
    '''
	collision_time = pos_x_debris_after_rotation / (abs_vel_net * math.cos(rotated_aim_1)) # the time from fireing to the collision in seconds

    # get the angle in the non-roted system
	aim_1 = math.degrees(rotated_aim_1) + degrees_of_rotation
    #aim_2 = math.degrees(rotated_aim_2) + degrees_of_rotation

	return aim_1%360, collision_time #, aim_2



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

	degrees_of_rotation = angle_goal - satallite.angle
	if  degrees_of_rotation < 180 and degrees_of_rotation >= 0:
		aa = F/(m*r) #rad/s^2 angular acceleration
	elif degrees_of_rotation > 180:
		degrees_of_rotation = -(degrees_of_rotation - 180)
		aa = -(F/(m*r))
	elif degrees_of_rotation < -180:
		degrees_of_rotation = -(degrees_of_rotation + 180)
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
	return math.sqrt((bx-ax)**2 + (by-ay)**2)

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
	# start by figuring out of satellite are moving clockwise or counter clockwise
	satelite_rot_dir = orbit_direction(satellite.orbit)
	result_list = list()
	index_count = -1
	for debris in debris_list: # we are testing for all debris
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
		print(nr_of_crossings)


		for i in range(0, int(nr_of_crossings.value)*2): #searching at different samples - nyquist style amount of samples
			start_time = (seek_time_unit/(int(nr_of_crossings.value)*2)) * i 
			time_increments = seek_time_unit/(int(nr_of_crossings.value)*2)

			gradient_decent = True
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

			if dist1 < 2000000 and slobe < 0: # see if this is a valid point already
				satellite_temp = satellite
				debris_temp = debris

				satellite_temp.orbit = sat_orbit_copy
				debris_temp.orbit = deb_orbit_copy

				satellite_temp.set_vel_vector()
				debris_temp.set_vel_vector()

				aim_angle, time_to_collision = get_net_angle_immediate(math.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

				debris_temp.orbit.propagate(time_to_collision * u.s)
				debris_temp.set_vel_vector()

				debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

				collision_angle = aim_angle - debris_angle_flip

				if collision_angle < -180:
					collision_angle = -(collision_angle%360)
				
				if collision_angle < 15 and collision_angle > -15:
					print("succes")
					result_list[index_count].append((gradiant_time, collision_angle, aim_angle))
				else:
					print("failed")
					result_list[index_count].append(False)
				#SAVE SUCCES
				continue

			#We test edge case for this samples search interval
			if slobe <= 0: 
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

				if dist1 < 2000000 and slobe < 0:
					satellite_temp = satellite
					debris_temp = debris

					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(math.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit.propagate(time_to_collision * u.s)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					collision_angle = aim_angle - debris_angle_flip

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < 15 and collision_angle > -15:
						print("succes")
						result_list[index_count].append((gradiant_time + gradiant_diff, collision_angle, aim_angle))
						#SAVE SUCCES
					else:
						print("failed")
						result_list[index_count].append(False)
					continue

				if slobe <= 0:
					print("Not here")
					continue
					#decide what to do with failure
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

				if dist1 < 2000000 and slobe < 0:
					satellite_temp = satellite
					debris_temp = debris

					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(math.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit.propagate(time_to_collision * u.s)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					collision_angle = aim_angle - debris_angle_flip

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < 15 and collision_angle > -15:
						print("succes")
						result_list[index_count].append((gradiant_time - gradiant_diff, collision_angle, aim_angle))
					#SAVE SUCCES
						
					else:
						print("failed")
						result_list[index_count].append(False)
					continue

				if slobe > 0:
					print("Not here")
					continue
					#decide what to do with failure
				else:# sets intial parameters for gradient decent - binary search style
					#Gradiant decent in negative direction
					gradiant_diff *= 0.5
					gradiant_time -= gradiant_diff

			
			while gradient_decent: #Binary search style
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

				if dist1 < 2000000 and slobe < 0:
					satellite_temp = satellite
					debris_temp = debris

					satellite_temp.orbit = sat_orbit_copy
					debris_temp.orbit = deb_orbit_copy

					satellite_temp.set_vel_vector()
					debris_temp.set_vel_vector()

					aim_angle, time_to_collision = get_net_angle_immediate(math.sqrt(satellite_temp.vel_vector[0]**2 + satellite_temp.vel_vector[1]**2), satellite_temp, debris_temp)

					debris_temp.orbit.propagate(time_to_collision * u.s)
					debris_temp.set_vel_vector()

					debris_angle_flip = (vec_to_angle_2d(debris_temp.vel_vector[0], debris_temp.vel_vector[1]) + 180)%360

					collision_angle = aim_angle - debris_angle_flip

					if collision_angle < -180:
						collision_angle = -(collision_angle%360)
				
					if collision_angle < 15 and collision_angle > -15:
						print("succes")
						result_list[index_count].append((gradiant_time, collision_angle, aim_angle))
					#SAVE SUCCES
						break
					else:
						print("failed")
						result_list[index_count].append(False)
						break
				
				if slobe <= 0: # we move in positive time direction
					gradiant_diff *= 0.5
					gradiant_time += gradiant_diff

				else: # we move in negative time direction
					gradiant_diff *= 0.5
					gradiant_time -= gradiant_diff

	print(result_list)
	
