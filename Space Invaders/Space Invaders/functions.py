import math
import random

import numpy as np

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.maneuver import Maneuver

from poliastro.plotting import plot

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

    # get the angle in the non-roted system
	aim_1 = math.degrees(rotated_aim_1) + degrees_of_rotation
    #aim_2 = math.degrees(rotated_aim_2) + degrees_of_rotation

	return aim_1 #, aim_2



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
	r_length = random.uniform(160.0 + Er, 2000.0 + Er) # km
	angle = random.uniform(0.0, math.pi)

	r = [math.cos(angle)*r_length, math.sin(angle)*r_length, 0.0]

	G = 6.67408*10**-11 #Gravitational constant m^3*kg^-1*s^-2
	M_e = 5.9722*10**24 #mass of earth kg 
	cirular_abs_vel = math.sqrt((G*M_e) / (r_length*10**3))


	ran_add_angle = random.uniform(-0.5, 0.5)
	ran_seed = random.uniform(-1.0, 1.0)

	if(ran_seed < 0):
		vel_angle = angle + ((math.pi/2) * -1) + ran_add_angle

	else:
		vel_angle = angle + (math.pi/2) + ran_add_angle

	vel = [math.cos(vel_angle) * cirular_abs_vel, math.sin(vel_angle) * cirular_abs_vel, 0]

	return r, vel

def rotate_satellite(satallite, angle_goal, start_time):
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
	t = math.sqrt((math.radians(degrees_of_rotation))/aa)

	return [True, t, degrees_of_rotation, aa, start_time, satallite.angle]

def update_satellite_rotation(satellite, current_time):
	t = satellite.rotation_info[1]
	degrees_of_rotation = satellite.rotation_info[2]
	aa = satellite.rotation_info[3]
	start_time = satellite.rotation_info[4]
	start_angle = satellite.rotation_info[5]

	if current_time > 2*t+start_time:
		return (start_angle + degrees_of_rotation) % 360, [False, 0, 0, 0, 0, 0]

	time_passed = current_time - start_time
	
	if time_passed <= t:
		rotated = (aa/2)*(time_passed**2)
		new_angle = math.degrees(rotated) + start_angle
		
	else:
		rotated =  math.radians(degrees_of_rotation)/2 + aa*t*(time_passed-t) + (-aa/2)*((time_passed-t)**2) 
		new_angle = math.degrees(rotated) + start_angle

	new_angle = new_angle % 360

	return new_angle, [True, t, degrees_of_rotation, aa, start_time, start_angle]