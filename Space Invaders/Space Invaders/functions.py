import math


import numpy as np

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

from poliastro.maneuver import Maneuver

from poliastro.plotting import plot


# 2000kms is the distance of the further debris 
# since our satelite is in the middle of the screen [300,300]
# 2000km = 300 pixels
scaling_factor = 3 / 20 # 0.15

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
	a = orbit.state.a
	e = orbit.state.ecc
	theta = orbit.state.nu / u.rad



	#Equation from https://en.wikipedia.org/wiki/Kepler_orbit#Johannes_Kepler
	#The radius equation
	upper_part = a*(1-e*e)
	lower_part = 1+e*math.cos(theta )

	distance = upper_part/lower_part

	x,y = angle_to_vec_2d(math.degrees(theta))

	x *= distance
	y *= distance

	
	#print("x:{0:08.2f} y:{1:08.2f}".format(x,y))

	return x.value,y.value

def orbit_impulse(orbit, vector):
	

	dv = vector * u.m / u.s
	man = Maneuver.impulse(dv)

	orbit = orbit.apply_maneuver(man)


def angle_to_vec_2d(angle):
	x = math.cos(math.radians(angle))
	y = math.sin(math.radians(angle))
	return x, y

def metric_to_canvas(kms):
	# currently not used anywhere
	return kms*scaling_factor

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