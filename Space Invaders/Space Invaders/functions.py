import math

def collision(projectile, debris):

	if abs(projectile.center_x - debris.center_x) < 5:
		if abs(projectile.center_y - debris.center_y) < 5:
			return True

	return False

def vec_to_angle_2d(x,y):
	
	radians = math.atan2(y,x)
	degrees = math.degrees(radians)
	
	return degrees

def angle_to_vec_2d(angle):
	x = math.cos(math.radians(angle))
	y = math.sin(math.radians(angle))
	return x, y

def get_net_angle_immediate(abs_vel_net, pos_x_satalite, pos_y_satalite, vel_x_debris, vel_y_debris, pos_x_debris, pos_y_debris):
    #First we move the satalite to the center of coordinate system
    pos_x_debris_center = pos_x_debris - pos_x_satalite
    pos_y_debris_center = pos_y_debris - pos_x_satalite
    
    # Now we rotate the system such that the debris is traveling parallel to the y axis
    abs_vel_debris = math.sqrt(vel_x_debris**2 + vel_y_debris**2)
    angle_before_rotate = vec_to_angle_2d(vel_x_debris,vel_y_debris)
    degrees_of_rotation = 0

    if angle_before_rotate >= 0:
        degrees_of_rotation = 90 - angle_before_rotate
    else:
        abs_vel_debris *= -1 # the debris will be moving the negative y direction
        degrees_of_rotation = -1*(90 + angle_before_rotate)

    pos_x_debris_after_rotation = pos_y_debris_center*math.sin(math.radians(degrees_of_rotation)) + pos_x_debris_center*math.cos(math.radians(degrees_of_rotation))
    pos_y_debris_after_rotation = pos_y_debris_center*math.cos(math.radians(degrees_of_rotation)) - pos_x_debris_center*math.sin(math.radians(degrees_of_rotation))

    #calculate the two possiple angles to aim for in rotated system
    rotated_aim_1 = math.atan2(((pos_x_debris_after_rotation**2) * abs_vel_debris + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2),
                              (-abs_vel_debris * pos_y_debris_after_rotation + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2))
    
    rotated_aim_2 = math.atan2(((pos_x_debris_after_rotation**2) * abs_vel_debris - math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2) * pos_y_debris_after_rotation)/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2), 
                               -((abs_vel_debris * pos_y_debris_after_rotation + math.sqrt((pos_x_debris_after_rotation**2+pos_y_debris_after_rotation**2) * abs_vel_net**2 - (pos_x_debris_after_rotation**2) * abs_vel_debris**2))*pos_x_debris_after_rotation/(pos_x_debris_after_rotation**2 + pos_y_debris_after_rotation**2)))
    
    # get the angle in the non-roted system
    aim_1 = math.degrees(rotated_aim_1) - degrees_of_rotation
    aim_2 = math.degrees(rotated_aim_2) - degrees_of_rotation

    return aim_1, aim_2