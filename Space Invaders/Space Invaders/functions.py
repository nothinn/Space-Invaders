import math

def vec_to_angle_2d(x,y):


	
	radians = math.atan2(y,x)
	degrees = math.degrees(radians)
	
	return degrees

def angle_to_vec_2d(angle):
	x = math.cos(math.radians(angle))
	y = math.sin(math.radians(angle))
	return x, y
