import arcade
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

class netted_debris():
	def __init__(self, projectile, debris):
		self.projectile = projectile
		self.debris = debris

		new_vel_vector = debris.vel_vector

		new_vel_vector[0] += projectile.vel_vector[0]
		new_vel_vector[1] += projectile.vel_vector[1]


		projectile.vel_vector = new_vel_vector
		debris.vel_vector = new_vel_vector

	def update(self, delta_time, canvas_info):
		self.projectile.update(delta_time, canvas_info)
		self.debris.update(delta_time, canvas_info)

	def kill(self):
		self.debris.kill()
		self.projectile.kill()



class debris(arcade.Sprite):

	def __init__(self, x, y, vel_vector, path, canvas_info, filename = "Images/debris.png"):
		super().__init__(filename=path,scale = 0.2)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, canvas_info) #center is the postion on canvas
		self.vel_vector = vel_vector


		#Create an orbit object, which a debris is.
		r = [-6045, -3490, 2500]
		v = [-3457, 6618, 2533]
		rtest, vtest = functions.get_random_ellipse_orbit()
		self.orbit = Orbit.from_vectors(Earth, rtest * u.km, vtest * u.m / u.s)


	def set_vel_vector(self):
		vel_from_orbit = self.orbit.state.v.to(u.m / u.s)
		self.vel_vector[0] = vel_from_orbit[0].value
		self.vel_vector[1] = vel_from_orbit[1].value


	def give_impulse(self):
		
		print(functions.get_vector_orbit(self.orbit))
		self.orbit = functions.orbit_impulse(self.orbit,[100,0,0])


	def update(self, delta_time, canvas_info):
		#self.model_x += self.vel_vector[0]*delta_time
		#self.model_y += self.vel_vector[1]*delta_time


		#We move the debris a certain time.
		self.orbit = self.orbit.propagate(delta_time * u.s)
		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)
		self.set_vel_vector()
		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)






class projectile(arcade.Sprite):
	def __init__(self, scale, x, y, vel_vector, canvas_info, filename = "Images/net.png"):
		super().__init__(filename, scale)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, canvas_info) #center is the postion on canvas
		self.vel_vector = vel_vector

		self.angle = functions.vec_to_angle_2d(vel_vector[0],vel_vector[1]) 


	def update(self, delta_time, canvas_info):
		self.model_x += self.vel_vector[0]*delta_time
		self.model_y += self.vel_vector[1]*delta_time

		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)


class earth(arcade.Sprite):
	def __init__(self, filename, scale):
		super().__init__(filename,scale)

		self.model_x = 0 
		self.model_y = 0
	
	def update(self, delta_time, canvas_info):
		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)





class satellite(arcade.Sprite):
	def __init__(self, x, y, filename, scale):
		super().__init__(filename,scale)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x = 300 #center is the postion on canvas
		self.center_y = 300
		self.vel_vector = [0, 0]
		self.angle = 0
		self.mass = 440 # kg - from new horizons
		self.radius_to_thruster = 1.2 #m - distance from center of mass og new horizons to is rotation thruster
		self.rotate_thrust_force = 0.9 #N - from new horizons orientations thrusters

		self.rotation_info = [False, 0, 0, 0, 0, 0] # [0: Rotation in motion, 1: rotation half time, 2: degrees of rotation, 3: angular accelearion, 4: start time, 5:start angle]
		self.rotation_start = [False, 0] 

		
		#Setup values needed for updating the satelllite
		self.has_objective = False
		self.angle_goal = 0
		self.time_to_shoot = float("inf")


		#Create an orbit object, which a debris is.
		r = [-6045, -3490, 2500]
		v = [-3457, 6618, 2533]
		rtest, vtest = functions.get_random_ellipse_orbit()
		self.orbit = Orbit.from_vectors(Earth, rtest * u.km, vtest * u.m / u.s)
		self.set_vel_vector()

	def set_vel_vector(self):
		vel_from_orbit = self.orbit.state.v.to(u.m / u.s)
		self.vel_vector[0] = vel_from_orbit[0].value
		self.vel_vector[1] = vel_from_orbit[1].value

	def __eq__(self, other):
		# equality metho for comapring satellite instance
		return self.__dict__ == other.__dict__

	def update_angle(self, debris):

		real_x = debris.center_x - self.center_x
		real_y = debris.center_y - self.center_y

		self.angle = functions.vec_to_angle_2d(real_x,real_y)


	def update(self, delta_time, canvas_info, total_time):

		if self.rotation_start[0]:
			self.rotation_info = functions.rotate_satellite(self, self.rotation_start[1], total_time)
			self.rotation_start[0] = False

		if self.rotation_info[0]:
			self.angle, self.rotation_info = functions.update_satellite_rotation(self, total_time)

		
		#We move the debris a certain time.
		self.orbit = self.orbit.propagate(delta_time * u.s)
		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)
		self.set_vel_vector()

		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)
	
	def start_rotation(self, angle_goal):
		self.rotation_start = [True, angle_goal]

	def rotate(self, delta_time, direction):
		self.angle += delta_time * direction*50


		if(self.has_objective):
			#Move the satellite towards the right angle
			if self.angle != self.angle_goal:
				movement = 1

				difference = abs(self.angle - self.angle_goal)

				if difference < 10:
					movement = difference/10



				if self.angle_goal > self.angle:
					self.rotate(delta_time,movement)
				else:
					self.rotate(delta_time,-movement)

			#Update the time to shoot:
			self.time_to_shoot -= delta_time




	def has_objective(self):
		return self.has_objective

	def give_objective(self, debris):
		self.has_objective = True

		wait_time, aim_angle = functions.get_time_to_shoot(1.0, self, debris)

		#Insert calculation for finding the angle the satellite should have
		self.angle_goal = aim_angle
		
		#Insert calculation for finding the time to shoot
		self.time_to_shoot = wait_time


	def get_projectile(self):
		shot = projectile(0.5,self.center_x, self.center_y, functions.angle_to_vec_2d(self.angle), "Images/net.png")

		return shot