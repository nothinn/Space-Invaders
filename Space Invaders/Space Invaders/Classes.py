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
		
		self.debris.orbit = functions.orbit_impulse(debris.orbit, functions.velocity_change(projectile, debris))

		print(self.debris.orbit.r_p - 6371*u.km)

	def update(self, delta_time, canvas_info):
		self.debris.update(delta_time, canvas_info)
		self.projectile.update(delta_time, canvas_info, self.debris)

	def kill(self):
		self.debris.kill()
		self.projectile.kill()


class debris(arcade.Sprite):

	def __init__(self, x, y, vel_vector, path, canvas_info,mass = 42*u.kg, filename = "Images/debris.png"):
		super().__init__(filename=path,scale = 0.2)

		self.model_x = x * u.km #model is the postion in the background model
		self.model_y = y * u.km
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, canvas_info) #center is the postion on canvas
		self.vel_vector = vel_vector
		self.speed = 0 
		self.altitude = 0

		#Create an orbit object, which a debris is.
		r = [-6045, -3490, 2500]
		v = [-3457, 6618, 2533]
		rtest, vtest = functions.get_random_ellipse_orbit()
		self.orbit = Orbit.from_vectors(Earth, rtest * u.km, vtest * u.m / u.s)

		self.update(0.00001,[0,0,0,0])
		self.mass = mass

	def set_vel_vector(self):
		vel_from_orbit = self.orbit.state.v.to(u.m / u.s)
		speed_vector = self.orbit.state.v
		self.speed = round(math.sqrt((speed_vector[0].value*1000)**2 + (speed_vector[1].value*1000)**2),3) #m/s
		self.vel_vector[0] = vel_from_orbit[0]
		self.vel_vector[1] = vel_from_orbit[1]

		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)

	def give_impulse(self):
		
		self.orbit = functions.orbit_impulse(self.orbit,[100,0,0])

	def update(self, delta_time, canvas_info):

		self.altitude = round(math.sqrt((self.model_x.value)**2 + (self.model_y.value)**2) - 6371, 3)

		#We move the debris a certain time.
		self.orbit = self.orbit.propagate(delta_time * u.s)
		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)
		self.set_vel_vector()
		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)


class projectile(arcade.Sprite):
	def __init__(self, scale, x, y, vel_vector, canvas_info, mass = 42*u.kg, filename = "Images/net.png"):
		super().__init__(filename, scale)


		if type(x) != type(1*u.m):
			x *= u.km
			y *= u.km
			print("Projectile called without units")

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, canvas_info) #center is the postion on canvas

		self.vel_vector = vel_vector

		if type(self.vel_vector) != type([0,0,0]*u.m/u.s):
			self.vel_vector = self.vel_vector *u.m/u.s

		self.angle = functions.vec_to_angle_2d(self.vel_vector[0].value,self.vel_vector[1].value) 
		self.mass = mass


	def update(self, delta_time, canvas_info, debris = None):
		if debris == None:

			self.model_x += self.vel_vector[0]*delta_time * u.s
			self.model_y += self.vel_vector[1]*delta_time * u.s
		else:
			self.model_x = debris.model_x
			self.model_y = debris.model_y
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
		self.speed = 0
		self.angle = 0
		self.mass = 440 # kg - from new horizons
		self.radius_to_thruster = 1.2 #m - distance from center of mass og new horizons to is rotation thruster
		self.rotate_thrust_force = 0.9 #N - from new horizons orientations thrusters

		self.rotation_info = [False, 0, 0, 0, 0, 0] # [0: Rotation in motion, 1: rotation half time, 2: degrees of rotation, 3: angular accelearion, 4: start time, 5:start angle]
		self.rotation_start = [False, 0] 
		
		#Setup values needed for updating the satelllite
		self.has_objective = False
		self.objective_list = list()
		self.objective_iterator = 0
		self.angle_goal = 0
		self.time_to_shoot = float("inf")*u.s
		self.time_to_hit = float("inf")*u.s


		r = [0, 6371+1100, 0]
		v = [7304.048234, 0, 0]
		rtest, vtest = functions.get_random_ellipse_orbit()
		self.orbit = Orbit.from_vectors(Earth, r * u.km, v * u.m / u.s)
		self.set_vel_vector()

	def set_vel_vector(self):
		vel_from_orbit = self.orbit.state.v.to(u.m / u.s)
		speed_vector = self.orbit.state.v
		self.speed = round(math.sqrt((speed_vector[0].value*1000)**2 + (speed_vector[1].value*1000)**2),3) #m/s
		self.vel_vector[0] = vel_from_orbit[0]
		self.vel_vector[1] = vel_from_orbit[1]
		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)

	def __eq__(self, other):
		# equality metho for comapring satellite instance
		return self.__dict__ == other.__dict__

	def update_angle(self, debris):

		real_x = debris.center_x - self.center_x
		real_y = debris.center_y - self.center_y

		self.angle = functions.vec_to_angle_2d(real_x,real_y)


	def update(self, delta_time, canvas_info, total_time):

		self.canvas_info = canvas_info

		if self.rotation_start[0]:
			self.rotation_info = functions.rotate_satellite(self, self.rotation_start[1], total_time)
			self.rotation_start[0] = False

		if self.rotation_info[0]:
			self.angle, self.rotation_info = functions.update_satellite_rotation(self, total_time)
		
		
		#We move the debris a certain time.
		self.orbit = self.orbit.propagate(delta_time * u.s)
		self.model_x, self.model_y = functions.orbit_to_position(self.orbit)
		self.set_vel_vector()
		#Update the time to shoot:
		self.time_to_shoot -= delta_time*u.s
		self.time_to_hit -= delta_time*u.s


		self.center_x, self.center_y = functions.get_canvas_pos(self.model_x, self.model_y, canvas_info)
	
	def start_rotation(self, angle_goal): #sets infromation to start update on next update
		self.rotation_start = [True, angle_goal]

	def rotate(self, delta_time, direction): #Manuel rotation - Not used any more
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


	def has_objective(self): #are the satellite about to shoot debris
		return self.has_objective

	def give_objective(self, possibilities): # give a debris a obejctive
		self.has_objective = True

		first_to_shoot = functions.get_lightest_shot(possibilities)

		if first_to_shoot == False:
			
			return
		wait_time = first_to_shoot[0]

		aim_angle = first_to_shoot[2]


		self.time_to_hit = first_to_shoot[3] + wait_time - 5*u.s

		self.weight_objective = first_to_shoot[4]

		#Insert calculation for finding the angle the satellite should have
		self.start_rotation(aim_angle.value)
		
		self.time_to_shoot = wait_time

	def give_objectives(self, possibilities): # Give multiples debis as obejctives
		self.has_objective = True

		objective_list = functions.get_ordered_target_list(possibilities)

		if len(objective_list) == 0:
			self.has_objective = False
			return

		self.objective_list = objective_list
		self.give_next_objetive()

	def give_next_objetive(self): # when a debris has been shot the next debris are being given as obejtive
		if self.has_objective:
			if self.objective_iterator >= len(self.objective_list):
				self.has_objective = False
				self.objective_iterator = 0
				self.objective_list.clear()
				self.time_to_shoot = float("inf") * u.s
			else:
				if self.objective_iterator == 0:
					self.weight_objective = self.objective_list[self.objective_iterator][4]
					self.time_to_shoot = self.objective_list[self.objective_iterator][0] + 0*u.s
					self.time_to_hit = self.objective_list[self.objective_iterator][3] + self.objective_list[self.objective_iterator][0] - 5*u.s
					self.start_rotation(self.objective_list[self.objective_iterator][2].value)
					self.objective_iterator += 1
					print(self.time_to_shoot)
					print(self.time_to_hit)
				else:
					old_hit_time = self.objective_list[self.objective_iterator - 1][3] + self.objective_list[self.objective_iterator - 1][0]

					self.weight_objective = self.objective_list[self.objective_iterator][4]
					self.time_to_shoot = self.objective_list[self.objective_iterator][0] - old_hit_time
					self.time_to_hit = self.objective_list[self.objective_iterator][3] + self.objective_list[self.objective_iterator][0] - old_hit_time - 5*u.s
					self.start_rotation(self.objective_list[self.objective_iterator][2].value)
					self.objective_iterator += 1
					print(self.time_to_shoot)
					print(self.time_to_hit)


	def get_projectile(self, mass = None): # shoots projectile

		if mass == None:
			mass = self.weight_objective

		#We change the angle to be in the angle of the satellite
		angle = self.angle
		#We use pythagoras to find the velocity

		vel_vector = functions.get_vector_orbit(self.orbit)
		velocity = np.sqrt(vel_vector[0]**2 + vel_vector[1]**2)

		#We convert the angle to vector:
		new_vector = functions.angle_to_vec_2d(angle)

		#And multiply by the velocity
		new_vector *= velocity


		shot = projectile(0.5,self.model_x, self.model_y,vel_vector = new_vector,mass = mass, canvas_info=self.canvas_info)

		return shot