import arcade
import functions

import math


class netted_debris():
	def __init__(self, projectile, debris):
		self.projectile = projectile
		self.debris = debris

		new_vel_vector = debris.vel_vector

		new_vel_vector[0] += projectile.vel_vector[0]
		new_vel_vector[1] += projectile.vel_vector[1]




		projectile.vel_vector = new_vel_vector
		debris.vel_vector = new_vel_vector

	def update(self, delta_time):
		self.projectile.update(delta_time)
		self.debris.update(delta_time)

	def kill(self):
		self.debris.kill()
		self.projectile.kill()



class debris(arcade.Sprite):

	def __init__(self, x, y, vel_vector,path):
		super().__init__(filename=path,scale = 0.2)
		self.center_x = x
		self.center_y = y
		self.vel_vector = vel_vector

	def update(self, delta_time):
		self.center_x += self.vel_vector[0]*delta_time
		self.center_y += self.vel_vector[1]*delta_time


class projectile(arcade.Sprite):
	def __init__(self, scale, x, y, vel_vector, filename = "Images/net.png"):
		super().__init__(filename, scale)

		self.center_x = x
		self.center_y = y
		self.vel_vector = vel_vector

		self.angle = functions.vec_to_angle_2d(vel_vector[0],vel_vector[1]) 


	def update(self, delta_time):
		self.center_x += self.vel_vector[0]*delta_time
		self.center_y += self.vel_vector[1]*delta_time






class satellite(arcade.Sprite):
	def __init__(self, filename, scale):
		super().__init__(filename,scale)

		self.center_x = 300
		self.center_y = 300
		self.angle = 0


		
		#Setup values needed for updating the satelllite
		self.has_objective = False
		self.angle_goal = 0
		self.time_to_shoot = float("inf")


	def update_angle(self, debris):

		real_x = debris.center_x - self.center_x
		real_y = debris.center_y - self.center_y

		self.angle = functions.vec_to_angle_2d(real_x,real_y)

	def rotate(self, delta_time, direction):
		self.angle += delta_time * direction*50


	def update(self, delta_time):
		if(self.has_objective):
			#Move the satellite towards the right angle
			if self.angle != self.angle_goal:
				if self.angle_goal > self.angle:
					self.rotate(delta_time,1)
				else:
					self.rotate(delta_time,-1)

			#Update the time to shoot:
			self.time_to_shoot -= delta_time




	def has_objective(self):
		return self.has_objective

	def give_objective(self, debris):
		self.has_objective = True

		#Insert calculation for finding the angle the satellite should have
		self.angle_goal = 42
		
		#Insert calculation for finding the time to shoot
		self.time_to_shoot = 42


	def get_projectile(self):
		shot = projectile(0.5,self.center_x, self.center_y, functions.angle_to_vec_2d(self.angle), "Images/net.png")

		return shot