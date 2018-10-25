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



	def update_angle(self, debris):

		real_x = debris.center_x - self.center_x
		real_y = debris.center_y - self.center_y

		self.angle = functions.vec_to_angle_2d(real_x,real_y)

	def rotate(self, delta_time, direction):
		self.angle += delta_time * direction*50
