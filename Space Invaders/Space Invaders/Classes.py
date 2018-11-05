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

	def update(self, delta_time, ref_x, ref_y):
		self.projectile.update(delta_time, ref_x, ref_y)
		self.debris.update(delta_time, ref_x, ref_y)

	def kill(self):
		self.debris.kill()
		self.projectile.kill()



class debris(arcade.Sprite):

	def __init__(self, x, y, vel_vector,path, ref_x, ref_y, filename = "Images/debris.png"):
		super().__init__(filename=path,scale = 0.2)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, ref_x, ref_y) #center is the postion on canvas
		self.vel_vector = vel_vector

	def update(self, delta_time, ref_x, ref_y):
		self.model_x += self.vel_vector[0]*delta_time
		self.model_y += self.vel_vector[1]*delta_time

		self.center_x, self.center_x = functions.get_canvas_pos(self.model_x, self.model_y, ref_x, ref_y)



class projectile(arcade.Sprite):
	def __init__(self, scale, x, y, vel_vector, ref_x, ref_y, filename = "Images/net.png"):
		super().__init__(filename, scale)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x, self.center_y = functions.get_canvas_pos(x, y, ref_x, ref_y) #center is the postion on canvas
		self.vel_vector = vel_vector

		self.angle = functions.vec_to_angle_2d(vel_vector[0],vel_vector[1]) 


	def update(self, delta_time, ref_x, ref_y):
		self.model_x += self.vel_vector[0]*delta_time
		self.model_y += self.vel_vector[1]*delta_time

		self.center_x, self.center_x = functions.get_canvas_pos(self.model_x, self.model_y, ref_x, ref_y)




class satellite(arcade.Sprite):
	def __init__(self, x, y, filename, scale):
		super().__init__(filename,scale)

		self.model_x = x #model is the postion in the background model
		self.model_y = y
		self.center_x = 300 #center is the postion on canvas
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

		wait_time, aim_angle = functions.get_time_to_shoot(1, self, debris)

		#Insert calculation for finding the angle the satellite should have
		self.angle_goal = aim_angle
		
		#Insert calculation for finding the time to shoot
		self.time_to_shoot = wait_time


	def get_projectile(self):
		shot = projectile(0.5,self.center_x, self.center_y, functions.angle_to_vec_2d(self.angle), "Images/net.png")

		return shot