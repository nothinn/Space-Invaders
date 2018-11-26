import arcade
import os
import math
import random
import numpy as np

import functions
import Classes
import test


from astropy import units as u

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600



class MyGame(arcade.Window):
	def __init__(self, width, height):
		super().__init__(width,height)

		self.TIME_MULTIPLIER = 200

		arcade.set_background_color(arcade.color.AMAZON)

		#Set path directory
		self.file_path = os.path.dirname(os.path.abspath(__file__))
		os.chdir(self.file_path)


		#Set lists for sprites
		self.sprites_list = None
		self.debris_list = None
		self.projectile_list = None
		self.netted_debris_list = None
		self.earth = None
		self.satellite = None

		self.right = 0
		self.left = 0
		self.total_time = 0.0

		# for canvas center and zoom factor - relative to model coordinate system
		self.canvas_info = [0, 0, 1/32768, 1]

		self.center_option = False # False satellite is center, True Earth is center

		#self.canvas_center_x = 0 # equal to satalite real in (set update)
		#self.canvas_center_y = 0
		#self.zoom_mult = 1 

		self.display_karman_line = True
		self.display_coordinates = False
		self.update = True
		#Used to check if time was slowed last
		self.slowed = False

		#Used when doing large calculations. The next update will be skipped to avoid time jumping.
		self.skip_update = False


	def setup(self):

		#Prepare lists to contain all sprites
		self.debris_list = arcade.SpriteList()
		self.sprites_list =	arcade.SpriteList()
		self.projectile_list = arcade.SpriteList()
		self.netted_debris_list = list()
		self.total_time = 0.0


		#Generate satellite sprite
		self.satellite = Classes.satellite(0, 705*10**3, "Images/satellite.png",0.2)
		self.sprites_list.append(self.satellite)
	
		#Generate earth
		self.earth = Classes.earth("Images/earth.png",0.82)
		self.sprites_list.append(self.earth)


		#Generate debris sprite
		debris_x = random.uniform(self.satellite.model_x.value - (SCREEN_WIDTH/2), self.satellite.model_x.value + (SCREEN_WIDTH/2))
		debris_y = random.uniform(self.satellite.model_y.value-(SCREEN_WIDTH/2),self.satellite.model_y.value + (SCREEN_WIDTH/2))
		debris = Classes.debris(debris_x, debris_y, [0,0],"Images/debris.png", [0, 705*10**3, 1])
		#debris = Classes.debris(300,500,[0,0], "Images/debris.png")


		self.debris_list.append(debris)
		self.sprites_list.append(debris)

		self.satellite.update_angle(self.debris_list[0])



	def on_key_press(self, symbol, modifiers):
		#Shoot projectile straight from satellite
		if symbol == arcade.key.ENTER:

			vec_x, vec_y = functions.angle_to_vec_2d(self.satellite.angle)
			projectile = self.satellite.get_projectile() #Classes.projectile(0.5,self.satellite.center_x,self.satellite.center_y,[vec_x,vec_y], self.canvas_info)

			self.projectile_list.append(projectile)
			self.sprites_list.append(projectile)
			
		#Make debris and shoot at it immediately
		#elif symbol == arcade.key.M:

		#	debris_vel = 0.15
		#	debris_x = random.uniform(self.canvas_info[0] - 300, self.canvas_info[0] + 300)
		#	debris_y = random.uniform(self.canvas_info[1] - 300, self.canvas_info[1] + 300)
		#	debris_vel_vec = [math.cos(random.uniform(-1*math.pi,math.pi))*debris_vel, math.sin(random.uniform(-1*math.pi,math.pi))*debris_vel]
		#	new_debris = Classes.debris(debris_x, debris_y, debris_vel_vec, "Images/debris.png",  self.canvas_info)
		#	self.debris_list.append(new_debris)
		#	self.sprites_list.append(new_debris)

		#	pro_angle_1 = functions.get_net_angle_immediate(0.4, self.satellite, new_debris)

		#	projetile_vel = 0.4
		#	projetile_vel_vec = [math.cos(math.radians(pro_angle_1))*projetile_vel, math.sin(math.radians(pro_angle_1))*projetile_vel]
		#	new_projectile = Classes.projectile(0.5, self.satellite.model_x, self.satellite.model_y, projetile_vel_vec, self.canvas_info)
		#	self.projectile_list.append(new_projectile)
		#	self.sprites_list.append(new_projectile)


		elif symbol == arcade.key.N:
			possibilities = functions.find_crossing_times(self.satellite,self.debris_list,86400)

			self.skip_update = True

			self.satellite.give_objective(possibilities)

		#Display coordinates of objects
		elif symbol == arcade.key.D:
			if(self.display_coordinates):
				self.display_coordinates = False
			else:
				self.display_coordinates = True

		#Show or hide the karman line
		elif symbol == arcade.key.K:
			if(self.display_karman_line):
				self.display_karman_line = False
			else:
				self.display_karman_line = True

		#zoom in
		elif symbol == arcade.key.I: 
			self.canvas_info[2] *= 2
			self.canvas_info[3] += 1
			self.earth.kill()

			if(self.canvas_info[3] < 0):
				self.earth = Classes.earth("Images/earth.png",0.82 * ( 1 / (-4 * self.canvas_info[3])))
			elif(self.canvas_info[3] == 0):
				self.earth = Classes.earth("Images/earth.png",0.82 * 0.5)
			else:
				self.earth = Classes.earth("Images/earth.png",0.82 * self.canvas_info[3])

			self.sprites_list.append(self.earth)
			
		#zoom out
		elif symbol == arcade.key.O: 
			self.canvas_info[2] *= 0.5
			self.canvas_info[3] -= 1
			self.earth.kill()

			if(self.canvas_info[3] < 0):
				self.earth = Classes.earth("Images/earth.png",0.82 * ( 1 / (-4 * self.canvas_info[3])))
			elif(self.canvas_info[3] == 0):
				self.earth = Classes.earth("Images/earth.png",0.82 * 0.5)
			else:
				self.earth = Classes.earth("Images/earth.png",0.82 * self.canvas_info[3])

			self.sprites_list.append(self.earth)

		#Create random debris
		elif symbol == arcade.key.C:
			new_debris = Classes.debris(0, 0, [0,0], "Images/debris.png", self.canvas_info)
			self.debris_list.append(new_debris)
			self.sprites_list.append(new_debris)

		#Apply impulse to debris 0
		elif symbol == arcade.key.L: 
			debris = self.debris_list[0]
			debris.give_impulse()


		#Pause
		elif symbol == arcade.key.P: 
			if(self.update):
				self.update = False
				#The following was used for following was used for distance measurements
				sx, sy = functions.orbit_to_position(self.satellite.orbit)
				dx, dy = functions.orbit_to_position(self.debris_list[-1].orbit)
				print(functions.distance_distance_two_objects(sx, sy, dx, dy))
			else:
				self.update = True

		#Rotate	satellite to random angle
		elif symbol == arcade.key.R:
			self.satellite.start_rotation(random.uniform(0, 360))

		#Toggle center of screen to be between earth and satellite
		elif symbol == arcade.key.Y:
			self.center_option = not self.center_option

		# Runs the position/time seek function
		elif symbol == arcade.key.T:
			self.skip_update = True
			search_list = functions.find_crossing_times(self.satellite, self.debris_list, 86400)
			
			first_target = functions.get_first_shoot(search_list)
			print(first_target)


		#Speed up the simulation
		elif symbol == arcade.key.NUM_ADD:
			self.TIME_MULTIPLIER += 10

		#Slow down the simulation
		elif symbol == arcade.key.NUM_SUBTRACT:
			self.TIME_MULTIPLIER -= 10


		#Slow test
		elif symbol == arcade.key.Q:
			self.skip_update = True
			import time
			time.sleep(10)


	
	def on_key_release(self, symbol, modifiers):
		if symbol == arcade.key.LEFT:
			self.left = 0
		elif symbol == arcade.key.RIGHT:
			self.right = 0

	def get_max(self, max_value, value):
		if(max_value < value):
			max_value = value
		return max_value

	def get_min(self, min_value, value):
		if(min_value > value):
			min_value = value
		return min_value



	def on_draw(self):
		arcade.start_render()

		# Calculate time
		hours = int(self.total_time) // 3600
		minutes = int(self.total_time) // 60 % 60
		seconds = int(self.total_time) % 60 # using a modulus (remainder)
        # Figure out our output
		time_list = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

		#Write text on the screen in the top left corner
		arcade.draw_text("Number of debris: {}\nNumber of nets: {}\nNumber of netted debris: {}\nTime: {}\nTime multiplier:{}\nTimeToShoot:{}".format(len(self.debris_list),len(self.projectile_list),len(self.netted_debris_list),(time_list),self.TIME_MULTIPLIER,self.satellite.time_to_shoot),
                         10, SCREEN_HEIGHT -10, arcade.color.BLACK, 12, anchor_x="left", anchor_y="top")

		
		#Update scale when zooming
		middle_scale = 75/self.canvas_info[2]
		right_scale = 150/self.canvas_info[2]
		str_left_scale = 0
		str_middle_scale = 0
		str_right_scale = 0

		if(self.canvas_info[2]):
			if(middle_scale >= 1000):
				middle_scale = middle_scale/1000
				right_scale = right_scale/1000
				str_left_scale = "0" + "km"
				str_middle_scale = f"{middle_scale}" + "km"
				str_right_scale = f"{right_scale}" + "km"
			else:
				str_left_scale = "0" + "m"
				str_middle_scale = f"{middle_scale}" + "m"
				str_right_scale = f"{right_scale}" + "m"

		#Draw scale in bottom left corner
		arcade.draw_line(10, 10, 160, 10, arcade.color.BLACK, 1) # line
		arcade.draw_line(10, 10, 10, 16, arcade.color.BLACK, 1) # left dash
		arcade.draw_line(85, 10, 85, 16, arcade.color.BLACK, 1) # middle dash
		arcade.draw_line(160, 10, 160, 16, arcade.color.BLACK, 1) # right dash
		arcade.draw_text(str_left_scale, 8, 18, arcade.color.BLACK, 8)
		arcade.draw_text(str_middle_scale, 76, 18, arcade.color.BLACK, 8)
		arcade.draw_text(str_right_scale, 148, 18, arcade.color.BLACK, 8)



		#Draw all sprites
		self.sprites_list.draw()

		#Draw coordinates for nets & debris
		if(self.display_coordinates):
			for x in self.debris_list:
				if(x != self.satellite):
					arcade.draw_text("    [{}]\n    [{}]".format(round(x.center_x,3), round(x.center_y,3)), x.center_x, x.center_y, arcade.color.BLACK, 8)

		if(self.display_karman_line):
			arcade.draw_circle_outline(self.earth.center_x, self.earth.center_y, (120 + 6371) * 1000 * self.canvas_info[2] , arcade.color.CYAN, 1)
			arcade.draw_circle_outline(self.earth.center_x, self.earth.center_y, (2000 + 6371) * 1000 * self.canvas_info[2] , arcade.color.CYBER_GRAPE, 1)


	def update(self, delta_time):

		if(self.update == True and self.skip_update == False):

			

			if self.slowed:
				self.TIME_MULTIPLIER = self.old_TIME_MULTIPLIER
				self.slowed = False

			#See if we move too fast for achieving the satellites objective
			if(self.satellite.time_to_shoot < delta_time * self.TIME_MULTIPLIER * u.s):
				self.old_TIME_MULTIPLIER = self.TIME_MULTIPLIER
				self.TIME_MULTIPLIER = (self.satellite.time_to_shoot / delta_time*u.s).value
				self.slowed = True

			if(self.satellite.time_to_hit < delta_time * self.TIME_MULTIPLIER * u.s):
				self.old_TIME_MULTIPLIER = 1# Slow way down when nearing something to hit self.TIME_MULTIPLIER
				self.TIME_MULTIPLIER = (self.satellite.time_to_hit / delta_time*u.s).value
				self.satellite.time_to_hit = float("inf")*u.s
				self.slowed = True
			
				

			self.total_time += delta_time * self.TIME_MULTIPLIER
			if self.center_option:
				self.canvas_info[0] = self.satellite.model_x.value * 1000
				self.canvas_info[1] = self.satellite.model_y.value * 1000
			else:
				self.canvas_info[0] = SCREEN_WIDTH/2
				self.canvas_info[1] = SCREEN_HEIGHT/2

			if self.satellite.time_to_shoot <= 0:
				shot = self.satellite.get_projectile()
				self.projectile_list.append(shot)
				self.sprites_list.append(shot)
				self.satellite.time_to_shoot = float("inf")

			self.satellite.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info, self.total_time)
			for member in self.projectile_list:
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)

			for member in self.debris_list:
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)
			for member in self.netted_debris_list:
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)

				#Check if debris is within the karman line plus the radius of earth
				if np.sqrt(member.debris.model_x**2 + member.debris.model_y**2) < 6371*u.km + 120*u.km:
					member.kill()
					print("Killed netted debris")
					self.netted_debris_list.remove(member)
			self.earth.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)

			#Test for collisions:
			for idx, projectile in enumerate(self.projectile_list):
				for idy, debris in enumerate(self.debris_list):
					if functions.collision(projectile, debris):
						print("Collision detected at X:{}, Y:{}".format(projectile.center_x, projectile.center_y))

						self.netted_debris_list.append(Classes.netted_debris(projectile,debris))
						self.projectile_list.remove(projectile)
						self.debris_list.remove(debris)
		
		self.skip_update = False


def main():

	if False:
		test.main()
	else:


		game = MyGame(SCREEN_WIDTH,SCREEN_HEIGHT)
		game.setup()

		arcade.run()


		print("Running")

if __name__ == "__main__":
	main()