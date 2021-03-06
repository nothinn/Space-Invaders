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
		self.canvas_info = [0, 0, 1/32768]

		self.center_option = False # False satellite is center, True Earth is center

		self.display_karman_line = True
		self.display_coordinates = False
		self.update = True
		#Used to check if time was slowed last
		self.slowed = False

		#Used when doing large calculations. The next update will be skipped to avoid time jumping.
		self.skip_update = False

		# Used for setting up shooting que:
		self.has_obejtive_que = False


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
		self.earth = Classes.earth("Images/earth.png",0)
		self.sprites_list.append(self.earth)


		#Generate debris sprite
		debris_x = random.uniform(self.satellite.model_x.value - (SCREEN_WIDTH/2), self.satellite.model_x.value + (SCREEN_WIDTH/2))
		debris_y = random.uniform(self.satellite.model_y.value-(SCREEN_WIDTH/2),self.satellite.model_y.value + (SCREEN_WIDTH/2))
		debris = Classes.debris(debris_x, debris_y, [0,0],"Images/debris.png", [0, 705*10**3, 1])
		#debris = Classes.debris(300,500,[0,0], "Images/debris.png")


		self.debris_list.append(debris)
		self.sprites_list.append(debris)

		self.satellite.update_angle(self.debris_list[0])



	def on_key_press(self, symbol, modifiers):#Buttons functions
		
		#Shoot projectile straight from satellite
		if symbol == arcade.key.ENTER:

			vec_x, vec_y = functions.angle_to_vec_2d(self.satellite.angle)
			projectile = self.satellite.get_projectile() #Classes.projectile(0.5,self.satellite.center_x,self.satellite.center_y,[vec_x,vec_y], self.canvas_info)

			self.projectile_list.append(projectile)
			self.sprites_list.append(projectile)
			

		elif symbol == arcade.key.N: #Search shoot possibilities and start mission to shoot it
			search_time = float(input("Write for how many hours the seach space should be: "))
			possibilities = functions.find_crossing_times(self.satellite,self.debris_list, ((search_time*u.h).to(u.s)).value)

			self.skip_update = True
			#functions.print_best_shots(possibilities)
			ordered_list = functions.get_ordered_target_list(possibilities)
			print(ordered_list)
			functions.plot_percentage(len(self.debris_list), ordered_list)
			#self.satellite.give_objective(possibilities)
			self.satellite.give_objectives(possibilities)


			functions.plot_result(possibilities)

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
			
		#zoom out
		elif symbol == arcade.key.O: 
			self.canvas_info[2] *= 0.5

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
			vinkel = input("write angle")
			vinkelt_int = int(vinkel)
			self.satellite.start_rotation(vinkelt_int)

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

		#Print information about debris
		elif symbol == arcade.key.W:
			functions.print_debris(self.debris_list)


	def get_max(self, max_value, value):
		if(max_value < value):
			max_value = value
		return max_value

	def get_min(self, min_value, value):
		if(min_value > value):
			min_value = value
		return min_value


	def on_draw(self): #Draws the canvas
		arcade.start_render()

		#Draw Earth and Karman line
		arcade.draw_circle_filled(self.earth.center_x, self.earth.center_y, 6371 * 1000 * self.canvas_info[2] , arcade.color.BLUE)
		if(self.display_karman_line):
			arcade.draw_circle_outline(self.earth.center_x, self.earth.center_y, (100 + 6371) * 1000 * self.canvas_info[2] , arcade.color.CYAN, 1)
			arcade.draw_circle_outline(self.earth.center_x, self.earth.center_y, (2000 + 6371) * 1000 * self.canvas_info[2] , arcade.color.CYBER_GRAPE, 1)

		#Draw all sprites
		self.sprites_list.draw()

		#Draw coordinates for nets & debris
		if(self.display_coordinates):
			for x in self.debris_list:
				if(x != self.satellite):
					arcade.draw_text("    [{}]km\n    [{}]m/s".format(x.altitude, x.speed), x.center_x, x.center_y, arcade.color.BLACK, 8)

		# Calculate time
		hours = int(self.total_time) // 3600
		minutes = int(self.total_time) // 60 % 60
		seconds = int(self.total_time) % 60 # using a modulus (remainder)
        # Figure out our output
		time_list = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

		#Write text on the screen in the top left corner
		arcade.draw_text("Number of debris: {}\nNumber of nets: {}\nNumber of netted debris: {}\nTime: {}\nTime multiplier: {}\nTimeToShoot: {:4.0f}\nSatellite weight: {}kg\nSatellite angle: {}°\nSatellite speed: {}m/s".
				   format(
					len(self.debris_list),
					len(self.projectile_list),
					len(self.netted_debris_list),
					(time_list),
					self.TIME_MULTIPLIER,
					self.satellite.time_to_shoot,
					self.satellite.mass,
					round(self.satellite.angle, 2),
					self.satellite.speed
					),10, SCREEN_HEIGHT -10, arcade.color.BLACK, 12, anchor_x="left", anchor_y="top")

		#Update scale when zooming
		middle_scale = 75/self.canvas_info[2]
		right_scale = 150/self.canvas_info[2]
		str_left_scale = 0
		str_middle_scale = 0
		str_right_scale = 0

		if(self.canvas_info[2]): #Sets string for scale in buttom right corner
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




	def update(self, delta_time): #Continuesly updates values

		if(self.update == True and self.skip_update == False): # we skip a frame when seaching for shoot possibilites as this takes a long time
			
			if self.slowed:
				self.TIME_MULTIPLIER = self.old_TIME_MULTIPLIER
				self.slowed = False
			else:
				if self.satellite.time_to_shoot.value != float("inf"):
					self.TIME_MULTIPLIER = int(self.satellite.time_to_shoot.value/10)*2
					if self.TIME_MULTIPLIER < 5:
						self.TIME_MULTIPLIER = 5

			#See if we move too fast for achieving the satellites objective
			if(self.satellite.time_to_shoot < delta_time * self.TIME_MULTIPLIER * u.s):
				self.old_TIME_MULTIPLIER = self.TIME_MULTIPLIER
				self.TIME_MULTIPLIER = (self.satellite.time_to_shoot / delta_time*u.s).value
				self.slowed = True

			if(self.satellite.time_to_hit < delta_time * self.TIME_MULTIPLIER * u.s):
				self.old_TIME_MULTIPLIER = 1# Slow way down when nearing something to hit self.TIME_MULTIPLIER
				self.TIME_MULTIPLIER = (self.satellite.time_to_hit / delta_time*u.s).value
				self.satellite.time_to_hit = float("inf")*u.s
				#self.satellite.give_next_objetive()
				self.slowed = True 
			
				
			
			self.total_time += delta_time * self.TIME_MULTIPLIER #Update total time elapsed
			if self.center_option: # when satellite is center of canvas
				self.canvas_info[0] = self.satellite.model_x.value * 1000
				self.canvas_info[1] = self.satellite.model_y.value * 1000
			else: #when earth is center of canvas
				self.canvas_info[0] = SCREEN_WIDTH/2
				self.canvas_info[1] = SCREEN_HEIGHT/2

			if self.satellite.time_to_shoot <= 0: #shoots projectile
				shot = self.satellite.get_projectile()
				self.projectile_list.append(shot)
				self.sprites_list.append(shot)
				self.satellite.time_to_shoot = float("inf") * u.s
				
			self.satellite.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info, self.total_time) #update satellite varibles

			for member in self.projectile_list: #update projectiles
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)

			for member in self.debris_list: #Update debris
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)
			for member in self.netted_debris_list: #Update nettet debris
				member.update(delta_time*self.TIME_MULTIPLIER, self.canvas_info)

				#Check if debris is within the karman line plus the radius of earth
				if np.sqrt(member.debris.model_x**2 + member.debris.model_y**2) < 6371*u.km + 100*u.km:
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

						self.satellite.give_next_objetive()
		
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