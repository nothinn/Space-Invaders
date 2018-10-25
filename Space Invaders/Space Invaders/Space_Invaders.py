import arcade
import os

import random
import functions
import Classes

import math


SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600



class MyGame(arcade.Window):
	def __init__(self, width, height):
		super().__init__(width,height)
		arcade.set_background_color(arcade.color.AMAZON)
		
		#Set path directory
		self.file_path = os.path.dirname(os.path.abspath(__file__))
		os.chdir(self.file_path)


		#Set lists for sprites
		self.sprites_list = None
		self.debris_list = None
		self.projectile_list = None
		self.netted_debris_list = None

		self.right = 0
		self.left = 0


	def setup(self):

		#Prepare lists to contain all sprites
		self.debris_list = arcade.SpriteList()
		self.sprites_list =	arcade.SpriteList()
		self.projectile_list = arcade.SpriteList()
		self.netted_debris_list = list()


		#Generate satellite sprite
		self.satellite = Classes.satellite("Images/satellite.png",0.2)
		self.sprites_list.append(self.satellite)


		#Generate debris sprite
		debris = Classes.debris(random.uniform(0,600),random.uniform(0,600),[0,0], "Images/debris.png")
		#debris = Classes.debris(300,500,[0,0], "Images/debris.png")


		self.debris_list.append(debris)
		self.sprites_list.append(debris)

		self.satellite.update_angle(self.debris_list[0])


		



		

	def on_key_press(self, symbol, modifiers):
		if symbol == arcade.key.ENTER:

			vec_x, vec_y = functions.angle_to_vec_2d(self.satellite.angle)
			projectile = Classes.projectile(0.5,self.satellite.center_x,self.satellite.center_y,[vec_x,vec_y])

			self.projectile_list.append(projectile)
			self.sprites_list.append(projectile)
		
		elif symbol == arcade.key.M:

			new_debris = Classes.debris(random.uniform(0,600), random.uniform(0,600), [math.cos(random.uniform(-1*math.pi,math.pi))*0.15, math.sin(random.uniform(-1*math.pi,math.pi))*0.15], "Images/debris.png")
			self.debris_list.append(new_debris)
			self.sprites_list.append(new_debris)

			pro_angle_1 = functions.get_net_angle_immediate(0.2, self.satellite.center_x, self.satellite.center_y, new_debris.vel_vector[0], new_debris.vel_vector[1], new_debris.center_x, new_debris.center_y)

			new_projectile = Classes.projectile(0.5, self.satellite.center_x, self.satellite.center_y,[math.cos(math.radians(pro_angle_1))*0.2, math.sin(math.radians(pro_angle_1))*0.2])
			self.projectile_list.append(new_projectile)
			self.sprites_list.append(new_projectile)

		elif symbol == arcade.key.LEFT:
			self.left = 1
		elif symbol == arcade.key.RIGHT:
			self.right = 1

	def on_key_release(self, symbol, modifiers):
		if symbol == arcade.key.LEFT:
			self.left = 0
		elif symbol == arcade.key.RIGHT:
			self.right = 0


	def on_draw(self):
		arcade.start_render()


		#Write text on the screen in the top left corner
		arcade.draw_text("Number of debris:{}\nNumber of nets:{}\nNumber of netted debris:{}".format(len(self.debris_list),len(self.projectile_list),len(self.netted_debris_list)),
                         10, SCREEN_HEIGHT -10, arcade.color.BLACK, 12, anchor_x="left", anchor_y="top")

		#Draw all sprites
		self.sprites_list.draw()



	def update(self, delta_time):
		self.satellite.rotate(delta_time,self.left - self.right)
		
		for member in self.projectile_list:
			member.update(delta_time*100)

		for member in self.debris_list:
			member.update(delta_time*100)
		for member in self.netted_debris_list:
			member.update(delta_time*100)
			if abs(member.debris.center_x) > 1000 or abs(member.debris.center_y) > 1000:
				member.kill()
				print("Killed netted debris")
				self.netted_debris_list.remove(member)


	


		#Test for collisions:
		for idx, projectile in enumerate(self.projectile_list):
			for idy, debris in enumerate(self.debris_list):
				if functions.collision(projectile, debris):
					print("Collision detected at X:{}, Y:{}".format(projectile.center_x, projectile.center_y))

					self.netted_debris_list.append(Classes.netted_debris(projectile,debris))
					self.projectile_list.remove(projectile)
					self.debris_list.remove(debris)
					
					
					


def main():
	game = MyGame(SCREEN_WIDTH,SCREEN_HEIGHT)
	game.setup()




	arcade.run()


	print("Running")

if __name__ == "__main__":
	main()