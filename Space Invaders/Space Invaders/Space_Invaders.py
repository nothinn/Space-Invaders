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


	def setup(self):

		#Prepare lists to contain all sprites
		self.debris_list = arcade.SpriteList()
		self.sprites_list =	arcade.SpriteList()
		self.projectile_list = arcade.SpriteList()


		#Generate satellite sprite
		self.satellite = Classes.satellite("Images/satellite.png",0.2)
		self.sprites_list.append(self.satellite)


		#Generate debris sprite
		debris = Classes.debris(random.uniform(0,600),random.uniform(0,600),[0,0], "Images/debris.png")
		#debris = Classes.debris(300,500,[0,0], "Images/debris.png")


		self.debris_list.append(debris)
		self.sprites_list.append(debris)


		

	def on_key_press(self, symbol, modifiers):
		if symbol == arcade.key.ENTER:

			vec_x, vec_y = functions.angle_to_vec_2d(self.satellite.angle)
			projectile = Classes.projectile(0.5,self.satellite.center_x,self.satellite.center_y,[vec_x,vec_y])

			self.projectile_list.append(projectile)
			self.sprites_list.append(projectile)


	def on_draw(self):
		arcade.start_render()

		#Draw all sprites
		self.sprites_list.draw()



	def update(self, delta_time):
		self.satellite.update_angle(self.debris_list[0])
		
		for member in self.projectile_list:
			member.update(delta_time*100)





def main():
	game = MyGame(SCREEN_WIDTH,SCREEN_HEIGHT)
	game.setup()




	arcade.run()


	print("Running")

if __name__ == "__main__":
	main()