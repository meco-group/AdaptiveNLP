import os
from PIL import Image

def create_gif(png_folder, gif_file, delay):
    # Get the list of PNG files in the folder and sort them based on frame number
    png_files = sorted([file for file in os.listdir(png_folder) if file.endswith('.png')], key=lambda x: int(x.split('_')[-1].split('.')[0]))

    images = []
    for png_file in png_files:
        image_path = os.path.join(png_folder, png_file)
        image = Image.open(image_path)
        images.append(image)

    # Save as an animated GIF
    images[0].save(gif_file, save_all=True, append_images=images[1:], optimize=False, duration=delay, loop=0)

import pathlib
png_folder = str(pathlib.Path(__file__).parent.resolve()) # 'figures/animation'
frame_format = '.png'
gif_file = png_folder + "/animation.gif" #'figures/animation/animation.gif'
delay = 300  # Delay between frames in milliseconds

create_gif(png_folder, gif_file, delay)