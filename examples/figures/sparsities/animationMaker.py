import os
from PIL import Image

def create_gif(png_folder, gif_file, delay, start_str):
    # Get the list of PNG files in the folder and sort them based on frame number
    png_files = sorted([file for file in os.listdir(png_folder) if file.startswith(start_str) and file.endswith('.png')], key=lambda x: int(x.split('_')[-1].split('.')[0]))

    images = []
    for png_file in png_files:
        image_path = os.path.join(png_folder, png_file)
        image = Image.open(image_path)
        images.append(image)

    # Save as an animated GIF
    images[0].save(gif_file, save_all=True, append_images=images[1:], optimize=False, duration=delay, loop=0)

import pathlib
png_folder = str(pathlib.Path(__file__).parent.resolve()) # 'figures/sparsities'
gif_folder = png_folder + "/"

delay = 300  # Delay between frames in milliseconds

create_gif(png_folder, gif_folder+'animation_jac_adaptive.gif', delay, 
           'jac_adaptive')
create_gif(png_folder, gif_folder+'animation_jac_casadi.gif', delay, 
           'jac_casadi')
create_gif(png_folder, gif_folder+'animation_hess_adaptive.gif', delay, 
           'hess_adaptive')
create_gif(png_folder, gif_folder+'animation_hess_casadi.gif', delay, 
           'hess_casadi')