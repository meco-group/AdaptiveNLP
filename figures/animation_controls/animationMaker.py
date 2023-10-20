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

# Example usage
png_folder = 'figures/animation_controls'
gif_file = 'figures/animation_controls/animation_controls.gif'
delay = 100  # Delay between frames in milliseconds

create_gif(png_folder, gif_file, delay)