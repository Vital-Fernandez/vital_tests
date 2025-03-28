import imageio.v2 as imageio
import os
from PIL import Image

# Define the folder containing your image frames
frames_folder = '/home/vital/Astrodata/CAPERS/CAPERS_UDS_V0.1/galaxy_plots/'
output_gif = "animation.gif"
frame_duration = 0.2  # seconds per frame

# Load, resize, and convert all frames to same size/mode
frame_files = sorted([
    os.path.join(frames_folder, f)
    for f in os.listdir(frames_folder)
    if f.endswith(".png")
])

# Open the first image to get reference size/mode
first_image = Image.open(frame_files[0])
ref_size = first_image.size
ref_mode = first_image.mode

# Collect processed frames
frames = []
for f in frame_files:
    img = Image.open(f).convert(ref_mode).resize(ref_size)
    frames.append(img)

# Save to GIF
frames[0].save(
    output_gif,
    save_all=True,
    append_images=frames[1:],
    duration=int(frame_duration * 1000),
    loop=0
)

print(f"GIF saved as {output_gif}")