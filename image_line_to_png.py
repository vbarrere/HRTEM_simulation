#!/usr/bin/env python3
import os

import numpy as np
from PIL import Image

image_file = "images.dat"
#image_file = "/home/victor/Data/HRTEM_data/AgCo/Dataset1/images_full.dat"
out_dir = "images_png"
#out_dir="/home/victor/Data/HRTEM_data/AgCo/Dataset1/images_png"
nx = 96

os.makedirs(out_dir, exist_ok=True)

with open(image_file) as f:
    for line_no, line in enumerate(f, start=1):
        line = line.strip()
        if not line:
            continue
        id_sim, _, pixels = line.partition(" ")
        img = np.array(pixels.split(), dtype=np.int16)
        img = (img + 128).astype(np.uint8).reshape(nx, nx)
        #out_png = os.path.join(out_dir, f"image_{line_no:05d}.png")
        out_png = os.path.join(out_dir, f"{id_sim}.png")
        Image.fromarray(img, mode="L").save(out_png)
        print(f"wrote {out_png}")