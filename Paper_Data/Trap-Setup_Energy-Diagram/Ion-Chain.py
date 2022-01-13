# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:14:12 2021

@author: brend
"""

from PIL import Image
import numpy as np


im = Image.open("Ion-Image_04_17_Ba138-Chain.bmp")
width, height = im.size
left = 10
right = 47
top = 18
bottom = 26
im = im.crop((left, top, right, bottom))
width, height = im.size
# im = im.convert("RGB")

image_array = np.zeros((height, width))
for i in range(width):
    for j in range(height):
        image_array[j, i] = im.getpixel((i, j))
        
# image_max = np.max(image_array)
    
# source = im.split()
# R, G, B = 0, 1, 2

# gval = 100
# bval = 253

# outR = source[R].point(lambda i: 0)
# outG = source[G].point(lambda i: i*gval/image_max)
# outB = source[B].point(lambda i: i*bval/image_max)

# sourceout = (outR, outG, outB)
# im = Image.merge("RGB", sourceout)

# image_array = np.zeros((height, width))
# for i in range(width):
#     for j in range(height):
#         image_array[j, i] = im.getpixel((i, j))[2]
        
im.save('Chain_Ba138.bmp')