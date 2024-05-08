# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 12:17:20 2024

@author: Vicenç Fernàndez Fernàndez
"""

"""
This code implement color dithering using Bayer matrices 
and error difusion algorithm
"""


import skimage
import matplotlib.pyplot as plt
import numpy as np


im = skimage.data.astronaut()

"""
BAYER DITHERING
"""

# STEP 1. Create Threshold array

B = np.array([[0,8,2,10],[12,4,14,6],[3,11,1,9],[15,7,13,5]])

T_prima=np.zeros(np.shape(im[:,:,1]))

for i in range(0,512,4):
    for j in range(0,512,4):
        T_prima[j:j+4,i:i+4]=B
        
T = T_prima - 0.5 # Binarization array

# Step 2: Image quantization
im_quan= im/16
im_quan_int= np.uint8(im_quan)

# Step 3. Binarization
im_R= im_quan[:,:,0]
im_G= im_quan[:,:,1]
im_B= im_quan[:,:,2]


# DITHERING

red=255*np.uint8(im_R>=T)
green=255*np.uint8(im_G>=T)
blue= 255*np.uint8(im_B>=T)

RGB1= np.zeros((512,512,3))
RGB1[:,:,0]= red
RGB1[:,:,1]= green
RGB1[:,:,2]= blue


# Plot
plt.figure()
plt.subplot(2,2,1)
plt.imshow(im)
plt.title('Original')
plt.axis('off')

plt.subplot(2,2,2)
plt.imshow(im_quan_int*16)
plt.title('Quantization')
plt.axis('off')


plt.subplot(2,2,3)
plt.imshow(np.uint8(RGB1))
plt.title('Binarization')
plt.axis('off')


"""
ERROR DIFFUSION ALGORITHM STUCKI MATRIX
"""
def im_difusion(im,th):
    for i in range(im.shape[0]-2):
        for j in range(im.shape[1]-2):
            if im[i][j] < th:
                error = im[i,j]
            if im[i][j] > th:
                error = im[i,j] - 255
                
            im[i,j+1] +=  error*8/42
            im[i,j+2] += error*4/42
            
            im[i+1,j-2]+= error*2/42
            im[i+1,j-1] += error*4/42
            im[i+1,j] += error*8/42
            im[i+1,j+1] += error*4/42
            im[i+1,j+2] += error*2/42
            
            im[i+2,j-1] += error*1/42
            im[i+2,j-2] += error*2/42
            im[i+2,j] += error*4/42
            im[i+2,j+1] += error*2/42
            im[i+2,j+2] += error*1/42
            
    return im



lut_bina_ = np.zeros(256,dtype=np.uint8)
for i in range(256):
    if i > 120:
        lut_bina_[i]=255

     
im_R= im[:,:,0]
im_G= im[:,:,1]
im_B= im[:,:,2]

im_dif_r= im_difusion(im_R,120)
im_dif_g= im_difusion(im_G,120)
im_dif_b= im_difusion(im_B,120)

ima_bina_r = lut_bina_[im_dif_r]
ima_bina_g = lut_bina_[im_dif_g]
ima_bina_b = lut_bina_[im_dif_b]



RGB= np.zeros((512,512,3))
RGB[:,:,0]= ima_bina_r
RGB[:,:,1]= ima_bina_g
RGB[:,:,2]= ima_bina_b

plt.subplot(2,2,4)
plt.imshow(np.uint8(RGB))
plt.title('Binarization')
plt.axis('off')
plt.show()



