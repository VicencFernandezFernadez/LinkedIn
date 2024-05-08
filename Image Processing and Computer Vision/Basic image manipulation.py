# Vicenç Fernàndez Fernàndez

import matplotlib.pyplot as plt
import skimage
import numpy as np

# Imports retina image with skiimage
im=skimage.data.retina()

plt.figure()

# Show the image
plt.imshow(im)
plt.show()

# Show the size of the image
print(im.size)
print(len(im))
print(im.shape)

# Set an arry with 0
red_im= np.zeros(im.shape, dtype=np.uint8)
green_im = np.zeros(im.shape, dtype=np.uint8)
blue_im = np.zeros(im.shape, dtype=np.uint8)


# CHANNEL EXTRACTION
# Red, Blue, green, and grises
red_im[:,:,0] = im[:,:,0]
green_im[:,:,1] = im[:,:,1]
blue_im[:,:,2] = im[:,:,2]

grayred=im[:,:,0]
graygreen=im[:,:,1]
grayblue=im[:,:,2]



#PLOTS THE RESULT

plt.figure(figsize=(10,5))

plt.subplot(2,3,1)
plt.imshow(grayred, cmap='gray')

plt.subplot(2,3,2)
plt.imshow(graygreen, cmap='gray')

plt.subplot(2,3,3)
plt.imshow(grayblue, cmap='gray')

plt.subplot(2,3,4)
plt.imshow(red_im)

plt.subplot(2,3,5)
plt.imshow(green_im)

plt.subplot(2, 3,6)
plt.imshow(blue_im)
plt.show()




# Lightness, colormaps and false color

# a) avarage

im=np.double(skimage.data.retina()) #implement double to prevent pasing 255

# Matrix that content the information for red blue and green channel
red=im[:,:,0]
green=im[:,:,1]
blue=im[:,:,2]

result= (red + green +blue)/3
print(result)

plt.imshow(result/result.max(), cmap='gray')
plt.show()

# b)  Luma s L=0.299*R+0.587*G+0.114*B

L = 0.299 * red +0.587*green + 0.114*blue
print(L)
plt.imshow(L/L.max(), cmap='gray')
plt.show()

####  c)  imm = im^g. Consider (i) 0<g<1 and (ii) g>1

cell= np.double(skimage.data.cell())

# Show the image of the cell
plt.imshow(cell)
plt.show()



imm_05= cell**0.5
imm_2= cell**2

plt.subplot(1,2,1)
plt.imshow(imm_05,cmap='gray')

plt.subplot(1,2,2)
plt.imshow(imm_2,cmap='gray')
plt.show()



### d) Aplications of LUTS

# 
cell_i = skimage.data.cell()

# Create all the luts with an array with 0

lut_contrast=np.zeros(256, dtype=np.uint8)
lut_linear = np.zeros(256, dtype=np.uint8)
lut_logistic = np.zeros(256, dtype=np.uint8)
lut_binarization = np.zeros(256, dtype=np.uint8)

# Calculate the luts
for i in range(256):
    lut_contrast[i]= 255 - i
    lut_linear[i] = min(2*i + 120, 255)
    lut_logistic[i] = 255/(1+np.exp(-0.2*(i-70)))
    
for i in range(128):
    lut_binarization[128+i]=255
    

# Aplication of the luts
cell_contrast=lut_contrast[cell_i]
cell_linear=lut_linear[cell_i]
cell_logistic=lut_logistic[cell_i]
cell_binarization=lut_binarization[cell_i]

# Show original
plt.subplot(2,4,(2,3))
plt.imshow(cell_i, cmap='gray')
plt.title('ORIGINAL')
plt.axis('off')

# Show contrast
plt.subplot(2,4,5)
plt.imshow(cell_contrast, cmap='gray')
plt.title('CONTRSAST')
plt.axis('off')


# Show linear lut
plt.subplot(2,4,6)
plt.imshow(cell_linear, cmap='gray')
plt.title('LINEAR')
plt.axis('off')

# Show logisitic
plt.subplot(2,4,7)
plt.imshow(cell_logistic, cmap='gray')
plt.title('LOGISTIC')
plt.axis('off')

# Show binarization
plt.subplot(2,4,8)
plt.imshow(cell_binarization, cmap='gray')
plt.title('BINARIZATION')
plt.axis('off')

plt.show()


# Comparation

MR= skimage.metrics.mean_squared_error(red,L)
MG= skimage.metrics.mean_squared_error(green,L)
MB= skimage.metrics.mean_squared_error(blue,L)

MSE= (MR+MG+MB)/3
print(MSE)




