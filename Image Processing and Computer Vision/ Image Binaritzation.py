'''
Practica 2: Image Binaritzation
Vicenç Fernandez Fernandez
'''

import skimage
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sp
from tqdm import tqdm

# Creation of an bidimiensional array that contains retina info
im= skimage.data.retina()

red = im[:,:,0]
green = im[:,:,1]
blue = im[:,:,2]

L = 0.299*red +0.587*green + 0.114*blue


# Creation of binarization lut with a threshold 

lut_bina = np.zeros(256,dtype=np.uint8)
for i in range(256):
    if i > 105:
        lut_bina[i]=255
    
    
    
L_int= np.uint8(L)    
retina_bina = lut_bina[L_int]


med = sp.medfilt(L_int, kernel_size=7)
median = 255*np.uint8(L_int>=med)



order= sp.order_filter(L_int, domain=np.ones((3,3)), rank=4)
order_=255*np.uint8(L_int>=order)

plt.figure()

plt.subplot(2,2,1)
plt.imshow(L,cmap='gray')
plt.title('ORIGINAL')
plt.axis('off')

plt.subplot(2,2,2)
plt.imshow(retina_bina, cmap='gray')
plt.title('LUT 105')
plt.axis('off')


plt.subplot(2,2,3)
plt.imshow(median, cmap='gray')
plt.title('MEDFILT')
plt.axis('off')

plt.subplot(2,2,4)
plt.imshow(order_, cmap='gray')
plt.title('ORDER')
plt.axis('off')
plt.show()


# DITHERING
th = 128
cell= skimage.data.cell()
new_cell = np.copy(cell)

for i in tqdm(range(cell.shape[0]-1)):
    for j in range(cell.shape[1]-1):
        if new_cell[i][j] < th:
            error = new_cell[i,j]
        if new_cell[i][j] > th:
            error = new_cell[i,j] - 255
            
        new_cell[i,j+1] +=  error*7/16
        new_cell[i+1,j] += error*5/16
        new_cell[i+1,j+1]+= error*1/16
        new_cell[i+1,j-1] += error*3/16
        

lut_bina_ = np.zeros(256,dtype=np.uint8)
for i in range(256):
    if i > th:
        lut_bina_[i]=255
  
new_cell_bina = lut_bina_[new_cell]

plt.figure()
plt.subplot(1,2,1)
plt.imshow(cell, cmap='gray')
plt.title('ORIGINAL')
plt.axis('off')

plt.subplot(1,2,2)
plt.imshow(new_cell_bina, cmap='gray')
plt.title('NEW')
plt.axis('off')


plt.show()
            

