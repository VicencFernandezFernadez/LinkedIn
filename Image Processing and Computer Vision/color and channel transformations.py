# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 12:50:20 2024

@author: Vicenç Fernandez Fernandez
"""

import skimage
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from numpy.linalg import inv

vector=np.load('RGB.npy')


# 4.1 RGB coordinates 
lamda = vector[:,0]
x_lam = vector[:,1]
y_lam = vector[:,2]
z_lam = vector[:,3]
E_lam = vector[:,4]

x= sp.integrate.simpson(E_lam*x_lam,lamda) 
y= sp.integrate.simpson(E_lam*y_lam,lamda)
z= sp.integrate.simpson(E_lam*z_lam,lamda)

X= x/(x+y+z)
Y= y/(x+y+z)
Z= z/(x+y+z)

R = np.uint8(255*(3.240479*X -1.537150*Y -0.498535*Z))
G = np.uint8(255*( -0.969256*X+1.875992*Y+ 0.041556*Z)) 
B = np.uint8(255*(0.055648*X -0.204043*Y + 1.057311*Z))

RGB_array= np.array([R,G,B])
N=300
RGB= np.zeros((N,N,3),dtype=np.uint8)
RGB[:,:,0]= RGB_array[0]
RGB[:,:,1]= RGB_array[1]
RGB[:,:,2]= RGB_array[2]

plt.figure()
plt.imshow(RGB,cmap='gray')
plt.show()

# YCbCr

Y = 0.299*R+ 0.587*G +0.114*B
Cb = 0.168736*R+ 0.331264*G+ 0.5*B + 128
Cr =0.5*R - 0.418688*G -0.081312*B + 128

astro= skimage.data.astronaut()


# STEP 1
Y_astro = 0.299*astro[:,:,0]+ 0.587*astro[:,:,1] +0.114*astro[:,:,2]
Cb_astro = 0.168736*astro[:,:,0]+ 0.331264*astro[:,:,1]+ 0.5*astro[:,:,2] + 128
Cr_astro =0.5*astro[:,:,0]- 0.418688*astro[:,:,1] -0.081312*astro[:,:,1] + 128

plt.figure()
plt.subplot(3,3,1)
plt.imshow(Y_astro,cmap='gray')
plt.title('Y')
plt.axis('off')

plt.subplot(3,3,2)
plt.imshow(Cb_astro, cmap='viridis')
plt.title('Cb')
plt.axis('off')


plt.subplot(3,3,3)
plt.imshow(Cr_astro, cmap='viridis')
plt.title('Cr')
plt.axis('off')

plt.subplot(3,3,4)
plt.imshow(Y_astro[0:200,120:320],cmap='gray')
plt.title('Y')
plt.axis('off')

plt.subplot(3,3,5)
plt.imshow(Cb_astro[0:200,120:320], cmap='viridis')
plt.title('Cb')
plt.axis('off')


plt.subplot(3,3,6)
plt.imshow(Cr_astro[0:200,120:320], cmap='viridis')
plt.title('Cr')
plt.axis('off')




# STEP 2. Pixelated version of Cb an Cr
Cb_pixel = np.zeros((512,512))
Cr_pixel = np.zeros((512,512))

for i in range(0,512,4):
    for j in range(0,512,4):
        Cb_pixel[i:i+4, j:j+4] = Cb_astro[i,j]
        Cr_pixel[i:i+4, j:j+4] = Cr_astro[i,j]




plt.subplot(3,3,8)
plt.imshow(Cb_pixel[0:200,120:320],cmap='viridis')
plt.title('Cb pixel')
plt.axis('off')

plt.subplot(3,3,9)
plt.imshow(Cr_pixel[0:200,120:320], cmap='viridis')
plt.title('Cr pixel')
plt.axis('off')
plt.show()

# STEP 4: convert YCb’Cr’ to new image R’G’B’
m=np.array([[0.299,0.587,0.114],[-0.168736,-0.331264, 0.5],
           [0.5,-0.418688,-0.081312]])

inverse = inv(m)

YCbCr= np.array([Y_astro,Cb_pixel,Cr_pixel])

R_new = (inverse[0,0]*Y_astro + inverse[0,1]*(Cb_pixel-128)
        + inverse[0,2]*(Cr_pixel-128))

G_new =(inverse[1,0]*Y_astro + inverse[1,1]*(Cb_pixel-128)
        + inverse[1,2]*(Cr_pixel-128)) 

B_new = (inverse[2,0]*Y_astro + inverse[2,1]*(Cb_pixel-128)
        + inverse[2,2]*(Cr_pixel-128)) 

def check(original_matrix):
    matrix=np.copy(original_matrix)
    for i in range(512):
        for j in range(512):
            if matrix[i,j]>255:
                matrix[i,j]=255
            if matrix[i,j]<0:
                matrix[i,j]=0
    return matrix

R_checked = check(R_new)
G_checked= check(G_new)
B_checked= check(B_new)



plt.figure()
plt.subplot(2,3,1)
plt.imshow(astro[:,:,0])
plt.title('R')
plt.axis('off')

plt.subplot(2,3,2)
plt.imshow(astro[:,:,1])
plt.title('G')
plt.axis('off')


plt.subplot(2,3,3)
plt.imshow(astro[:,:,2])
plt.title('B')
plt.axis('off')

plt.subplot(2,3,4)
plt.imshow(R_checked)
plt.title("R'")
plt.axis('off')

plt.subplot(2,3,5)
plt.imshow(G_checked)
plt.title("G'")
plt.axis('off')


plt.subplot(2,3,6)
plt.imshow(B_checked)
plt.title("B'")
plt.axis('off')
plt.show()

astro_prima=np.zeros((512,512,3), dtype=np.uint8)
astro_prima[:,:,0]= R_checked
astro_prima[:,:,1]= G_checked
astro_prima[:,:,2]= B_checked

plt.figure()
plt.imshow(astro_prima)
plt.title("New astro")
plt.axis('off')
plt.show()

for i in range(0,3):
    smss=skimage.metrics.structural_similarity(astro_prima[:,:,i],astro[:,:,i])
    print(smss)

# 4.3. The HSV (hue, saturation, value) color space. Use in image fusion

def RGBtoHSV(image_RGB):
    imageRGB=np.copy(image_RGB)
    imageRGB=imageRGB/imageRGB.max()
    imageHSV = np.zeros(imageRGB.shape)
    
    for i in range(imageRGB.shape[0]):
        for j in range(imageRGB.shape[1]):
            R=imageRGB[i,j,0]
            G=imageRGB[i,j,1]
            B=imageRGB[i,j,2]
            V=max(R,G,B)
            m=min(R,G,B)
            C=V-m
            if C==0:
                H=0
            elif V==R:
                H=((G-B)/C)%6
            elif V==G:
                H=(B-R)/C+2
            elif V==B:
                H=(R-G)/C+4
            H_t=H/6
            
            if V==0:
                S=0
            else:
                S=C/V
            imageHSV[i,j,0] = np.float64(H_t)
            imageHSV[i,j,1] = np.float64(S)
            imageHSV[i,j,2] = np.float64(V)
    
    return imageHSV

def f(n,H,V,S):
    #global H, V, S
    k=(n+6*H)%6
    return V-V*S*max(0,min(k,4-k,1))

def HSVtoRGB(image_HSV):
    imageHSV= np.copy(image_HSV)
    imageHSV= imageHSV/imageHSV.max()
    imageRGB = np.zeros(imageHSV.shape)
    
    for i in range(imageHSV.shape[0]):
        for j in range(imageHSV.shape[1]):
            H=imageHSV[i,j,0]
            S=imageHSV[i,j,1]
            V=imageHSV[i,j,2]
            
            imageRGB[i,j,0] = f(5,H,V,S)
            imageRGB[i,j,1] = f(3,H,V,S)
            imageRGB[i,j,2] = f(1,H,V,S)
            
    return imageRGB
    


imCafe = skimage.data.coffee()
imCafe_HSV = RGBtoHSV(imCafe)
imCafe_HSV2 = skimage.color.rgb2hsv(imCafe)
imCafe_RGB = HSVtoRGB(imCafe_HSV)


plt.figure()

plt.subplot(2, 3, 1)
plt.imshow(imCafe_HSV, cmap='viridis')
plt.title('Cafe HSV')
plt.axis('off')

plt.subplot(2, 3, 2)
plt.imshow(imCafe_RGB, cmap='gray')
plt.title('Cafe RGB Converted')
plt.axis('off')

plt.subplot(2, 3, 3)
plt.imshow(imCafe, cmap='gray')
plt.title('Cafe RGB Original')
plt.axis('off')

plt.subplot(2, 3, 4)
plt.imshow(imCafe_HSV[:,:,0], cmap='viridis')
plt.title('Cafe H')
plt.axis('off')

plt.subplot(2, 3, 5)
plt.imshow(imCafe_HSV[:,:,1], cmap='viridis')
plt.title('Cafe S')
plt.axis('off')

plt.subplot(2, 3, 6)
plt.imshow(imCafe_HSV[:,:,2], cmap='viridis')
plt.title('Cafe V')
plt.axis('off')


astro_HSV_teoric= skimage.color.rgb2hsv(astro)


astro_HSVtoRGB = skimage.color.hsv2rgb(astro_HSV_teoric)




# Histogram equalization. Image entropy.
#%%
#imBody = plt.imread('imatges/body.png')
imBody = skimage.data.cell()


 # old way
 
# Histogram and cumulative histogram of the original image
bins = 256
h = sp.ndimage.histogram(imBody, 0, 255, bins)
ch = h.cumsum()


# Equalized image histogram and cumulative histogram
M = imBody.shape[0]
N = imBody.shape[1]

imBody_eq = np.uint8(255 * ch[imBody] / (M * N))


h_eq = sp.ndimage.histogram(imBody_eq, 0, 255, bins)
ch_eq = h_eq.cumsum()


# fast way
def equalization(image):
    bins = 256
    M, N = image.shape[0], image.shape[1]
    
    h = sp.ndimage.histogram(image, 0, 255, bins)
    ch=h.cumsum()
    
    image_eq = np.uint8(255*ch[image]/(M * N))
    
    h_eq = sp.ndimage.histogram(image_eq, 0, 255, bins)
    ch_eq = h_eq.cumsum()
    
    return image_eq, h, h_eq, ch, ch_eq

imBody_eq, h, h_eq, ch, ch_eq = equalization(imBody)


# Subplot
plt.figure(figsize=(20, 10))
plt.subplots_adjust(hspace=0.4)

x = np.arange(bins)

plt.subplot(2, 2, 1)
plt.plot(x, h, color='red')
plt.title('Histogram Body')
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')

plt.subplot(2, 2, 2)
plt.plot(x, ch, color='red')
plt.title('CumSum Histogram Body')
plt.xlabel('Pixel Value')
plt.ylabel('Cumulative Frequency')

plt.subplot(2, 2, 3)
plt.plot(x, h_eq, color='blue')
plt.title('Histogram Eq Body')
plt.xlabel('Pixel Value')
plt.ylabel('Frequency')

plt.subplot(2, 2, 4)
plt.plot(x, ch_eq, color='blue')
plt.title('CumSum Histogram Eq Body')
plt.xlabel('Pixel Value')
plt.ylabel('Cumulative Frequency')



plt.figure()


plt.subplot(1, 2, 1)
plt.imshow(imBody, cmap='gray')
plt.title('Body Original')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(imBody_eq, cmap='gray')
plt.title('Body Equalized')
plt.axis('off')



