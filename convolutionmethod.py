import numpy as np
import cv2
import math
import holopy as hp
import tifffile 

import matplotlib.pyplot as plt

#load image
image_path="/home/alinsi/Desktop/inlinehologram/apirl20laser310x.bmp"
#image_path2="/home/owner/hologram/samples/dec22_3.tif"
#image_path3="/home/owner/hologram/samples/dec22_2.tif"

h = cv2.imread(image_path,0)
plt.imshow(h,cmap='Greys_r')
#k = cv2.imread(image_path2,0)
#l = cv2.imread (image_path3,0)

h1 = np.array(h).astype(float)

#plt.imshow(k,cmap='Greys_r')
#plt.imshow(h1,cmap='Greys_r')
#
# #load background
# back = cv2.imread('/home/owner/hologram/samples/background2.tif',0)
# background = back.astype(float)
# # #subtract background from image
# h1=h1-background



(Nx,Ny)= h1.shape[:2]
minN = min(Nx,Ny)
h1 = h1[:minN,:minN]
#hp.save('/home/alinsi/Desktop/lasersquare',h1)

(Nx2,Ny2) = h1.shape[:2]
#plt.imshow(h1,cmap='Greys_r')
#plt.imsave("/home/alinsi/Desktop/inlinehologram/processed/may1microbeads5original",h1,cmap='Greys_r')
lambda0 = 0.000488
delx=5.32/1024
dely=6.66/1280
i=complex(0,1)
pi=math.pi
maxd=110#looping of distance in mm from object to CCD , maximum is 22cm, 220mm,minmum is 60mm6
mind=50
steps=5#distanced looping step

xmax=Nx2*delx/2
ymax=Ny2*dely/2

nx = np.arange (-Nx2/2,Nx2/2,1)
ny = np.arange (-Ny2/2,Ny2/2,1)

X=nx*delx
Y=ny*dely

[XX,YY]=np.meshgrid(X,Y)

# X2=nx/(Nx2*delx)
# Y2=ny/(Ny2*dely)
#
# [XX2,YY2]=np.meshgrid(X2,Y2)
#gy maginification factor for the reconstruction
Gy = 10
k=2*pi/lambda0



#mag0 = np.arange[0,10,5]
#mag = 15
imageslices=(maxd-mind)/steps
threeD=np.zeros([minN,minN,imageslices])
imaginary=np.zeros([minN,minN,imageslices])
original=np.zeros([minN,minN,imageslices])


dp = np.arange(mind, maxd, steps)

for d in dp:

    #d0= dp*mag
    #f= (1/dp + 1/d0)**(-1)
    #L = np.exp(i*pi/lambda0/f*(XX**2,YY**2)) #reference wave
    #P = #factor , correction of reconstructed wave field caused by lens(phase aberrations)

    num = np.exp(-i*2*pi/lambda0*np.sqrt((d**2+XX**2+YY**2)))
    den = np.sqrt(d**2+XX**2+YY**2)
    g=i/lambda0*num/den#g is the impulse response of propagation
    #############################################
    # SPHEREICAL Waves##
    # zi = -Gy * d;
    # zc = 1 / (1 / d + 1 / zi)
    # sphere = np.exp(i * k/2/zc*(XX**2+YY**2))

    #################################################

    G=np.fft.fft2(g)

    H=np.fft.fft2(h1)#plane wave is assumed to be a constant amplitude

    Rec_image=np.fft.fftshift(np.fft.ifft2(G*H))
    #amp_rec_image = np.abs(Rec_image)

    #contrast=np.arctan(np.imag(Rec_image)/np.real(Rec_image))`` ##this is phase contrast

    #height=lambda0/4/pi*contrast



    ind=(d-mind)/steps
    #

    threeD[:,:,ind]=np.abs(Rec_image)
    threeD=threeD.astype(np.float32)


  #  imaginary[:, :, ind] = np.imag(Rec_image)
  #  imaginary= imaginary.astype(np.float32)

    #
   # original[:, :, ind] = np.angle(Rec_image)
   # original= original.astype(np.float32)




     #to close windows plt.close('all')
      #plt.savefig("/home/owner/hologram/samples/1/convolution/" + str(d) + ".tif")

# #
# tifffile.imsave("/home/owner/hologram/samples/1/convolution/dec15micro2/50220angle.tif",original.transpose())
# tifffile.imsave("/home/owner/hologram/samples/1/convolution/dec15micro2/50220imaginary.tif",imaginary.transpose())
# tifffile.imsave("/home/owner/hologram/samples/1/convolution/dec15micro2/50220amp.tif",threeD.transpose())
# # #
#hp.show(original)
#hp.show(imaginary)

hp.show(threeD)

#tifffile.imsave("/home/alinsi/Desktop/inlinehologram/processed/may1Two10Xmicrobeads5",threeD.transpose())
