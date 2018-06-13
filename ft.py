#!/bin/env python

#---------------------------------------------
# using the fourier transformation equation to
# demenstrate FFT
#----------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

n=256
dt=0.02
x=np.arange(n)*dt
y=np.sin(2*3.141592768*x)+0.5*np.sin(2*3.141592768*5*x)

m=int(np.floor(n/2)+1)

a=np.zeros(m,dtype=float)
b=np.zeros(m,dtype=float)
c=np.zeros(m,dtype=float)

cx=np.arange(m)/(n*dt)
#--------------------------------
# do fft
#--------------------------------
for k in np.arange(m):
    for ii in np.arange(n):
        a[k]=a[k]+2./n*y[ii]*np.cos(2*3.141592768*k*ii/n)
        b[k]=b[k]+2./n*y[ii]*np.sin(2*3.141592768*k*ii/n)

    c[k]=np.sqrt(a[k]**2+b[k]**2)

#---------------------------------
# Do ifft
#---------------------------------
yy=np.zeros(n,dtype=float)

for ii in np.arange(n):
    yy[ii]=a[0]/2.
    for k in np.arange(m):
        yy[ii]=yy[ii]+a[k]*np.cos(2*3.141592768*k*ii/n)+b[k]*np.sin(2*3.141592768*k*ii/n)

#---------------------------------
# Begin plot 
#---------------------------------
plt.subplot(3,1,1)
plt.plot(x,y)
plt.subplot(3,1,2)
plt.plot(cx,c)
plt.subplot(3,1,3)
plt.plot(x,yy)
plt.show()
