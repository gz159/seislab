#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 09:45:48 2017
Geophys. J. Int. (2014),1719-1735
Receiver function deconvolution using transdimensional hierarchical Bayesian inference
@author: seis
"""
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


import numpy as np
import matplotlib.pylab as plt

lamda=0.2
ommiga=4.4
size=50
rij1=np.empty((size,size))
rij2=np.empty((size,size))
rij3=np.empty((size,size))
for i in range(size):
    for j in range(size):
        rij1[i,j]=np.exp(-1.*lamda*np.abs(j-i))
        rij2[i,j]=np.exp(-1.*lamda*lamda*np.abs(j-i))
        rij3[i,j]=np.exp(-1.*lamda*np.abs(j-i))*np.cos(lamda*ommiga*np.abs(j-i))
        
#plt.figure(figsize=[12,3])
fig,ax=plt.subplots(1,3,figsize=(12,3))
#plt.gca().invert_yaxis()
#plt.subplot(1,3,1)
ax[0].set_title(r'$R_{ij}=e^{-\lambda|t_j-t_i|}$')
ax[0].set_ylim(size,0)
ax[0].set_xlim(0,size)
ax[0].set_ylabel(r'$Time_i(s)$')
ax[0].set_xlabel(r'$Time_j(s)$')
pcm=ax[0].pcolor(rij1,cmap='gist_rainbow',vmin=-1,vmax=1)
#colorbar(cb)
fig.colorbar(pcm,ax=ax[0],extend='both')

#plt.subplot(1,3,2)
ax[1].set_title(r'$R_{ij}=e^{-\lambda^2|t_j-t_i|}$')
ax[1].set_ylim(size,0)
ax[1].set_xlim(0,size)
pcm=ax[1].pcolor(rij2,cmap='gist_rainbow',vmin=-1,vmax=1)
#colorbar(cb)
fig.colorbar(pcm,ax=ax[1],extend='both')

#plt.subplot(1,3,3)
ax[2].set_title(r'$R_{ij}=e^{-\lambda|t_j-t_i|}\cos(\lambda\omega_0|t_j-t_i|)$')
ax[2].set_ylim(size,0)
ax[2].set_xlim(0,size)
pcm=ax[2].pcolor(rij3,cmap='gist_rainbow',vmin=-1,vmax=1)
#colorbar(cb)
fig.colorbar(pcm,ax=ax[2],extend='both')
plt.show()
