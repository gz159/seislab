import matplotlib.pyplot as plt
import numpy as np


# Generate a triangle source time function of time shift ts and half duration of hd
def cstf(ts,hd,dt):
    
    # sometime the ts less than hd, so we need to choose the larger one as the length of stf
    stflen = 2*max(ts,hd)
    
    stlen  = int(stflen/float(dt))
    
    # for the convience of ploting stf, we define t
    t     = np.linspace(-stflen,stflen,2*stlen)
    stf   = np.zeros(2*stlen)
    
    tst   = ts/float(dt)
    hdt   = hd/float(dt)
     
    for j in np.arange(2*stlen):
        if (tst-hdt+stlen) <= j <= (tst+stlen) or (tst+stlen) <j< (tst+hdt+stlen):
            stf[j] = 1-abs((j-tst-stlen)/(hdt))
    stf /= hd
    
    # Trim the STF, only use the part start with amplitude greater than zero
    # and following
    if (ts-hd) < 0:
        stf = stf[-(stlen-int((ts-hd)/float(dt))):-1]
        t   = t[-(stlen-int((ts-hd)/float(dt))):-1]
    else:
        stf = stf[stlen:-1]
        t   = t[stlen:-1]
    #plt.plot(t,stf)
    #plt.show()
    return t, stf

#----------------------------------------------------
# here we use the full mode of convolve, which means
# the total lenght of convolution result is 
# (lenght of syn + length stf - 1)
# syn and stf need to have same dt, samplingrate
#-----------------------------------------------------
def conv(syn,stf,dt):

    lensyn = len(syn)
    lenstf = len(stf)
    
    tr = np.zeros(lensyn + lenstf -1)
    print(len(tr))
    for m in np.arange(lensyn):
        for n in np.arange(lenstf):
            tr[m+n] = tr[m+n] + syn[m]*stf[n]

    data = tr*float(dt)

    t = np.linspace(0,lensyn + lenstf -1,lensyn + lenstf -1)
    
    return t, data
##--------------------------------------------------

ts=10
hd=10
dt=0.1

t,stf=cstf(ts,hd,dt)
plt.plot(t,stf)
plt.show()

t,data=conv(stf,stf,dt)
plt.plot(t,data)
plt.show()
    
