from numpy import zeros,pi,sin,cos,round,arange,empty
from subprocess import call
from sys import stdout, stderr
from obspy import read, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np

def conv_isosc(trace, hd, ts):
    x   = trace.data
    t2  = trace.copy()
    y   = t2.data
    dt  = trace.stats['delta']
    n   = len(x)
    nh  = int(round(hd/dt))
    ns  = int(round(ts/dt))
    
    
    # Note the stf divided by hd and 
    # tr multiplied by dt are used to
    # make sure the area of STF is unit
    
    # Generate source time function
    stf = np.zeros(2*nh+1)
    for k in range(2*nh+1):
        if (0 <= k <= nh) or (nh <k< 2*nh):
            stf[k]=1-abs(float(k-nh)/nh)
    stf/=hd
    
    # tr keeps the convolution of syn with STF
    tr=np.zeros(n+len(stf)-1)
    
    # start convolute syn with STF
    for i in range(n):
        for j in range(len(stf)):
            tr[i+j] += x[i]*stf[j]
            
    tr *= dt
    
    
    tri =  arange(float(nh+1))/float(nh)**2
    for i in range(n):
        y[i] = tri[-1]*x[i]
        for j in range(1,min(nh,i+1)): y[i] += tri[-1-j]*x[i-j]  # left
        for j in range(1,min(nh,n-i)): y[i] += tri[-1-j]*x[i+j]  # right
    t2.stats['starttime'] += ts
    plt.plot(arange(n+len(stf)-1),tr,color='r')
    plt.plot(arange(n)+ts,y,color='b')
    plt.show()
    return t2

def cstf(ts,td,dt):
    
    # sometime the ts less than td, so we need to choose the larger one as the length of stf
    stflen = 2*max(ts,td)
    
    stlen  = int(stflen/float(dt))
    
    # for the convience of ploting stf, we define t
    t     = np.linspace(-stflen,stflen,2*stlen)
    stf   = np.zeros(2*stlen)
    
    tst   = ts/float(dt)
    tdt   = td/float(dt)
     
    for j in np.arange(2*stlen):
        if (tst-tdt+stlen) <= j <= (tst+stlen) or (tst+stlen) <j< (tst+tdt+stlen):
            stf[j] = 1-abs((j-tst-stlen)/(tdt))
    stf /= td
    
    # Trim the STF, only use the part start with amplitude greater than zero
    # and following
    if (ts-td) < 0:
        stf = stf[-(stlen-int((ts-td)/float(dt))):-1]
        t   = t[-(stlen-int((ts-td)/float(dt))):-1]
    else:
        stf = stf[stlen:-1]
        t   = t[stlen:-1]
    
    plt.plot(t,stf)
    plt.show()
    
    return t, stf

def sdr2mom(strike, dip, rake, M0=1.):
    phi    = strike*pi/180.
    delta  =    dip*pi/180.
    lambd  =   rake*pi/180.
    n    = zeros(3)
    n[0] = -sin(delta)*sin(phi)
    n[1] =  sin(delta)*cos(phi)
    n[2] = -cos(delta)
    s = zeros(3)
    s[0] =             cos(lambd)*cos(phi) + cos(delta)*sin(lambd)*sin(phi)
    s[1] =             cos(lambd)*sin(phi) - cos(delta)*sin(lambd)*cos(phi)
    s[2] = -sin(delta)*sin(lambd)
    M = zeros((3,3))
    for j in range(3):
        M[j][j] = 2.*n[j]*s[j]
        for k in range(j+1, 3):
            M[j][k] = n[j]*s[k]+n[k]*s[j]
            M[k][j] = M[j][k]
    H = zeros(6)
    H[0] =  M0*M[2][2]
    H[1] =  M0*M[0][0]
    H[2] =  M0*M[1][1]
    H[3] =  M0*M[0][2]
    H[4] = -M0*M[1][2]
    H[5] = -M0*M[0][1]
    return H

def bring_GF(model,evla,evlo,evdpmts,knetwk,kstnm,stla,stlo):
    url = 'http://service.iris.edu/irisws/syngine/1/query'
    labels = ['rr','tt','pp','rt','rp','tp']
    formsac = '%s.%s.%s.%s.sac'
    for j,label in enumerate(labels):
        mom = zeros(6)
        mom[j]=1.e21
        pref = 'M%s'%label
        req = pref+'.req'
        zip = pref+'.zip'
        fd = open(req,'w')
        fd.write('model=%s\n'%model)
        fd.write('sourcelatitude=%010.4f\n'%evla)
        fd.write('sourcelongitude=%010.4f\n'%evlo)
        fd.write('sourcedepthinmeters=%d\n'%int(evdpmts))
        fd.write('sourcemomenttensor=')
        for j in range(6):
            fd.write('%+7.2e'%mom[j])
            if j != 5: fd.write(',')
        fd.write('\n')
        line = '%10.5f %10.5f STACODE=%s NETCODE=%s LOCCODE=%s'%(stla,stlo,kstnm,knetwk,label)
        line = line.strip()
        fd.write(line+'\n')
        fd.close()
        cmd = 'wget --post-file=%s -O %s %s'%(req,zip,url)
        call(cmd, shell=True, stdout=stdout, stderr=stderr)
        call('unzip -o %s'%zip, shell=True, stdout=stdout, stderr=stderr)
        for c in 'ZNE':
            sac = formsac%(knetwk,kstnm,label,'BX%c'%c)
            gf  = read(sac)[0]
            gf.decimate(5)
            gf.write(sac,'sac')
    return

def calc_syn(knetwk,kstnm,strike, dip, rake, mw, hd, ts, OT):
    labels = ['rr','tt','pp','rt','rp','tp']
    M0nm = pow(10., 1.5*mw + 9.1)
    H    = sdr2mom(strike, dip, rake, M0nm)
    formsac = '%s.%s.%s.%s.sac'
    t0 = read(formsac%(knetwk,kstnm,'rr','BXZ'))[0]
    t0.data *= 0.
    for c in 'ZNE':
        kcmpnm = 'BX%c'%c
        t = t0.copy()
        for j in range(6):
            gf = read(formsac%(knetwk,kstnm,labels[j],kcmpnm))[0]
            t.data += H[j]*gf.data/1.e21
        t.stats['starttime'] = OT
        del t.stats['sac'] 
        #t.data *= 0.; t.data[int(len(t.data)/2)] = 1./t.stats['delta']
        #t.write('before_%s.%s.%s.sac'%(knetwk,kstnm,kcmpnm),'sac')
        t2 = conv_isosc(t, hd, ts)
        t2.write('%s.%s.%s.sac'%(knetwk,kstnm,kcmpnm),'sac')
    return
    
# main
model = 'ak135f_2s'
# Source
OT = UTCDateTime(1990,12,16,12,5,48,0)
evla, evlo, evdpmts = 36.65,105.,10000.
strike, dip, rake   = 130., 85., 5.
mw, hd, ts          = 7.9, 10., 20.

# Receiver
knetwk = 'XX'
stats = (('JENA', 50.9519, 11.5833),
         ('DBN',  52.100,   5.1833))

cstf(10,10,0.1)

for stat in stats:
    kstnm,stla,stlo = stat 
    #bring_GF(model,evla,evlo,evdpmts,knetwk,kstnm,stla,stlo)
    calc_syn(knetwk,kstnm,strike,dip,rake,mw,hd,ts,OT)
