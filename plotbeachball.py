import matplotlib.pyplot as plt
import numpy as np

# beachball plot of earthquake focal mechanism
# the beachball is circle of radius 1


### calculate moment tensor
# use formula from SRL(1989),60(2) paper "A student guide and reviewer of moment tensors"
def dislocationtomt(strike,dip,rake,m0):

   m = np.empty((3,3))
   phi      = (strike*np.pi)/180.
   twophi   = phi*2.
   delta    = (dip*np.pi)/180.
   twodelta = delta*2.
   lamda    = (rake*np.pi)/180.

   mxx = -(np.sin(delta)*np.cos(lamda)*np.sin(twophi)
           +np.sin(twodelta)*np.sin(lamda)*np.sin(phi)*np.sin(phi))

   myy =  (np.sin(delta)*np.cos(lamda)*np.sin(twophi)
           -np.sin(twodelta)*np.sin(lamda)*np.cos(phi)*np.cos(phi))

   mzz =  (np.sin(twodelta)*np.sin(lamda))

   mxy =  (np.sin(delta)*np.cos(lamda)*np.cos(twophi)
           +0.5*np.sin(twodelta)*np.sin(lamda)*np.sin(twophi))

   mxz = -(np.cos(delta)*np.cos(lamda)*np.cos(phi)
           +np.cos(twodelta)*np.sin(lamda)*np.sin(phi))

   myz = -(np.cos(delta)*np.cos(lamda)*np.sin(phi)
           -np.cos(twodelta)*np.sin(lamda)*np.cos(phi))

   mrr =  mzz*m0
   mrt =  mxz*m0
   mrp = -myz*m0
   mtt =  mxx*m0
   mtp = -mxy*m0
   mpp =  myy*m0

   m[0][0] = mxx
   m[0][1] = mxy
   m[0][2] = mxz
   m[1][1] = myy
   m[1][2] = myz
   m[2][2] = mzz
   m[1][0] = m[0][1]
   m[2][0] = m[0][2]
   m[2][1] = m[1][2]

   return m
   #return mrr,mtt,mpp,mrt,mrp,mtp

   # using focal mechansim parameter strike dip rake to
   # get its stereographic equal area projction plane
def pnodal(strike,dip,rake):
    strike = (strike*np.pi)/180.
    dip    = (dip*np.pi)/180.
    rake   = (rake*np.pi)/180.
    px = []
    py = []
    #px.append(np.cos(strike))
    #py.append(np.sin(strike))
    for i in range(0,180):
        ii = i*np.pi/180
        x = np.cos(strike)*np.cos(ii)-np.sin(strike)*np.sin(ii)*np.cos(dip)
        y = np.sin(strike)*np.cos(ii)+np.cos(strike)*np.sin(ii)*np.cos(dip)
        z = np.sin(ii)*np.sin(dip)
        azimuth = np.atan2(y,x)
        ain = np.atan2(np.sqrt(1-z**2),z)
        r = np.sqrt(2)*np.sin(ain/2.0)
        px.append(r*np.sin(azimuth))
        py.append(r*np.cos(azimuth))

    return px,py

# GETAUX returns auxilary fault plane strike, dip & rake,
# given strike,dip,rake of main fault plane.
def getaux(strike1, dip1, rake1):

   degrad = 180./3.1415927
   s1 = strike1/degrad
   d1 = dip1/degrad
   r1 = rake1/degrad

   d2 = np.arccos(np.sin(r1)*np.sin(d1))

   sr2 = np.cos(d1)/np.sin(d2)
   cr2 = -np.sin(d1)*np.cos(r1)/np.sin(d2)
   r2 = np.arctan2(sr2, cr2)

   s12 = np.cos(r1)/np.sin(d2)
   c12 = -1./(np.tan(d1)*np.tan(d2))
   s2 = s1 - np.arctan2(s12, c12)

   strike2 = s2*degrad
   dip2 = d2*degrad
   rake2 = r2*degrad

   if (dip2 > 90.):
      strike2 = strike2 + 180.
      dip2 = 180. - dip2
      rake2 = 360. - rake2
   if (strike2 > 360.): strike2 = strike2 - 360.

   return strike2, dip2, rake2

######################## main program ####################

# Focal mechanism represents with strike dpi, rake, m0
strike, dip, rake = 203.0, 10.0, 88.0
strike1, dip1, rake1 = getaux(strike, dip, rake)

mt = dislocationtomt(strike,dip,rake,1)

plt.figure(figsize=(6, 6))
r = [0.0,0.0,0.0]
aP = mt[0][0]
aT = mt[0][0]
jxP, jyP = 0.0, 0.0
jxT, jyT = 0.0, 0.0
for i in range(-20,20):
    x = i/20.0
    ymax = np.sqrt(1-x**2)
    m = int(20.0*ymax)
    #print(m)
    if m > 0:
        #print(m)
        for j in range(-m,m):
            y    = j/20.0
            rad  = np.sqrt(x**2+y**2)
            azi  = np.atan2(x,y)
            ain  = 2.0*np.asin(rad/np.sqrt(2.0))
            r[0] = np.sin(ain)*np.cos(azi)
            r[1] = np.sin(ain)*np.sin(azi)
            r[2] = np.cos(ain)
            rpp  = 0.0
            for ii in range(3):
                for jj in range(3):
                    rpp += r[ii]*mt[ii][jj]*r[jj]
            if rpp > 0.0:
                plt.scatter(x, y, s=rpp*55, c='red', edgecolors='white', linewidths=1)
                #plt.text(x, y, str(rpp), fontsize=8, ha='center', va='center')
            if rpp < 0.0:
                plt.scatter(x, y, s=abs(rpp*55), c='blue', edgecolors='white', linewidths=1)
                #plt.text(x, y, str(rpp), fontsize=8, ha='center', va='center')
            if rpp > -0.033 and rpp < 0.033:
                plt.scatter(x, y, s=45, c='black')
            if rpp < aP:
                    aP = rpp
                    jxP = x
                    jyP = y
            if rpp > aT:
                    aT = rpp
                    jxT = x
                    jyT = y

plt.scatter(jxP, jyP, s=abs(aP*55), marker='^', c='blue')
plt.text(jxP, jyP, str(aP), fontsize=8, ha='center', va='center')
plt.scatter(jxT, jyT, s=abs(aT*55), marker='^', c='green')
plt.text(jxT, jyT, str(aT), fontsize=8, ha='center', va='center')

circle = plt.Circle((0,0),1,color='black',fill=False)
plt.gca().add_patch(circle)

x,y=pnodal(strike,dip,rake)
plt.plot(x,y,color='green',linewidth=2)

x,y=pnodal(strike1,dip1,rake1)
plt.plot(x,y,color='purple',linewidth=2)

plt.xlim(-1, 1)
plt.ylim(-1, 1)
title='strike=%d,dip=%d,rake=%d'%(strike,dip,rake)
plt.title(title)

# 显示网格
#plt.grid(True)

plt.show()
