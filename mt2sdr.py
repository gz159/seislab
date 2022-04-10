import numpy as np
import math

#----------------------------------------------------------------------------------
# According to moment tensor theory, the moment tensor Mij can be represent by it
# eigen value mij and eigen vector ei
#
#             |m11  m12  m13|       |e1|       |e|        | |
#       Mij = |m21  m22  m23| = m11*|e2| + m22*| | + m33* | |
#             |m31  m32  m33|       |e3|       | |        | |
#
#     M = m11e1 + m22e2 + m33e3  and here we defined m11 >= m22 >= m33
#
# In this order, e1 vector repreent T axis, 
#                e2 vector reprent null axis 
#                e3 vector repreent P axis
#-----------------------------------------------------------------------------------

#         mtrep - moment tensor representation
#                 = 0 for spherical coordinates (r=upper,t=south,p=east)
#                   ( Up, South, East )
#                 = 1 for cartesian coordinates (x=north,y=east,z=down)
#                   ( North, East, Down )
#                 = 2 for f1,...,f6 notation
#
#         mf or mt - components of moment tensor
#                 | mrr mrt mrp |    |  mxx -mxy  mxz |    | f1 f2 f3 |
#                 | mtr mtt mtp | == | -myx  myy -myz | or | f2 f4 f5 |
#                 | mpr mpt mpp |    |  mzx -mzy  mzz |    | f3 f5 f6 |
#
#              NOTE:
#                   f1 = mrr =  mxx
#                   f2 = mtt =  myy
#                   f3 = mpp =  mzz
#                   f4 = mrt = -mxy
#                   f5 = mrp =  mxz
#                   f6 = mtp = -myz



# find fault plane's strike, dip, rake angles using fault normal (vn) and slip vector (vs)
def vn2sdr(vn,vs):

	epsi = 0.001

	if vn[0] < 0.: # Upwards normal
		vn = -1.*vn
		vs = -1.*vs

	if vn[0] > 1.0 - epsi: # Horizontal plane
		strike = 0
		dip = 0
		rake = np.rad2deg(np.arctan2(-vs[2],-vs[1]))
	elif vn[0] < epsi: # Vertical plane
		strike = np.rad2deg(np.arctan2(vn[1],vn[2]))
		dip = np.rad2deg(np.pi/2.)
		rake =  np.rad2deg(np.arctan2(vs[0],-vs[1]*vn[2]+vs[2]*vn[1]))
	else:  # Oblique plane
		strike = np.rad2deg(np.arctan2(vn[1],vn[2]))
		dip = np.rad2deg(np.arccos(vn[0]))
		rake = np.rad2deg(np.arctan2((-vs[1]*vn[1]-vs[2]*vn[2]),(-vs[1]*vn[2]+vs[2]*vn[1])*vn[0]))
    
	if strike < 0:
		strike = strike + 360.

	if rake < -180:
		rake = rake + 360.
	if rake > 180:
		rake = rake - 360.

	return strike,dip,rake

#################------------------------------------------------------
# using global CMT result as an input example
#        Mrr     Mtt     Mpp     Mrt     Mrp     Mtp
# CMT   2.180  -2.170  -0.018   1.260  -0.044   0.985
# Fault plane:  strike=303    dip=33   slip=108
# Fault plane:  strike=101    dip=59   slip=79
# Eigenvector:  eigenvalue:  2.54   plunge: 73   azimuth: 342
# Eigenvector:  eigenvalue:  0.29   plunge: 10   azimuth: 107
# Eigenvector:  eigenvalue: -2.84   plunge: 13   azimuth: 199
#
#          | Mrr Mrt Mrp |
#      m = | Mtr Mtt Mtp |
#          | Mpr Mpt Mpp |
#
#################------------------------------------------------------


cmt=[2.180,-2.170,-0.018,1.260,-0.044,0.985]
#      Mrr    Mtt   Mpp   Mrt    Mrp    Mtp
# cmt=[-0.917,-0.028,0.946,0.069,-0.004,-0.025]


mt  = np.array([
	[cmt[0],cmt[3],cmt[4]],
	[cmt[3],cmt[1],cmt[5]],
	[cmt[4],cmt[5],cmt[2]]])

m,v = np.linalg.eig(mt)
#idx = m.argsort()[::-1]  # sort m in descending order 
idx = np.argsort(m) # sort eigen value m in ascending order
m = m[idx] # sort eigen value m in ascending order (m1<m2<m3)
# After sort m in ascening order,
# the first element of m correspond to P axis, 
# the second element correspond to null axis, 
# the third element correspond to T axis
v = v[:,idx] # sort eigen vector according to eigen value

# Now first colum of matrix v correspond to vector of P axis, the second colum of matrix v
# correspond to vector of null axis and the third colum of v correspond to the vector of T axis
# calculate P, null and T axis plunge and azimuth angle using eigen vector
for i in range(3):
	azm = np.rad2deg(np.arctan2(v[2][i],-v[1][i]))
	scale = v[1][i]*v[1][i]+v[2][i]*v[2][i]
	plg = np.rad2deg(np.arctan2(-v[0][i],np.sqrt(scale)))
	if plg < 0.0:
		plg = -1.0 * plg
		azm = azm + 180.

	azm=math.fmod(azm,360.)

	if azm < 0:
		azm = azm + 360.

	print("eigval={0:.2f} plunge={1:.1f} azimuth={2:.1f}".format(m[i],plg,azm))



# construct eigen vector correspond P, null and T axis
p = v[:,0]
n = v[:,1]
t = v[:,2]

# calculate fault normal and slip vector
v1 = (t+p)/np.sqrt(2.)
v2 = (t-p)/np.sqrt(2.)


# focal parameter for fault plane I
strike,dip,rake=vn2sdr(v1,v2)

print("strike={0:.1f} dip={1:.1f} rake={2:.1f}".format(strike,dip,rake))


# focal parameter for fault plane II
strike,dip,rake=vn2sdr(v2,v1)

print("strike={0:.1f} dip={1:.1f} rake={2:.1f}".format(strike,dip,rake))
