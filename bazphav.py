#!/usr/bin/env python

import matplotlib.pyplot as plt

from obspy.taup.taup_geo import calc_dist_azi

file = open('paths10','r')
phav = []
baz = []
i = 0

for line in file:
    if line.strip():
        par = line.split()
        
        if not par[0].startswith("#"):
            phav.append('0.0e00')
            baz.append('0.0e00')
        
file.close()

file = open('paths10','r')
for line in file:
    if line.strip():
        par = line.split()
        
        if not par[0].startswith("#"):
            stlat = float(par[0])
            stlon = float(par[1])
            evlat = float(par[2])
            evlon = float(par[3])
            
            phav[i] = float(par[4])
            
            result=calc_dist_azi(source_latitude_in_deg=stlat,source_longitude_in_deg=stlon,receiver_latitude_in_deg=evlat,receiver_longitude_in_deg=evlon,radius_of_planet_in_km=6378.137,flattening_of_planet=0.0033528106647474805)
            
            baz[i] = float(result[2])
            
            i = i+1


plt.plot(baz,phav,'bo')
plt.show()
