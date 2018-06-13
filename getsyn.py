#!/usr/bin/python

from obspy.clients.syngine import Client
from obspy import read
from obspy.taup import TauPyModel
from obspy.core import UTCDateTime
from subprocess import call
from sys import stdout, stderr
import numpy as np


client = Client()

##------------------------------------------------------------------------
# we use uniformly distributed random earthquake focal mechanisms to generate
# synthetic seismogram
# The IRIS syngine Green function database have been used
# http://service.iris.edu/irisws/syngine/1/
##-------------------------------------------------------------------------
for i in range(1,5):
    # the source depth range from 10km to 70km
    evdpinmeter=int(np.random.uniform(10,70)*1000)
    strike=int(np.random.uniform(0.1,1)*360)
    rake=int(np.random.uniform(0.1,1)*180)
    dip=int(np.random.uniform(0.1,1)*90)
    # epicenter distance from 30 degree to 90 degree
    dist=np.random.uniform(30,90)
    # the moment magnitude from 5.5 to 7.5
    mw=np.random.uniform(5.5,7.5)
    m0=np.power(10.,1.5*mw+9.1)
    model = TauPyModel(model="ak135")
    arrt=model.get_travel_times(source_depth_in_km=evdpinmeter/1000.,distance_in_degree=dist,phase_list=["P"])
    at=divmod(arrt[0].time,60.)
    orgtime=UTCDateTime(1970, 1, 1, 1, int(at[0]), int(at[1]))
    #print(arrt)
    st = client.get_waveforms(model="ak135f_1s", origintime=orgtime, receiverlatitude=0.0, receiverlongitude=dist, sourcelatitude=0.0, sourcelongitude=0.0, sourcedepthinmeters=evdpinmeter, sourcedoublecouple=[strike, dip, rake, m0], starttime="P-30", endtime="P+150", components="Z", units="displacement", dt=0.01)
    fn="rdSyn_"+str(i)+".bhz.sac"
    print("Generated synthetic seismogram: "+str(i))
    st[0].write(fn,format='sac')
    
    cmd='''
    sac<<FIN
    read %s
    ch o -%s
    ch t0 %s
    ch gcarc %s
    ch stla %s
    ch stlo %s
    ch evla %s
    ch evlo %s
    w over
    quit
    FIN'''%(fn,arrt[0].time,arrt[0].time,dist,0.0,dist,0.0,0.0)
    call(cmd,shell=True,stdout=stdout,stderr=stderr)

