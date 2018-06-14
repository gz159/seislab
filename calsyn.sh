#!/bin/sh

# this script compute and compare synthetic seismogram generate
# by using Herrman CPS, Lupei Zhu and IRIS Axism program


cd herrmann
/bin/rm -rf 0450 B00101Z00.sac  B00102N00.sac B00103E00.sac date.end date.start dfile file96 hspec96.dat hspec96.grn
DOIT
hpulse96 -t -V -l 1 | fmech96 -MW 7.9 -D 45 -S 315 -R 90 -A 0 -B 180 | f96tosac

cd ../zhulupei
/bin/rm -rf IRISak135_lupei_45 junk.p junk.s S30.r S30.t S30.z
fk.pl -MIRISak135_lupei/45 -N4096/1 -S2 9990
syn -M7.9/315/45/90 -D0.5 -A0 -OS30.z -GIRISak135_lupei_45/9990.grn.0

cd ..
wget -O synM7.9.zip "http://service.iris.edu/irisws/syngine/1/query?model=ak135f_2s&format=saczip&components=ZNE&units=velocity&dt=0.2&receiverlatitude=90&receiverlongitude=0&sourcelatitude=0&sourcelongitude=0&sourcedepthinmeters=45000&sourcedoublecouple=315,45,90,8.9e20"

unzip synM7.9.zip
