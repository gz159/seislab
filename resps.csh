#!/bin/csh -f
sac<<FIN
fg impu n 16000
trans from none to polezero s $1
fft
wsp displacement_response
divomega
wsp velocity_response
r displacement_response.am velocity_response.am
xlim .001 .5
grid on
qdp off
loglog
ppk
q
FIN
