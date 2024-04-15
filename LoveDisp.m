# This matlab script use to determine love wave phase and group
# velocity by solve Haskell propagation matrix 

#c     Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
#c     mode computations, in  Methods in Computational Physics, 
#c         Volume 11,
#c     Seismology: Surface Waves and Earth Oscillations,  
#c         B. A. Bolt (ed),
#c     Academic Press, New York
#c
#c     Love Wave Equations  (1)-(8), pp 90-92

#------------------------------------------------------

function f = love(c,omga,dm,betam,roum)

mium = roum * betam^2;
n = length(dm);

k = omga / c;

for ii = 1 : n
   if (c > betam(ii))
   	   rbetam(ii)=sqrt((c / betam(ii))^2 - 1);
   else
   	   rbetam(ii)=-i*sqrt(1 - (c / betam(ii))^2);
end

Qm = k * rbetam * dm;

A = [1,0;0,1];

for ii = n-1 : -1 : 1
    am = [cos(Qm(ii)),i*sin(Qm(ii))/rbetam(ii)/mium(ii);
    	  i*mium(ii)*rbetam(ii)*sin(Qm(ii)),cos(Qm(ii))];
    A = A * am;
end

f = (A(2,1) + mium(n) * rbetam(n) * A(1,1));

return


#---------------------------------------------------------

function rt = findroot(h,prea,omga,dm,betam,roum)

x = prea;

F0 = love(x,omga,dm,betam,roum);

x = x + h;

F1 = love(x,omga,dm,betam,roum);

while(F1 / F0 > 0.0)
   F0 = F1;
   if(x - 10.0 >= 0)
   	   'Can not find root';
   else
       x = x + h;
       F1 = love(x,omga,dm,betam,roum);
   end
end

i = 0;

t1 = x - h;
t2 = x;

while(i < 100)
   y = t2 - love(t2,omga,dm,betam,roum) / (love(t2,omga,dm,betam,roum) - love(t1,omga,dm,betam,roum)) * (t2 - t1);
   if(abs(y - t2) > 10^(-6))
       t1 = t2;
       t2 = y;
   else
   	   break;
   end
   i = i + 1;
end

rt = y;

return

#-----------------------------------------------------

dm = [35,1000];
betam = [3.5,4.5];
roum = [2.7,3.3];

T = [5.0,6.0,7.0,8.0];

n = length(T);
vph = zeros(1,n);
vgrp = vph;
prea = min(betam);

for ii = 1 : n
	omga = 2 * pi / T(ii);
    vph(ii) = findroot(0.05,prea,omga,dm,betam,roum);
    vpha = findroot(0.05,prea,2 * pi(T(ii) - 0.01),dm,betam,roum);
    vphb = findroot(0.05,prea,2 * pi(T(ii) + 0.01),dm,betam,roum);
    dcdt = (vphb - vpha) / 0.02;
    vgrp(ii) = vph(ii) / (1.0 + T(ii) / vph(ii)*dcdt);
end

[T',vph', vgrp'];

plot(T,vph,T,vgrp,':');
