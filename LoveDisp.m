function f = love(c,omga,dm,betam,roum)

	mium = roum. * betam^2;
	n = length(dm);
	k = omga / c;
	for ii = 1:n
		if(c>betam(ii))
			rbetam(ii) = sqrt((c/betam(ii))^2-1);
		else
			rbetam(ii) = -i * sqrt(1-(c/betam(ii))^2);
		end

	Qm = k * rbetam * dm;
	A = [1,0;0,1]
	for ii = n-1:-1:1
		am = [cos(Qm(ii)),i*sin(Qm(ii))/rbetam(ii)/mium(ii);i*mium(ii)*rbetam(ii)*sin(Qm(ii)),cos(Qm(ii))];
	    A = A * am;
	end

	f = (A(2,1) + mium(n) * rbetam(n) * A(1,1));

return

function rt = findroot(h,prea,omga,dm,betam,roum)

	x = prea;
	F0 = love(x,omga,dm,betam,roum);
	x = x + h;

	while(F1/F0>0.0)
		F0 = F1;
		if( x - 10.0 >= 0.0)
			'Cannot find the answer';
		else 
			x = x + h;
			F1 = love(x,omga,dm,betam,roum);
		end
	end

	i = 0;
	t1 = x - h;
	t2 = x;

	while (i <= 100)

		y = t2 - love(t2,omga,dm,betam,roum)/(love(t2,omga,dm,betam,roum)-love(t1,omga,dm,betam,roum))*(t2-t1);

		if abs(y - t2) > 10^(-6)
			t1 = t2;
			t2 = y;
		else 
			break;
		end

		i = i + 1;
	end
	
	rt = y;
return

