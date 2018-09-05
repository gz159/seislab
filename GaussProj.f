!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!//////6度带宽 54年北京坐标系 
!//高斯投影由经纬度(Unit:DD)反算大地坐标(含带号，Unit:Metres) 


      subroutine GaussProjCal(longitude,latitude,longitude00,X,Y) 

	  real*4 longitude,latitude,longitude00,X,Y
	  integer ProjNo
	  integer ZoneWide !////带宽 
	  real*4 longitude1,latitude1, longitude0,latitude0, X0,Y0, xval,yval
	  real*4 a,f, e2,ee, NN, T,C,A1, M, iPI
	  ProjNo=0
	  iPI = 0.0174532925199433  !////3.1415926535898/180.0; 
	  ZoneWide = 6              !////6度带宽 
	  a=6378245.0
	  f=1.0/298.3               !//54年北京坐标系参数 
                              !////a=6378140.0; f=1/298.257; //80年西安坐标系参数 
!	  write(*,*)'f=',f
	  ProjNo = int(longitude / ZoneWide)  
!	  write(*,*)'ProjNo=',ProjNo
!	  longitude0 = ProjNo * ZoneWide + ZoneWide / 2; 
	  longitude0 =longitude00
	  longitude0 = longitude0 * iPI 
	  latitude0=0 
	  longitude1 = longitude * iPI  !//经度转换为弧度
	  latitude1 = latitude * iPI    !//纬度转换为弧度
	  e2=2*f-f*f;
	  ee=e2*(1.0-e2);
	  NN=a/sqrt(1.0-e2*sin(latitude1)*sin(latitude1));
	  T=tan(latitude1)*tan(latitude1);
	  C=ee*cos(latitude1)*cos(latitude1);
	  A1=(longitude1-longitude0)*cos(latitude1);
	  M=a*((1-e2/4-3*e2*e2/64-5*e2*e2*e2/256)*latitude1-(3*e2/8+3*e2*e2/32+45*e2*e2*e2/1024)*sin(2*latitude1)+(15*e2*e2/256+45*e2*e2*e2/1024)*sin(4*latitude1)-(35*e2*e2*e2/3072)*sin(6*latitude1))
	  xval = NN*(A1+(1-T+C)*A1*A1*A1/6+(5-18*T+T*T+72*C-58*ee)*A1*A1*A1*A1*A1/120);
	  yval = M+NN*tan(latitude1)*(A1*A1/2+(5-T+9*C+4*C*C)*A1*A1*A1*A1/24+(61-58*T+T*T+600*C-330*ee)*A1*A1*A1*A1*A1*A1/720)
	  X0 = 0. !1000000. *(ProjNo+1)+500000. ; 
	  Y0 = 0; 
!	  write(*,*)x0*1000.,y0*1000.
	  xval = xval+X0; yval = yval+Y0; 
	  X = xval*0.001;
	  Y = yval*0.001;
	  return
	  end


!//高斯投影由大地坐标(Unit:Metres)反算经纬度(Unit:DD)
      subroutine GaussProjInvCal(X,Y,longitude0,longitude,latitude) 

	  real*4 longitude,latitude,X,Y
	  integer ProjNo
	  integer ZoneWide            !////带宽 
	  real*4  longitude1,latitude1, longitude0,latitude0, X0,Y0, xval,yval
	  real*4  e1,e2,f,a, ee, NN, T,C, M, D,R,u,fai, iPI;
	  iPI = 0.0174532925199433    !////3.1415926535898/180.0; 
	  a = 6378245.0
	  f  = 1.0/298.3               !//54年北京坐标系参数 
                                !////a=6378140.0; f=1/298.257; //80年西安坐标系参数 
	  ZoneWide = 6                !////6度带宽 
	  ProjNo = int(X/1000000.)   !//查找带号
!	  longitude0 = (ProjNo-1) * ZoneWide + ZoneWide / 2 
!	  longitude0 =125.
	  longitude0 = longitude0 * iPI !//中央经线
	  X0 = 0 !ProjNo*1000000+500000 
	  Y0 = 0; 
	  X=X*1000.;Y=Y*1000.
	  xval = X-X0; yval = Y-Y0;   !//带内大地坐标
	  e2 = 2*f-f*f;
	  e1 = (1.0-sqrt(1-e2))/(1.0+sqrt(1-e2));
	  ee = e2/(1-e2);
	  M = yval;
	  u = M/(a*(1-e2/4-3*e2*e2/64-5*e2*e2*e2/256));
	  fai = u+(3*e1/2-27*e1*e1*e1/32)*sin(2*u)+(21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*u)+(151*e1*e1*e1/96)*sin(6*u)+(1097*e1*e1*e1*e1/512)*sin(8*u)
	  C = ee*cos(fai)*cos(fai);
	  T = tan(fai)*tan(fai);
	  NN = a/sqrt(1.0-e2*sin(fai)*sin(fai));
	  R = a*(1-e2)/sqrt((1-e2*sin(fai)*sin(fai))*(1-e2*sin(fai)*sin(fai))*(1-e2*sin(fai)*sin(fai)))
	  D = xval/NN;
!//计算经度(Longitude) 纬度(Latitude)
	  longitude1 = longitude0+(D-(1+2*T+C)*D*D*D/6+(5-2*C+28*T-3*C*C+8*ee+24*T*T)*D*D*D*D*D/120)/cos(fai);
	  latitude1 = fai -(NN*tan(fai)/R)*(D*D/2-(5+3*T+10*C-4*C*C-9*ee)*D*D*D*D/24+(61+90*T+298*C+45*T*T-256*ee-3*C*C)*D*D*D*D*D*D/720); 
!//转换为度 DD
	  longitude = longitude1 / iPI; 
	  latitude = latitude1 / iPI;
	  return
	  end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
