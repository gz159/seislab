/* 
 *  Strike, dip, rake of the fault plane computed from P and T axis.
 *  Angles are in degrees.
 *
 *  Genevieve Patau, 16 june 1997
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define EPSIL 0.0001
double computed_rake2(str1, dip1, str2, dip2, fault)
/* 
Compute rake in the second nodal plane when strike and dip
for first and second nodal plane are given with a double
characterizing the fault :
+1. inverse fault
-1. normal fault.
Angles are in degrees.
*/
double str1, dip1, str2, dip2, fault;
{
	double rake2, sinrake2;
	double D2R = M_PI / 180.;
	str1 *= D2R; dip1 *= D2R;
	str2 *= D2R; dip2 *= D2R;	

	if(fabs(dip2 - M_PI_2) < EPSIL)
		sinrake2 = fault * cos(dip1);
	else
		sinrake2 = -fault * sin(dip1) * cos(str1 - str2) / cos(dip2);

	rake2 = atan2(sinrake2, - fault * sin(dip1) * sin(str1 - str2));
	rake2 /= D2R;
	return(rake2);
}

void axe2dc(pp,dp,pt,dt,p1,d1,g1,p2,d2,g2)
double pp, dp, pt, dt;
double *p1, *d1, *g1;
double *p2, *d2, *g2;
{
	double D2R = M_PI / 180.;
	double PII = M_PI * 2.;
	double cdp, sdp, cdt, sdt;
	double cpt, spt, cpp, spp;
	double amz, amy, amx;
	double im;

	pp *= D2R; dp *= D2R; pt *= D2R; dt *= D2R;
	cdp = cos(dp); sdp = sin(dp);
	cdt = cos(dt); sdt = sin(dt);
	cpt = cos(pt)*cdt; spt = sin(pt)*cdt;
	cpp = cos(pp)*cdp; spp = sin(pp)*cdp;
	amz = sdt + sdp; amx = spt + spp; amy = cpt + cpp; 
	*d1 = atan2(sqrt(amx*amx + amy*amy), amz);
	*p1 = atan2(amy, -amx);

	if(*d1 > M_PI_2) {
		*d1 = M_PI - *d1;
		*p1 += M_PI;	
		if(*p1 > PII) *p1 -= PII;
	}
	if(*p1 < 0.) *p1 += PII;
	amz = sdt - sdp; amx = spt - spp; amy = cpt - cpp; 
	*d2 = atan2(sqrt(amx*amx + amy*amy), amz);
	*p2 = atan2(amy, -amx);
	if(*d2 > M_PI_2) {
		*d2 = M_PI - *d2;
		*p2 += M_PI;
	if(*p2 > PII) *p2 -= PII;
	}

	if(*p2 < 0.) *p2 += PII;
	*d1 /= D2R; *p1 /= D2R; *d2 /= D2R; *p2 /= D2R;
	im = 1;
	if(dp > dt) im = -1;
	*g1 = computed_rake2(*p2,*d2,*p1,*d1,im);
	*g2 = computed_rake2(*p1,*d1,*p2,*d2,im);
}

/*********************************************************************/
main(argc, argv)
int argc;
char *argv[];
{
	double pp, dp, pt, dt;
	double p1, d1, g1;
	double p2, d2, g2;
	if(argc == 1) {
		fprintf(stderr, "Enter strike plungement of P and T axis in degrees: ");
		scanf("%lf %lf %lf %lf", &pp, &dp, &pt, &dt);
	} else {
		pp = atof(argv[1]); dp = atof(argv[2]);
		pt = atof(argv[3]); dt = atof(argv[4]);
	}

	axe2dc(pp,dp,pt,dt,&p1,&d1,&g1,&p2,&d2,&g2);
//	fprintf(stderr,"The output is strike1 dip1 rake1 strike2 dip2 rake2\n");
	fprintf(stdout,"%.0lf %.0lf %.0lf %.0lf %.0lf %.0lf\n",p1,d1,g1,p2,d2,g2);
}

