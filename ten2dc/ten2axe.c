/*
 *  Compute Eigen values,strike,and plungement of T,B,P axis from moment tensor.
 *
 *  The component of the tensor are in the following order :
 *  	mrr mtt mff mrt mrf mtf
 *
 *  The following variables are defined
 *  	pclvd = 1 or 2 for a pure CLVD tensor
 *  	piso = 1 for a pure implosion-explosion
 *  	pdc = 1 for a pure double couple
 *  	value, strike (in degrees), plungement (in degrees) for axes T, N, P
 *  	are respectively tv, ta, td, nv, na, nd, pv, pa, pd
 *  	i characterizes isotropy (varies from 0.0 to 1.0)
 *  	f characterizes clvd (varies from 0.0 to 0.5)
 *
 *  	Genevieve Patau, 5 june 1997
 *  	translated into english 26 january 1998
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define EPSIL 0.000001

void ten2axe(mrr, mtt, mff, mrt, mrf, mtf, pclvd, piso, pdc, tv, ta, td, nv, na, nd, pv, pa, pd)
/* Input */
double mrr, mtt, mff, mrt, mrf, mtf;
/* Output */
double *tv, *ta, *td;
double *nv, *na, *nd;
double *pv, *pa, *pd;
int *pclvd, *piso, *pdc;
{
	int i;
	double b, c, d;
	double delta, x1, x2, xm, y, x0;
	double q, r;
	double m[3], az[3], pl[3];
	double z, n, e, norm;
	double R2D = 180. / M_PI;

	*piso = 0; *pdc = 0; *pclvd = 0;
	/* det(M-xI) = 0. */
	/* a = -1 */
	b = mrr + mtt + mff;
	c = mrt*mrt + mtf*mtf + mrf*mrf - mrr*mtt - mtt*mff - mff*mrr;
	d = mrr*mtt*mff + 2.*mrt*mrf*mtf - mrr*mtf*mtf - mtt*mrf*mrf - mff*mrt*mrt;

	/* solving a*x**3 + b*x**2 + c*x *d = 0 */
	/* derivee nulle => extrema secondaires */
	delta = sqrt(b*b + 3.*c);
	x1 = (b - delta) / 3.;
	x2 = (b + delta) / 3.;
	if (fabs(x1 - x2) < EPSIL) {
		*piso = 1; /* pure implosion-explosion */
		for(i=0; i<=2; i++) { m[i] = x1; }
		az[0] = 0.; az[1] = 0.; az[2] = M_PI_2;
		pl[0] = M_PI_2; pl[1] = 0.; pl[2] = 0.; 
	} else if(fabs(- x1*x1*x1 + b*x1*x1 + c*x1 +d) < EPSIL) {
		*pclvd = 1; /* pur clvd */
		/*p = a = -1. */
		/* a*x**3 + b*x**2 + c*x *d = (x - x1)*(x - x1)*(p*x + q) */
		q = b - 2.*x1;
		m[0] = x1;
		m[1] = x1;
		m[2] = q;
	} else if(fabs(- x2*x2*x2 + b*x2*x2 + c*x2 +d) < EPSIL) {
		*pclvd = 2; /* pur clvd */
		/*p = a = -1. */
		/* a*x**3 + b*x**2 + c*x *d = (x - x2)*(x - x2)*(p*x + q) */
		q = b - 2.*x2;
		m[0] = q;
		m[1] = x2;
		m[2] = x2;
	} else {
		xm = (x1 + x2) / 2.;
		y = - xm*xm*xm + b*xm*xm + c*xm + d;
		while(fabs(y) - EPSIL > 0.) {
			if(y < 0.) { x1 = (x1 + x2) / 2.; }
			else { x2 = (x1 + x2) / 2.; }
			xm = (x1 + x2) / 2.;
			y = - xm*xm*xm + b*xm*xm + c*xm + d;
		}
		/* solution intermediaire */
		x0 = xm;
		if(fabs(x0) < EPSIL) *pdc = 1; /* pure double couple */
		/* a*x**3 + b*x**2 + c*x *d = (x - x0)*(p*x**2 + q*x +r) */
		/*p = a = -1. */
		q = - x0 + b;
		r = q*x0 + c;
		delta = sqrt(q*q + 4.*r);
	
		m[0] = (q - delta) / 2.;
		m[1] = x0;
		m[2] = (q + delta) / 2.;
		*pv = m[0];
		*nv = m[1];
		*tv = m[2];
		/* eigen vectors */
		/* z = -r; n = -t; */
		for(i=0; i<=2; i=i+1+*pdc) {
			delta = (mrr-m[i])*(mtt-m[i]) - mrt*mrt;
			if(fabs(delta) > EPSIL) {
				z = (- mrt*mtf + mrf*(mtt-m[i])) / delta;
				n = (- mrt*mrf + mtf*(mrr-m[i])) / delta;
				e = 1.;
			} else {
				delta = (mrr-m[i])*(mff-m[i]) - mrf*mrf;
				if(fabs(delta) > EPSIL) {
					z = (- mrf*mtf + mrt*(mff-m[i])) / delta;
					e = (mrt*mrf - mtf*(mrr-m[i])) / delta;
					n = 1.;
				} else {
					delta = (mtt-m[i])*(mff-m[i]) - mtf*mtf;
					n = (- mrt*mtf + mrt*(mff-m[i])) / delta;
					e = (mrt*mtf - mrt*(mtt-m[i])) / delta;
					z = 1.;
				}
			}		
			norm = sqrt(z*z + n*n + e*e);
			z /= norm;
			n /= norm;
			e /= norm;
			if(fabs(n) < EPSIL && fabs(e) < EPSIL) pl[i] = 0.;
			else pl[i] = asin(z);
			az[i] = atan2(e,n);
			if(pl[i] < 0.) {pl[i] = - pl[i]; az[i] += M_PI;}
			if(az[i] < 0.) az[i] += M_PI * 2.;
		}

		if(*pdc) {
			e = cos(pl[0])*cos(az[0])*sin(pl[2]) - cos(pl[2])*cos(az[2])*sin(pl[0]);
			n = sin(pl[0])*cos(pl[2])*sin(az[2]) - sin(pl[2])*cos(pl[0])*sin(az[0]);
			z = cos(pl[0])*sin(az[0])*cos(pl[2])*cos(az[2]) - cos(pl[2])*sin(az[2])*cos(pl[0])*cos(az[0]);
			norm = sqrt(z*z + n*n + e*e);	
			z /= norm;
			n /= norm;
			e /= norm;
			pl[1] = asin(z);
			az[1] = atan2(e,n);
			if(pl[1] < 0.) {pl[1] = - pl[1]; az[1] += M_PI;}
			if(az[1] < 0.) az[1] += M_PI * 2.;
		}
	}

	if(*pclvd > 0) {
		i = 2*(2-*pclvd);
		delta = (mrr-m[i])*(mtt-m[i]) - mrt*mrt;
		if(fabs(delta) > EPSIL) {
			z = (- mrt*mtf + mrf*(mtt-m[i])) / delta;
			n = (- mrt*mrf + mtf*(mrr-m[i])) / delta;
			e = 1.;
		} else {
			delta = (mrr-m[i])*(mff-m[i]) - mrf*mrf;
			if(fabs(delta) > EPSIL) {
				z = (- mrf*mtf + mrt*(mff-m[i])) / delta;
				e = (mrt*mrf - mtf*(mrr-m[i])) / delta;
				n = 1.;
			} else {
				delta = (mtt-m[i])*(mff-m[i]) - mtf*mtf;
				n = (- mrt*mtf + mrt*(mff-m[i])) / delta;
				e = (mrt*mtf - mrt*(mtt-m[i])) / delta;
				z = 1.;
			}
		}
	
		norm = sqrt(z*z + n*n + e*e);
		z /= norm;
		n /= norm;
		e /= norm;
		if(z < 0.) {z = -z; n = -n; e = -e;}
		pl[i] = asin(z);
		az[i] = atan2(e,n);
		if(az[i] < 0.) az[i] += M_PI * 2.;
	
		if(fabs(z) > EPSIL) {
			pl[1] = asin(- n);
			az[1] = atan2(0,z);
			if(pl[1] < 0.) {pl[1] = - pl[1]; az[1] += M_PI;}
			if(az[1] < 0.) az[1] += M_PI * 2.;
			pl[2-i] = asin(- e);
			az[2-i] = atan2(z,0);
			if(pl[2-i] < 0.) {pl[2-i] = - pl[2-i]; az[2-i] += M_PI;}
			if(az[2-i] < 0.) az[2-i] += M_PI * 2.;
		} else {
			az[1] = 0.;
			pl[1] = M_PI_2;
			az[2-i] = az[i] + M_PI_2;
			if(az[2-i] > M_PI * 2.) az[2-i] -= M_PI * 2.;
			pl[2-i] = 0.;
		}
	}

	*pv = m[0]; *pa = az[0] * R2D; *pd = pl[0] * R2D;
	*nv = m[1]; *na = az[1] * R2D; *nd = pl[1] * R2D;
	*tv = m[2]; *ta = az[2] * R2D; *td = pl[2] * R2D;
}


main(argc, argv)
int argc;
char *argv[];
{
	double mrr, mtt, mff, mrt, mrf, mtf;
	double tv, nv, pv;
	double ta, na, pa;
	double td, nd, pd;
	double v_iso, tv1, nv1, pv1;
	double iso, clvd;
	double modc; /* scalar moment */
	int pclvd, piso, pdc;

	if(argc == 1) {
		fprintf(stderr, "Enter mrr, mtt, mff, mrt, mrf, mtf : ");
		fprintf(stderr, "(the exponent is not taken into account)");
		scanf("%lf %lf %lf %lf %lf %lf", &mrr, &mtt, &mff, &mrt, &mrf, &mtf);
	} else {
		mrr = atof(argv[1]); mtt = atof(argv[2]); mff = atof(argv[3]);
		mrt = atof(argv[4]); mrf = atof(argv[5]); mtf = atof(argv[6]); 
	}

	ten2axe(mrr, mtt, mff, mrt, mrf, mtf, &pclvd, &piso, &pdc, &tv, &ta, &td, &nv, &na, &nd, &pv, &pa, &pd);
	v_iso = (tv + nv + pv) / 3.;
	tv1 = tv - v_iso; nv1 = nv - v_iso; pv1 = pv - v_iso;
	if(tv1 > - pv1) { clvd = - nv1 / tv1; iso = v_iso / tv1; }
	else { clvd = - nv1 / pv1; iso = v_iso / pv1; }
	modc = tv + nv / 2.;

//	fprintf(stderr, "pclvd piso pdc tval tazim tdip nval nazim ndip pval pazim pdip iso clvd mo_scalaire\n");
	fprintf(stdout, "%d %d %d %.2lf %.0lf %.0lf %.2lf %.0lf %.0lf %.2lf %.0lf %.0lf %.2lf %.2lf %.2e \n", pclvd, piso, pdc, tv, ta, td, nv, na, nd, pv, pa, pd, iso, clvd, modc);
}
