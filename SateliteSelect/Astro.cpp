
#include <cmath> 
#include "Sgpsdp.h" 
#include "vector.h" 


VECTOR CSGP4_SDP4::CalculateSolarPosition(double time)	//, solar_vector : vector);
{
	static VECTOR solar_vector;
	double DeltaET_year;
	double mjd,year,T,M,L,e,C,O,Lsa,nu,R,eps;
	mjd  = time - 2415020.0;
	year = 1900.0 + mjd/365.25;
	DeltaET_year = 26.465 + 0.747622*(year - 1950.0) + 1.886913*sin(2.0 * PI * (year - 1975.0)/33.0);

	T = (mjd + DeltaET_year/secday)/36525.0;
	M = (Modulus(358.47583 + Modulus(35999.04975*T,360.0) - (0.000150 + 0.0000033*T)*sqr(T),360.0)) * PI/180.0;
	L = (Modulus(279.69668 + Modulus(36000.76892*T,360.0) + 0.0003025*sqr(T),360.0)) * PI/180.0;
	e = 0.01675104 - (0.0000418 + 0.000000126*T)*T;
	C = ((1.919460 - (0.004789 + 0.000014*T)*T)*sin(M) + (0.020094 - 0.000100*T)*sin(2*M) + 0.000293*sin(3*M)) * PI/180.0;
	O = (Modulus(259.18 - 1934.142*T,360.0)) * PI/180.0;
	Lsa = Modulus(L + C - (0.00569 - 0.00479*sin(O)) * PI/180.0, PI * 2.0);
	nu  = Modulus(M + C,PI * 2.0);
	R   = 1.0000002*(1 - sqr(e))/(1 + e*cos(nu));
	eps = (23.452294 - (0.0130125 + (0.00000164 - 0.000000503*T)*T)*T + 0.00256*cos(O)) * PI/180.0;
	R   = AU*R;
	solar_vector.x = R*cos(Lsa);
	solar_vector.y = R*sin(Lsa)*cos(eps);
	solar_vector.z = R*sin(Lsa)*sin(eps);
	solar_vector.w = R;
	// end; {Procedure Calculate_Solar_Position}
	return solar_vector;
}

double CSGP4_SDP4::DepthOfEclipse(double time, VECTOR r1)
{
	double r1_r1,r1_r2,r2_r2,k,d,ds;
	VECTOR r2;
	Magnitude(&r1);
	r2 = CalculateSolarPosition(time);
//	solar_pos = r2;
	r1_r1 = sqr(r1.w);
	r1_r2 = -Dot(r1,r2);
	r2_r2 = sqr(r2.w);
	k = r1_r2/r2_r2;
	// {Calculate perpendicular distance from anti-solar vector}
	d = sqrt(r1_r1 - sqr(r1_r2)/r2_r2);
	// {Calculate shadow distance ds}
	ds = xkmper + k * (sr - xkmper);
	// {If d < ds, then satellite is in eclipse}
	if ( (k > 0.0) && (d < ds) )
		m_bEclipsed = true;
	else
		m_bEclipsed = false;
	return d - ds;
	// Depth_of_Eclipse := d - ds
	// end; {Function Depth_of_Eclipse}
}

bool CSGP4_SDP4::GetEclipsed()
{
	return m_bEclipsed;
}

void CSGP4_SDP4::CalculateLatLonAlt(double jdTime)
{
	m_vLLA = CalculateLatLonAlt(m_vPOS, jdTime);
	m_bLatLonAlt = true;
}

VECTOR CSGP4_SDP4::CalculateLatLonAlt(VECTOR vPOS, double time)
{
// Reference:  The 1992 Astronomical Almanac, page K12. 
	static VECTOR vLLA;
	double lat,lon,alt;
	double theta,r,e2,phi,c;
	double arg1, arg2;

	vLLA.x = vLLA.y = vLLA.z = vLLA.w = 0.0;
	lat = lon = alt = 0.0;
	theta = r = e2 = phi = c = 0.0;

//	theta = atan2(vPOS.y,vPOS.x);
	theta = AcTan(vPOS.y,vPOS.x);
	
	arg1 = ThetaG(time);
	arg1 = theta - arg1;
	arg2 = 2.0* PI;

//	lon = Modulus(theta - ThetaG(time),2.0*PI);
	lon = Modulus(arg1, arg2);

	r = sqrt(sqr(vPOS.x) + sqr(vPOS.y));
	e2 = f*(2.0 - f);
	lat = AcTan(vPOS.z,r);
	do	{
		phi = lat;
		c = 1.0/sqrt(1.0 - e2*sqr(sin(phi)));
		lat = AcTan( vPOS.z + xkmper*c*e2*sin(phi),r);
	}	while (fabs(lat - phi) > 1E-10);//1E-7); For speeding up calculation 7 digit
										//is exact enough (123.45
	alt = r/cos(lat) - xkmper*c;

	vLLA.x = lat*180.0/PI;   // radians
	vLLA.y = lon*180.0/PI;   // radians
	vLLA.z = alt;			// kilometers
	vLLA.w = theta*180.0/PI; // radians
	return vLLA;
}

/*----------------------------------------------------------------------*/
//int CSGP4_SDP4::Eclipse( double *r, double *sun, double *moon, int *em, int *ee, char **which )
//int CSGP4_SDP4::Eclipse( double *r, VECTOR *vSun, VECTOR *vMoon, int *em, int *ee, char **which )
//{
//  /* function to compute if object at location r is eclipsed by
//     either the earth || moon.  ECI coordinates for all vectors.
//     Inputs:
//	r is the location of s/c
//	sun is location of sun
//	moon is location of moon
//     Outputs
//	which points to  message about result of the SatteliteModel1
//	em is zero if not eclipsed by moon
//	ee is zero if not eclipsed by earth
//     Returns 0 if no eclipse.
//  */
//	double us[3], um[3], ue[3];
//	double alpha, beta, beta_sun;
//	double mdist, edist, sdist;
//	double x[3];
//	double *sun = (double *)&vSun;
//	double *moon= (double *)&vMoon; 
//  *em = 0;               /* assume no eclipse */
//  *ee = 0;
//  /* find distances && directions to sun, moon && earth */
//  VecDiff( sun, r, x, 3 ); sdist = VecMag( x, 3 );
//  UnitVec( x, us, 3 );   /* direction towards the sun */            
//  VecDiff( moon, r, x, 3 ); mdist = VecMag( x, 3 );
//  UnitVec( x, um, 3 );   /* direction towards the moon */            
//  VecScale( -1.0, r, ue, 3 ); /* vector to earth center */
//  UnitVec( ue, ue, 3 );  /* direction to earth center */
//  edist = VecMag( r, 3 );
//  /* cannot have eclipse if sun is closest */
//  *which = "No Eclipse";
//  if (( sdist <= edist ) && (sdist <= mdist )) return 0;
//  beta_sun = asin( RSUN / sdist );    /* half angle of sun */
//  /* look for eclipse by earth */
//  beta = asin( REarth / edist );      /* half angle of earth */
//  alpha = acos( VecDot( ue, us, 3 )); /* angle from earth to sun */
//  if ( alpha < (beta + beta_sun)) {   /* some kind of eclipse */
//    *ee = 1;
//    if (( beta >= beta_sun ) && ( alpha <= (beta - beta_sun))){ /* total */
//      *which = "Total Eclipse by Earth"; return 1; }
//    else {   /* partial */
//      *which = "Partial Eclipse by Earth"; return 1; }
//  }
//  /* look for eclipse by moon */
//  beta = asin( RMOON / mdist );      /* half angle of moon */
//  alpha = acos( VecDot( um, us, 3 ));/* angle from moon to sun */
//  if ( alpha < (beta + beta_sun)) {  /* some kind of eclipse */
//    *em = 1;
//    if (( beta >= beta_sun ) && ( alpha <= (beta - beta_sun))){ /* total */
//      *which = "Total Eclipse by Moon"; return 1; }
//    else {   /* partial */
//      *which = "Partial Eclipse by Moon"; return 1; }
//  }
//  return 0;
//}
/*----------------------------------------------------------------------*/