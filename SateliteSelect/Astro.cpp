#include "Sgpsdp.h" 
#include "myVector.h" 

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