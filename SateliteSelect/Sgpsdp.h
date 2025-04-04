#ifndef SGPSDP_H
#define SGPSDP_H

#include <iostream>
#include<cstring>
#include <string>
#include <stdio.h>
#include <ctime>
#include <time.h>
#include <chrono>
#include <cmath> 

//* Already defined in  >stdafx.h<
typedef struct tagVECTOR {
    double x, y, z, w;
}VECTOR;

typedef struct tagSATELLITE {
    char    cSatelliteName[23];
    int     iSecondMeanMotion;
    int     iSatelliteNumber;
    int     iLaunchYear;
    int     iLaunchNumber;
    char    cLaunchPiece[3];
    int     iEpochYear;
    double  fEpochDay;
    int     iEpochDay;
    double  fEpochFraction;
    double  fBalisticCoefficient;
    double  fSecondMeanMotion;
    double  fRadiationCoefficient;
    char    cEmphemeristType[2];
    int     iElementNumber;

    double  fInclination;
    double  fRightAscending;
    double  fEccentricity;
    double  fPeregee;
    double  fMeanAnomaly;
    double  fMeanMotion;
    int     iRevAtEpoch;
    double  fJulianEpoch;
}SATELLITE;

#ifndef PI
#define PI 3.141592653589793
#endif

class CSGP4_SDP4 {
protected:
    char        m_cLine0[23], m_cLine1[70], m_cLine2[70];
    SATELLITE   m_Sat;
    VECTOR      m_vLLA;
    VECTOR      m_vPOS;
    VECTOR      m_vVEL;
    VECTOR      m_vUPos;
    VECTOR      m_vUVel;
    VECTOR      m_vObs;
    VECTOR      m_vRad;
    VECTOR      m_vSolarPosition;
    double      m_fTime;
protected:
    bool        m_bLatLonAlt;
    double      ae;
    double      tothrd;
    double      xkmper;
    double      f;
    double      ge;
    double      J2;
    double      J3;
    double      J4;
    double      ck2;
    double      ck4;
    double      xj3;
    double      qo;
    double      s;
    double      e6a;
    double      dpinit;
    double      dpsec;
    double      dpper;
    long        secday;
    double      omega_E;
    double      sr;
    double      AU;

    int         iflag, ideep;
    double      xke, xmnpda;

    double      eqsq, siniq, cosiq, rteqsq, ao, cosq2, sinomo, cosomo;
    double      bsq, xlldot, omgdt, xnodot, xnodp;
    double      xll, omgasm, xnodes, _em, xinc, xn, t;
    double      qoms2t;
    bool        m_bEclipsed;
public:
    double      Modulus(double a1, double a2);
    void        ConvertData();
    void        Init();
    long        round(double arg);
    double      Fmod2p(double arg);
    double      AcTan(double sinx, double cosx);
    double      sqr(double arg);
    void        Magnitude(VECTOR* pVector);
    double      Dot(VECTOR v1, VECTOR v2);

protected:
    bool SGP4(double tsince, int* iflag, VECTOR* pos, VECTOR* vel);
    bool SDP4(double tsince, int* iflag, VECTOR* pos, VECTOR* vel);
    void Call_dpper(double* e, double* xincc, double* omgadf, double* xnode, double* xmam);
    void Call_dpsec(double* xmdf, double* omgadf, double* xnode, double* emm, double* xincc,
        double* xnn, double* tsince);
    void Call_dpinit(double* eosq, double* sinio, double* cosio, double* betao, double* aodp,
        double* theta2, double* sing, double* cosg, double* betao2, double* xmdot,
        double* omgdot, double* xnodott, double* xnodpp);
    bool Deep(int ideep);
public:
    CSGP4_SDP4(char* m_cLine0, char* m_cLine1, char* m_cLine2);
    CSGP4_SDP4(SATELLITE* pSat);
    CSGP4_SDP4();
    ~CSGP4_SDP4();
    void SetSatellite(char* m_cLine0, char* m_cLine1, char* m_cLine2, bool bInitSatellite = true);
    void SetSatellite(SATELLITE* pSat, bool bInitSatellite = true);
    void InitSatellite();

    bool SGP(double time);
    void ConvertSatState(VECTOR* pos, VECTOR* vel);
    void CalculateLatLonAlt(double time);
    VECTOR CalculateLatLonAlt(VECTOR vPos, double time);
    void CalculateUserPosVel(VECTOR* geodetic, double time);
    bool CalculateObs(VECTOR pos, VECTOR vel, VECTOR geodetic, double time);
    void CalculateRADec(VECTOR pos, VECTOR vel, VECTOR geodetic, double time);

    double      GetFloat(int iStart, int iEnd, char* cLine);
    long        GetInt(int iStart, int iEnd, char* cLine);
    char* GetString(int iStart, int iEnd, char* cLine);
    double      RadToDeg(double arg);
    double      DegToRad(double arg);
    double      JulianDate(tm st);
    double      JulianDate(double st);
    double      JulianDateOfYear(int yr);
    int         EpocheYear(int iYear);
    int         DayOfYear(int, int, int);
    double      FractionOfDay(int, int, int);
    double      ThetaG(double jd);
    tm          CalendarDate(double dJulian);
    double      SideralTime(double jd);
    VECTOR      GetPos();
    VECTOR      GetVel();
    double      GetLat();
    double      GetLon();
    double      GetAlt();
    double      GetTime();
    VECTOR      GetUserPos();
    VECTOR      GetUserVel();
    VECTOR      GetObserver();
    VECTOR      GetRADec();
    void* GetSatellite();
    int         GetNORAD();

    double      VecDot(double* X, double* Y, int N);
    void        VecCross(double* X, double* Y, double* Z, int N);
    double      VecMag(double* X, int N);
    void        UnitVec(double* X, double* Y, int N);
    void        VecDiff(double* X, double* Y, double* Z, int N);
    void        VecSum(double* X, double* Y, double* Z, int N);
    void        VecScale(double u, double* X, double* Y, int N);
};


#endif //SGPSDP_H