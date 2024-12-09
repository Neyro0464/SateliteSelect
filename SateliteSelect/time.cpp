#include <cmath>
#include "Sgpsdp.h"
#include <ctime>
#include <chrono>

using namespace std::chrono;

double CSGP4_SDP4::JulianDate(tm st)
{
    double fDay, fDate;
    int iYear = st.tm_year + 1900;
    int iMonth = st.tm_mon + 1;
    int iDay = st.tm_mday;
    int iHour = st.tm_hour;
    int iMin = st.tm_min;
    int iSec = st.tm_sec;

    int iEpocheYear = EpocheYear(iYear) * 1000; // 1999 to 99000 && 2032 to 32000.
    int iDayOfYear = DayOfYear(iYear, iMonth, iDay);
    double fFraction = FractionOfDay(iHour, iMin, iSec);

    fDate = iEpocheYear + iDayOfYear + fFraction;
    fDay = JulianDate(fDate);
    return fDay; // The JulianDate of the CTime Date
}

double CSGP4_SDP4::JulianDate(double st)
// Astronomical Formulae for Calculators, Jean Meeus, pages 23-25
// Calculate Julian Date of 0.0 Jan year
{
    // input like : 99200.27104438 => 1999, 200.27104438 days.
    double fYear, fDay;
    int iYear = (int)(st / 1000.0);
    fDay = st - iYear * 1000.0; // = 99200,271 - 99000.0 = 200,271
    if (iYear < 57) iYear += 2000;
    else iYear += 1900;
    fYear = JulianDateOfYear(iYear);
    return fYear + fDay;
}

double CSGP4_SDP4::JulianDateOfYear(int yr)
{
    long A, B;
    int iYear = yr - 1;
    double fYear;
    A = (long)(iYear / 100.0);
    B = 2 - A + (long)(A / 4.0);
    fYear = (long)(365.25 * iYear) + 428.0 + 1720994.5 + B;
    return fYear;
}

int CSGP4_SDP4::EpocheYear(int iYear)
{
    // returns 99 for 1999 && 23 for 2023
    int iEpocheYear;
    if (iYear > 2000)
        iEpocheYear = iYear - 2000;
    else
        iEpocheYear = iYear - 1900;
    return iEpocheYear;
}

int CSGP4_SDP4::DayOfYear(int yr, int mo, int dy)
{
    // yr=1999 / mo=12 / dy=24 => 12/24/1999=christmas 1999
    // December is never used since the last Month of a year
    // is added by getting the days from that month.
    int days[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    int i, iDay;

    iDay = 0;
    for (i = 0; i < mo - 1; i++) // mo-2 because zero based && just add until last month
        iDay += days[i];
    iDay += dy;

    if (((yr % 4) == 0) && // every four years
        (((yr % 100) != 0) || ((yr % 400) == 0)) && // but not at centuries other then every 400 years
        (mo > 2)) iDay += 1; // Schaltjahr ???
    return iDay;
}

double CSGP4_SDP4::FractionOfDay(int hr, int min, int sec)
{
    double fFraction;
    fFraction = (double)(60.0 * 60.0 * hr +
        60.0 * min +
        (double)sec) / (24.0 * 60.0 * 60.0);
    return fFraction;
}

tm CSGP4_SDP4::CalendarDate(double dJulian)
{
    //  { Astronomical Formulae for Calculators, Jean Meeus, pages 26-27 }
    int month, alpha;
    double Z, A, B, C, D, E, F;
    double day, year;
    double hour, minute, second; //, milliseconds;

    static tm ret;
    F = modf(dJulian + 0.5, &Z);
    if (Z < 2299161)
        A = Z;
    else {
        alpha = (int)((Z - 1867216.25) / 36524.25);
        A = Z + 1.0 + alpha - (int)(alpha / 4);
    }
    B = A + 1524.0;
    modf((B - 122.1) / 365.25, &C);
    modf(365.25 * C, &D);
    modf((B - D) / 30.6001, &E);
    day = B - D - (int)(30.6001 * E) + F;

    if (E < 13.5) month = round(E - 1.0);
    else month = round(E - 13.0);
    if (month > 2.5) year = C - 4716.0;
    else year = C - 4715.0;

    A = modf(day, &Z); // A = fraction of day ...
    A = modf(A * 24.0, &hour);
    A = modf(A * 60.0, &minute);
    A = modf(A * 60.0, &second);

    ret.tm_year = (int)year - 1900;
    ret.tm_mon = (int)month - 1;
    ret.tm_mday = (int)day;
    ret.tm_hour = (int)hour;
    ret.tm_min = (int)minute;
    ret.tm_sec = (int)++second;
    ret.tm_wday = 0;
    ret.tm_yday = 0;
    ret.tm_isdst = 0;
    return ret;
}

double CSGP4_SDP4::SideralTime(double jd)
{
    // Reference:  The 1992 Astronomical Almanac, page B6.
    double modulo;
    double UT, TU, GMST;
    UT = modf(jd + 0.5, &TU);
    jd -= UT;
    TU = (double)(jd - 2451545.0) / 36525.0;
    GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6));
    GMST = (GMST + secday * omega_E * UT); //%secday;

    modf(GMST / secday, &modulo);
    modulo = GMST - modulo * secday;
    if (modulo < 0.0)
        modulo += secday;
    return PI * 2.0 * modulo / secday;
}

double CSGP4_SDP4::ThetaG(double jd)
{
    // Reference:  The 1992 Astronomical Almanac, page B6.
    double UT, TU, GMST;
    UT = modf(jd + 0.5, &jd);
    TU = (jd - 2451545.0) / 36525.0;
    GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6));
    GMST = Modulus(GMST + secday * omega_E * UT, secday);
    return (2.0 * PI * GMST / secday);
}
