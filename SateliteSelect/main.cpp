#include <iostream>
#include <fstream>
#include "Sgpsdp.h"
#include <ctime>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#pragma warning(disable : 4996)

using namespace std;
using namespace std::chrono;

struct LLAPos // latitude, logtitude and altitude
{
    double Lat = 0, Lon = 0, Alt = 0;
};

struct Satellite
{
    std::string  SatelliteName;
    int SatelliteNumber;
    double AzmDeg;
    double ElvDeg;
    double Lat;
    double Lon;
    double Alt;
    tm time;
    double velocity;
    bool direction;
    double timeStamp;
};

// Reading TLE file and make new class objects
bool readFileAndCallFunction(std::string filename, std::vector<CSGP4_SDP4>& SatelliteModel, std::streampos* here) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ERRO: File is not open " << std::endl;
        return false;
    }

    char m_cLine0[70];
    char m_cLine1[70];
    char m_cLine2[70];
    file.seekg(*here);
    if (file.getline(m_cLine0, sizeof(m_cLine0)) &&
        file.getline(m_cLine1, sizeof(m_cLine1)) &&
        file.getline(m_cLine2, sizeof(m_cLine2))) {
        SatelliteModel.push_back(CSGP4_SDP4(m_cLine0, m_cLine1, m_cLine2));
        *here = file.tellg();
        file.close();
        return true;
    }

    file.close();
    return false; const int N = 600;
}

// Define station possition and check constraints
bool defStationPos(LLAPos* station) {
    std::cout << "Enter the coordinates (latitude longitude altitude) of the station:";
    double Lat = 0.0, Lon = 0.0, Alt = 0.0;
    std::cin >> Lat >> Lon >> Alt;
    if ((Lat < -90 || Lat > 90) || (Lon < -180 || Lon > 180)) {
        std::cout << "ERROR: wrong parametrs(Lat || Lon)";
        return false;
    }
    station->Lat = Lat;
    station->Lon = Lon;
    station->Alt = Alt;
    return true;
}

bool defStationVision(double& minAzm, double& maxAzm, double& minElv, double& maxElv, int& timeMinObserveSec) {
    std::cout << "Enter vision parameters of station: ";
    std::cin >> minAzm >> maxAzm >> minElv >> maxElv >> timeMinObserveSec;
    if (minAzm < 0 || minAzm >= maxAzm || maxAzm > 360 || minElv < 10 || maxElv < minElv || maxElv > 90 || timeMinObserveSec < 1) { // break this condition on many to easier detection of the problem
        std::cout << "ERROR: wrong parameters of station";
        return false;
    }
    return true;
}

// Convertation LLA coordinates to Decart coordinates
// (0;0;0) is center of planet
// (6371;0;0) is lan = 0, lon = 0, alt = 0
VECTOR ConvertLLAtoCoordinate(LLAPos* Object) {
    double x = (6371.0 + Object->Alt) * cos(Object->Lat) * cos(Object->Lon);
    double y = (6371.0 + Object->Alt) * cos(Object->Lat) * sin(Object->Lon);
    double z = (6371.0 + Object->Alt) * sin(Object->Lat);
    return { x, y, z };
}

// Calculate satellite position with class function
// then convert to Decart coordinates 
VECTOR SatellitePos(CSGP4_SDP4 SGP, LLAPos* SatelliteLLA_1, VECTOR Satellite_1, double time) {
    SGP.CalculateLatLonAlt(time);
    SatelliteLLA_1->Lat = SGP.GetLat();
    SatelliteLLA_1->Lon = SGP.GetLon();
    SatelliteLLA_1->Alt = SGP.GetAlt();

    SatelliteLLA_1->Lat = SGP.DegToRad(SatelliteLLA_1->Lat);
    SatelliteLLA_1->Lon = SGP.DegToRad(SatelliteLLA_1->Lon);


    Satellite_1 = ConvertLLAtoCoordinate(SatelliteLLA_1);
    return Satellite_1;
}

// Transfer satellite coordinate system in station coodinate system
VECTOR Transferring(VECTOR object, LLAPos Station) {
    double x = object.x;
    double y = object.y;
    double z = object.z;

    double alpha = Station.Lon;
    double beta = Station.Lat;
    double R = 6371.0 + Station.Alt;

    x = x + R;

    double x0 = cos(beta) * (x * cos(alpha) - y * sin(alpha)) - z * sin(beta);
    double y0 = x * sin(alpha) + y * cos(alpha);
    double z0 = sin(beta) * (x * cos(alpha) - y * sin(alpha)) + z * cos(beta);

    return { x0, y0, z0 };
}

// Checking intersection: if satellite in view area of station
bool IntersectCheck(LLAPos StationLLA_1, VECTOR Satellite_1, double& minAzm, double& maxAzm, double& minElv, double& maxElv, double& theta, double& fi) {
    double x = Satellite_1.x;
    double y = Satellite_1.y;
    double z = Satellite_1.z;

    // longtitude and latitude of the station to rotate around the sphere(Earth) 
    double alpha = StationLLA_1.Lon;
    double beta = StationLLA_1.Lat;
    double R = 6371.0 + StationLLA_1.Alt;

    

    // Azm(gamma) and Elv(delta) rotate around the new coordinate system.
    // koef_x and koef_y using to transform vision cone (narrow it down).
    double gamma = (maxAzm + minAzm) / 2 * PI / 180;
    double delta = PI / 2 - (minElv + maxElv) / 2 * PI / 180;

    double koef_y = tan(PI / 180 * (maxAzm - minAzm)) / 2;
    double koef_x = tan(PI / 180 * (maxElv - minElv)) / 2;
    gamma = -gamma; //invert gamma to make Azimut work (increasing the angle -> go clockwise)

    double new_x = z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha));
    double new_y = y * cos(alpha) - x * sin(alpha);
    double new_z = (z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha))) - R;

    fi = atan(new_y / new_x);
    theta = atan(sqrt(pow(new_x, 2) + pow(new_y, 2)) / new_z);

    /* Previous variant with bigger cone of view 
    double prepare = z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha));

    if (prepare <= R)
        return false;

    double koef = tan(PI / 18.0); // �� 10 �� 90 �������� (10-170)
    double intersection = pow(z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha)), 2) + pow(y * cos(alpha) - x * sin(alpha), 2) - pow(z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha)) - R, 2) / koef;
    */
    

    double prepare = new_z * cos(delta) + sin(delta) * (-new_y * sin(gamma) + new_x * cos(gamma));
    if (prepare <= 0)
        return false;

    double equation = 1 / koef_x * pow(-new_z * sin(delta) + cos(delta) * (new_x * cos(gamma) - new_y * sin(gamma)), 2) + 1 / koef_y * pow(new_y * cos(gamma) + new_x * sin(gamma), 2) - pow(new_z * cos(delta) + sin(delta) * (-new_y * sin(gamma) + new_x * cos(gamma)), 2);
/*
    fi = atan((new_y * cos(gamma) + new_x * sin(gamma) / -new_z * sin(delta) + cos(delta) * (new_x * cos(gamma) - new_y * sin(gamma))));
    theta = atan((sqrt(pow(-new_z * sin(delta) + cos(delta) * (new_x * cos(gamma) - new_y * sin(gamma)), 2) + pow(new_y * cos(gamma) + new_x * sin(gamma), 2))) / new_z * cos(delta) + sin(delta) * (-new_y * sin(gamma) + new_x * cos(gamma)));
    */
    if (equation < 0)
        return true;
    else
        return false;
}


bool VelocityCheck(LLAPos Station, std::vector<CSGP4_SDP4>& SatelliteModel) {

    VECTOR Satellite_1 = { 0, 0, 0, 0 };
    VECTOR Satellite_2 = { 0, 0, 0, 0 };
    LLAPos SatelliteLLA;

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm Time2 = now_tm;
    Time2.tm_min += 1;

    double time = SatelliteModel[0].JulianDate(now_tm);
    double time2 = SatelliteModel[0].JulianDate(Time2);


    SatelliteModel.back().SGP(time);
    Satellite_1 = SatellitePos(SatelliteModel.back(), &SatelliteLLA, Satellite_1, time);
    Satellite_1 = Transferring(Satellite_1, Station);

    double theta = atan(sqrt(pow(Satellite_1.x, 2) + pow(Satellite_1.y, 2)) / Satellite_1.z);
    theta = SatelliteModel[0].RadToDeg(theta);

    SatelliteModel.back().SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel.back(), &SatelliteLLA, Satellite_2, time);
    Satellite_2 = Transferring(Satellite_2, Station);

    double theta2 = atan(sqrt(pow(Satellite_2.x, 2) + pow(Satellite_2.y, 2)) / Satellite_2.z);
    theta2 = SatelliteModel[0].RadToDeg(theta2);

    return (abs(theta2 - theta) >= 1.9);
}

//Function that save only satellites with suitable speed
void SatelliteFilter(std::string filename, LLAPos Station, std::vector<CSGP4_SDP4>& SatelliteModel) {
    
    std::streampos here = 0;
    while (readFileAndCallFunction(filename, SatelliteModel, &here)) {
        
        if (VelocityCheck(Station, SatelliteModel)) {
            SatelliteModel.pop_back();
            continue;
        }
    }
}



// checks the intersection of the satellite with the station's visibility cone (+30 sec to +30 min)
//  Filter satellites that would be in cone of visibility for timeMinObserve(sec)
double TimeIntersect(CSGP4_SDP4 SatelliteModel, LLAPos StationLLA_1, double& minAzm, double& maxAzm, double& minElv, double& maxElv, int& timeMinObserveSec, double &theta, double& fi) { //remake
    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm lowerTime = now_tm;  // lower bound of time (current time + 30 sec)
    lowerTime.tm_sec += 30;

    tm temp = now_tm;       // time that we shift to find required time (current time)

    tm upperTime = now_tm;  // upper bound of time (current time + 30 min)
    upperTime.tm_min += 30;


    double lowerTm = SatelliteModel.JulianDate(lowerTime); // lower bound of time interval
    double tempTm = SatelliteModel.JulianDate(temp); // time cursor that move between lower and upper bounds of time interval
    double upperTm = SatelliteModel.JulianDate(upperTime); // upper bound of time interval

    LLAPos SatelliteLLA_1;
    VECTOR Satellite_1 = { 0, 0, 0, 0 };
    bool check1 = 0;

    while ((check1 == false || tempTm <= lowerTm) && tempTm <= upperTm) { 

        tempTm = SatelliteModel.JulianDate(temp);
        SatelliteModel.SGP(tempTm);
        Satellite_1 = SatellitePos(SatelliteModel, &SatelliteLLA_1, Satellite_1, tempTm);
        check1 = IntersectCheck(StationLLA_1, Satellite_1, minAzm, maxAzm, minElv, maxElv, theta, fi);
        temp.tm_sec += 1;
    }

    if (tempTm > upperTm || tempTm < lowerTm) {
        return 0;
    }

    // block for time filter
    tm filterTime = temp;
    filterTime.tm_sec += timeMinObserveSec - 1;
    double checkTime = SatelliteModel.JulianDate(filterTime);
    Satellite_1 = SatellitePos(SatelliteModel, &SatelliteLLA_1, Satellite_1, checkTime);
    if (IntersectCheck(StationLLA_1, Satellite_1, minAzm, maxAzm, minElv, maxElv, theta, fi)) {
        return tempTm;
    }
    else
        return 0;
}

// Check direction by defining sign of the angle between station and satellite
bool getDiraction(LLAPos Station, CSGP4_SDP4 SatelliteModel) {
    VECTOR Satellite_1 = { 0, 0, 0, 0 };
    VECTOR Satellite_2 = { 0, 0, 0, 0 };
    LLAPos SatelliteLLA;

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm Time2 = now_tm;
    Time2.tm_min += 1;

    double time = SatelliteModel.JulianDate(now_tm);
    double time2 = SatelliteModel.JulianDate(Time2);

    SatelliteModel.SGP(time);
    Satellite_1 = SatellitePos(SatelliteModel, &SatelliteLLA, Satellite_1, time);
    Satellite_1 = Transferring(Satellite_1, Station);

    double theta = atan(sqrt(pow(Satellite_1.x, 2) + pow(Satellite_1.y, 2)) / Satellite_1.z);
    theta = SatelliteModel.RadToDeg(theta);/////////////////

    SatelliteModel.SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel, &SatelliteLLA, Satellite_2, time);
    Satellite_2 = Transferring(Satellite_2, Station);
    double theta2 = atan(sqrt(pow(Satellite_2.x, 2) + pow(Satellite_2.y, 2)) / Satellite_2.z);
    theta2 = SatelliteModel.RadToDeg(theta2);
    if (theta2 - theta > 0) { // If diff > 0 that means satellite go higher and closer
        return true;
    }
    else {                  // If diff <= 0 that means satellite go lower and further
        return false;
    }
}


int main() {
    std::string filename = "TLE.txt";
    std::vector<CSGP4_SDP4> SatelliteModel;
    std::vector<Satellite> Objects;
    

    double minAzm = 0, maxAzm = 0, minElv = 0, maxElv = 0; 
    int timeMinObserveSec = 0;

    LLAPos StationLLA_1;
    VECTOR Station_1{ 0, 0, 0, 0 };
    if (!defStationPos(&StationLLA_1))
        return 0;
    if (!defStationVision(minAzm, maxAzm, minElv, maxElv, timeMinObserveSec)) {
        return 0;
    }
    std::cout << minAzm << ' ' << maxAzm << ' ' << minElv << ' ' << maxElv << ' ' << timeMinObserveSec << std::endl;

    StationLLA_1.Lat = SatelliteModel[0].DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatelliteModel[0].DegToRad(StationLLA_1.Lon);


    SatelliteFilter(filename, StationLLA_1, SatelliteModel);

    int N = SatelliteModel.size();
    double time = 0;
    int satNumer = 0;
    for (int i = 0; i < N; i++) {
        satNumer = SatelliteModel[i].GetNORAD();
        
        SATELLITE* s = (SATELLITE*)SatelliteModel[i].GetSatellite(); 

        double theta = 0;
        double fi = 0;
        time = TimeIntersect(SatelliteModel[i], StationLLA_1, minAzm, maxAzm, minElv, maxElv, timeMinObserveSec, theta, fi);
        if (time != 0) {
            bool dir = getDiraction(StationLLA_1, SatelliteModel[i]); 
            SatelliteModel[i].CalculateLatLonAlt(time);
            double Lat = SatelliteModel[i].GetLat();
            double Lon = SatelliteModel[i].GetLon();
            double Alt = SatelliteModel[i].GetAlt();
            fi = SatelliteModel[i].RadToDeg(fi);
            theta = SatelliteModel[i].RadToDeg(theta);
            Objects.push_back({ s->cSatelliteName, satNumer, fi, theta, Lat, Lon, Alt, SatelliteModel[i].CalendarDate(time), SatelliteModel[i].GetVel().w, dir, time});
        }


    }

    std::sort(Objects.begin(), Objects.end(), [](const Satellite& a, const Satellite& b) {
        return a.timeStamp < b.timeStamp;
        });

    std::cout << "# Satellite Name " << " \t | " << "# Satellite Number" << " \t| " << "AZM   |" << "ELV   |" << "Lat   | " << "Lon    |" << "Alt    | " << "Velocity km/s\t | " << "Direction\t| " << "DD/MM/YYYY\thh:mm:ss" << std::endl;

    for (const auto& item : Objects) {
        if (item.direction) {
            std::cout << "# " << item.SatelliteName << " | " << item.SatelliteNumber << " \t\t| " << item.AzmDeg << "| " << item.ElvDeg << "| " << item.Lat << " |" << item.Lon << " |" << item.Alt << " |" << item.velocity << "\t\t | \t + \t | " << item.time.tm_mday << " / " << item.time.tm_mon + 1 << " / " << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << " : " << item.time.tm_sec << " " << std::endl;
        }
        else {
            std::cout << "# " << item.SatelliteName << " | " << item.SatelliteNumber << " \t\t| " << item.AzmDeg << "| " << item.ElvDeg << "| " << item.Lat << " |" << item.Lon << " |" << item.Alt << " |" << item.velocity << "\t\t | \t - \t | " << item.time.tm_mday << " / " << item.time.tm_mon + 1 << " / " << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << " : " << item.time.tm_sec << " " << std::endl;
        }
    }
    return 0;
}