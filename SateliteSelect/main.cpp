#include <iostream>
#include <iomanip>
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
    int SatelliteNumber;    // NORAD number of the satellite
    double AzmDeg;          // Azimut in degrees relative to the station
    double ElvDeg;          // Elevation in degree relative to the station
    double Lat;             // Latitude 
    double Lon;             // Longtitude
    double Alt;             // Altitude
    tm time;                // time to filter
    double velocity;        // velocity km/s
    bool direction;         // direction (to/from) station
    double timeStamp;       // time when the satellite will be in the cone of view
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
    return false;
}

// Define station possition and check constraints
bool defStationPos(LLAPos* station) {
    std::cout << "Enter the coordinates of the station (latitude longitude altitude):";
    double Lat = 0.0, Lon = 0.0, Alt = 0.0;
    std::cin >> Lat >> Lon >> Alt;
    if (Lat < -90 || Lat > 90) {
        std::cout << "ERROR: wrong parametrs (-90 < Lat < 90)";
        return false;
    }
    if (Lon < -180 || Lon > 180) {
        std::cout << "ERROR: wrong parametrs (-180 < Lon < 180)";
        return false;
    }
    station->Lat = Lat;
    station->Lon = Lon;
    station->Alt = Alt;
    return true;
}


//Function to enter station's vision parametres
bool defStationVision(double& minAzm, double& maxAzm, double& minElv, double& maxElv, int& timeMinObserveSec) {
    std::cout << "Enter vision parameters of the station (minAzm maxAzm minElv maxElv time): ";
    std::cin >> minAzm >> maxAzm ;
    if (minAzm < 0 || minAzm > maxAzm || maxAzm >= 360) { 
        std::cout << "ERROR: wrong parameters of station - Azimut";
        return false;
    }
    std::cin >> minElv >> maxElv;
    if (minElv < 10 || maxElv < minElv || maxElv > 180){     
        std::cout << "ERROR: wrong parameters of station - Elevation";
        return false;
    }
    std::cin >> timeMinObserveSec;
    if (timeMinObserveSec < 1) {                           
        std::cout << "ERROR: wrong parameters of station - time";
        return false;
    }
    return true;
}

// Convertation LLA coordinates to Decart coordinates
// previously converted to radians!!!!
// (0;0;0) is center of planet
// (6371;0;0) is lan = 0, lon = 0, alt = 0
VECTOR ConvertLLAtoCoordinate(LLAPos& Object) {
    double x = (6371.0 + Object.Alt) * cos(Object.Lat) * cos(Object.Lon);
    double y = (6371.0 + Object.Alt) * cos(Object.Lat) * sin(Object.Lon);
    double z = (6371.0 + Object.Alt) * sin(Object.Lat);
    return {x, y, z, 0.0};
}

// Calculate satellite position with class function
// then convert to Decart coordinates 
VECTOR SatellitePos(CSGP4_SDP4 SGP, LLAPos& SatelliteLLA_1, VECTOR Satellite_1, double time) {
    SGP.CalculateLatLonAlt(time);
    SatelliteLLA_1.Lat = SGP.GetLat();
    SatelliteLLA_1.Lon = SGP.GetLon();
    SatelliteLLA_1.Alt = SGP.GetAlt();

    SatelliteLLA_1.Lat = SGP.DegToRad(SatelliteLLA_1.Lat);
    SatelliteLLA_1.Lon = SGP.DegToRad(SatelliteLLA_1.Lon);

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

    return { x0, y0, z0, 0.0 };
}


// Checking if velocity of satellite does not exceed 1.9 deg/min
bool VelocityCheck(LLAPos StationLLA, VECTOR Station, std::vector<CSGP4_SDP4>& SatelliteModel) {

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
    Satellite_1 = SatellitePos(SatelliteModel.back(), SatelliteLLA, Satellite_1, time);
    Satellite_1 = Transferring(Satellite_1, StationLLA);
    
    double theta = atan(sqrt(pow(Satellite_1.x + 6371.0 - Station.x, 2) + pow(Satellite_1.y - Station.y, 2)) / (Satellite_1.z - Station.z));
    theta = SatelliteModel[0].RadToDeg(theta);

    SatelliteModel.back().SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel.back(), SatelliteLLA, Satellite_2, time);
    Satellite_2 = Transferring(Satellite_2, StationLLA);

    double theta2 = atan(sqrt(pow(Satellite_2.x + 6371.0 - Station.x, 2) + pow(Satellite_2.y - Station.y, 2)) / (Satellite_2.z - Station.z));
    theta2 = SatelliteModel[0].RadToDeg(theta2);

    return (abs(theta2 - theta) >= 1.9);
}

//Function that save only satellites with suitable speed
void SatelliteFilter(std::string filename, LLAPos StationLLA, VECTOR Station, std::vector<CSGP4_SDP4>& SatelliteModel) {

    std::streampos here = 0;
    while (readFileAndCallFunction(filename, SatelliteModel, &here)) {

        if (VelocityCheck(StationLLA, Station, SatelliteModel)) {
            SatelliteModel.pop_back();
        }
    }
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


    double new_x = z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha));
    double new_y = y * cos(alpha) - x * sin(alpha);
    double new_z = (z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha))) - R;

    if (new_z <= 0) {
        return false;
    }

    fi = atan(new_y / new_x) + PI; // +PI to normalize lower and upper bounds
    theta = PI / 2 - acos(new_z / sqrt((pow(new_x, 2) + pow(new_y, 2) + pow(new_z, 2)))); // PI/2 to start from the bottom to the top
    
    if (fi <= maxAzm && fi >= minAzm) {
        if (theta <= maxElv && theta >= minElv) {
            return true;
        }
    }
    return false;
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

    //searching if Satellite will be in vision cone sometime in the time interval
    while ((check1 == false || tempTm <= lowerTm) && tempTm <= upperTm) { 
        tempTm = SatelliteModel.JulianDate(temp);
        SatelliteModel.SGP(tempTm);
        Satellite_1 = SatellitePos(SatelliteModel, SatelliteLLA_1, Satellite_1, tempTm);
        check1 = IntersectCheck(StationLLA_1, Satellite_1, minAzm, maxAzm, minElv, maxElv, theta, fi);
        temp.tm_sec += 1;
    }

    if (tempTm > upperTm || tempTm < lowerTm) {
        return 0;
    }
    
    double tmpfi = fi;
    double tmtheta = theta;
    
    // block for time filter
    // checking if satellite is in vision cone after minTimeObserve
    tm filterTime = temp;
    filterTime.tm_sec += timeMinObserveSec - 1;
    double checkTime = SatelliteModel.JulianDate(filterTime);
    SatelliteModel.SGP(checkTime);
    Satellite_1 = SatellitePos(SatelliteModel, SatelliteLLA_1, Satellite_1, checkTime);
    if (IntersectCheck(StationLLA_1, Satellite_1, minAzm, maxAzm, minElv, maxElv, tmtheta, tmpfi)) {
        return tempTm;
    }
    else
        return 0;
}

// Check direction by defining sign of the angle between station and satellite
bool getDiraction(VECTOR Station, CSGP4_SDP4 SatelliteModel) {
    VECTOR Satellite_1 = { 0, 0, 0, 0 }; // coordinates at the first moment
    VECTOR Satellite_2 = { 0, 0, 0, 0 }; // coordinates at the second moment
    LLAPos SatelliteLLA;

    double len = 0.0, len2 = 0.0; // to calculate radius vector

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm Time2 = now_tm;
    Time2.tm_min += 1;

    double time = SatelliteModel.JulianDate(now_tm);
    double time2 = SatelliteModel.JulianDate(Time2);

    SatelliteModel.SGP(time);
    Satellite_1 = SatellitePos(SatelliteModel, SatelliteLLA, Satellite_1, time);
    len = sqrt(pow(Satellite_1.x - Station.x, 2) + pow(Satellite_1.y - Station.y, 2) + pow(Satellite_1.z - Station.z, 2));

    SatelliteModel.SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel, SatelliteLLA, Satellite_2, time);
    len2 = sqrt(pow(Satellite_2.x - Station.x, 2) + pow(Satellite_2.y - Station.y, 2) + pow(Satellite_2.z - Station.z, 2));

    if (len2 - len < 0) {   // If diff < 0 that means satellite go closer
        return true;
    }
    else {                  // If diff > 0 that means satellite go further
        return false;
    }
}


int main() {
    std::string filename = "TLE.txt"; 
    std::vector<CSGP4_SDP4> SatelliteModel;
    std::vector<Satellite> Objects;
    

    double minAzm = 0, maxAzm = 0, minElv = 0, maxElv = 0;  // parameters of the station vision area
    int timeMinObserveSec = 0; // the time during which the satellite must be in the zone

    // station data preparation (start)
    LLAPos StationLLA_1;
    VECTOR Station_1{0.0, 0.0, 0.0, 0.0};

    if (!defStationPos(&StationLLA_1)) {
        return 0;
    }

    StationLLA_1.Lat = SatelliteModel[0].DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatelliteModel[0].DegToRad(StationLLA_1.Lon);
    Station_1 = ConvertLLAtoCoordinate(StationLLA_1);

    if (!defStationVision(minAzm, maxAzm, minElv, maxElv, timeMinObserveSec)) {
        return 0;
    }
    std::cout << minAzm << ' ' << maxAzm << ' ' << minElv << ' ' << maxElv << ' ' << timeMinObserveSec << std::endl;

    minAzm = SatelliteModel[0].DegToRad(minAzm);
    maxAzm = SatelliteModel[0].DegToRad(maxAzm);
    minElv = SatelliteModel[0].DegToRad(minElv);
    maxElv = SatelliteModel[0].DegToRad(maxElv);
    
    std::cout << minElv << " " << maxElv << "\t" << minAzm << " " << maxAzm << endl;
    // station data preparation (end)
    

    SatelliteFilter(filename, StationLLA_1, Station_1, SatelliteModel);    // pre-filter by satellite speed

    int N = SatelliteModel.size();  // number of satellites that passed verification by speed
    double time = 0;                // time when satellite will be in cone of view
    int satNumer = 0;               // NORAD satellite number
    double theta = 0;               // ELV
    double fi = 0;                  // AZM

    for (int i = 0; i < N; i++) {
        satNumer = SatelliteModel[i].GetNORAD();
        SATELLITE* s = (SATELLITE*)SatelliteModel[i].GetSatellite(); 
        time = TimeIntersect(SatelliteModel[i], StationLLA_1, minAzm, maxAzm, minElv, maxElv, timeMinObserveSec, theta, fi);
        if (time != 0) {
            bool dir = getDiraction(Station_1, SatelliteModel[i]); 
            SatelliteModel[i].SGP(time);
            SatelliteModel[i].CalculateLatLonAlt(time);
            double Lat = SatelliteModel[i].GetLat();
            double Lon = SatelliteModel[i].GetLon();
            double Alt = SatelliteModel[i].GetAlt();
            fi = SatelliteModel[i].RadToDeg(fi);
            theta = SatelliteModel[i].RadToDeg(theta);
            Objects.push_back({ s->cSatelliteName, satNumer, fi, theta, Lat, Lon, Alt, SatelliteModel[i].CalendarDate(time), SatelliteModel[i].GetVel().w, dir, time});
        }


    }

    std::sort(Objects.begin(), Objects.end(), [](const Satellite& a, const Satellite& b)
        {return a.timeStamp < b.timeStamp;}
    );

    std::cout << "Name" << "\t\t\t|" << "# Number" << "\t|" << "AZM   |" << "ELV  |" << "Lat   | " << "Lon    |" << "Alt    | " << "Velocity km/s\t | " << "Direction | " << "DD/MM/YYYY\thh:mm:ss" << std::endl;
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    for (const auto& item : Objects) {
        if (item.direction) {
            std::cout << item.SatelliteName << "\t|" << item.SatelliteNumber << " \t\t|" << item.AzmDeg << "|" << item.ElvDeg << "| " << item.Lat << " |" << item.Lon << " |" << item.Alt << " |" << item.velocity << "\t\t |\t+    | " << item.time.tm_mday << " / " << item.time.tm_mon + 1 << " / " << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << " : " << item.time.tm_sec << " " << std::endl;
        }
        else {
            std::cout << item.SatelliteName << "\t|" << item.SatelliteNumber << " \t\t|" << item.AzmDeg << "|" << item.ElvDeg << "| " << item.Lat << " |" << item.Lon << " |" << item.Alt << " |" << item.velocity << "\t\t |\t-    | " << item.time.tm_mday << " / " << item.time.tm_mon + 1 << " / " << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << " : " << item.time.tm_sec << " " << std::endl;
        }
    }

    return 0;
}
