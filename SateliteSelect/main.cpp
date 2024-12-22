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
    tm time;
    double velocity;
    double timeStamp;
    bool direction; 
};

// Reading TLE file and make new class objects
bool readFileAndCallFunction(std::string filename, std::vector<CSGP4_SDP4> &SatelliteModel, std::streampos* here) {
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
bool IntersectCheck(LLAPos StationLLA_1, VECTOR Satellite_1) {
    double x = Satellite_1.x;
    double y = Satellite_1.y;
    double z = Satellite_1.z;

    double alpha = StationLLA_1.Lon;
    double beta = StationLLA_1.Lat;
    double R = 6371.0 + StationLLA_1.Alt;

    double prepare = z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha));

    if (prepare <= R)
        return false;

    //coefficient for z
    //double koef = tan(PI / 18.0);
    //double intersection = pow(z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha)), 2) + pow(y * cos(alpha) - x * sin(alpha), 2) - pow(z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha)) - R, 2) / koef;
    
    double koef = tan(PI/18.0);
    double intersection = pow(z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha)), 2) + pow(y * cos(alpha) - x * sin(alpha), 2) - pow(z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha)) - R, 2) / koef;

    if (intersection < 0)
        return true;
    else
        return false;
}

void VelocityCheck(std::string filename, LLAPos Station, std::vector<CSGP4_SDP4> &SatelliteModel) {
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

    std::streampos here = 0;
    while (readFileAndCallFunction(filename, SatelliteModel, &here)) {
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

        if (abs(theta2 - theta) >= 1.9){
            SatelliteModel.pop_back();
        }
    }
}

double TimeIntersect(CSGP4_SDP4 SatelliteModel, LLAPos StationLLA_1) {
    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm temp = now_tm;
    tm Time2 = now_tm;
    Time2.tm_min += 30;
    temp.tm_sec += 30;

    double time = SatelliteModel.JulianDate(now_tm); // lower bound of time interval
    double tempTime = SatelliteModel.JulianDate(temp); // time cursor that move between lower and upper bounds of time interval
    double time2 = SatelliteModel.JulianDate(Time2); // upper bound of time interval

    LLAPos SatelliteLLA_1;
    VECTOR Satellite_1 = { 0, 0, 0, 0 };
    bool check1 = 0;

    while ((check1 == false || tempTime <= time) && tempTime <= time2) {
        tempTime = SatelliteModel.JulianDate(temp);
        SatelliteModel.SGP(tempTime);
        Satellite_1 = SatellitePos(SatelliteModel, &SatelliteLLA_1, Satellite_1, tempTime);

        check1 = IntersectCheck(StationLLA_1, Satellite_1);
        temp.tm_sec += 1;
    }

    if (tempTime > time2 || tempTime < time) {
        return 0;
    }
    return tempTime;
}

bool getDiraction(LLAPos Station, CSGP4_SDP4 SatelliteModel){
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
    theta = SatelliteModel.RadToDeg(theta);
    SatelliteModel.SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel, &SatelliteLLA, Satellite_2, time);
    Satellite_2 = Transferring(Satellite_2, Station);
    double theta2 = atan(sqrt(pow(Satellite_2.x, 2) + pow(Satellite_2.y, 2)) / Satellite_2.z);
    theta2 = SatelliteModel.RadToDeg(theta2);
    if (theta2 - theta > 0){ // If diff > 0 that means satellite go higher and closer
        return true;
    }   
    else {                  // If diff <= 0 that means satellite go lower and further
        return false;
    }
}


int main() {
    std::string filename = "TLE.txt"; // std::string
    std::vector<CSGP4_SDP4> SatelliteModel; // REMAKE in dinamic storage
    std::vector<Satellite> Objects;

    LLAPos StationLLA_1;
    VECTOR Station_1{ 0, 0, 0, 0 };
    if (!defStationPos(&StationLLA_1))
        return 0;

    StationLLA_1.Lat = SatelliteModel[0].DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatelliteModel[0].DegToRad(StationLLA_1.Lon);


    VelocityCheck(filename, StationLLA_1, SatelliteModel);

    int N = SatelliteModel.size();
    double time = 0;
    int test = 0;
    for (int i = 0; i < N; i++) {
        test = SatelliteModel[i].GetNORAD();
        SATELLITE* s = (SATELLITE*)SatelliteModel[i].GetSatellite(); // 3LINE
        
        time = TimeIntersect(SatelliteModel[i], StationLLA_1);
        if (time != 0) {
            Objects.push_back({ s->cSatelliteName, test, SatelliteModel[i].CalendarDate(time), SatelliteModel[i].GetVel().w, time, getDiraction(StationLLA_1, SatelliteModel[i])});
        }
        
        
    }

    std::sort(Objects.begin(), Objects.end(), [](const Satellite& a, const Satellite& b) {
        return a.timeStamp < b.timeStamp;
        });

    std::cout << "# Satellite Name " << " \t | "<< "# Satellite Number" << " \t | " << "Velocity km/s\t | " << "Direction\t| " << "DD/MM/YYYY\thh:mm:ss" << std::endl; // 3LINE

    for (const auto& item : Objects) {
        if(item.direction){
            std::cout << "# " << item.SatelliteName << " | " << item.SatelliteNumber << " \t\t | " << item.velocity << "\t\t |\t+\t| " << item.time.tm_mday << "/" << item.time.tm_mon + 1 << "/" << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << ":" << item.time.tm_sec << " " << std::endl;
        }
        else{
            std::cout << "# " << item.SatelliteName << " | " << item.SatelliteNumber << " \t\t | " << item.velocity << "\t\t |\t-\t| " << item.time.tm_mday << "/" << item.time.tm_mon + 1 << "/" << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << ":" << item.time.tm_sec << " " << std::endl;
        }
    }
    //system("pause");
    return 0;
}