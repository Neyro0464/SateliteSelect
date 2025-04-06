#include "Satselect.h"
#pragma warning(disable : 4996)

using namespace std;
using namespace std::chrono;
using namespace Sgpsdp;

SetelliteSelect::SetelliteSelect(const std::string& filename, Station& param)
{
    Satellites = {};
    this->filename = filename;
    param.geo.Lon = DegToRad(param.geo.Lon);
    param.geo.Lat = DegToRad(param.geo.Lat);
    param.lim.minAzm = DegToRad(param.lim.minAzm);
    param.lim.maxAzm = DegToRad(param.lim.maxAzm);
    param.lim.minElv = DegToRad(param.lim.minElv);
    param.lim.maxElv = DegToRad(param.lim.maxElv);
    station = param;
    dataPrepare();
}

double SetelliteSelect::DegToRad(double arg) {
    return (double)(arg / 360.0 * (2.0 * PI));
}
double SetelliteSelect::RadToDeg(double arg) {
    return (double)(arg / (2.0 * PI) * 360.0);
}
// Convertation LLA coordinates to Decart coordinates
// previously converted to radians!!!!
// (0;0;0) is center of planet
// (6371;0;0) is lan = 0, lon = 0, alt = 0
CoordDecart SetelliteSelect::ConvertGEOtoDecart(CoordGeodetic& Object) {
    double x = (6371.0 + Object.Alt) * cos(Object.Lat) * cos(Object.Lon);
    double y = (6371.0 + Object.Alt) * cos(Object.Lat) * sin(Object.Lon);
    double z = (6371.0 + Object.Alt) * sin(Object.Lat);
    Sgpsdp::CoordDecart result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}

// Transfer satellite coordinate system in station coodinate system
CoordDecart SetelliteSelect::Transferring(CoordDecart object, CoordGeodetic Station) {
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
    Sgpsdp::CoordDecart result;
    result.x = x0;
    result.y = y0;
    result.z = z0;
    return result;
}

// Calculate satellite position with class function
// then convert to Decart coordinates 
CoordDecart SetelliteSelect::SatellitePos(CSGP4_SDP4 SGP, CoordGeodetic& SatelliteLLA_1, CoordDecart Satellite_1, double time) {

    SGP.CalculateLatLonAlt(time);
    SatelliteLLA_1.Lat = SGP.GetLat();
    SatelliteLLA_1.Lon = SGP.GetLon();
    SatelliteLLA_1.Alt = SGP.GetAlt();

    SatelliteLLA_1.Lat = DegToRad(SatelliteLLA_1.Lat);
    SatelliteLLA_1.Lon = DegToRad(SatelliteLLA_1.Lon);

    Satellite_1 = ConvertGEOtoDecart(SatelliteLLA_1);
    return Satellite_1;
}

// Checking if velocity of satellite does not exceed 1.9 deg/min
bool SetelliteSelect::VelocityCheck(Station station, CSGP4_SDP4 SatelliteModel) {

    CoordDecart Satellite_1{};
    CoordDecart Satellite_2 {};
    CoordGeodetic SatelliteGEO{};

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Working in UTC time

    tm Time2 = now_tm;
    Time2.tm_min += 1;

    double time = SatelliteModel.JulianDate(now_tm);
    double time2 = SatelliteModel.JulianDate(Time2);

    SatelliteModel.SGP(time);
    Satellite_1 = SatellitePos(SatelliteModel, SatelliteGEO, Satellite_1, time);
    Satellite_1 = Transferring(Satellite_1, station.geo);

    double theta = atan(sqrt(pow(Satellite_1.x + 6371.0 - station.dec.x, 2) + pow(Satellite_1.y - station.dec.y, 2)) / (Satellite_1.z - station.dec.z));
    theta = RadToDeg(theta);

    SatelliteModel.SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel, SatelliteGEO, Satellite_2, time);
    Satellite_2 = Transferring(Satellite_2, station.geo);

    double theta2 = atan(sqrt(pow(Satellite_2.x + 6371.0 - station.dec.x, 2) + pow(Satellite_2.y - station.dec.y, 2)) / (Satellite_2.z - station.dec.z));
    theta2 = RadToDeg(theta2);

    return (abs(theta2 - theta) >= 1.9);
}

// Reading TLE file and make new class objects
//void SetelliteSelect::dataPrepare(std::vector<CSGP4_SDP4>& Satellites) {
//    std::ifstream file(filename);
//    if (!file.is_open()) {
//        std::cout << "ERRO: File is not open " << std::endl;
//        abort();
//    }
//
//    char m_cLine0[70];
//    char m_cLine1[70];
//    char m_cLine2[70];
//
//    while (file.getline(m_cLine0, sizeof(m_cLine0)) &&
//        file.getline(m_cLine1, sizeof(m_cLine1)) &&
//        file.getline(m_cLine2, sizeof(m_cLine2))) {
//        Satellites.push_back(CSGP4_SDP4(m_cLine0, m_cLine1, m_cLine2));
//        if (VelocityCheck(station, Satellites.back())) {
//            Satellites.pop_back();
//        }
//    }
//    file.close();
//}

void SetelliteSelect::dataPrepare() {
    auto Sat = Satellites;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "ERROR: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::cout << "File opened successfully, attempting to read" << std::endl;

    std::string m_cLine0, m_cLine1, m_cLine2;
    int count = 0;
    int totalLines = 0;

    while (true) {
        // ������ ��� ��������
        if (!std::getline(file, m_cLine0)) break;  // ����� ��� ���������� ����� ����� ��� ������
        if (m_cLine0.empty()) {
            std::cout << "Skipping empty line" << std::endl;
            continue;
        }

        // ������ ������ ������ TLE (������ ���������� � "1")
        if (!std::getline(file, m_cLine1) || m_cLine1.empty() || m_cLine1[0] != '1') {
            std::cout << "Invalid or missing Line 1 after: " << m_cLine0 << std::endl;
            break;
        }

        // ������ ������ ������ TLE (������ ���������� � "2")
        if (!std::getline(file, m_cLine2) || m_cLine2.empty() || m_cLine2[0] != '2') {
            std::cout << "Invalid or missing Line 2 after: " << m_cLine0 << " | " << m_cLine1 << std::endl;
            break;
        }
        totalLines++;

        // ����������� std::string � char[] ��� CSGP4_SDP4
        char line0[70], line1[70], line2[70];
        strncpy(line0, m_cLine0.c_str(), sizeof(line0) - 1); line0[sizeof(line0) - 1] = '\0';
        strncpy(line1, m_cLine1.c_str(), sizeof(line1) - 1); line1[sizeof(line1) - 1] = '\0';
        strncpy(line2, m_cLine2.c_str(), sizeof(line2) - 1); line2[sizeof(line2) - 1] = '\0';

        try {
            Sat.push_back(CSGP4_SDP4(line0, line1, line2));
            if (VelocityCheck(station, Sat.back())) {
                Sat.pop_back();
            }
            else {
                count++;
            }
        }
        catch (const std::exception& e) {
            std::cout << "Error creating CSGP4_SDP4: " << e.what() << std::endl;
        }
    }

    file.close();
    Satellites = Sat;
    std::cout << "Total TLE sets read: " << totalLines << std::endl;
    std::cout << "Loaded " << count << " satellites after velocity check" << std::endl;
}


// Define station possition and check constraints
bool SetelliteSelect::SetStationPos(Station& station) {
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
    station.geo.Lat = DegToRad(Lat);
    station.geo.Lon = DegToRad(Lon);
    station.geo.Alt = Alt;
    station.dec = ConvertGEOtoDecart(station.geo);
    return true;
}

//Enter station vision parametres
bool SetelliteSelect::SetFilter(Station& station)
{
    std::cout << "Enter vision parameters of the station (minAzm maxAzm minElv maxElv time): ";
    std::cin >> station.lim.minAzm >> station.lim.maxAzm;
    if (station.lim.minAzm < 0 || station.lim.minAzm > station.lim.maxAzm || station.lim.maxAzm >= 360) {
        std::cout << "ERROR: wrong parameters of station - Azimut";
        return false;
    }
    std::cin >> station.lim.minElv >> station.lim.maxElv;
    if (station.lim.minElv < 10 || station.lim.maxElv < station.lim.minElv || station.lim.maxElv > 180) {
        std::cout << "ERROR: wrong parameters of station - Elevation";
        return false;
    }
    std::cin >> station.lim.timeMinObserveSec;
    if (station.lim.timeMinObserveSec < 1) {
        std::cout << "ERROR: wrong parameters of station - time";
        return false;
    }
    std::cout << std::endl;
    station.lim.minAzm = DegToRad(station.lim.minAzm);
    station.lim.maxAzm = DegToRad(station.lim.maxAzm);
    station.lim.minElv = DegToRad(station.lim.minElv);
    station.lim.maxElv = DegToRad(station.lim.maxElv);
    return true;
}






// Checking intersection: if satellite in view area of station
bool SetelliteSelect::IntersectCheck(CoordDecart& Satellite_1, Station station, double& theta, double& fi) {
    double x = Satellite_1.x;
    double y = Satellite_1.y;
    double z = Satellite_1.z;

    // longtitude and latitude of the station to rotate around the sphere(Earth) 
    double alpha = station.geo.Lon;
    double beta = station.geo.Lat;
    double R = 6371.0 + station.geo.Alt;


    double new_x = z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha));
    double new_y = y * cos(alpha) - x * sin(alpha);
    double new_z = (z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha))) - R;

    if (new_z <= 0) {
        return false;
    }

    fi = atan(new_y / new_x) + PI; // +PI to normalize lower and upper bounds
    theta = PI / 2 - acos(new_z / sqrt((pow(new_x, 2) + pow(new_y, 2) + pow(new_z, 2)))); // PI/2 to start from the bottom to the top

    if (fi <= station.lim.maxAzm && fi >= station.lim.minAzm) {
        if (theta <= station.lim.maxElv && theta >= station.lim.minElv) {
            return true;
        }
    }
    return false;
}


// checks the intersection of the satellite with the station's visibility cone (+30 sec to +30 min)
//  Filter satellites that would be in cone of visibility for timeMinObserve(sec)
double SetelliteSelect::TimeIntersect(CSGP4_SDP4& SatelliteModel, Station& station, double& theta, double& fi) { //remake
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

    CoordGeodetic SatelliteLLA_1;
    CoordDecart Satellite_1{};
    bool check1 = 0;

    //searching if Satellite will be in vision cone sometime in the time interval
    while ((check1 == false || tempTm <= lowerTm) && tempTm <= upperTm) {
        tempTm = SatelliteModel.JulianDate(temp);
        SatelliteModel.SGP(tempTm);
        Satellite_1 = SatellitePos(SatelliteModel, SatelliteLLA_1, Satellite_1, tempTm);
        check1 = IntersectCheck(Satellite_1, station, theta, fi);
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
    filterTime.tm_sec += station.lim.timeMinObserveSec - 1;
    double checkTime = SatelliteModel.JulianDate(filterTime);
    SatelliteModel.SGP(checkTime);
    Satellite_1 = SatellitePos(SatelliteModel, SatelliteLLA_1, Satellite_1, checkTime);
    if (IntersectCheck(Satellite_1, station, tmtheta, tmpfi)) {
        return tempTm;
    }
    else
        return 0;
}

// Check direction by defining sign of the angle between station and satellite
bool SetelliteSelect::getDiraction(Station station, CSGP4_SDP4& SatelliteModel) {
    CoordDecart Satellite_1{}; // coordinates at the first moment
    CoordDecart Satellite_2{}; // coordinates at the second moment
    CoordGeodetic SatelliteLLA;

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
    len = sqrt(pow(Satellite_1.x - station.dec.x, 2) + pow(Satellite_1.y - station.dec.y, 2) + pow(Satellite_1.z - station.dec.z, 2));

    SatelliteModel.SGP(time2);
    Satellite_2 = SatellitePos(SatelliteModel, SatelliteLLA, Satellite_2, time);
    len2 = sqrt(pow(Satellite_2.x - station.dec.x, 2) + pow(Satellite_2.y - station.dec.y, 2) + pow(Satellite_2.z - station.dec.z, 2));

    if (len2 - len < 0) {   // If diff < 0 that means satellite go closer
        return true;
    }
    else {                  // If diff > 0 that means satellite go further
        return false;
    }
}

const std::vector<NORAD_DATA> SetelliteSelect::GetSatArray() {

    std::vector<NORAD_DATA> Objects;
    COORDS dummy;
    Station stationtmp = station;
    int N = Satellites.size();      // number of satellites that passed verification by speed
    double time = 0;                // time when satellite will be in cone of view
    unsigned int satNumer = 0;      // NORAD satellite number
    double theta = 0;               // ELV
    double fi = 0;                  // AZM
    std::stringstream onTime1;
    tm timeDate;
    bool dir = 0;
    for (int i = 0; i < N; i++) {
        time = TimeIntersect(Satellites[i], stationtmp, theta, fi);
        if (time != 0) {
            NORAD_DATA dummy{};
            dummy.noradNumber = Satellites[i].GetNORAD();
            SATELLITE* s = (SATELLITE*)Satellites[i].GetSatellite();
            dummy.dirPositive = getDiraction(stationtmp, Satellites[i]);
            Satellites[i].SGP(time);
            Satellites[i].CalculateLatLonAlt(time);
            dummy.coords.geo.Lat = Satellites[i].GetLat();
            dummy.coords.geo.Lon = Satellites[i].GetLon();
            dummy.coords.geo.Alt = Satellites[i].GetAlt();
            dummy.coords.topo.azm = Satellites[i].RadToDeg(fi);
            dummy.coords.topo.elv = Satellites[i].RadToDeg(theta);
            timeDate = Satellites[i].CalendarDate(time);
            dummy.name = s->cSatelliteName;
            onTime1 << std::to_string(timeDate.tm_mday)
                << "/" << std::to_string(timeDate.tm_mon + 1)
                << "/" << std::to_string(timeDate.tm_year + 1900)
                << "\t" << std::to_string(timeDate.tm_hour)
                << ":" << std::to_string(timeDate.tm_min)
                << ":" << std::to_string(timeDate.tm_sec);
            dummy.onTime = onTime1.str();
            dummy.time = time;

            Objects.push_back(dummy);
            onTime1.str("");
        }

    }

    std::sort(Objects.begin(), Objects.end(), [](const NORAD_DATA& a, const NORAD_DATA& b)
        {return a.time < b.time; }
    );

    return Objects;
}

void SetelliteSelect::showSat(const std::vector<NORAD_DATA>& SatellitesRes) {
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    for (const auto& item : SatellitesRes) {

        std::cout << item.name << "|" << item.noradNumber << " |" << item.coords.topo.azm << "|" << item.coords.topo.elv << "| "
            << item.coords.geo.Lat << " |" << item.coords.geo.Lon << " |" << item.coords.geo.Alt << " | "
            << item.dirPositive << " | " << item.onTime << std::endl;
    }
}


