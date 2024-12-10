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

const int N = 1000;

struct LLAPos // Структура содержащая геопозицию
{
    double Lat = 0, Lon = 0, Alt = 0;
};

struct Sattelite
{
    int SatteliteNumber;
    tm time;
    double velocity;
    double timeStamp;
};

// Чтение TLE файла
bool readFileAndCallFunction(const char* filename, CSGP4_SDP4* SatteliteModel, std::streampos* here, int i) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ERRO: File is not open " << std::endl;
        return false; // Возвращаем false, если файл не удалось открыть
    }

    char m_cLine0[70];
    char m_cLine1[70];
    char m_cLine2[70];
    file.seekg(*here);
    // Чтение строк из файла
    if (file.getline(m_cLine0, sizeof(m_cLine0)) &&
        file.getline(m_cLine1, sizeof(m_cLine1)) &&
        file.getline(m_cLine2, sizeof(m_cLine2))) {
        // Создание объекта класса
        SatteliteModel[i] = CSGP4_SDP4(m_cLine0, m_cLine1, m_cLine2);
        *here = file.tellg();
        file.close();
        return true; // Возвращаем true, если данные успешно считаны
    }

    file.close();
    return false; // Возвращаем false, если не удалось считать данные
}

// Определение геопозиции станции
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

VECTOR ConvertLLAtoCoordinate(LLAPos* Object) {
    double x = (6371.0 + Object->Alt) * cos(Object->Lat) * cos(Object->Lon);
    double y = (6371.0 + Object->Alt) * cos(Object->Lat) * sin(Object->Lon);
    double z = (6371.0 + Object->Alt) * sin(Object->Lat);
    return { x, y, z };
}

VECTOR SattelitePos(CSGP4_SDP4 SGP, LLAPos* SatteliteLLA_1, VECTOR Sattelite_1, double time) {
    // Расчет и получение широты, долготы и высоты Спутника
    SGP.CalculateLatLonAlt(time);
    SatteliteLLA_1->Lat = SGP.GetLat();
    SatteliteLLA_1->Lon = SGP.GetLon();
    SatteliteLLA_1->Alt = SGP.GetAlt();

    SatteliteLLA_1->Lat = SGP.DegToRad(SatteliteLLA_1->Lat);
    SatteliteLLA_1->Lon = SGP.DegToRad(SatteliteLLA_1->Lon);

    // Создание вектора с декартовыми координатами спутника
    Sattelite_1 = ConvertLLAtoCoordinate(SatteliteLLA_1);
    return Sattelite_1;
}

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

// Проверка вхождения спутника в область видимости станции связи
bool IntersectCheck(LLAPos StationLLA_1, VECTOR Sattelite_1) {
    double x = Sattelite_1.x;
    double y = Sattelite_1.y;
    double z = Sattelite_1.z;

    double alpha = StationLLA_1.Lon;
    double beta = StationLLA_1.Lat;
    double R = 6371.0 + StationLLA_1.Alt;

    double prepare = z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha));

    if (prepare <= R)
        return false;

    double koef = tan(PI / 18.0); // от 10 до 90 градусов (10-170)
    double intersection = pow(z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha)), 2) + pow(y * cos(alpha) - x * sin(alpha), 2) - pow(z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha)) - R, 2) / koef;

    if (intersection < 0)
        return true;
    else
        return false;
}

void VelocityCheck(const char* filename, LLAPos Station, CSGP4_SDP4* SatteliteModel) {
    VECTOR Sattelite_1 = { 0, 0, 0, 0 };
    VECTOR Sattelite_2 = { 0, 0, 0, 0 };
    LLAPos SatteliteLLA;
    VECTOR StationXYZ = ConvertLLAtoCoordinate(&Station);

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Используем gmtime для получения времени в UTC

    tm Time2 = now_tm;
    Time2.tm_min += 1;

    double time = SatteliteModel[0].JulianDate(now_tm);
    double time2 = SatteliteModel[0].JulianDate(Time2);

    std::streampos here = 0;
    int i = 0;
    while (readFileAndCallFunction(filename, SatteliteModel, &here, i)) {
        SatteliteModel[i].SGP(time);
        Sattelite_1 = SattelitePos(SatteliteModel[i], &SatteliteLLA, Sattelite_1, time);
        Sattelite_1 = Transferring(Sattelite_1, Station);

        double theta = atan(sqrt(pow(Sattelite_1.x, 2) + pow(Sattelite_1.y, 2)) / Sattelite_1.z);
        theta = SatteliteModel->RadToDeg(theta);

        SatteliteModel[i].SGP(time2);
        Sattelite_2 = SattelitePos(SatteliteModel[i], &SatteliteLLA, Sattelite_2, time);
        Sattelite_2 = Transferring(Sattelite_2, Station);

        double theta2 = atan(sqrt(pow(Sattelite_2.x, 2) + pow(Sattelite_2.y, 2)) / Sattelite_2.z);
        theta2 = SatteliteModel->RadToDeg(theta2);

        if (abs(theta2 - theta) < 1.9)
            i++;
    }
}

double TimeIntersect(CSGP4_SDP4 SatteliteModel, LLAPos StationLLA_1) {
    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_tm = *gmtime(&now_time_t); // Используем gmtime для получения времени в UTC

    tm temp = now_tm;
    tm Time2 = now_tm;
    Time2.tm_min += 30;
    temp.tm_sec += 30;

    double time = SatteliteModel.JulianDate(now_tm); // нижняя граница
    double tempTime = SatteliteModel.JulianDate(temp); // Временной промежуток
    double time2 = SatteliteModel.JulianDate(Time2); // верхняя граница

    LLAPos SatteliteLLA_1;
    VECTOR Sattelite_1 = { 0, 0, 0, 0 };
    bool check1 = 0;

    while ((check1 == false || tempTime <= time) && tempTime <= time2) {
        tempTime = SatteliteModel.JulianDate(temp);
        SatteliteModel.SGP(tempTime);
        Sattelite_1 = SattelitePos(SatteliteModel, &SatteliteLLA_1, Sattelite_1, tempTime);

        check1 = IntersectCheck(StationLLA_1, Sattelite_1);
        temp.tm_sec += 1;
    }

    if (tempTime > time2 || tempTime < time) {
        return 0;
    }
    return tempTime;
}

int main() {
    const char* filename = "TLE.txt";
    CSGP4_SDP4 SatteliteModel[N];
    std::vector<Sattelite> Objects;

    LLAPos StationLLA_1;
    VECTOR Station_1{ 0, 0, 0, 0 };
    if (!defStationPos(&StationLLA_1))
        return 0;

    StationLLA_1.Lat = SatteliteModel[0].DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatteliteModel[0].DegToRad(StationLLA_1.Lon);

    Station_1 = ConvertLLAtoCoordinate(&StationLLA_1);

    VelocityCheck(filename, StationLLA_1, SatteliteModel);

    double time = 0;
    int test = 0;
    for (int i = 0; i < N; i++) {
        test = SatteliteModel[i].GetNORAD();
        if (test != 0 && test != 32763) {
            time = TimeIntersect(SatteliteModel[i], StationLLA_1);
            if (time != 0) {
                Objects.push_back({ test, SatteliteModel[i].CalendarDate(time), SatteliteModel[i].GetVel().w, time });
            }
        }
    }

    std::sort(Objects.begin(), Objects.end(), [](const Sattelite& a, const Sattelite& b) {
        return a.timeStamp < b.timeStamp;
        });

    std::cout << "# Sattelite Number" << " \t| " << "velocity km/s\t | " << "Datetime\tDD/MM/YYYY\thh:mm:ss" << std::endl;

    for (const auto& item : Objects) {
        std::cout << "# " << item.SatteliteNumber << " \t\t| " << item.velocity << "\t | \t\t" << item.time.tm_mday << "/" << item.time.tm_mon + 1 << "/" << item.time.tm_year + 1900 << "\t" << item.time.tm_hour << ":" << item.time.tm_min << ":" << item.time.tm_sec << " " << std::endl;
    }
    //system("pause");
    return 0;
}
