#include<iostream>
#include <fstream>
#include <cstring>
#include "Sgpsdp.h"
#include <ctime>

#pragma warning(disable : 4996)

struct LLAPos // Структура содержащая геопозицию
{
    double Lat = 0, Lon = 0, Alt = 0; 
};

struct Sattelite
{
    std::string name;
    SYSTEMTIME time;
};

//Чтение TLE файла
void readFileAndCallFunction(const char* filename, CSGP4_SDP4 *SatteliteModel1) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Ошибка при открытии файла!" << std::endl;
        return;
    }

    char m_cLine0[256];
    char m_cLine1[256];
    char m_cLine2[256];
    
    // Чтение строк из файла
    file.getline(m_cLine0, sizeof(m_cLine0));
    file.getline(m_cLine1, sizeof(m_cLine1));
    file.getline(m_cLine2, sizeof(m_cLine2));
    // Закрытие файла
    file.close();

    // Создание объекта класса
    *SatteliteModel1 = CSGP4_SDP4(m_cLine0, m_cLine1, m_cLine2);
}

//Определение геопозиции станции
void defStationPos(LLAPos* station) {
    std::cout << "Enter the coordinates (latitude longitude altitude) of the station:";
    double Lat = 0.0, Lon = 0.0, Alt = 0.0;
    std::cin >> Lat >> Lon >> Alt;
    if ((Lat < -90 || Lat > 90) || (Lon < -180 || Lon > 180)) {
        std::cout << "ERROR: wrong parametrs(Lat or Lon)";
        return;
    }
    station->Lat = Lat;
    station->Lon = Lon;
    station->Alt = Alt;
    std::cout << "LLA: " << station->Lat << ':' << station->Lon << ':' << station->Alt << std::endl;
}

VECTOR ConvertLLAtoCoordinate(LLAPos *Object) {
    double x = (6371.0 + Object->Alt) * cos(Object->Lat) * cos(Object->Lon);
    double y = (6371.0 + Object->Alt) * cos(Object->Lat) * sin(Object->Lon);
    double z = (6371.0 + Object->Alt) * sin(Object->Lat);

    std::cout<< "Coordinates of sattelite: " << x << ' ' << y << ' ' << z << std::endl;

    return { x, y, z };
}
VECTOR SattelitePos(CSGP4_SDP4 *SGP, LLAPos *SatteliteLLA_1, VECTOR Sattelite_1, double time) {
    // Расчет и получение широты, долготы и высоты Спутника
    SGP->CalculateLatLonAlt(time);
    SatteliteLLA_1->Lat = SGP->GetLat();
    SatteliteLLA_1->Lon = SGP->GetLon();
    SatteliteLLA_1->Alt = SGP->GetAlt();
    std::cout << "Staelite geoposition(LLA): " << SatteliteLLA_1->Lat << "Lat; " << SatteliteLLA_1->Lon << "Lon; " << SatteliteLLA_1->Alt << "Alt" << std::endl;


    SatteliteLLA_1->Lat = SGP->DegToRad(SatteliteLLA_1->Lat);
    SatteliteLLA_1->Lon = SGP->DegToRad(SatteliteLLA_1->Lon);

    //Создание вектора с декартовыми координатами спутника
    Sattelite_1 = ConvertLLAtoCoordinate(SatteliteLLA_1);
    return Sattelite_1;
}

//Проверка вхождения спутника в область видимости станции связи
bool IntersectCheck(LLAPos StationLLA_1, VECTOR Sattelite_1){
    
    //Координаты спутника
    double x = Sattelite_1.x;
    double y = Sattelite_1.y;
    double z = Sattelite_1.z;

    double alpha = StationLLA_1.Lon;
    double beta = -StationLLA_1.Lat;
    double R = 6371.0 + StationLLA_1.Alt;
   
    double prepare = z * (-sin(beta)) + cos(beta) * (y * sin(alpha) + x * cos(alpha));

    if (prepare <= R)
        return false;
    //Коэффициент, который определяет угол видимости станции связи *10 градусов*
    double koef = tan(PI / 18.0); //от 10 до 90 градусов (10-170)
    //Уравнение конуса со смещением оси координат
    double intersection = pow(z*cos(beta)+sin(beta)*(y*sin(alpha)+x*cos(alpha)), 2) + pow(y*cos(alpha)-x*sin(alpha), 2) - pow(z*(-sin(beta))+cos(beta)*(y*sin(alpha)+x*cos(alpha)) - R, 2) / koef;
    //Проверка уравнения: Если меньше ноля, то спутник входит в область, Иначе не входит
    if(intersection < 0)
        return true;
    else
        return false;
};

bool VelosityCheck() {
    return 1;
}




int main() {
    //setlocale(LC_ALL, "");
    //Чтение файла TLE состоящего из 3-ч строк
    const char* filename = "TLE.txt"; // 
    CSGP4_SDP4 SatteliteModel1;
    readFileAndCallFunction(filename, &SatteliteModel1);
    
    bool check1 = 0; //Проверка вхождения спутника в область видимости станции через 30 секунд
    bool check2 = 0; //Проверка 


    //------------Определение времени-------------------
    //SYSTEMTIME currentTime;
    SYSTEMTIME TimeForVel;
    SYSTEMTIME Time;
    SYSTEMTIME Time2;

    //GetSystemTime(&currentTime); //Текущее время
    GetSystemTime(&Time); // Текущее время + 30 секунд
    GetSystemTime(&TimeForVel); // Текущее время + 90 сек
    GetSystemTime(&Time2); // Текущее время + 30 минут

    TimeForVel.wSecond = TimeForVel.wSecond + 30 + 60;
    Time.wSecond = Time.wSecond + 30;
    Time2.wMinute = Time2.wMinute + 30;

    double time = SatteliteModel1.JulianDate(Time); //нижняя граница
    double timePlusMin = SatteliteModel1.JulianDate(TimeForVel);
    double time2 = SatteliteModel1.JulianDate(Time2);//верхняя граница

    //------------Определение времени-------------------


    //Запуск модели
    SatteliteModel1.SGP(time);



    //--------------------Пространственное определение СПУТНИКА-------------------------
    LLAPos SatteliteLLA_1;
    VECTOR Sattelite_1 = {0, 0, 0, 0};
    Sattelite_1 = SattelitePos(&SatteliteModel1,&SatteliteLLA_1, Sattelite_1, time);
    //--------------------Пространственное определение СПУТНИКА-------------------------
    

    
    //--------------------Пространственное определение СТАНЦИИ-------------------------
    LLAPos StationLLA_1;
    VECTOR Station_1{0, 0, 0, 0};
    defStationPos(&StationLLA_1); //Ручной ввод и проверка на аддекватность вводимых данных

    //Перевод градусов в радианы для станции
    StationLLA_1.Lat = SatteliteModel1.DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatteliteModel1.DegToRad(StationLLA_1.Lon);

    //Создание вектора с декартовыми координатами станции
    Station_1 = ConvertLLAtoCoordinate(&StationLLA_1);
    //--------------------Пространственное определение СТАНЦИИ-------------------------


    //--------------------Проверка 1: вхождение в область видимости станции-------------------------
    check1 = IntersectCheck(StationLLA_1, Sattelite_1);
    std::cout<< "Intersection: " << check1 << std::endl;
    //--------------------Проверка 1: вхождение в область видимости станции-------------------------

    check2 = 0;
    if (check1) {
        
        //Запуск модели секунду после
        SatteliteModel1.SGP(timePlusMin);
        LLAPos SatteliteLLA_2;
        VECTOR Sattelite_2 = { 0, 0, 0, 0 };
        Sattelite_2 = SattelitePos(&SatteliteModel1, &SatteliteLLA_2, Sattelite_2, timePlusMin);


        //2 часть - работа со скоростями
        double DistStationTOSattelite1 = std::sqrt(pow(Station_1.x - Sattelite_1.x, 2) + pow(Station_1.y - Sattelite_1.y, 2) + pow(Station_1.z - Sattelite_1.z, 2));
        double DistStationTOSattelite2 = std::sqrt(pow(Station_1.x - Sattelite_2.x, 2) + pow(Station_1.y - Sattelite_2.y, 2) + pow(Station_1.z - Sattelite_2.z, 2));
        double DistSattelite1TOSattelite2 = std::sqrt(pow(Sattelite_1.x - Sattelite_2.x, 2) + pow(Sattelite_1.y - Sattelite_2.y, 2) + pow(Sattelite_1.z - Sattelite_2.z, 2));

        double alpha = acos((pow(DistStationTOSattelite1, 2) + pow(DistStationTOSattelite2, 2) - pow(DistSattelite1TOSattelite2, 2)) / (2 * DistStationTOSattelite1 * DistStationTOSattelite2));

        std::cout << "alpha radian = " << alpha << std::endl;
        alpha = SatteliteModel1.RadToDeg(alpha);
        std::cout << "alpha degree = " << alpha << std::endl;
        check2 = VelosityCheck();
    }
    return 0;
}