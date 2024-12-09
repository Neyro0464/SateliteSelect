#include<iostream>
#include <fstream>
#include <cstring>
#include "Sgpsdp.h"
#include <ctime>
#include <vector>
#include <algorithm>

#pragma warning(disable : 4996)

const int N = 1000;

struct LLAPos // Структура содержащая геопозицию
{
    double Lat = 0, Lon = 0, Alt = 0; 
};

struct Sattelite
{
    int SatteliteNumber;
    SYSTEMTIME time;
    double velocity;
    double timeStamp;
};

//void CapacityUP(CSGP4_SDP4* SatteliteModel) {
//    CSGP4_SDP4* SatteliteModelNEW;
//    SatteliteModelNEW = new CSGP4_SDP4[N + N];
//    std::copy(SatteliteModel, SatteliteModel + N, SatteliteModelNEW);
//    delete[] SatteliteModel;
//    SatteliteModel = SatteliteModelNEW;
//}

//Чтение TLE файла
bool readFileAndCallFunction(const char* filename, CSGP4_SDP4 *SatteliteModel, std::streampos *here, int i) {
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


//Определение геопозиции станции
bool defStationPos(LLAPos* station) {
    std::cout << "Enter the coordinates (latitude longitude altitude) of the station:";
    double Lat = 0.0, Lon = 0.0, Alt = 0.0;
    std::cin >> Lat >> Lon >> Alt;
    if ((Lat < -90 || Lat > 90) || (Lon < -180 || Lon > 180)) {
        std::cout << "ERROR: wrong parametrs(Lat or Lon)";
        return false;
    }
    station->Lat = Lat;
    station->Lon = Lon;
    station->Alt = Alt;
    //std::cout << "LLA: " << station->Lat << ':' << station->Lon << ':' << station->Alt << std::endl;
    return true;
}

VECTOR ConvertLLAtoCoordinate(LLAPos *Object) {
    double x = (6371.0 + Object->Alt) * cos(Object->Lat) * cos(Object->Lon);
    double y = (6371.0 + Object->Alt) * cos(Object->Lat) * sin(Object->Lon);
    double z = (6371.0 + Object->Alt) * sin(Object->Lat);
    return { x, y, z };
}
VECTOR SattelitePos(CSGP4_SDP4 SGP, LLAPos *SatteliteLLA_1, VECTOR Sattelite_1, double time) {
    // Расчет и получение широты, долготы и высоты Спутника
    SGP.CalculateLatLonAlt(time);
    SatteliteLLA_1->Lat = SGP.GetLat();
    SatteliteLLA_1->Lon = SGP.GetLon();
    SatteliteLLA_1->Alt = SGP.GetAlt();
    //std::cout << "Staelite geoposition(LLA): " << SatteliteLLA_1->Lat << "Lat; " << SatteliteLLA_1->Lon << "Lon; " << SatteliteLLA_1->Alt << "Alt" << std::endl;


    SatteliteLLA_1->Lat = SGP.DegToRad(SatteliteLLA_1->Lat);
    SatteliteLLA_1->Lon = SGP.DegToRad(SatteliteLLA_1->Lon);

    //Создание вектора с декартовыми координатами спутника
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
    //std::cout << x << ", " << y << ", " << z << std::endl;

    /*double z0 = z * cos(beta) - sin(beta) * (y * sin(alpha) + x * cos(alpha));
    double y0 = y * cos(alpha) - x * sin(alpha);
    double x0 = z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha));*/
    
    double x0 = cos(beta) * (x * cos(alpha) - y * sin(alpha)) - z * sin(beta);
    double y0 = x * sin(alpha) + y * cos(alpha);
    double z0 = sin(beta) * (x * cos(alpha) - y * sin(alpha)) + z * cos(beta);
    
    return { x0, y0, z0 };

}

//Проверка вхождения спутника в область видимости станции связи
bool IntersectCheck(LLAPos StationLLA_1, VECTOR Sattelite_1){
    
    //Координаты спутника
    double x = Sattelite_1.x;
    double y = Sattelite_1.y;
    double z = Sattelite_1.z;

    double alpha = StationLLA_1.Lon;
    double beta = StationLLA_1.Lat;
    double R = 6371.0 + StationLLA_1.Alt;
   
    double prepare = z * sin(beta) + cos(beta) * (y * sin(alpha) + x * cos(alpha));

    if (prepare <= R)
        return false;
    //Коэффициент, который определяет угол видимости станции связи *10 градусов*
    double koef = tan(PI / 18.0); //от 10 до 90 градусов (10-170)
    //Уравнение конуса со смещением оси координат
    double intersection = pow(z*cos(beta)-sin(beta)*(y*sin(alpha)+x*cos(alpha)), 2) + pow(y*cos(alpha)-x*sin(alpha), 2) - pow(z*sin(beta)+cos(beta)*(y*sin(alpha)+x*cos(alpha)) - R, 2) / koef;
    //Проверка уравнения: Если меньше ноля, то спутник входит в область, Иначе не входит
    if(intersection < 0)
        return true;
    else
        return false;
};

void VelocityCheck(const char* filename, LLAPos Station, CSGP4_SDP4* SatteliteModel) {

    VECTOR Sattelite_1 = { 0, 0, 0, 0 };
    VECTOR Sattelite_2 = { 0, 0, 0, 0 };
    LLAPos SatteliteLLA;
    VECTOR StationXYZ = ConvertLLAtoCoordinate(&Station);

    SYSTEMTIME Time;
    SYSTEMTIME Time2;
    GetSystemTime(&Time);
    Time2 = Time;
    Time2.wMinute += 1;
    double time = SatteliteModel[0].JulianDate(Time);
    double time2 = SatteliteModel[0].JulianDate(Time2);


    std::streampos here = 0;
    int i = 0;
    while (readFileAndCallFunction(filename, SatteliteModel, &here, i)) {
      /*  if (i % N == 0) {
            CapacityUP(SatteliteModel);
        }*/
        SatteliteModel[i].SGP(time);
        Sattelite_1 = SattelitePos(SatteliteModel[i], &SatteliteLLA, Sattelite_1, time);
        //std::cout << Sattelite_1.x << ", " << Sattelite_1.y << ", " << Sattelite_1.z << std::endl;
        Sattelite_1 = Transferring(Sattelite_1, Station);
        //std::cout << Sattelite_1.x << ", " << Sattelite_1.y << ", " << Sattelite_1.z << std::endl;

        double theta = atan(sqrt(pow(Sattelite_1.x, 2) + pow(Sattelite_1.y, 2)) / Sattelite_1.z);
        theta = SatteliteModel->RadToDeg(theta);

        SatteliteModel[i].SGP(time2);
        Sattelite_2 = SattelitePos(SatteliteModel[i], &SatteliteLLA, Sattelite_2, time);
        //std::cout << Sattelite_1.x << ", " << Sattelite_1.y << ", " << Sattelite_1.z << std::endl;
        Sattelite_2 = Transferring(Sattelite_2, Station);
        //std::cout << Sattelite_2.x << ", " << Sattelite_2.y << ", " << Sattelite_2.z << std::endl;
        
        double theta2 = atan(sqrt(pow(Sattelite_2.x, 2) + pow(Sattelite_2.y, 2)) / Sattelite_2.z);
        theta2 = SatteliteModel->RadToDeg(theta2);

        if (abs(theta2 - theta) < 1.9)
            i++;
    }
    //std::cout << i << "-------------" << std::endl;
}

double TimeIntersect(CSGP4_SDP4 SatteliteModel, LLAPos StationLLA_1) {
    
    SYSTEMTIME Time;
    SYSTEMTIME temp;
    SYSTEMTIME Time2;
    GetSystemTime(&Time); // Текущее время + 30 секунд
    GetSystemTime(&Time2); // Текущее время + 30 минут

    temp = Time;
    Time.wSecond = Time.wSecond + 30;
    Time2.wMinute = Time2.wMinute + 30;

    double time = SatteliteModel.JulianDate(Time); //нижняя граница
    double tempTime = SatteliteModel.JulianDate(temp); //Временной промежуток
    double time2 = SatteliteModel.JulianDate(Time2);//верхняя граница

    ////--------------------Пространственное определение СПУТНИКА-------------------------
    //SatteliteModel.SGP(tempTime);
    LLAPos SatteliteLLA_1;
    VECTOR Sattelite_1 = { 0, 0, 0, 0 };
    bool check1 = 0;
    //Sattelite_1 = SattelitePos(&SatteliteModel, &SatteliteLLA_1, Sattelite_1, time);
    ////--------------------Пространственное определение СПУТНИКА-------------------------

    ////--------------------Проверка 1: вхождение в область видимости станции-------------------------
    //bool check1 = IntersectCheck(StationLLA_1, Sattelite_1);
    //std::cout << "Intersection: " << check1 << std::endl;
    ////--------------------Проверка 1: вхождение в область видимости станции-------------------------


    while ((check1 == false || tempTime <= time) && tempTime <= time2 ) {
        
        tempTime = SatteliteModel.JulianDate(temp);
        //std::cout << temp.wDay << "d " << temp.wHour << "h " << temp.wMinute << "m " << temp.wSecond << "s " << std::endl;
        //--------------------Пространственное определение СПУТНИКА-------------------------
        SatteliteModel.SGP(tempTime);
        Sattelite_1 = SattelitePos(SatteliteModel, &SatteliteLLA_1, Sattelite_1, tempTime);


        //--------------------Пространственное определение СПУТНИКА-------------------------

        //--------------------Проверка 1: вхождение в область видимости станции-------------------------
        check1 = IntersectCheck(StationLLA_1, Sattelite_1);
        //std::cout << "Intersection: " << check1 << std::endl;
        //--------------------Проверка 1: вхождение в область видимости станции-------------------------
        temp.wSecond += 1;
    }

    if (tempTime > time2 || tempTime < time) {
        return 0;
    }
    return tempTime;

}


int main() {
    //setlocale(LC_ALL, "");
    //Чтение файла TLE состоящего из 3-ч строк
    const char* filename = "TLE.txt"; // 
    CSGP4_SDP4 SatteliteModel[N];
    std::vector <Sattelite> Objects;

    //--------------------Пространственное определение СТАНЦИИ-------------------------
    LLAPos StationLLA_1;
    VECTOR Station_1{ 0, 0, 0, 0 };
    if (!defStationPos(&StationLLA_1)) //Ручной ввод и проверка на аддекватность вводимых данных
        return 0;
    //Перевод градусов в радианы для станции
    StationLLA_1.Lat = SatteliteModel[0].DegToRad(StationLLA_1.Lat);
    StationLLA_1.Lon = SatteliteModel[0].DegToRad(StationLLA_1.Lon);
    //--------------------Пространственное определение СТАНЦИИ-------------------------
    Station_1 = ConvertLLAtoCoordinate(&StationLLA_1);
    //std::cout << Station_1.x << ' ' << Station_1.y << ' ' << Station_1.z << std::endl;

    VelocityCheck(filename, StationLLA_1, SatteliteModel);
   
    double time = 0;
    int test = 0;
    for (int i = 0; i < N; i++)
    {
        if (SatteliteModel[i].GetNORAD() != -858993460) {
            time = TimeIntersect(SatteliteModel[i], StationLLA_1);
            test = SatteliteModel[i].GetNORAD();
            //std::cout << test << std::endl;
            if (time != 0) {
                Objects.push_back({ SatteliteModel[i].GetNORAD(), SatteliteModel[i].CalendarDate(time), SatteliteModel[i].GetVel().w, time });
            }
        }
    }

    std::sort(Objects.begin(), Objects.end(), [](const Sattelite& a, const Sattelite& b) {
        return a.timeStamp < b.timeStamp; // Сортировка по значению
        });

    std::cout << "# Sattelite Number" << " \t| " << "velocity km/s\t | " << "Datetime\tDD/MM/YYYY\thh:mm:ss" << std::endl;


    for (const auto& item : Objects) {
        std::cout << "# " << item.SatteliteNumber << " \t\t| " << item.velocity << "    \t | \t\t" << item.time.wDay << "/" << item.time.wMonth << "/" << item.time.wYear << "\t" << item.time.wHour << ":" << item.time.wMinute << ":" << item.time.wSecond << " " << std::endl;
    }

    return 0;
}