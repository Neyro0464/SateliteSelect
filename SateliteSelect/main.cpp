#include<iostream>
#include <fstream>
#include <cstring>
#include "Sgpsdp.h"
#include <ctime>
#include<chrono>

#pragma warning(disable : 4996)

void readFileAndCallFunction(const char* filename, CSGP4_SDP4 *test) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "������ ��� �������� �����!" << std::endl;
        return;
    }

    char m_cLine0[256];
    char m_cLine1[256];
    char m_cLine2[256];
    
    // ������ ����� �� �����
    file.getline(m_cLine0, sizeof(m_cLine0));
    file.getline(m_cLine1, sizeof(m_cLine1));
    file.getline(m_cLine2, sizeof(m_cLine2));
    // �������� �����
    file.close();

    // �������� ������� ������
    *test = CSGP4_SDP4(m_cLine0, m_cLine1, m_cLine2);
}

int main() {
    // ������ ����� TLE ���������� �� 3-� �����
    const char* filename = "TLE.txt"; //
    CSGP4_SDP4 test;
    readFileAndCallFunction(filename, &test);
    
    // ����������� �������� ������� ��� ������� ������ 
    SYSTEMTIME t2;
    GetSystemTime(&t2);
    std::cout << t2.wMonth << ' ' << t2.wDay << ' ' << t2.wHour << ' ' << t2.wMinute << ' ' << t2.wSecond << std::endl;
    double time = test.JulianDate(t2);

    //������ ������
    test.SGP(time);

    //��������� ������� ��������� ������
    VECTOR possition = test.GetVel();
    std::cout << possition.x << ' ' << possition.y << ' ' << possition.z << ' ' << possition.w << std::endl;

    // ������ � ��������� ������, ������� � ������
    test.CalculateLatLonAlt(time);
    std::cout << test.GetLat() << ' ' << test.GetLon() << ' ' << test.GetAlt() << std::endl;

    return 0;
}