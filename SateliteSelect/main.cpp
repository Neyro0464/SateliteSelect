#include<iostream>
#include"SGP4.h"

#include <fstream>
#include <cstring>

using namespace SGP4Funcs;

const int MAX_TLE_LINES = 100; // ������������ ���������� TLE �������
const int MAX_LINE_LENGTH = 130; // ������������ ����� ������ TLE

void readTLEFile(const std::string& filename, char longstr1[MAX_TLE_LINES * MAX_LINE_LENGTH], char longstr2[MAX_TLE_LINES * MAX_LINE_LENGTH], int& count) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "������ �������� �����: " << filename << std::endl;
        return;
    }

    count = 0;
    while (count < MAX_TLE_LINES) {
        // ������ ������ ������ TLE
        if (file.getline(&longstr1[count * MAX_LINE_LENGTH], MAX_LINE_LENGTH)) {
            // ������ ������ ������ TLE
            if (!file.getline(&longstr2[count * MAX_LINE_LENGTH], MAX_LINE_LENGTH)) {
                std::cerr << "������������ ����� � ����� ��� ��������� ������." << std::endl;
                break;
            }
            count++;
        }
        else {
            break;
        }
    }

    file.close();
}

int main() {
    // ������������� ������� ������
    setlocale(LC_ALL, "Russian");

    char longstr1[MAX_TLE_LINES * MAX_LINE_LENGTH];
    char longstr2[MAX_TLE_LINES * MAX_LINE_LENGTH];
    int count = 0;

    readTLEFile("TLE.txt", longstr1, longstr2, count);

    // ����� ����������� TLE ������
    std::cout << "����������� TLE ������:\n" << std::endl;
    for (int i = 0; i < count; ++i) {
        std::cout << "TLE " << i + 1 << ":\n";
        std::cout << "������ ������: " << &longstr1[i * MAX_LINE_LENGTH] << "\n";
        std::cout << "������ ������: " << &longstr2[i * MAX_LINE_LENGTH] << "\n" << std::endl;
    }
    gravconsttype prototype;
    twoline2rv(longstr1, longstr2, 'm', 'm', 'a', prototype, );
}