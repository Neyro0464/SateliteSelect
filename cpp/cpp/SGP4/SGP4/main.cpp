#include<iostream>
#include "SGP4.h"
#include <iostream>
#include <fstream>
#include <cstring>

const int MAX_TLE_LINES = 100; // Максимальное количество TLE записей
const int MAX_LINE_LENGTH = 100;


using namespace SGP4Funcs;

void readTLEFile(const std::string& filename, char tle1[MAX_TLE_LINES][MAX_LINE_LENGTH], char tle2[MAX_TLE_LINES][MAX_LINE_LENGTH], int& count) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
        return;
    }

    count = 0;
    while (count < MAX_TLE_LINES && file.getline(tle1[count], MAX_LINE_LENGTH)) {
        if (file.getline(tle2[count], MAX_LINE_LENGTH)) {
            count++;
        }
        else {
            std::cerr << "Недостаточно строк в файле для последней записи." << std::endl;
            break;
        }
    }

    file.close();
}

int main() {
    char tle1[MAX_TLE_LINES][MAX_LINE_LENGTH];
    char tle2[MAX_TLE_LINES][MAX_LINE_LENGTH];
    int count = 0;

    readTLEFile("TLE.txt", tle1, tle2, count);

    // Вывод прочитанных TLE данных
    for (int i = 0; i < count; ++i) {
        std::cout << "TLE " << i + 1 << ":\n";
        std::cout << tle1[i] << "\n" << tle2[i] << "\n" << std::endl;
    }

    return 0;
}


