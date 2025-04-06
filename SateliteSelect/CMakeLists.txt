cmake_minimum_required(VERSION 3.10)
project(SGPSDP)

# Установка стандарта C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Сборка статической библиотеки из исходных файлов
add_library(SGPSDP STATIC
    Astro.cpp
    Math.cpp
    observer.cpp
    SGP4.cpp
    Sgpsdp.cpp
    time.cpp
    Utils.cpp
)

# Добавление директории с заголовочными файлами
target_include_directories(SGPSDP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Подавление предупреждений компилятора для GCC/Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(SGPSDP PRIVATE -w)
endif()

# Отключение специфичных для MSVC предупреждений (если нужно)
if(MSVC)
    target_compile_options(SGPSDP PRIVATE /wd4996)
endif()