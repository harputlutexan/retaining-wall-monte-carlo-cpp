cmake_minimum_required(VERSION 3.16)
project(retaining_wall_monte_carlo LANGUAGES CXX)

add_executable(retaining_wall main.cpp)
target_compile_features(retaining_wall PRIVATE cxx_std_17)

if (MSVC)
    target_compile_options(retaining_wall PRIVATE /W4 /utf-8 /O2)
else()
    target_compile_options(retaining_wall PRIVATE -Wall -Wextra -Wpedantic -O2)
endif()

find_package(Threads REQUIRED)
target_link_libraries(retaining_wall PRIVATE Threads::Threads)
