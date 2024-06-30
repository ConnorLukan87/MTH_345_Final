cmake_minimum_required(VERSION 3.25)
project(RungeKuttaMTH345)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)

add_executable(RungeKuttaMTH345
        particle_simulation.cpp)

target_link_libraries(RungeKuttaMTH345 PUBLIC OpenMP::OpenMP_CXX ${OpenCV_LIBS})
