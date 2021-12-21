cmake_minimum_required(VERSION 3.5)

project(MyProject VERSION 1.0 DESCRIPTION "my project" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED HINTS $ENV{mkBoostPrefix})

find_package(Boost 1.40.0 COMPONENTS filesystem system iostreams REQUIRED)

find_package(Eigen3 REQUIRED HINTS $ENV{mkEigenPrefix})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../externals/gnuplot-iostream" ../external/include SYSTEM ${EIGEN3_INCLUDE_DIRS} ../include)


add_executable(main src/main.cpp ../src/output.cpp ../src/matrix.cpp ../src/parameters.cpp ../src/muparser_fun.cpp ../src/functions.cpp ../src/transport_decay.cpp)


target_include_directories(main PUBLIC src)


target_link_libraries(main ${Boost_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/../external/lib/libmuparser.so)




