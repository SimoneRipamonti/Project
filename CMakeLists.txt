cmake_minimum_required(VERSION 3.5)

project(MyProject VERSION 1.0 DESCRIPTION "my project" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED HINTS $ENV{mkBoostPrefix})

find_package(Boost 1.40.0 COMPONENTS filesystem system iostreams REQUIRED)

find_package(Eigen3 REQUIRED HINTS $ENV{mkEigenPrefix})

#find_package(muparser HINTS /home/simoripa96/Scrivania/pacs-examples/Examples/lib)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/externals/gnuplot-iostream" /home/simoripa96/Scrivania/pacs-examples/Examples/include SYSTEM ${EIGEN3_INCLUDE_DIRS})

add_executable(my_exec main.cpp src/matrix.cpp src/muparser_fun.cpp src/output.cpp src/parameters.cpp src/systems.cpp src/data.pot)

target_include_directories(my_exec PUBLIC src)

target_link_libraries(my_exec ${Boost_LIBRARIES} /home/simoripa96/Scrivania/pacs-examples/Examples/lib/libmuparser.so)

#target_link_libraries(my_exec /home/simoripa96/Scrivania/pacs-examples/Examples/lib/libmuparser.so)







