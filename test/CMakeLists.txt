cmake_minimum_required(VERSION 3.5)

project(MyProject VERSION 1.0 DESCRIPTION "my project" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED HINTS $ENV{mkBoostPrefix})

find_package(Boost 1.40.0 COMPONENTS filesystem system iostreams REQUIRED)

find_package(Eigen3 REQUIRED HINTS $ENV{mkEigenPrefix})

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/externals/gnuplot-iostream" ${CMAKE_CURRENT_SOURCE_DIR}/../pacs-examples/Examples/include SYSTEM ${EIGEN3_INCLUDE_DIRS} include)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/externals/gnuplot-iostream" external/include SYSTEM ${EIGEN3_INCLUDE_DIRS} include)




#add_executable(my_exec main.cpp solo_trasporto_sparso/output.cpp solo_trasporto_sparso/matrix.cpp solo_trasporto_sparso/parameters.cpp solo_trasporto_sparso/muparser_fun.cpp solo_trasporto_sparso/systems.cpp)

add_executable(my_exec main.cpp solo_trasporto_sparso/output.cpp src/matrix.cpp solo_trasporto_sparso/parameters.cpp src/muparser_fun.cpp solo_trasporto_sparso/systems.cpp)


target_include_directories(my_exec PUBLIC solo_trasporto_sparso)


#target_link_libraries(my_exec ${Boost_LIBRARIES} ./external/lib/libmuparser.so)
#target_link_libraries(my_exec ${Boost_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/../pacs-examples/Examples/lib/libmuparser.so)
#target_link_libraries(my_exec ${Boost_LIBRARIES} /home/simoripa96/Scrivania/my_project/external/lib/libmuparser.so)
target_link_libraries(my_exec ${Boost_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/libmuparser.so)






