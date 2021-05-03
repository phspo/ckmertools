cmake_minimum_required(VERSION 3.14)
project(ckmertools)

set(CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS program_options log log_setup REQUIRED)
link_libraries(jsoncpp)
find_package(jsoncpp REQUIRED)
#get_target_property(JSON_INC_PATH jsoncpp_lib_static INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${JSON_INC_PATH})
include_directories(${Boost_INCLUDE_DIR})

add_executable(c_kmertools src/main.cpp src/probabilistic.cpp src/probabilistic.h src/logsystem.cpp src/logsystem.h src/parsing.cpp src/parsing.h src/KmersWrapper.cpp)

# removed static?
target_link_libraries (c_kmertools ${Boost_LIBRARIES} jsoncpp_lib cnpy z)