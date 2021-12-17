cmake_minimum_required(VERSION 3.13.0)

# create project
project(moe-mpc VERSION 0.1.0 LANGUAGES CXX)

# defines conventional GNU installation directories
include(GNUInstallDirs)

# set compiler flags
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -03") # all warnings
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3") # warning level 4
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # multicore build
endif()

# Use fetch content to get libraries that moe is dependent on
include(FetchContent)

# MAHI DAQ
FetchContent_Declare(moe GIT_REPOSITORY https://github.com/mahilab/moe.git
                         GIT_TAG        origin/main)
FetchContent_MakeAvailable(moe)

# MAHI MPC
FetchContent_Declare(mahi-mpc GIT_REPOSITORY https://github.com/mahilab/mahi-mpc.git
                              GIT_TAG        origin/parallelize-update)
FetchContent_MakeAvailable(mahi-mpc)

# add definitions
add_definitions(-D_CRT_SECURE_NO_WARNINGS -DNOMINMAX -D_WINSOCK_DEPRECATED_NO_WARNINGS)

add_executable(rom_demo_mpc src/rom_demo_mpc.cpp)
target_link_libraries(rom_demo_mpc moe mahi::mpc)

add_executable(rom_demo_mpc_one_dof src/rom_demo_mpc_one_dof.cpp)
target_link_libraries(rom_demo_mpc_one_dof moe mahi::mpc)

add_executable(moe_model_generate src/moe_model_generate.cpp)
target_link_libraries(moe_model_generate moe mahi::mpc)

add_executable(moe_model_generate_one_dof src/moe_model_generate_one_dof.cpp)
target_link_libraries(moe_model_generate_one_dof moe mahi::mpc)