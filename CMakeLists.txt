#cmake_minimum_required(VERSION 3.12)
#project(MonteCarloEngine)

#set(CMAKE_CXX_STANDARD 14)


#add_executable(MonteCarloEngine main.cpp Option.cpp )

cmake_minimum_required(VERSION 3.12)
project(MonteCarloEngine)

# 1. Define some well-known folders:
# - Place for compile binaries

set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR}/build)

# - Root for include files
set(INCLUDE_HOME ${PROJECT_SOURCE_DIR}/src)

# - In-house cmake functions
set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-lib)

# - set boost home
SET (BOOST_ROOT /usr/include/boost)


# Print some stuff
message("--> CMake extra functions: " ${CMAKE_MODULE_PATH})
message("--> Binaries: " ${CMAKE_BINARY_DIR})
message("--> Include base path: " ${INCLUDE_HOME})

# 2. Load in-house cmake modules
message("--> Loading CMake modules: ")
include(${CMAKE_MODULE_PATH}/CMakeBoostTest.cmake)
include(${CMAKE_MODULE_PATH}/CMakeAuxFunctions.cmake)

# 3. Set CMake environment
# - Set folders
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# - C++ file filter for adding files in cmake's projects
set(CXX_FILE_EXTENSIONS "*.h" "*.cpp")

# 4. Set Compilation environment
include_directories(${INCLUDE_HOME})
include_directories(/usr/local/lib)
include_directories(/usr/local/include/)

add_definitions(-std=gnu++14)

# 5. Add subdirs
add_subdirectory(src)
