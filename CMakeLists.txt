cmake_minimum_required(VERSION 3.0.2)
project(adaptive_nlp)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## System dependencies are found with CMake's conventions
#############################
## Find CasADi
#############################
## -----------
## Option 1
## -----------
# set(CASADI_DIR "/usr/local")
# set(CMAKE_PREFIX_PATH ${CASADI_DIR} ${CMAKE_PREFIX_PATH})
# find_package(casadi 
#   REQUIRED 
#   HINTS ${CASADI_DIR} ${CASADI_DIR}/lib
# )
# set(CASADI_LIBRARY_DIRS ${CASADI_DIR}/lib)
# set(CASADI_INCLUDEDIR ${CASADI_DIR}/include)
# set(CASADI_LIBRARIES casadi)

## -----------
## Option 2
## -----------
## Test if pkg-config is able to find your Ipopt installation by using 'pkg-config --libs ipopt'
set(CASADI_DIR "/usr/local")
set(CMAKE_PREFIX_PATH ${CASADI_DIR} ${CMAKE_PREFIX_PATH})
set(PKG_CONFIG_PATH ${CASADI_DIR}/lib/pkgconfig ${PKG_CONFIG_PATH})
find_package(PkgConfig REQUIRED)
pkg_search_module(CASADI casadi REQUIRED)

# ## -----------

# # if(NOT CASADI_LIBRARIES)
# #   message(SEND_ERROR "Could not find casadi.pc - did you set PKG_CONFIG_PATH?")
# # endif()

# message("CASADI_LIBRARY_DIRS: ${CASADI_LIBRARY_DIRS}")
# message("CASADI_LIBRARIES: ${CASADI_LIBRARIES}")
# message("CASADI_INCLUDEDIR: ${CASADI_INCLUDEDIR}")

#############################
## Find IPOPT
#############################
## Test if pkg-config is able to find your Ipopt installation by using 'pkg-config --libs ipopt'
# set(PKG_CONFIG_PATH "/usr/local/include/coin-or/")
# set(PKG_CONFIG_PATH "/usr/local/lib/pkgconfig/")
find_package(PkgConfig)
pkg_search_module(IPOPT ipopt REQUIRED)
set(IPOPT_INC "/usr/local/include/coin-or/")

#############################
## Find Eigen
#############################
set(EIGEN_INC "/usr/include/eigen3")

###########
## Build ##
###########
link_directories(${CASADI_LIBRARY_DIRS})

add_subdirectory(src)

add_subdirectory(examples)