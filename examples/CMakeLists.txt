include_directories(${CASADI_INCLUDEDIR})
include_directories(${EIGEN_INC})
include_directories(${IPOPT_INC})
include_directories(${CMAKE_SOURCE_DIR})

set(COMMON_EXAMPLE_SOURCES)

LIST(APPEND COMMON_EXAMPLE_SOURCES
    makeBuildingBlocks.cpp
    makeBuildingBlocks.hpp
    plotter.cpp 
    plotter.hpp
)

#####################
## MINIMAL EXAMPLE ##
#####################

add_executable(minimalExample
    minimalExample/minimalExample.cpp
    ${COMMON_EXAMPLE_SOURCES}
)

target_include_directories(minimalExample PRIVATE 
    ${CMAKE_SOURCE_DIR}/examples
)

target_link_libraries(minimalExample
    ${PROJECT_NAME}
    ${CASADI_LIBRARIES}
    ipopt
)

#############################################
## ADAPTIVE CORRIDOR EXAMPLE (MPC-EXAMPLE) ##
#############################################

add_executable(adaptiveCorridorExample
    adaptiveCorridorExample/adaptiveCorridorExample.cpp
    adaptiveCorridorExample/adaptiveCorridorExampleHelpers.cpp
    adaptiveCorridorExample/adaptiveCorridorExampleHelpers.hpp
    ${COMMON_EXAMPLE_SOURCES}
)

target_include_directories(adaptiveCorridorExample PRIVATE 
    ${CMAKE_SOURCE_DIR}/examples
)

target_link_libraries(adaptiveCorridorExample
    ${PROJECT_NAME}
    ${CASADI_LIBRARIES}
    ipopt
)

########################
## MOONLANDER EXAMPLE ##
########################

add_executable(moonlanderExample
    moonlanderExample/moonlanderExample.cpp
    moonlanderExample/moonlanderExampleHelper.cpp
    moonlanderExample/moonlanderExampleHelper.hpp
    moonlanderExample/errorEstimator.cpp
    moonlanderExample/errorEstimator.hpp
    ${COMMON_EXAMPLE_SOURCES}
)

target_include_directories(moonlanderExample PRIVATE 
    ${CMAKE_SOURCE_DIR}/examples
)

target_link_libraries(moonlanderExample
    ${PROJECT_NAME}
    ${CASADI_LIBRARIES}
    ipopt
)

############################
## MANY OBSTACLES EXAMPLE ##
############################

add_executable(manyObstaclesExample
    manyObstaclesExample/manyObstaclesExample.cpp
    manyObstaclesExample/manyObstaclesHelper.cpp
    manyObstaclesExample/manyObstaclesHelper.hpp
    ${COMMON_EXAMPLE_SOURCES}
)

target_include_directories(manyObstaclesExample PRIVATE 
    ${CMAKE_SOURCE_DIR}/examples
)

target_link_libraries(manyObstaclesExample
    ${PROJECT_NAME}
    ${CASADI_LIBRARIES}
    ipopt
)