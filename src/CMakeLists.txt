## Create AdaptiveNLP library
set(${PROJECT_NAME}_SOURCES)

LIST(APPEND ${PROJECT_NAME}_SOURCES 
  ${CMAKE_SOURCE_DIR}/src/adaptiveNLP.cpp
  ${CMAKE_SOURCE_DIR}/src/bookkeeper.cpp
  ${CMAKE_SOURCE_DIR}/src/buildingBlocks.cpp
  ${CMAKE_SOURCE_DIR}/src/interfaceTester.cpp
  ${CMAKE_SOURCE_DIR}/src/NLPInterface.cpp
)

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CASADI_INCLUDEDIR} ${EIGEN_INC} ${IPOPT_INC})