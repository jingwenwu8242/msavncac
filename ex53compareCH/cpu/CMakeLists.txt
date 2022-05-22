cmake_minimum_required(VERSION 3.0.0)
project(TwoPhase VERSION 0.1.0)

add_subdirectory(3componentCACcpu)
add_subdirectory(3componentCHcpu)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
