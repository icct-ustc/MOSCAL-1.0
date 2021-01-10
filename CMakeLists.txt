cmake_minimum_required (VERSION 2.8)
project (DEOM)
set (DEOM_VERSION_MAJOR 0)
set (DEOM_VERSION_MINOR 0)

set (CMAKE_VERBOSE_MAKEFILE OFF)

file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/share)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

set (CMAKE_CXX_FLAGS "-std=c++11 -g -Wall -Ofast -march=native")

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/src/blockdeom 
                    ${PROJECT_SOURCE_DIR}/src/blockdeom2 
                    ${PROJECT_SOURCE_DIR}/src/deom
                    ${PROJECT_SOURCE_DIR}/src/ideom
                    ${PROJECT_SOURCE_DIR}/src/deom2
                    ${PROJECT_SOURCE_DIR}/thirdparty/trie
                   )

link_directories(${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/blockdeom)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/blockdeom2)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/deom)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ideom)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/deom2)

FILE(GLOB apps_src ${PROJECT_SOURCE_DIR}/src/apps/*.cpp)

foreach (p ${apps_src})
    get_filename_component(p_bin ${p} NAME_WE)
    add_executable(${p_bin} ${p})
    target_link_libraries(${p_bin} json11 openblas lapack armadillo deom blockdeom deom2 ideom blockdeom2)
endforeach ()

#enable_language(Fortran)
#include(ExternalProject)

#if (EXISTS ${PROJECT_SOURCE_DIR}/lib/libjson11.a)
#else()
#  ExternalProject_Add(json11
#    GIT_REPOSITORY https://github.com/dropbox/json11.git
#    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
#    INSTALL_DIR ${PROJECT_SOURCE_DIR}
#    INSTALL_COMMAND cp <SOURCE_DIR>/json11.hpp <INSTALL_DIR>/include && cp <BINARY_DIR>/libjson11.a <INSTALL_DIR>/lib
#  )
#endif()
##
#if (EXISTS ${PROJECT_SOURCE_DIR}/include/armadillo)
#else()
#  ExternalProject_Add(armadillo
#    URL http://sourceforge.net/projects/arma/files/armadillo-7.400.2.tar.xz
#    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
#    INSTALL_DIR ${PROJECT_SOURCE_DIR}
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#  )
#endif()
