cmake_minimum_required(VERSION 3.15)
project(Denigma)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)


add_executable(Denigma
        src/ArrayUtilities.cpp
        src/CypherUtilities.cpp
        src/decipher.cpp
        src/decipher.h
        src/MCMC.cpp
        src/Randomize.cpp
        src/FineSearch.cpp)
