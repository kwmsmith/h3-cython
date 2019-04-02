cmake_minimum_required(VERSION 3.5.0)

project(h3cy)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)

include(CMakeDependentOption)
include(CheckIncludeFile)

include(CheckAlloca)
include(CheckVLA)

check_alloca(have_alloca)
check_vla(have_vla)

add_subdirectory(h3cy)
