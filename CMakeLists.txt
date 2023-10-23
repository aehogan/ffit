cmake_minimum_required(VERSION 2.8)

project(ffit)

set(SRC
src/main.cpp
src/es.cpp
src/lj.cpp
src/pol.cpp
src/system.cpp
src/phahst.cpp
src/axilrod.cpp
)

include_directories(src)

add_compile_options(-g -O2)

add_executable(ffit ${SRC})

