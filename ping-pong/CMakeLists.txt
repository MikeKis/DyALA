cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(ping-pong)

set(CMAKE_CXX_STANDARD 20)

include_directories(../../ArNI)
include_directories(.)
link_directories(/usr/local/lib) #it is - to find boost!
link_directories(~/E/ArNI/lib) #it is - to find boost!
link_directories(../../ArNI/lib)

add_subdirectory(environment)

if (DEFINED ENV{GRAPHIC_INTERFACE})
	add_subdirectory(ping-pong-visualizer)
endif()

