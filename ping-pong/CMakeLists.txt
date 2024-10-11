cmake_minimum_required(VERSION 3.14)

project(ping-pong)

set(CMAKE_CXX_STANDARD 14)

include_directories(../../ArNI)
include_directories(.)
link_directories(/usr/local/lib) #it is - to find boost!

add_subdirectory(environment)

if (DEFINED ENV{GRAPHIC_INTERFACE})
	add_subdirectory(ping-pong-visualizer)
endif()

