cmake_minimum_required(VERSION 3.23)
include(FetchContent)

message(STATUS "cmake source dir: ${CMAKE_SOURCE_DIR}")
message(STATUS "library postfix: ${CMAKE_DEBUG_POSTFIX}")

project(solver2d LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED YES)
set(CMAKE_C_EXTENSIONS YES)
set(CMAKE_CXX_STANDARD 17)

if(MSVC)
	# enable hot reloading
	add_compile_options("$<$<CONFIG:Debug>:/ZI>")
	add_link_options("$<$<CONFIG:Debug>:/INCREMENTAL>")
endif()

add_subdirectory(src)
add_subdirectory(samples)

# default startup project for Visual Studio
if (MSVC)
	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT samples)
	set_property(TARGET samples PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()	


# https://stackoverflow.com/questions/72635402/how-to-set-working-directory-in-cmake-visual-studio-2022
# launch.vs.json
# "currentDir": "${workspaceRoot}"