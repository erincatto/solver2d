# solver2d samples app

# glad for OpenGL API
set(GLAD_DIR ${CMAKE_SOURCE_DIR}/extern/glad)

add_library(
	glad STATIC
	${GLAD_DIR}/src/glad.c
	${GLAD_DIR}/include/glad/glad.h
	${GLAD_DIR}/include/KHR/khrplatform.h
)
target_include_directories(glad PUBLIC ${GLAD_DIR}/include)

# glfw for windowing and input
SET(GLFW_BUILD_DOCS OFF CACHE BOOL "GLFW Docs")
SET(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "GLFW Examples")
SET(GLFW_BUILD_TESTS OFF CACHE BOOL "GLFW Tests")
SET(GLFW_INSTALL OFF CACHE BOOL "GLFW Install")

FetchContent_Declare(
	glfw
	GIT_REPOSITORY https://github.com/glfw/glfw.git
	GIT_TAG master
	GIT_SHALLOW TRUE
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(glfw)

# imgui and glfw backend for GUI
# https://gist.github.com/jeffamstutz/992723dfabac4e3ffff265eb71a24cd9
FetchContent_Populate(imgui
  URL https://github.com/ocornut/imgui/archive/docking.zip
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/build/imgui
)

set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/build/imgui)

add_library(imgui STATIC
	${IMGUI_DIR}/imgui.cpp
	${IMGUI_DIR}/imgui_draw.cpp
	${IMGUI_DIR}/imgui_demo.cpp
	${IMGUI_DIR}/imgui_tables.cpp
	${IMGUI_DIR}/imgui_widgets.cpp

	${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
	${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)
target_link_libraries(imgui PUBLIC glfw glad)
target_include_directories(imgui PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)

# jsmn for json
set(JSMN_DIR ${CMAKE_SOURCE_DIR}/extern/jsmn)

add_library(jsmn INTERFACE ${JSMN_DIR}/jsmn.h)
target_include_directories(jsmn INTERFACE ${JSMN_DIR})

set(S2D_SAMPLES
	draw.cpp
	draw.h
	main.cpp
	sample.cpp
	sample.h
	settings.h
	settings.cpp

	collection/sample_contact.cpp
)

add_executable(samples ${S2D_SAMPLES})

set_target_properties(samples PROPERTIES
	CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_include_directories(samples PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${JSMN_DIR})
target_link_libraries(samples PUBLIC solver2d imgui glfw glad)

# message(STATUS "runtime = ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
# message(STATUS "binary = ${CMAKE_CURRENT_BINARY_DIR}")

# Copy font files, etc
add_custom_command(
	TARGET samples POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
			${CMAKE_CURRENT_SOURCE_DIR}/data/
			${CMAKE_CURRENT_BINARY_DIR}/data/)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${S2D_SAMPLES})
