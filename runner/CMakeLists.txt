# Headless runner. Reuses the sample scenes without a window and dumps CSV per solver.

set(SAMPLES_DIR ${CMAKE_SOURCE_DIR}/samples)

set(S2D_RUNNER
	draw_null.cpp
	main.cpp
	probe.c
	probe.h

	${SAMPLES_DIR}/sample.cpp
	${SAMPLES_DIR}/collection/human.cpp
	${SAMPLES_DIR}/collection/sample_contact.cpp
	${SAMPLES_DIR}/collection/sample_joints.cpp
	${SAMPLES_DIR}/collection/sample_far.cpp
)

add_executable(runner ${S2D_RUNNER})

set_target_properties(runner PROPERTIES
	C_STANDARD 17
	C_STANDARD_REQUIRED YES
	CXX_STANDARD 17
	CXX_STANDARD_REQUIRED YES
	CXX_EXTENSIONS NO
)

# The probe reads solver internals
target_include_directories(runner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${SAMPLES_DIR} ${CMAKE_SOURCE_DIR}/src)

# Only for GLFW headers pulled in by the sample code. No window is created.
target_link_libraries(runner PRIVATE solver2d glfw)

if (MSVC)
	set_property(TARGET runner PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()
