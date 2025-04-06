cmake_minimum_required(VERSION 3.16)
project(FluidSimulator LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Find SFML package
find_package(SFML 2.6.2 COMPONENTS graphics window system audio network REQUIRED)

# Automatically find and add all source files in the src directory
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.hpp" "src/*.h" "src/*.c" "src/fluid/*")

# Add source files
add_executable(${PROJECT_NAME} ${SOURCES})

# Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()

# Link SFML libraries
target_link_libraries(${PROJECT_NAME} PRIVATE sfml-graphics)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)