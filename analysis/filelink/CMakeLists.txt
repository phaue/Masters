cmake_minimum_required(VERSION 3.18)
project(filelink)

# Find ROOT
find_package(ROOT REQUIRED)


# Find all headers
file(GLOB HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/cmake/*.h")
# Set the include directories for ROOT
include_directories(${ROOT_INCLUDE_DIRS})

add_executable(linkfiles filelink.cpp)

target_link_libraries(linkfiles PRIVATE ROOT::Core ROOT::TreePlayer)

target_include_directories(linkfiles PRIVATE ${CMAKE_SOURCE_DIR}/../../cmake)
