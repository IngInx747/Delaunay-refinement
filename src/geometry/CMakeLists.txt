project(libgeometry)

file(GLOB SRCS
  "*.h" "*.hh" "*.hpp"
  "*.c" "*.cc" "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} OpenMeshCore)
add_dependencies     (${PROJECT_NAME} OpenMeshCore)

target_include_directories(${PROJECT_NAME} PUBLIC .)
