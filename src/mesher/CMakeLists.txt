project(libmesher)

file(GLOB SRCS
  "*.h" "*.hh" "*.hpp"
  "*.c" "*.cc" "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} predicates)
add_dependencies     (${PROJECT_NAME} predicates)

target_link_libraries(${PROJECT_NAME} libmesh)
add_dependencies     (${PROJECT_NAME} libmesh)

target_include_directories(${PROJECT_NAME} PUBLIC .)
