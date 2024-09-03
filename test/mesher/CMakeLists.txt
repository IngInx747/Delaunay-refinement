project(mesher)

file(GLOB SRCS "*.h" "*.hh" "*.hpp" "*.c" "*.cc" "*.cpp")

add_executable(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} libmesher)
add_dependencies     (${PROJECT_NAME} libmesher)