project (predicates)
file(GLOB SRCS "*.h" "*.c")
add_library(${PROJECT_NAME} ${SRCS})
target_include_directories(${PROJECT_NAME} PUBLIC .)