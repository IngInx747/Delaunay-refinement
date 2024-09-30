project(libgeometry)

add_library(${PROJECT_NAME} INTERFACE)

add_dependencies(${PROJECT_NAME} OpenMeshCore)

target_include_directories(${PROJECT_NAME} INTERFACE .)
