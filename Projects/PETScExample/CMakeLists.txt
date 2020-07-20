file(GLOB header *.h)

set(sources main.cpp tutorials.cpp)
add_executable(PETScExample ${sources} ${header})
target_link_libraries(PETScExample pz)