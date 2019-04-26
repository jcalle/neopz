file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS main.cpp ${sources} ${headers})

add_executable(SBFemTests ${main_SRCS} ${header})

target_link_libraries(SBFemTests pz)