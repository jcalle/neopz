########## next target ##############

file(GLOB headers *.h 
  ../NewModalAnalysis/SPZModalAnalysisData.h
  ../NewModalAnalysis/SPZModalAnalysisDataReader.h)
file(GLOB sources *.cpp ../NewModalAnalysis/SPZModalAnalysisDataReader.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(StepFiberAnalysisPzMesh ${main_SRCS})
target_link_libraries(StepFiberAnalysisPzMesh pz)
