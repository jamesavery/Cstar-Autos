cmake_minimum_required(VERSION 3.2.0)
project(shiftspace_autos)


set(CMAKE_CXX_STANDARD 11)

include_directories(AFTER SYSTEM src/)
add_subdirectory(src/algebra)
add_subdirectory(src/wordgraph)
add_subdirectory(src/automaton)
add_subdirectory(src/output)
add_subdirectory(src/shiftspace_auto)

add_library(shiftspace_autos
  $<TARGET_OBJECTS:algebra>
  $<TARGET_OBJECTS:wordgraph>
  $<TARGET_OBJECTS:automaton>
  $<TARGET_OBJECTS:output>
  $<TARGET_OBJECTS:shiftspace_auto>
)


add_subdirectory(programs)


