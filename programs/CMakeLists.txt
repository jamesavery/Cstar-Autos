set(PROGRAMS generate_autos)

foreach(program ${PROGRAMS})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} shiftspace_autos)
endforeach(program ${PROGRAMS})
