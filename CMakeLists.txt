cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
#set(CMAKE_CXX_FLAGS "-g -Wall -fsanitize=thread")
set(CMAKE_CXX_FLAGS "-O3 -Wall")

option(openmp "Compile with OpenMP for shared memory parallelism" ON)

include_directories(include vendor/pcg vendor/docopt)

project(cmmct)

add_executable(cmmct src/main.cpp
                     src/cross_sections.cpp
                     src/tallies.cpp
                     src/utils.cpp
                     src/transporter.cpp
                     src/substeper.cpp
                     src/direct_sampler.cpp
                     src/independent_trials.cpp
                     src/fixed_source_all.cpp
                     src/settings.cpp
                     vendor/docopt/docopt.cpp
              )

if(openmp)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(cmmct PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()
