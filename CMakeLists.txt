cmake_minimum_required(VERSION 3.9)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

project(1DEuler)

include_directories(
  /home/nguyenly/local/1DEuler/1DEuler
  /home/nguyenly/local/eigen/3.2.9/eigen
  /usr/local/include
)

link_directories(

)

add_executable(
1DEuler 
main.cpp
GridOps/grid1d.cpp GridOps/grid1d.h
EOS/genericeos.h EOS/genericeos.cpp
EOS/idealgaseos.cpp EOS/idealgaseos.h
Solvers/genericsolver.cpp Solvers/genericsolver.h
#toml++/toml.hpp
Solvers/Interface.h
EOS/pengrobinsoneos.cpp EOS/pengrobinsoneos.h
)

target_link_libraries(
1DEuler
pugixml
)

target_compile_definitions(
1DEuler PUBLIC
NASAPOLYPATH="${CMAKE_SOURCE_DIR}/NASAPOLY/nasa_poly.xml"
)
