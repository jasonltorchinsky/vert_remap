# Vertical remapping algorithm.
# Modifying Author: Jason Torchisnky
# Date: Summer 2021

# Define path variables used throuhgout the project
set(ENV{SOURCE_DIR} ${CMAKE_BINARY_DIR}/src)
set(ENV{OUTPUT_MOD_DIR} $ENV{SOURCE_DIR}/output_mod)
set(ENV{UTILS_MOD_DIR} $ENV{SOURCE_DIR}/utils_mod)
set(ENV{VERTREMAP_MOD_DIR} $ENV{SOURCE_DIR}/vertremap_mod)

add_subdirectory(output_mod)
add_subdirectory(utils_mod)
add_subdirectory(vertremap_mod)

add_executable(vert_remap
               main.F90
              )

target_link_libraries(vert_remap
                      PUBLIC output_mod
                      PUBLIC utils_mod
                      PUBLIC vertremap_mod
                     )

target_include_directories(vert_remap
                           PUBLIC $ENV{OUTPUT_MOD_DIR}
                           PUBLIC $ENV{UTILS_MOD_DIR}
                           PUBLIC $ENV{VERTREMAP_MOD_DIR}
                          )
                           
