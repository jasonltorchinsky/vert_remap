# Vertical remapping algorithm.
# Modifying Author: Jason Torchisnky
# Date: Summer 2021

# Add options for compiling on and off Cori-Knl
option(OnCori "OnCori" OFF)

# Define path variables used throuhgout the project
set(ENV{SOURCE_DIR} ${CMAKE_BINARY_DIR}/src)
set(ENV{CONV_COMB_MOD_DIR} $ENV{SOURCE_DIR}/conv_comb_mod)
set(ENV{GET_INSECS_MOD_DIR} $ENV{SOURCE_DIR}/get_insecs_mod)
set(ENV{MASS_BORROW_MOD_DIR} $ENV{SOURCE_DIR}/mass_borrow_mod)
set(ENV{OUTPUT_MOD_DIR} $ENV{SOURCE_DIR}/output_mod)
set(ENV{UTILS_MOD_DIR} $ENV{SOURCE_DIR}/utils_mod)
set(ENV{VERTREMAP_MOD_DIR} $ENV{SOURCE_DIR}/vertremap_mod)

add_subdirectory(conv_comb_mod)
add_subdirectory(get_insecs_mod)
add_subdirectory(mass_borrow_mod)
add_subdirectory(output_mod)
add_subdirectory(utils_mod)
add_subdirectory(vertremap_mod)

add_executable(vert_remap
               main.F90
              )

if(OnCori)
  target_link_libraries(vert_remap
                        PUBLIC conv_comb_mod
                        PUBLIC get_insecs_mod
                        PUBLIC mass_borrow_mod
                        PUBLIC output_mod
                        PUBLIC utils_mod
                        PUBLIC vertremap_mod
                       )

  target_include_directories(vert_remap
                             PUBLIC $ENV{CONV_COMB_MOD_DIR}
                             PUBLIC $ENV{GET_INSECS_MOD_DIR}
                             PUBLIC $ENV{MASS_BORROW_MOD_DIR}
                             PUBLIC $ENV{OUTPUT_MOD_DIR}
                             PUBLIC $ENV{UTILS_MOD_DIR}
                             PUBLIC $ENV{VERTREMAP_MOD_DIR}
                            )
else()
  find_package(netCDF-Fortran REQUIRED)
  target_link_libraries(vert_remap
                        PRIVATE netCDF::netcdff
                        PUBLIC conv_comb_mod
                        PUBLIC get_insecs_mod
			PUBLIC mass_borrow_mod
                        PUBLIC output_mod
                        PUBLIC utils_mod
                        PUBLIC vertremap_mod
                       )

  target_include_directories(vert_remap
                             PUBLIC $ENV{CONV_COMB_MOD_DIR}
                             PUBLIC $ENV{GET_INSECS_MOD_DIR}
                             PUBLIC $ENV{MASS_BORROW_MOD_DIR}
                             PUBLIC $ENV{OUTPUT_MOD_DIR}
                             PUBLIC $ENV{UTILS_MOD_DIR}
                             PUBLIC $ENV{VERTREMAP_MOD_DIR}
                            )
endif()
                           
