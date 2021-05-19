add_library(output_mod
            output_mod.F90
	   )

# Find needed packages
if(OnCori)
  target_link_libraries(output_mod
                        PUBLIC utils_mod
                       )

  target_include_directories(output_mod
                             PUBLIC $ENV{UTILS_MOD_DIR}
                            )

else()
  find_package(netCDF-Fortran REQUIRED)

  target_link_libraries(output_mod
                        PRIVATE netCDF::netcdff
                        PUBLIC utils_mod
                       )

  target_include_directories(output_mod
                             PUBLIC $ENV{UTILS_MOD_DIR}
                            )
endif()


