# Defines utils_mod as a library
add_library(utils_mod
            utils_mod.F90
            netcdf_utils_mod.F90
            sys_utils_mod.F90
           )

# Find needed packages
#find_package(netCDF-Fortran REQUIRED)

#target_link_libraries(utils_mod
#                      PRIVATE netCDF::netcdff
#                     )
