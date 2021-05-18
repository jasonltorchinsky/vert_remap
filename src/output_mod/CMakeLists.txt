add_library(output_mod
            output_mod.F90
	   )

target_link_libraries(output_mod
                      PUBLIC utils_mod
                     )

target_include_directories(output_mod
                           PUBLIC $ENV{UTILS_MOD_DIR}
                          )
