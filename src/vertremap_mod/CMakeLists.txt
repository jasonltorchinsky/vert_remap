add_library(vertremap_mod
            vertremap_base.F90
	    vertremap_redux.F90
           )

target_link_libraries(vertremap_mod
                      PUBLIC lapack
		      PUBLIC blas
		     )
