!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for output subroutines.
! TO-DO: Add file meta data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module output_mod

  use netcdf
  use utils_mod

  private

  public :: netcdf_add_metadata
  public :: netcdf_write_output

contains


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Adds metadata to the output file.
  ! TO-DO: - Add subroutine meta data
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine netcdf_add_metadata(outputFile)

    use netcdf
    use utils_mod

    implicit none

    type(outfile), intent(in) :: outputFile

    ! Go back into define mode.
    call netcdf_err_check( nf90_redef(outputFile%id) )

    ! Add metadata for each variable.
    !! grid1
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(1), 'description', &
         & 'original grid levels') )

    !! dp1
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(2), 'description', &
         & 'original grid spacing') )

    !! grid1_stg
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(3), 'description', &
         & 'original grid cell centers') )

    !! Q1
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(4), 'description', &
         & 'densities for original grid centers') )
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(4), 'units', &
         & 'kg m^(-1)') )

    !! Qdp1
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(5), 'description', &
         & 'masses for original grid centers') )
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(5), 'units', &
         & 'kg') )

    !! grid2
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(6), 'description', &
         & 'new grid levels') )

    !! dp2
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(7), 'description', &
         & 'new grid spacing') )

    !! grid2_stg
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(8), 'description', &
         & 'new grid cell centers') )

    !! Q2
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(9), 'description', &
         & 'densities for new grid centers') )
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(9), 'units', &
         & 'kg m^(-1)') )

    !! Qdp2
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(10), 'description', &
         & 'masses for original grid centers') )
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(10), 'units', &
         & 'kg') )

    !! QTrue
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(11), 'description', &
         & 'true densities for new grid centers') )
    call netcdf_err_check( nf90_put_att(outputFile%id, &
         & outputFile%var_ids(11), 'units', &
         & 'kg m^(-1)') )

    ! Exit define mode, i.e., tell netCDF we are done defining meta-data.
    call netcdf_err_check( nf90_enddef(outputFile%id) )

  end subroutine netcdf_add_metadata

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Write output to file.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine netcdf_write_output(nlev, ncell, grid1, dp1, grid1_stg, Q1, &
       & Qdp1, grid2, dp2, grid2_stg, Q2, Qdp2, QTrue, outputFile)

    use iso_fortran_env, only: int32, real64
    use netcdf
    use utils_mod

    implicit none

    integer(int32), intent(in) :: nlev, ncell
    real(real64), intent(in) :: grid1(nlev), grid2(nlev)
    real(real64), intent(in) :: dp1(ncell), grid1_stg(ncell), Q1(ncell), &
         & Qdp1(ncell), dp2(ncell), grid2_stg(ncell), Q2(ncell), &
         & Qdp2(ncell), QTrue(ncell)
    type(outfile), intent(in) :: outputFile

    ! Write each variable to file.
    !! grid1
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(1), grid1) )

    !! dp1
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(2), dp1) )

    !! grid1_stg
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(3), grid1_stg) )

    !! Q1
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(4), Q1) )

    !! Qdp1
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(5), Qdp1) )

    !! grid2
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(6), grid2) )

    !! dp2
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(7), dp2) )

    !! grid2_stg
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(8), grid2_stg) )

    !! Q2
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(9), Q2) )

    !! Qdp2
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(10), Qdp2) )

    !! QTrue
    call netcdf_err_check( nf90_put_var(outputFile%id, &
         & outputFile%var_ids(11), QTrue) )


  end subroutine netcdf_write_output

end module output_mod
