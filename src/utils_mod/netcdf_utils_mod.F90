!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for netCDF subroutines.
! TO-DO: Add file meta data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (utils_mod) netcdf_utils_mod

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Check netCDF calls for errors.
  ! TO-DO: Add subroutine metadata.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  module subroutine netcdf_err_check_sbr(statusFlag)

    use iso_fortran_env, only: int32
    use netcdf

    implicit none

    integer(int32), intent(in) :: statusFlag

    if (statusFlag .ne. nf90_noerr) then
       print *, 'NetCDF error code: ', statusFlag
       print *, trim(nf90_strerror(statusFlag))
       error stop '  Stopped due to netCDF error.'
    end if

  end subroutine netcdf_err_check_sbr

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Initializes output file (pre-defined format), returns handle
  ! TO-DO: - Add subroutine meta data
  !        - var_dims only accepts one dimension per variable. This is okay
  !          for the 1-D problem we are starting with, but may become a problem
  !        - in the future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  module subroutine netcdf_init_outfile_sbr(ndims, dim_names, dims, nvars, &
       & var_names, var_types, var_dims, outDir, fname, outputFile)

    use iso_fortran_env, only: int32, real64
    use netcdf

    implicit none

    integer(int32), intent(in) :: ndims, nvars
    character(len=15), intent(in) :: dim_names(ndims)
    character(len=15), intent(in) :: var_names(nvars)
    integer(int32), intent(in) :: dims(ndims)
    integer(int32), intent(in) :: var_types(nvars)
    integer(int32), intent(in) :: var_dims(nvars)
    character(len=50), intent(in) :: outDir, fname
    type(outfile), intent(inout) :: outputFile
    integer(int32) :: ii

    ! Create the output file.
100 format (A, A, A)
    write(outputFile%fname, 100) trim(adjustl(outDir)), '/', &
         & trim(adjustl(fname))

    call netcdf_err_check_sbr( nf90_create(trim(adjustl(outputFile%fname)), &
         & nf90_clobber, outputFile%id) )

    ! Define dimensions for file:
    outputFile%ndims = ndims
    allocate(outputFile%dim_ids(ndims))
    do ii = 1, ndims
       call netcdf_err_check_sbr( nf90_def_dim(outputFile%id, &
            & trim(adjustl(dim_names(ii))), dims(ii), &
            & outputFile%dim_ids(ii)) )
    end do

    ! Define variables for the file
    outputFile%nVars = nvars
    allocate(outputFile%var_ids(nvars))
    do ii = 1, nvars
       call netcdf_err_check_sbr( nf90_def_var(outputFile%id, &
            & trim(adjustl(var_names(ii))), var_types(ii), &
            & outputFile%dim_ids(var_dims(ii)), outputFile%var_ids(ii)) )
    end do

    ! Exit define mode, i.e., tell netCDF we are done defining meta-data.
    call netcdf_err_check_sbr( nf90_enddef(outputFile%id) )

  end subroutine netcdf_init_outfile_sbr

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Close the output file.
  ! TO-DO: - Add subroutine meta-data.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  module subroutine netcdf_close_outfile_sbr(outputFile)

    use netcdf

    implicit none

    type(outfile), intent(inout) :: outputFile

    call netcdf_err_check_sbr( nf90_close(outputFile%id) )

    deallocate(outputFile%dim_ids)
    deallocate(outputFile%var_ids)
    outputFile%fname = ' '
    outputFile%id = 0
    outputFile%ndims = 0
    outputFile%nvars = 0

  end subroutine netcdf_close_outfile_sbr

end submodule netcdf_utils_mod
