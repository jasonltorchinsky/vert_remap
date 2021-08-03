!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! General utilities module.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module utils_mod

  use iso_fortran_env, only: int32, real64
  use netcdf

  implicit none

  private

  ! Struct to hold information about the output file
  public :: outfile
  type :: outfile
     character(len=50) :: fname ! File name
     integer(int32) :: id ! NetCDF ID 
     integer(int32) :: ndims ! Number of dimensions
     integer(int32), allocatable :: dim_ids(:) ! Dimension netCDF IDs
     integer(int32) :: nvars ! Number of variables
     integer(int32), allocatable :: var_ids(:) ! Variable netCDF Ids
  end type outfile

  ! Visible subroutines
  !! netcdf_utils_mod
  public :: netcdf_err_check
  interface netcdf_err_check
     module subroutine netcdf_err_check_sbr(statusFlag)
       use iso_fortran_env, only: int32
       use netcdf
       implicit none
       integer(int32), intent(in) :: statusFlag
     end subroutine netcdf_err_check_sbr
  end interface netcdf_err_check

  public :: netcdf_init_outfile
  interface netcdf_init_outfile
     module subroutine netcdf_init_outfile_sbr(ndims, dim_names, &
          & dims, nvars, var_names, var_types, var_dims, &
          & outDir, fname, outputFile)
       use iso_fortran_env, only: int32
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
     end subroutine netcdf_init_outfile_sbr
  end interface netcdf_init_outfile

  public :: netcdf_close_outfile
  interface netcdf_close_outfile
     module subroutine netcdf_close_outfile_sbr(outputFile)
       use netcdf
       implicit none
       type(outfile), intent(inout) :: outputFile
     end subroutine netcdf_close_outfile_sbr
  end interface netcdf_close_outfile

  !! remap_check_utils_mod
  public :: check_global_bounded
  interface check_global_bounded
     module subroutine check_global_bounded_sbr(ncell, Q, ubnd, lbnd, flag)
       use iso_fortran_env, only: int32, real64
       implicit none
       integer(int32), intent(in) :: ncell
       real(real64), intent(in) :: Q(ncell)
       real(real64), intent(in) :: ubnd, lbnd
       integer(int32), intent(inout) :: flag
     end subroutine check_global_bounded_sbr
  end interface check_global_bounded

  !! sys_utils_mod
  public :: mkdir
  interface mkdir
     module subroutine mkdir_sbr(dirName)
       implicit none
       character(len=50), intent(in) :: dirName
     end subroutine mkdir_sbr
  end interface mkdir

end module utils_mod
