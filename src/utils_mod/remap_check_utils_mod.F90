!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for checking the remapped solution
! TO-DO: Add file meta data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (utils_mod) remap_check_utils_mod

contains

  ! Checks the global bounds for a remapped solution.
  ! 0 - Not globally bounded
  ! 1 - Globally bounded
  module subroutine check_global_bounded_sbr(ncell, Q, ubnd, lbnd, flag)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell
    real(real64), intent(in) :: Q(ncell)
    real(real64), intent(in) :: ubnd, lbnd
    integer(int32), intent(inout) :: flag
    real(real64) :: maxQ, minQ
    real(real64), parameter :: tol = 3.5e-14

    flag = 1_int32
    maxQ = maxval(Q)
    minQ = minval(Q)

    if ((maxQ - ubnd .gt. tol) .or. (lbnd - minQ .gt. tol)) then
       flag = 0_int32
    end if
       
    
  end subroutine check_global_bounded_sbr

end submodule remap_check_utils_mod
