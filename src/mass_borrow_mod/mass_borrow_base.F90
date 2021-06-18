module mass_borrow_mod

  implicit none

  private

  public :: borrow_mass

contains

  subroutine borrow_mass(ncell, Qdp, dp2, ubnd, lbnd)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells
    real(real64), intent(inout) :: Qdp(ncell) ! Mass of each cell
    real(real64), intent(in) :: dp2(ncell) ! Width of each cell
    real(real64), intent(in) :: ubnd, lbnd ! Upper and lower bound for each
    ! cell density, that needs to be carried over.
    real(real64) :: avgdens(ncell) ! Average density of each cell.
    real(real64) :: densdiff ! Difference between average density of cell
    ! and bound
    integer(int32) :: ii ! Counter for DO loops

    avgdens = Qdp / dp2

    do ii = 1, (ncell-1)/2
       if (avgdens(ii) .gt. ubnd) then
          densdiff = avgdens(ii) - ubnd
          avgdens(ii) = ubnd
          avgdens(ii+1) = avgdens(ii+1) + densdiff
       end if
       if (avgdens(ncell-ii+1) .gt. ubnd) then
          densdiff = avgdens(ncell-ii+1) - ubnd
          avgdens(ncell-ii+1) = ubnd
          avgdens(ncell-ii) = avgdens(ncell-ii) + densdiff
       end if
    end do

    Qdp = avgdens * dp2

  end subroutine borrow_mass

end module mass_borrow_mod
