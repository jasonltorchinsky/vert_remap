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
    real(real64) :: massdiff ! Mass difference associated with the density
    ! difference.
    real(real64) :: dens ! Average density of the cell we are interested in
    integer(int32) :: ii ! Counter for DO loops
    integer(int32) :: l_cnt, r_cnt
    real(real64), parameter :: tol = 3.5e-14
    real(real64) :: max11, min11

    avgdens = Qdp / dp2
    l_cnt = 0
    r_cnt = 0

    max11 = maxval(avgdens)
    min11 = minval(avgdens)

    if ((max11 - ubnd .lt. tol) .and. (lbnd - min11 .lt. tol)) then
       ! q_alg = 11 algorithm is good
    else
       do ii = 1, (ncell-1)/2
          ! Ensure left->right
          dens = Qdp(ii) / dp2(ii)
          if (dens - ubnd .gt. tol) then
             ! Too big, give mass to the right
             densdiff = dens - ubnd
             massdiff = densdiff * dp2(ii)
             Qdp(ii) = Qdp(ii) - massdiff
             Qdp(ii+1) = Qdp(ii+1) + massdiff
             l_cnt = l_cnt + 1
          else if (lbnd - dens .gt. tol) then
             ! Too small, take mass from the right
             densdiff = lbnd - dens
             massdiff = densdiff * dp2(ii)
             Qdp(ii) = Qdp(ii) + massdiff
             Qdp(ii+1) = Qdp(ii+1) - massdiff
             l_cnt = l_cnt + 1
          end if
          
          ! Ensure right->left
          dens = Qdp(ncell-ii+1)/dp2(ncell-ii+1)
          if (dens - ubnd .gt. tol) then
             ! Too big, give mass to the left
             densdiff = dens - ubnd
             massdiff = densdiff * dp2(ncell-ii+1)
             Qdp(ncell-ii+1) = Qdp(ncell-ii+1) - massdiff
             Qdp(ncell-ii) = Qdp(ncell-ii) + massdiff
             r_cnt = r_cnt + 1
          else if (lbnd - dens .gt. tol) then
             ! Too small, take mass from left
             densdiff = lbnd - dens
             massdiff = densdiff * dp2(ncell-ii+1)
             Qdp(ncell-ii+1) = Qdp(ncell-ii+1) + massdiff
             Qdp(ncell-ii) = Qdp(ncell-ii) - massdiff
             r_cnt = r_cnt + 1
          end if
       end do

       print *, ' ~~ Cells adjusted:'
       print *, '    l_cnt: ', l_cnt
       print *, '    r_cnt: ', r_cnt
       
    end if

  end subroutine borrow_mass

end module mass_borrow_mod
