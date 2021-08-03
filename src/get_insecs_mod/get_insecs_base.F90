module get_insecs_mod

  implicit none

  private

  public :: get_insecs

contains

  subroutine get_insecs(ncell, grid1, grid2, cellinsecs)

    use iso_fortran_env, only : int32, real64

    implicit none

    integer(int32), intent(in) :: ncell
    real(real64), intent(in) :: grid1(0:ncell), grid2(0:ncell)
    integer(int32), intent(out) :: cellinsecs(2,ncell)
    integer(int32) :: jj, kk, kidx
    
    ! We're going to get a list of pairs to see in what original cells the new cell centers are in.
    cellinsecs(1, 1) = 1_int32 ! Left bound of first new cell is in first old cell
    kidx = 1_int32 ! Cell of old grid were are checking the right boundary of
    do jj = 1, ncell-1 ! Loop through all cells of new grid
       do kk = kidx, ncell ! Loop through all cell boundaries of new grid
          if (grid2(jj) .le. grid1(kk)) then ! In the kkth cell
             cellinsecs(2, jj) = kk
             cellinsecs(1, jj+1) = kk ! Left bound of next cell is same as right boundary of current cell
             kidx = kk
             exit
          end if
       end do
    end do
    cellinsecs(2, ncell) = ncell ! Right bound of last new cell in last old cell.

  end subroutine get_insecs

end module get_insecs_mod
