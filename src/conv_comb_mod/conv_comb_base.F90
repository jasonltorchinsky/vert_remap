module conv_comb_mod

  implicit none

  private

  public :: conv_comb

contains

  subroutine conv_comb(ncell, Qdp10, Qdp11, dp2, ubnd, lbnd)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells
    real(real64), intent(in) :: Qdp10(ncell) ! Mass of each cell for q_alg = 10
    real(real64), intent(inout) :: Qdp11(ncell) ! Mass of each
    ! cell for q_alg = 11, will be re-written.
    real(real64), intent(in) :: dp2(ncell) ! Width of each cell
    real(real64), intent(in) :: ubnd, lbnd ! Upper and lower bound for each
    ! cell density, that needs to be carried over.
    real(real64) :: avgdens10(ncell), avgdens11(ncell) ! Average density of
    ! each cell for q_alg = 10, q_alg = 11.
    real(real64) :: ext10, max11, min11 ! Extremum for q_alg = 10,
    ! maximum and minimum average densities for Q_alg = 11
    real(real64) :: theta ! Convex scalar for combo
    real(real64), parameter :: tol = 3.5e-14
    integer(int32) :: ii ! Counter for DO loops

    avgdens10 = Qdp10 / dp2
    avgdens11 = Qdp11 / dp2

    max11 = maxval(avgdens11)
    min11 = minval(avgdens11)

    if ((max11 - ubnd .lt. tol) .and. (lbnd - min11 .lt. tol)) then
       ! q_alg = 11 algorithm is good
       
    else ! q_alg 11 algorithm is not good
       if (max11 - ubnd .gt. lbnd - min11) then
          ! q_alg violates the upper bound more
          ext10 = maxval(avgdens10)
          if (abs(max11 - ext10) .lt. tol) then
             ! The extrema are equal, go with limited solution
             theta = 0.0_real64
          else
             theta = (ubnd - ext10) / (max11 - ext10)
          end if
       else
          ! q_alg violates lower bound more
          ext10 = minval(avgdens10)
          if (abs(ext10 - min11) .lt. tol) then
             ! The extrema are equal, go with limited solution
             theta = 0.0_real64
          else
             theta = (lbnd - ext10)/(min11 - ext10)
          end if
       end if

       if (theta .ne. 1.0_real64) then
          print *, ' ~~ Convex combination performed:'
          print *, '    theta: ', theta
       end if

       Qdp11 = theta * Qdp11 + (1.0_real64 - theta) * Qdp10

    end if

  end subroutine conv_comb

end module conv_comb_mod
