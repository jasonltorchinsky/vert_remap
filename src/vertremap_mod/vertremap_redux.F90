!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! New vertical remapping code, should take in generally the same arguments
! as the original, just simplified.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (vertremap_mod) vertremap_redux

contains

  module subroutine new_remap_sbr(qdp, ncell, dp1, dp2, remap_alg, verbosity)

    use iso_fortran_env, only: int32, real32, real64

    implicit none

    integer(int32), intent(in)  :: ncell ! Number of cells in the grid
    real(real64), intent(inout) :: qdp(ncell) ! Mass of each cell in original grid
    real(real64), intent(in)    :: dp1(ncell), dp2(ncell) ! Width of cells for
    ! original, new grids
    integer(int32), intent(in)  :: remap_alg ! Algorithm flag to use for
    ! remapping
    integer(int32), intent(in)  :: verbosity ! Print debug messages (1) or not (0)
    real(real64)                :: grid1(0:ncell) ! Reconstruction of the
    ! original grid.
    real(real64)                :: ccells1(ncell) ! Reconstruction of the
    ! original cell centers.
    real(real64)                :: avgdens1(ncell) ! Average density of each
    ! original cell. (<a>_j)
    real(real64)                :: grid2(0:ncell) ! Reconstruction of the
    ! new grid.
    real(real64)                :: ccells2(ncell) ! Reconstruction of the
    ! new grid.
    real(real64)                :: totmass(0:ncell) ! Cumulative mass array
    integer(int32), parameter   :: interpord = 5 ! Interpolation order for
    ! cumulative mass function
    integer(int32)              :: K(interpord, 0:ncell) ! Indices for each
    ! interpolation polynomial, for getting edge values.
    real(real64)                :: A(interpord,interpord) ! Coefficient matrix for
    ! getting polynomial interpolant. Ax=b
    real(real64)                :: x(interpord) ! Holds coefficients for polynomial
    ! interpolant. Ax=b
    real(real64)                :: b(interpord) ! Holds total mass at the grid
    ! points for the interpolant. Ax=b
    real(real64)                :: edgevals(0:ncell) ! Estimated edge values
    integer(int32)              :: limreq(0:ncell) ! Array for holding whether
    ! an edge estimate needs to be limited.
    real(real64), parameter     :: C = 1.25 ! Multiple used in the limiter.
    real(real64)                :: parabvals(3, ncell) ! Holds aj,+, aj,-, and
    ! a6,j for the parabolic interpolant in each cell.
    integer(int32)              :: cellinsecs(2, ncell) ! Holds the cell number
    ! of the original grid that the cell boundaries of the new grid intersect.
    real(real64)                :: qdp2(ncell) ! Mass of each cell in new grid
    real(real64)                :: temps_dp(4) ! Temporary real values used in
    ! calculations.
    integer(int32)              :: temps_int(2) ! Temporary integer values used
    ! in calculations.
    integer(int32)              :: ipiv(interpord) ! Array for dsgesv call.
    real(real64)                :: work(interpord, 1) ! Work array for dsgesv call.
    real(real32)                :: swork(interpord, interpord+1) ! Swork array for dsgev call.
    integer(int32)              :: jj, kk, kidx, mm, midx, ll, lidx
    ! Counters for do loops


    ! Reconstruct the original grid from 0 to end
    grid1(0) = 0.0_real64
    do jj = 1, ncell
       grid1(jj) = grid1(jj-1) + dp1(jj)
    end do
    ! Get the cell centers explicitly
    do jj = 1, ncell
       ccells1(jj) = (grid1(jj) + grid1(jj-1))/2.0_real64
    end do
    ! Get the average density of each cell
    do jj = 1, ncell
       avgdens1(jj) = qdp(jj)/dp1(jj)
    end do

    ! Get the cumulative mass function
    totmass(0) = 0
    do jj = 1, ncell
       totmass(jj) = totmass(jj-1) + qdp(jj)
    end do

    ! Set the arrays for hold the indices of the points to use for
    ! interpolation
    do jj = 0, interpord/2_int32
       K(:, jj) = [(kk, kk = 0, interpord-1)]
    end do
    do jj = interpord/2_int32 + 1, ncell - (interpord/2_int32)
       K(:, jj) = [(kk, kk = jj - interpord/2_int32, jj + interpord/2_int32)]
    end do
    do jj = ncell-(interpord/2_int32) + 1, ncell
       K(:, jj) = [(kk, kk = ncell-(interpord-1), ncell)]
    end do

    ! Get the primary estimates for edge values. We'll do the system
    ! of multiplications from Wikipedia to get the coefficients
    ! and calculate the edge value from that.
    edgevals = 0.0_real64
    do jj = 0, ncell
       edgevals(jj) = get_edgeval(ncell, interpord, K(:,jj), grid1, totmass, jj)
       ! Global limiter for exponential
!!$       if (edgevals(jj) .gt. exp(1.0_real64)) then
!!$          edgevals(jj) = exp(1.0_real64)
!!$       else if (edgevals(jj) .lt. 1.0_real64) then
!!$          edgevals(jj) = 1.0_real64
!!$       end if
    end do

    ! Now we figure out which edge value estimates need to be limited
    ! Limiter at left boundary
    if ((avgdens1(2) - avgdens1(1)) * (avgdens1(1) - edgevals(0)) .gt. 0.0_real64) then
       limreq(0) = 0
    else
       limreq(0) = 1
    end if
    ! Limiter on interior
    do jj = 1, ncell-1
       if ((avgdens1(jj+1) - edgevals(jj))*(edgevals(jj) - avgdens1(jj)) .gt. 0) then
          ! Edge value is between
          limreq(jj) = 0
       else ! Edge value estimate is outside
          limreq(jj) = 1
       end if
    end do
    ! Limiter at right boundary
    if ((edgevals(ncell) - avgdens1(ncell)) * (avgdens1(ncell) - avgdens1(ncell-1)) .gt. 0) then
       limreq(ncell) = 0
    else
       limreq(ncell) = 1
    end if

    ! /Interpolating face values./
    ! Use a modified version of the limiter from Colella08.
    ! It's fine on the interior and we don't limit on the boundary, but
    ! near the boundary we don't have quite enough information.
    do jj = 0, ncell
       if (limreq(jj) .eq. 1) then ! Limiting is required, do calculations
          edgevals(jj) = correct_edgeval(ncell, dp1, avgdens1, edgevals(jj), C, jj)
       end if
    end do

    ! /Constructing the parabolic interpolant./
    ! First we fill in the preliminary parabolic interpolant parameters for
    ! each cell. In the same loop, we correct the values appropriately.
    ! This loop can definitely be streamlined, I just wanted it to match the write up for easier debugging.

    do jj = 1, ncell
       parabvals(1, jj) = edgevals(jj-1) ! aj,-
       parabvals(2, jj) = edgevals(jj)   ! aj,+
       parabvals(3, jj) = 6.0_real64 * avgdens1(jj) - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj))
#if 1
       ! Try old parabolic correction algorithm.
       call correct_parabvals(ncell, dp1, avgdens1, parabvals(:, jj), C, jj)
#endif
#if 0
       ! Try new parabolic correction algorithm.
       call correct_parabvals_2(ncell, avgdens1, parabvals(:, jj), jj)
#endif
       
    end do

    ! Check that each parabolic piece is monotone
    if (verbosity .eq. 1_int32) then
       temps_int(1) = 0_int32
       do jj = 1, ncell
          if (parab_piece_monotone_check(ncell, parabvals(:, jj), jj) .eq. 0_int32) then ! If not monotone, trip flag.
             temps_int(1) = 1_int32
             !print *, '  ~~ !WARNING! Parabolic piece is not monotone:', jj
          end if
       end do

       if (temps_int(1) .eq. 1_int32) then
          print *, '  ~~ !WARNING! At least one parabolic piece is not monotone.'
       end if
    end if

    ! Check that each parabolic piece is bounds-preserving
    if (verbosity .eq. 1_int32) then
       temps_int(1) = 0_int32
       do jj = 1, ncell
          if (parab_piece_local_bnd_preserve_check(ncell, parabvals(:, jj), avgdens1, jj) .eq. 0_int32) then ! If not
             ! bounds-preserving, trip flag.
             temps_int(1) = 1_int32
             print *, '  ~~ !WARNING! Parabolic piece is not bounds-preserving:', jj
          end if
       end do

       if (temps_int(1) .eq. 1_int32) then
          !print *, '  ~~ !WARNING! At least one parabolic piece is not bounds-preserving.'
       end if
    end if


    ! Now I have the piecewise parabolic reconstruction, now I have to get the mass in each of the new cells.
    ! I want a neat formula for integrating the piecewise parabolic reeconstruction.
    ! I have the neat formula, I just need to find which cells then new cell boundaries are in.
    ! This will require reconstructing the new grid.

    ! To get qdp on the new grid, we must first reconstruct the new grid.
    grid2(0) = 0.0_real64
    do jj = 1, ncell
       grid2(jj) = grid2(jj-1) + dp2(jj)
    end do
    ! And cell centers
    do jj = 1, ncell
       ccells2(jj) = (grid2(jj) + grid2(jj-1))/2.0_real64
    end do


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

    ! Now we integrate the parabolic pieces to get the new masses
    do jj = 1, ncell
       temps_dp = 0.0_real64

       if (cellinsecs(2, jj) .gt. cellinsecs(1, jj) + 1) then ! New cell contains at least one entire old cell
          ! Left part of new cell is in part of an old cell
          temps_dp(1) = integrate_parab_piece(grid2(jj-1), grid1(cellinsecs(1, jj)), &
               & grid1(cellinsecs(1, jj)-1), dp1(cellinsecs(1, jj)), parabvals(:, cellinsecs(1, jj)))
          ! Middle part of new cell spans at least one cell
          do kk = cellinsecs(1, jj) + 1, cellinsecs(2, jj) - 1
             temps_dp(2) = temps_dp(2) + qdp(kk)
          end do
          ! Right part of new cell is in part of an old cell
          temps_dp(3) = integrate_parab_piece(grid1(cellinsecs(2, jj)-1), grid2(jj), &
               & grid1(cellinsecs(2, jj)-1), dp1(cellinsecs(2, jj)), parabvals(:, cellinsecs(2, jj)))
       else if (cellinsecs(2, jj) .eq. cellinsecs(1, jj) + 1) then ! New cell intersects two old cells
          ! Left part of new cell is in part of an old cell
          temps_dp(1) = integrate_parab_piece(grid2(jj-1), grid1(cellinsecs(1, jj)), &
               & grid1(cellinsecs(1, jj)-1), dp1(cellinsecs(1, jj)), parabvals(:, cellinsecs(1, jj)))
          ! Right part of new cell is in part of an old cell
          temps_dp(3) = integrate_parab_piece(grid1(cellinsecs(2, jj)-1), grid2(jj), &
               & grid1(cellinsecs(2, jj)-1), dp1(cellinsecs(2, jj)), parabvals(:, cellinsecs(2, jj)))
       else if (cellinsecs(2, jj) .eq. cellinsecs(1, jj)) then ! New cell in one old cell
          ! New cell is in part of an old cell
          temps_dp(1) = integrate_parab_piece(grid2(jj-1), grid2(jj), &
               & grid1(cellinsecs(1, jj)-1), dp1(cellinsecs(1, jj)), parabvals(:, cellinsecs(1, jj)))
       end if

       qdp2(jj) = sum(temps_dp(1:3))
    end do

    qdp = qdp2

  end subroutine new_remap_sbr

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Gets the edge value approximation.
  ! In particular, gets the derivtaive of a interpolated polynomial for
  ! the cumulative mass function

  function get_edgeval(ncell, interpord, K, grid1, totmass, cell) result(edgeval)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells
    integer(int32), intent(in) :: interpord ! Order of polynomial interpolation
    integer(int32), intent(in) :: K(interpord) ! Points for polynomial intrpolation
    real(real64), intent(in) :: grid1(0:ncell) ! Grid points of original grid
    real(real64), intent(in) :: totmass(0:ncell) ! Cumulative mass function
    integer(int32), intent(in) :: cell ! Cell we are working with
    real(real64) :: edgeval !Approximated edge value
    integer(int32) :: kk, kidx, ll, lidx ! For iteration
    real(real64) :: temps_dp(3) ! Array for storing working variables (real)

    temps_dp = 0.0_real64
    edgeval = 0.0_real64
    do kk = 1, interpord
       kidx = K(kk)
       if (kidx .ne. cell) then
          temps_dp(1) = 1.0_real64
          do ll = 1, interpord
             lidx = K(ll)
             if ((lidx .ne. cell) &
                  & .and. (lidx .ne.kidx) ) then
                temps_dp(1) = temps_dp(1) * (grid1(cell) - grid1(lidx)) &
                     & / (grid1(kidx) - grid1(lidx))
             end if
          end do
          edgeval = edgeval + totmass(kidx) / (grid1(kidx) - grid1(cell)) * temps_dp(1)

       else if (kidx .eq. cell) then
          temps_dp(1) = 0.0_real64
          do ll = 1, interpord
             lidx = K(ll)
             if (lidx .ne. cell) then
                temps_dp(1) = temps_dp(1) + (1.0_real64 / (grid1(cell) - grid1(lidx)))
             end if
          end do
          edgeval = edgeval + totmass(kidx) * temps_dp(1)
       end if
    end do

  end function get_edgeval

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! First step of the limiter, make sure edge values are monotone with neighboring cell averages.
  function correct_edgeval(ncell, dp1, avgdens, edgeval, C, cell) result(corr_edgeval)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells in the grid
    real(real64), intent(in) :: dp1(ncell) ! Width of cells on original grid
    real(real64), intent(in) :: avgdens(ncell) ! Average density of original cells
    real(real64), intent(in) :: edgeval ! Original edge value
    real(real64), intent(in) :: C ! Scale factor for second derivative approximations
    integer(int32), intent(in) :: cell ! Current cell we are correcting the edge value for
    real(real64) :: corr_edgeval ! Corrected edge value
    real(real64) :: temps_dp(4) ! Holds working real variables

    temps_dp = 0.0_real64
    if (cell .eq. 0) then

       ! Left domain boundary
       temps_dp(1) = 8.0_real64 * &
            & ( edgeval / (dp1(1) * (2.0_real64 * dp1(1) + dp1(2))) &
            &   - avgdens(1) / (dp1(1) * (dp1(1) + dp1(2))) &
            &   + avgdens(2) / ((dp1(1) + dp1(2)) * (2.0_real64 * dp1(1) + dp1(2))) )
       temps_dp(2) = 8.0_real64 * &
            & ( avgdens(1) / ((dp1(1) + dp1(2)) * (dp1(1) + 2.0_real64 * dp1(2) + dp1(3))) &
            &   - avgdens(2) / ((dp1(1) + dp1(2)) * (dp1(2) + dp1(3))) &
            &   + avgdens(3) / ((dp1(2) + dp1(3)) * (dp1(1) + 2.0_real64 * dp1(2) + dp1(3))) )

       ! Check is approxs have same sign
       if (temps_dp(1) * temps_dp(2) .gt. 0) then
          temps_dp(1) = sign(min(C * abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1))
       else
          temps_dp(1) = 0.0_real64
       end if

       ! Set temps_dp(2) to grid spacing for correction
       temps_dp(2) = (dp1(1) + dp1(2)) / 2.0_real64

       corr_edgeval = avgdens(1) + (temps_dp(2))**2/4.0_real64 * temps_dp(1)

    else if (cell .eq. 1) then
       ! Right boundary of leftmost cell
       temps_dp(1) = 8.0_real64 * &
            & ( avgdens(1) / (dp1(1) * (dp1(2) + dp1(1))) &
            &   - edgeval / (dp1(1) * dp1(2)) &
            &   + avgdens(2) / (dp1(2) * (dp1(2) + dp1(1))) )
       temps_dp(2) = 8.0_real64 * &
            & ( avgdens(1) / (dp1(2) * (dp1(3) + dp1(2))) &
            &   - avgdens(2) / (dp1(2) * dp1(3)) &
            &   + avgdens(3) / (dp1(3) * (dp1(3) + dp1(2))) )

       ! Check is approxs have same sign
       if (temps_dp(1) * temps_dp(2) .gt. 0) then
          temps_dp(1) = sign(min(C * abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1))
       else
          temps_dp(1) = 0.0_real64
       end if

       ! Set temps_dp(2) to grid spacing for correction
       temps_dp(2) = (dp1(1) + dp1(2)) / 2.0_real64

       corr_edgeval = (avgdens(1) + avgdens(2))/2.0_real64 - (temps_dp(2))**2/4.0_real64 * temps_dp(1)

    else if (cell .eq. ncell-1) then
       ! Left boundary of rightmost cell
       temps_dp(1) = 8.0_real64 * &
            & ( avgdens(ncell-1) / (dp1(ncell-1) * (dp1(ncell) + dp1(ncell-1))) &
            &   - edgeval / (dp1(ncell-1) * dp1(ncell)) &
            &   + avgdens(ncell) / (dp1(ncell) * (dp1(ncell) + dp1(ncell-1))) )
       temps_dp(2) = 8.0_real64 * &
            & ( avgdens(ncell-2) / (dp1(ncell-1) * (dp1(ncell) + dp1(ncell-1))) &
            &   - avgdens(ncell-1) / (dp1(ncell-1) * dp1(ncell)) &
            &   + avgdens(ncell) / (dp1(ncell) * (dp1(ncell) + dp1(ncell-1))) )

       ! Check is approxs have same sign
       if (temps_dp(1) * temps_dp(2) .gt. 0) then
          temps_dp(1) = sign(min(C * abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1))
       else
          temps_dp(1) = 0.0_real64
       end if

       ! Set temps_dp(2) to grid spacing for correction
       temps_dp(2) = (dp1(ncell-1) + dp1(ncell)) / 2.0_real64

       corr_edgeval = (avgdens(ncell-1) + avgdens(ncell))/2.0_real64 - (temps_dp(2))**2/4.0_real64 * temps_dp(1)   

    else if (cell .eq. ncell) then ! Right domain boundary
       temps_dp(1) = 8.0_real64 * &
            & ( avgdens(ncell-1) / ((dp1(ncell) + dp1(ncell-1) * (2.0_real64 * dp1(ncell) + dp1(ncell-1)))) &
            &   - avgdens(ncell) / (dp1(ncell) * (dp1(ncell) + dp1(ncell-1))) &
            &   + edgeval / (dp1(ncell) * (2.0_real64 * dp1(ncell) + dp1(ncell-1))) )
       temps_dp(2) = 8.0_real64 * &
            & ( avgdens(ncell-2) / ((dp1(ncell-1) + dp1(ncell-2)) &
            &                         * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) &
            &   - avgdens(ncell-1) / ((dp1(ncell) + dp1(ncell-1)) * (dp1(ncell-1) + dp1(ncell-2))) &
            &   + avgdens(ncell) / ((dp1(ncell) + dp1(ncell-1)) * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) )

       ! Check is approxs have same sign
       if (temps_dp(1) * temps_dp(2) .gt. 0) then
          temps_dp(1) = sign(min(C * abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1))
       else
          temps_dp(1) = 0.0_real64
       end if

       ! Set temps_dp(2) to grid spacing for correction
       temps_dp(2) = (dp1(1) + dp1(2)) / 2.0_real64

       corr_edgeval = avgdens(ncell) + (temps_dp(2))**2/4.0_real64 * temps_dp(1)

    else ! Boundaries of interior cells
       temps_dp(1) = 8.0_real64 * &
            & ( avgdens(cell) / (dp1(cell) * (dp1(cell+1) + dp1(cell))) &
            &   - edgeval / (dp1(cell) * dp1(cell+1)) &
            &   + avgdens(cell+1) / (dp1(cell+1) * (dp1(cell+1) + dp1(cell))) )
       temps_dp(2) = 8.0_real64 * &
            & ( avgdens(cell-1) / (dp1(cell-1) * (dp1(cell+1) + dp1(cell) + dp1(cell-1))) &
            &   - avgdens(cell) / (dp1(cell-1) * (dp1(cell+1) + dp1(cell))) &
            &   + avgdens(cell+1) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+1) + dp1(cell) + dp1(cell-1))) )
       temps_dp(3) = 8.0_real64 * &
            & ( avgdens(cell) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) &
            &   - avgdens(cell+1) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+2) + dp1(cell+1))) &
            &   + avgdens(cell+2) / ((dp1(cell+2) + dp1(cell+1)) * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) )

       ! Check if approxs have same sign
       if ((temps_dp(1) * temps_dp(2) .gt. 0) &
            & .and. (temps_dp(2) * temps_dp(3) .gt. 0)) then
          temps_dp(1) = sign(min(abs(temps_dp(1)), C * abs(temps_dp(2)), C * abs(temps_dp(3))), temps_dp(1))
       else
          temps_dp(1) = 0.0_real64
       end if

       ! Set temps_dp(2) to grid spacing for correction
       temps_dp(2) = (dp1(cell) + dp1(cell+1)) / 2.0_real64

       corr_edgeval = (avgdens(cell) + avgdens(cell+1))/2.0_real64 + (temps_dp(2))**2/4.0_real64 * temps_dp(1)

    end if

  end function correct_edgeval

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine correct_parabvals(ncell, dp1, avgdens, parabvals, C, cell)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells in the grid
    real(real64), intent(in) :: dp1(ncell) ! Width of cells on original grid
    real(real64), intent(in) :: avgdens(ncell) ! Average density of original cells
    real(real64), intent(inout) :: parabvals(3) ! Original parabola values
    real(real64), intent(in) :: C ! Scale factor for second derivative approximations
    integer(int32), intent(in) :: cell ! Current cell we are correcting the edge value for
    real(real64) :: corr_edgeval ! Corrected edge value
    real(real64) :: temps_dp(4) ! Holds working real variables
    integer(int32) :: temps_int(1)

    ! The boundary cells must be handled differently than Colella08.
    if (cell .eq. 1) then ! Leftmost cell
       if ((parabvals(2) - avgdens(1)) * (avgdens(1) - parabvals(1)) .le. 0.0_real64) then ! At local extremum
          temps_dp(1) = 4.0_real64 / (dp1(1)**2) * (parabvals(1) - 2.0_real64 * avgdens(1) + parabvals(2))
          temps_dp(2) = 8.0_real64 * &
               & (avgdens(1) / ((dp1(2) + dp1(1)) * (dp1(3) + 2.0_real64 * dp1(2) + dp1(1))) &
               &  - avgdens(2) / ((dp1(2) + dp1(1)) * (dp1(3) + dp1(2))) &
               &  + avgdens(3) / ((dp1(3) + dp1(2)) * (dp1(3) + 2.0_real64 * dp1(2) + dp1(1))) )

          ! Check if second derivative signs match. Put D^2a_1,lim in temps_dp(2)
          if (temps_dp(1) * temps_dp(2) .gt. 0.0_real64) then
             temps_dp(2) = sign( min(abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1) )
          else
             temps_dp(2) = 0.0_real64
          end if

          ! Correct parabolic piece values
          if (temps_dp(1) .ne. 0) then
             parabvals(1) = avgdens(1) + (parabvals(1) - avgdens(1)) * temps_dp(2) / temps_dp(1)
             parabvals(2) = avgdens(1) + (parabvals(2) - avgdens(1)) * temps_dp(2) / temps_dp(1)
          else
             parabvals(1) = avgdens(1)
             parabvals(2) = avgdens(1)
          end if

       else ! Not at local extremum
          ! Set s in temps_int(1)
          if (avgdens(2) - avgdens(1) .gt. 0.0_real64) then
             temps_int(1) = 1_int32
          else
             temps_int(1) = -1_int32
          end if

          temps_dp(1) = parabvals(1) - avgdens(1) ! alpha_{1,-}
          temps_dp(2) = parabvals(2) - avgdens(1) ! alpha_{1,+}

          if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(2)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(2) - avgdens(1) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(2) = avgdens(1) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(1)))
             end if
          end if

       end if

       parabvals(3) = 6.0_real64 * avgdens(1) - 3.0_real64 * (parabvals(1) + parabvals(2))

    else if (cell .eq. 2) then ! Second leftmost cell
       if (((parabvals(2) - avgdens(2)) * (avgdens(2) - parabvals(1)) .le. 0.0_real64) &
            & .or. ((avgdens(3) - avgdens(2)) * (avgdens(2) - avgdens(1)) .le. 0.0_real64)) then ! At local extremum
          temps_dp(1) = 4.0_real64 / (dp1(2)**2) * (parabvals(1) - 2.0_real64 * avgdens(2) + parabvals(2))
          temps_dp(2) = 8.0_real64 * &
               & (avgdens(1) / ((dp1(2) + dp1(1)) * (dp1(3) + 2.0_real64 * dp1(2) + dp1(1))) &
               &  - avgdens(2) / ((dp1(2) + dp1(1)) * (dp1(3) + dp1(2))) &
               &  + avgdens(3) / ((dp1(3) + dp1(2)) * (dp1(3) + 2.0_real64 * dp1(2) + dp1(1))) )
          temps_dp(3) = 8.0_real64 * &
               & (avgdens(2) / ((dp1(3) + dp1(2)) * (dp1(4) + 2.0_real64 * dp1(3) + dp1(2))) &
               &  - avgdens(3) / ((dp1(3) + dp1(2)) * (dp1(4) + dp1(3))) &
               &  + avgdens(4) / ((dp1(4) + dp1(3)) * (dp1(4) + 2.0_real64 * dp1(3) + dp1(2))) )

          ! Check if second derivative signs match. Put D^2a_1,lim in temps_dp(2)
          if ((temps_dp(1) * temps_dp(2) .gt. 0.0_real64) &
               & .and. (temps_dp(2) * temps_dp(3) .gt. 0.0_real64)) then
             temps_dp(2) = sign( min(abs(temps_dp(1)), C * abs(temps_dp(2)), C * abs(temps_dp(3))), temps_dp(1) )
          else
             temps_dp(2) = 0.0_real64
          end if

          ! Correct parabolic piece values
          if (temps_dp(1) .ne. 0) then
             parabvals(1) = avgdens(2) + (parabvals(1) - avgdens(2)) * temps_dp(2) / temps_dp(1)
             parabvals(2) = avgdens(2) + (parabvals(2) - avgdens(2)) * temps_dp(2) / temps_dp(1)
          else
             parabvals(1) = avgdens(2)
             parabvals(2) = avgdens(2)
          end if

       else ! Not at local extremum
          ! Set s in temps_int(1)
          if (avgdens(3) - avgdens(1) .gt. 0.0_real64) then
             temps_int(1) = 1_int32
          else
             temps_int(1) = -1_int32
          end if

          temps_dp(1) = parabvals(1) - avgdens(2) ! alpha_{2,-}
          temps_dp(2) = parabvals(2) - avgdens(2) ! alpha_{2,+}

          if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(2)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(3) - avgdens(2) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(2) = avgdens(2) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(1)))

             end if

          else if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(1)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(1) - avgdens(2) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(1) = avgdens(2) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(2)))

             end if
          end if
       end if

       parabvals(3) = 6.0_real64 * avgdens(2) - 3.0_real64 * (parabvals(1) + parabvals(2))

    else if (cell .eq. ncell - 1) then ! Second rightmost cell
       if (((parabvals(2) - avgdens(ncell-1)) * (avgdens(ncell-1) - parabvals(1)) .le. 0.0_real64) &
            & .or. ((avgdens(ncell) - avgdens(ncell-1)) * (avgdens(ncell-1) - avgdens(ncell-2)) .le. 0.0_real64)) then ! At local extremum
          temps_dp(1) = 4.0_real64 / (dp1(ncell-1)**2) &
               & * (parabvals(1) - 2.0_real64 * avgdens(ncell-1) + parabvals(2))
          temps_dp(2) = 8.0_real64 * &
               & (avgdens(ncell-3) / ((dp1(ncell-2) + dp1(ncell-3)) &
               &                        * (dp1(ncell-1) + 2.0_real64 * dp1(ncell-2) + dp1(ncell-3))) &
               &  - avgdens(ncell-2) / ((dp1(ncell-2) + dp1(ncell-3)) &
               &                          * (dp1(ncell-1) + dp1(ncell-2))) &
               &  + avgdens(ncell-1) / ((dp1(ncell-1) + dp1(ncell-2)) &
               &                         * (dp1(ncell-1) + 2.0_real64 * dp1(ncell-2) + dp1(ncell-3))) )
          temps_dp(3) = 8.0_real64 * &
               & (avgdens(ncell-2) / ((dp1(ncell-1) + dp1(ncell-2)) &
               &                        * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) &
               &  - avgdens(ncell-1) / ((dp1(ncell-1) + dp1(ncell-2)) &
               &                          * (dp1(ncell) + dp1(ncell-1))) &
               &  + avgdens(ncell) / ((dp1(ncell) + dp1(ncell-1)) &
               &                        * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) )

          ! Check if second derivative signs match. Put D^2a_n-1,lim in temps_dp(2)
          if ((temps_dp(1) * temps_dp(2) .gt. 0.0_real64) &
               & .and. (temps_dp(2) * temps_dp(3) .gt. 0.0_real64)) then
             temps_dp(2) = sign( min(abs(temps_dp(1)), C * abs(temps_dp(2)), C * abs(temps_dp(3))), temps_dp(1) )
          else
             temps_dp(2) = 0.0_real64
          end if

          ! Correct parabolic piece values
          if (temps_dp(1) .ne. 0) then
             parabvals(1) = avgdens(ncell-1) + (parabvals(1) - avgdens(ncell-1)) * temps_dp(2) / temps_dp(1)
             parabvals(2) = avgdens(ncell-1) + (parabvals(2) - avgdens(ncell-1)) * temps_dp(2) / temps_dp(1)
          else
             parabvals(1) = avgdens(ncell-1)
             parabvals(2) = avgdens(ncell-1)
          end if

       else ! Not at local extremum
          ! Set s in temps_int(1)
          if (avgdens(ncell) - avgdens(ncell-2) .gt. 0.0_real64) then
             temps_int(1) = 1_int32
          else
             temps_int(1) = -1_int32
          end if

          temps_dp(1) = parabvals(1) - avgdens(ncell-1) ! alpha_{n-1,-}
          temps_dp(2) = parabvals(2) - avgdens(ncell-1) ! alpha_{n-1,+}

          if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(2)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(ncell) - avgdens(ncell-1) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(2) = avgdens(ncell-1) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(1)))

             end if

          else if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(1)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(ncell-2) - avgdens(ncell-1) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(1) = avgdens(ncell-1) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(2)))

             end if
          end if
       end if

       parabvals(3) = 6.0_real64 * avgdens(ncell-1) - 3.0_real64 * (parabvals(1) + parabvals(2))

    else if (cell .eq. ncell) then ! Rightmost cell
       if ((parabvals(2) - avgdens(ncell)) * (avgdens(ncell) - parabvals(1)) .le. 0.0_real64) then ! At local extremum
          temps_dp(1) = 4.0_real64 / (dp1(ncell)**2) * (parabvals(1) - 2.0_real64 * avgdens(ncell) + parabvals(2))
          temps_dp(2) = 8.0_real64 * &
               & (avgdens(ncell-2) / ((dp1(ncell-1) + dp1(ncell-2)) * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) &
               &  - avgdens(ncell-1) / ((dp1(ncell-1) + dp1(ncell-2)) * (dp1(ncell) + dp1(ncell-1))) &
               &  + avgdens(ncell) / ((dp1(ncell) + dp1(ncell-1)) * (dp1(ncell) + 2.0_real64 * dp1(ncell-1) + dp1(ncell-2))) )

          ! Check if second derivative signs match. Put D^2a_n,lim in temps_dp(2)
          if (temps_dp(1) * temps_dp(2) .gt. 0.0_real64) then
             temps_dp(2) = sign( min(abs(temps_dp(1)), C * abs(temps_dp(2))), temps_dp(1) )
          else
             temps_dp(2) = 0.0_real64
          end if

          ! Correct parabolic piece values
          if (temps_dp(1) .ne. 0) then
             parabvals(1) = avgdens(ncell) + (parabvals(1) - avgdens(ncell)) * temps_dp(2) / temps_dp(1)
             parabvals(2) = avgdens(ncell) + (parabvals(2) - avgdens(ncell)) * temps_dp(2) / temps_dp(1)
          else
             parabvals(1) = avgdens(ncell)
             parabvals(2) = avgdens(ncell)
          end if

       else ! Not at local extremum
          ! Set s in temps_int(1)
          if (avgdens(ncell) - avgdens(ncell-1) .gt. 0.0_real64) then
             temps_int(1) = 1_int32
          else
             temps_int(1) = -1_int32
          end if

          temps_dp(1) = parabvals(1) - avgdens(ncell) ! alpha_{n,-}
          temps_dp(2) = parabvals(2) - avgdens(ncell) ! alpha_{n,+}

          if (abs(temps_dp(1)) .gt. 2.0_real64 * abs(temps_dp(2))) then
             temps_dp(3) = -1.0_real64 * temps_dp(1)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(ncell-1) - avgdens(ncell) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(1) = avgdens(ncell) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(2)))
             end if
          end if

       end if

       parabvals(3) = 6.0_real64 * avgdens(ncell) - 3.0_real64 * (parabvals(1) + parabvals(2))

    else ! Interior cell
       if (((parabvals(2) - avgdens(cell)) * (avgdens(cell) - parabvals(1)) .le. 0.0_real64) &
            & .or. ((avgdens(cell+1) - avgdens(cell)) * (avgdens(cell) - avgdens(cell-1)) .le. 0.0_real64)) then ! At local extremum
          temps_dp(1) = 4.0_real64 / (dp1(cell)**2) * (parabvals(1) - 2.0_real64 * avgdens(cell) + parabvals(2))
          temps_dp(2) = 8.0_real64 * &
               & (avgdens(cell-2) / ((dp1(cell-1) + dp1(cell-2)) * (dp1(cell) + 2.0_real64 * dp1(cell-1) + dp1(cell-2))) &
               &  - avgdens(cell-1) / ((dp1(cell-1) + dp1(cell-2)) * (dp1(cell) + dp1(cell-1))) &
               &  + avgdens(cell) / ((dp1(cell) + dp1(cell-1)) * (dp1(cell) + 2.0_real64 * dp1(cell-1) + dp1(cell-2))) )
          temps_dp(3) = 8.0_real64 * &
               & (avgdens(cell-1) / ((dp1(cell) + dp1(cell-1)) * (dp1(cell+1) + 2.0_real64 * dp1(cell) + dp1(cell-1))) &
               &  - avgdens(cell) / ((dp1(cell) + dp1(cell-1)) * (dp1(cell+1) + dp1(cell))) &
               &  + avgdens(cell+1) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+1) + 2.0_real64 * dp1(cell) + dp1(cell-1))) )
          temps_dp(4) = 8.0_real64 * &
               & (avgdens(cell) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) &
               &  - avgdens(cell+1) / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+2) + dp1(cell+1))) &
               &  + avgdens(cell+2) / ((dp1(cell+2) + dp1(cell+1)) * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) )

          ! Check if second derivative signs match. Put D^2a_n,lim in temps_dp(2)
          if ((temps_dp(1) * temps_dp(2) .gt. 0.0_real64) &
               & .and. (temps_dp(2) * temps_dp(3) .gt. 0.0_real64) &
               & .and. (temps_dp(3) * temps_dp(4) .gt. 0.0_real64)) then
             temps_dp(2) = sign( min(abs(temps_dp(1)), C * abs(temps_dp(2)), &
                  &                  C * abs(temps_dp(3)), C * abs(temps_dp(4))), temps_dp(1) )
          else
             temps_dp(2) = 0.0_real64
          end if

          ! Correct parabolic piece values
          if (temps_dp(1) .ne. 0) then
             parabvals(1) = avgdens(cell) + (parabvals(1) - avgdens(cell)) * temps_dp(2) / temps_dp(1)
             parabvals(2) = avgdens(cell) + (parabvals(2) - avgdens(cell)) * temps_dp(2) / temps_dp(1)
          else
             parabvals(1) = avgdens(cell)
             parabvals(2) = avgdens(cell)
          end if

       else ! Not at local extremum
          ! Set s in temps_int(1)
          if (avgdens(cell+1) - avgdens(cell-1) .gt. 0.0_real64) then
             temps_int(1) = 1_int32
          else
             temps_int(1) = -1_int32
          end if

          temps_dp(1) = parabvals(1) - avgdens(cell) ! alpha_{j,-}
          temps_dp(2) = parabvals(2) - avgdens(cell) ! alpha_{j,+}

          if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(2)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(cell+1) - avgdens(cell) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(2) = avgdens(cell) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(1)))

             end if

          else if (abs(temps_dp(2)) .gt. 2.0_real64 * abs(temps_dp(1))) then
             temps_dp(3) = -1.0_real64 * temps_dp(1)**2 / (4.0_real64 * (temps_dp(2) + temps_dp(1))) ! delta I_ext
             temps_dp(4) = avgdens(cell-1) - avgdens(cell) ! delta a
             if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                parabvals(1) = avgdens(cell) &
                     & - 2.0_real64 * (temps_dp(4) &
                     &                 + temps_int(1) * sqrt(temps_dp(4)**2 - temps_dp(4) * parabvals(2)))

             end if
          end if
       end if

       parabvals(3) = 6.0_real64 * avgdens(cell) - 3.0_real64 * (parabvals(1) + parabvals(2))

    end if

  end subroutine correct_parabvals

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine correct_parabvals_2(ncell, avgdens, parabvals, cell)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells in the grid
    real(real64), intent(in) :: avgdens(ncell) ! Average density of original cells
    real(real64), intent(inout) :: parabvals(3) ! Original parabola values
    integer(int32), intent(in) :: cell ! Current cell we are correcting the edge value for
    real(real64) :: crit_pt ! Critical point for parabolic piece
    real(real64) :: epsilon ! Adjustment number for endpoints.

    ! Endpoints same value, must be constant.
    if (abs(parabvals(1) - parabvals(2)) .le. 1.0e-12_real64) then
       parabvals(1) = avgdens(cell)
       parabvals(2) = avgdens(cell)
       parabvals(3) = 0.0_real64
    end if

    if (abs(parabvals(3)) .ge. 1.0e-12_real64) then ! Is non-linear, might
       ! not be monotone.
       crit_pt = ((parabvals(2) - parabvals(1)) / parabvals(3) + 1) &
            & / 2.0_real64
       if ((crit_pt .lt. 1.0_real64) .and. (crit_pt .gt. 0.0_real64)) then
          ! Critical point is in interval of interest.
          if (((parabvals(3) .lt. 0) .and. (parabvals(1) .lt. parabvals(2))) &
               & .or. ((parabvals(3) .gt. 0) .and. (parabvals(1) .gt. parabvals(2)))) then

             epsilon = 2.0_real64 * parabvals(1) + parabvals(2) &
                  & - 3.0_real64 * avgdens(cell)

             parabvals(1) = parabvals(1) - (1.0_real64 / 3.0_real64) * epsilon
             parabvals(2) = parabvals(2) - (1.0_real64 / 3.0_real64) * epsilon
             
          else if (((parabvals(3) .gt. 0) .and. (parabvals(1) .lt. parabvals(2))) &
               & .or. ((parabvals(3) .lt. 0) .and. (parabvals(1) .gt. parabvals(2)))) then

             epsilon = 2.0_real64 * parabvals(2) + parabvals(1) &
                  & - 3.0_real64 * avgdens(cell)

             parabvals(1) = parabvals(1) - (1.0_real64 / 3.0_real64) * epsilon
             parabvals(2) = parabvals(2) - (1.0_real64 / 3.0_real64) * epsilon
             
          end if

          parabvals(3) = 6.0_real64 * avgdens(cell) &
               & - 3.0_real64 * (parabvals(1) + parabvals(2))
       end if
    end if

  end subroutine correct_parabvals_2
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Integrates a piece of the piecewise-parabolic reconstruction
  function integrate_parab_piece(lb, rb, lbcell, dcell, parabvals) result(res)

    use iso_fortran_env, only: int32, real64

    implicit none

    real(real64), intent(in)  :: lb, rb ! Bounds for integration.
    real(real64), intent(in)  :: lbcell, dcell ! Left boundary of cell and cell width
    real(real64), intent(in)  :: parabvals(3) ! aj,-, aj,+, and a6,j for parabolic piece
    real(real64)              :: res
    real(real64)              :: coeffs(3) ! Coefficients of a + bx + cx^2

    coeffs(1) = parabvals(1) &
         & - (parabvals(2) - parabvals(1) + parabvals(3)) * lbcell / dcell &
         & - parabvals(3) * lbcell**2 / dcell**2
    coeffs(2) = (parabvals(2) - parabvals(1) + parabvals(3)) / dcell &
         & + 2.0_real64 * parabvals(3) * lbcell / dcell**2
    coeffs(3) = - parabvals(3) / dcell**2

    res = coeffs(1) * (rb - lb) &
         & + coeffs(2)/2.0_real64 * (rb**2 - lb**2) &
         & + coeffs(3)/3.0_real64 * (rb**3 - lb**3)

  end function integrate_parab_piece

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Checks that a parabolic piece is monotone.
  function parab_piece_monotone_check(ncell, parabvals, cell) result(res)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell
    real(real64), intent(in) :: parabvals(3)
    integer(int32), intent(in) :: cell
    integer(int32) :: res
    real(real64) :: temp_dp

    temp_dp = parabvals(1) - parabvals(2) + parabvals(3)
    if (abs(parabvals(3)) .le. 1.0e-10_real64) then
       ! Is monotone (sufficient, not necessary).
       res = 1_int32
    else if ((temp_dp / (2.0_real64 * parabvals(3)) .gt. 0.0_real64) &
         & .and. (temp_dp / (2.0_real64 * parabvals(3)) .le. 1.0_real64)) then
       ! Is not monotone (necessary and sufficient).
       res = 0_int32
    end if

  end function parab_piece_monotone_check

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Checks that a parabolic piece is locally bounds-preserving.
  function parab_piece_local_bnd_preserve_check(ncell, parabvals, avgdens, cell) result(res)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell
    real(real64), intent(in) :: parabvals(3)
    real(real64), intent(in) :: avgdens(ncell)
    integer(int32), intent(in) :: cell
    integer(int32) :: res
    real(real64) :: lo_bnd,up_bnd ! Upper and lower bound for parabolic piece
    real(real64) :: in_ext ! Inner extremum of parabolic piece, if not monotone
    real(real64) :: min_parab, max_parab ! Minimun and maximum values of parabola in the interval
    real(real64) :: inf_pt ! Xi coordinate of infelction point
    integer(int32) :: have_extrema ! Flag to indicate the extrema have been calculated

    have_extrema = 0

    ! Check if parabolic piece has inflection point in [0, 1]
    if (abs(parabvals(3)) .gt. 1.0e-10_real64) then ! Not linear, get inflection point
       inf_pt = (parabvals(1) - parabvals(2) + parabvals(3)) &
            & / parabvals(3) ! Coordinate of inflection point
       if ((inf_pt .le. 1.0_real64) &
            & .and. (inf_pt .ge. 0.0_real64)) then ! If inflection point in the interval,
          ! Then get possible internal extremum, and calculate extrema
          in_ext = parabvals(1) &
               & + (parabvals(2) - parabvals(1) + parabvals(3)) * inf_pt &
               & - parabvals(3) * (inf_pt)**2
          min_parab = min(parabvals(1), parabvals(2), in_ext)
          max_parab = max(parabvals(1), parabvals(2), in_ext)
          have_extrema = 1_int32 ! Mark that extrema have been found
       end if
    end if

    if (have_extrema .eq. 0) then ! Have not gotten extrema, so piece is monotone
       min_parab = min(parabvals(1), parabvals(2))
       max_parab = max(parabvals(1), parabvals(2))
       have_extrema = 1_int32 ! Mark that extrema have been found
    end if

    ! Get the bounds the parabola must be between
    if (cell .eq. 1_int32) then ! Leftmost cell
       lo_bnd = min(avgdens(cell), avgdens(cell+1))
       up_bnd = max(avgdens(cell), avgdens(cell+1))
    else if (cell .eq. ncell) then ! Rightmost cell
       lo_bnd = min(avgdens(cell-1), avgdens(cell))
       up_bnd = max(avgdens(cell-1), avgdens(cell))
    else ! Interior cell
       lo_bnd = min(avgdens(cell-1), avgdens(cell), avgdens(cell+1))
       up_bnd = max(avgdens(cell-1), avgdens(cell), avgdens(cell+1))
    end if

    ! Compare max, min values
    if ((min_parab .le. lo_bnd) .or. (max_parab .ge. up_bnd)) then
       res = 0_int32
    else
       res = 1_int32
    end if

  end function parab_piece_local_bnd_preserve_check

end submodule vertremap_redux
