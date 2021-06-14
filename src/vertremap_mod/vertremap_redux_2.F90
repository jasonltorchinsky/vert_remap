!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! New vertical remapping code version 2, should take in generally the same
! arguments as the original, just simplified.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (vertremap_mod) vertremap_redux_2

contains

  module subroutine new_remap_2_sbr(qdp, ncell, dp1, dp2, remap_alg, verbosity)

    use iso_fortran_env, only: int32, real32, real64

    implicit none

    integer(int32), intent(in)  :: ncell ! Number of cells in the grid
    real(real64), intent(inout) :: qdp(ncell) ! Mass of each cell in original
    ! grid.
    real(real64), intent(in)    :: dp1(ncell), dp2(ncell) ! Width of cells for
    ! original, new grids.
    integer(int32), intent(in)  :: remap_alg ! Algorithm flag to use for
    ! remapping.
    integer(int32), intent(in)  :: verbosity ! Print debug messages (1) or
    ! not (0).
    integer(int32), parameter   :: interpcnt = 5 ! Number of points to
    ! extrapolate/interpolate over.
    real(real64), parameter     :: C = 1.25 ! C from Colella08.
    real(real64)                :: dp1ext(-1:ncell+2) ! Width of cells
    ! for original grid including ghost cells.
    real(real64)                :: grid1(-2:ncell+2) ! Reconstruction of the
    ! original grid.
    integer(int32)              :: K(interpcnt,0:ncell) ! Indices for
    ! extrapolation and interpolation.
    real(real64)                :: totmassfor(-2:ncell+2) ! Cumulative mass
    ! function on original grid, from left to right.
    real(real64)                :: totmassrev(-2:ncell+2) ! Cumulative mass
    ! function on original grid, from right to left. (The negative of.)
    real(real64)                :: avgdens(-1:ncell+2) ! Average density of
    ! all cells, including ghost cells.
    real(real64)                :: edgevals(0:ncell) ! Interpolated interface
    ! values.
    integer(int32)              :: limreq(0:ncell) ! Array to store if limiting
    ! is required (1) or not (0).
    real(real64)                :: parabvals(3, ncell) ! Parabola parameters
    ! for each cell.
    real(real64)                :: grid2(0:ncell) ! New grid interfaces
    real(real64)                :: ccells2(ncell) ! New grid cell centers
    integer(int32)              :: cellinsecs(2, ncell) ! Where new grid cells
    ! intersect old grid cells.
    real(real64)                :: qdp2(ncell) ! Mass in each new cell.
    real(real64)                :: temps_dp(4) ! Working real values
    integer(int32)              :: temps_int(1)
    integer(int32)              :: ii, jj, kk, kidx ! Counters for DO loops


    ! Get original cell widths including ghost cells.
    temps_dp(1) = minval(dp1) ! Hold smallest cell width.
    dp1ext(-1) = temps_dp(1)
    dp1ext(0)  = temps_dp(1)
    do ii = 1, ncell
       dp1ext(ii) = dp1(ii)
    end do
    dp1ext(ncell+1) = temps_dp(1)
    dp1ext(ncell+2) = temps_dp(1)

    ! Reconstruct the original grid (i.e., list of interface coordinates)
    ! including for ghost cells.
    grid1(-2) = -2.0_real64 * temps_dp(1) ! temps_dp(1) still holds smallest
    ! cell width
    grid1(-1) = -1.0_real64 * temps_dp(1)
    grid1(0)  = 0.0_real64
    do ii = 1, ncell+2
       grid1(ii) = grid1(ii-1) + dp1ext(ii)
    end do

    ! Set the arrays for hold the indices of the points to use for
    ! interpolation and extrapolation
    do ii = 0, interpcnt/2_int32
       K(:, ii) = [(kk, kk = 0, interpcnt-1)]
    end do
    do ii = interpcnt/2_int32 + 1, ncell - (interpcnt/2_int32)
       K(:, ii) = [(kk, kk = ii - interpcnt/2_int32, ii + interpcnt/2_int32)]
    end do
    do ii = ncell-(interpcnt/2_int32) + 1, ncell
       K(:, ii) = [(kk, kk = ncell-(interpcnt-1), ncell)]
    end do

    ! Extrapolate the total mass function to the ghost cell boundaries.
    totmassfor = 0.0_real64
    do ii = 1, ncell
       totmassfor(ii) = totmassfor(ii-1) + qdp(ii)
    end do
    totmassfor(-2) = extrap_totmass(ncell, interpcnt, totmassfor(0:ncell), &
         & grid1, grid1(-2), K(:, 0))
    totmassfor(-1) = extrap_totmass(ncell, interpcnt, totmassfor(0:ncell), &
         & grid1, grid1(-1), K(:, 0))
    totmassfor(ncell+1) = extrap_totmass(ncell, interpcnt, totmassfor(0:ncell), &
         & grid1, grid1(ncell+1), K(:, ncell))
    totmassfor(ncell+2) = extrap_totmass(ncell, interpcnt, totmassfor(0:ncell), &
         & grid1, grid1(ncell+2), K(:, ncell))

    totmassrev = 0.0_real64
    do ii = ncell-1, 0, -1
       totmassrev(ii) = totmassrev(ii+1) + qdp(ii+1)
    end do
    totmassrev(-2) = extrap_totmass(ncell, interpcnt, totmassrev(0:ncell), &
         & grid1, grid1(-2), K(:, 0))
    totmassrev(-1) = extrap_totmass(ncell, interpcnt, totmassrev(0:ncell), &
         & grid1, grid1(-1), K(:, 0))
    totmassrev(ncell+1) = extrap_totmass(ncell, interpcnt, totmassrev(0:ncell), &
         & grid1, grid1(ncell+1), K(:, ncell))
    totmassrev(ncell+2) = extrap_totmass(ncell, interpcnt, totmassrev(0:ncell), &
         & grid1, grid1(ncell+2), K(:, ncell))

    ! Get the cell averages. Use reverse cumulative mass for left ghost cells,
    ! forward for all of the rest.
    do ii = -1, 0
       avgdens(ii) = -(totmassrev(ii) - totmassrev(ii-1)) / dp1ext(ii)
    end do
    do ii = 1, ncell
       avgdens(ii) = qdp(ii) / dp1ext(ii)
    end do
    do ii = ncell+1, ncell+2
       avgdens(ii) = (totmassfor(ii) - totmassfor(ii-1)) / dp1ext(ii)
    end do

    ! Get the primary estimates for interface values of real ghost cells.
    ! But first get new Ks
    edgevals = 0.0_real64
    do ii = 0, ncell
       edgevals(ii) = get_edgeval(ncell, interpcnt, K(:,ii), grid1, &
            & totmassfor, totmassrev, ii)
    end do

    ! Now we figure out which edge value estimates need to be limited
    ! Limiter on interior
    do ii = 0, ncell
       if ((avgdens(ii+1) - edgevals(ii))*(edgevals(ii) - avgdens(ii)) &
            & .gt. 0) then
          ! Edge value is between
          limreq(ii) = 0_int32
       else ! Edge value estimate is outside
          limreq(ii) = 1_int32
       end if
    end do

    ! /Interpolating face values./
    ! Use a modified version of the limiter from Colella08.
    ! It's fine on the interior and we don't limit on the boundary, but
    ! near the boundary we don't have quite enough information.
    do ii = 0, ncell
       if (limreq(ii) .eq. 1) then ! Limiting is required, do calculations
          edgevals(ii) = correct_edgeval(ncell, dp1ext, avgdens, &
               & edgevals(ii), C, ii)
       end if
    end do

    ! /Constructing the parabolic interpolant./
    ! First we fill in the preliminary parabolic interpolant parameters for
    ! each cell. In the same loop, we correct the values appropriately..

    do ii = 1, ncell
       parabvals(1, ii) = edgevals(ii-1) ! ai,-
       parabvals(2, ii) = edgevals(ii)   ! ai,+
       parabvals(3, ii) = 6.0_real64 * avgdens(ii) - 3.0_real64 * (parabvals(1, ii) + parabvals(2, ii))
       call correct_parabvals(ncell, dp1ext, avgdens, parabvals(:, ii), C, ii)
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
          if (parab_piece_local_bnd_preserve_check(ncell, parabvals(:, jj), avgdens(0:ncell), jj) .eq. 0_int32) then ! If not
             ! bounds-preserving, trip flag.
             temps_int(1) = 1_int32
             !print *, '  ~~ !WARNING! Parabolic piece is not locally bounds-preserving:', jj
          end if
       end do

       if (temps_int(1) .eq. 1_int32) then
          print *, '  ~~ !WARNING! At least one parabolic piece is not locally bounds-preserving.'
       end if
    end if

    ! Check that each parabolic piece is globally bounds-preserving
    if (verbosity .eq. 1_int32) then
       temps_int(1) = 0_int32
       do jj = 1, ncell
          if (parab_piece_global_bnd_preserve_check(ncell, parabvals(:, jj), jj) .eq. 0_int32) then ! If not
             ! bounds-preserving, trip flag.
             temps_int(1) = 1_int32
             !print *, '  ~~ !WARNING! Parabolic piece is not globally bounds-preserving:', jj
          end if
       end do

       if (temps_int(1) .eq. 1_int32) then
          print *, '  ~~ !WARNING! At least one parabolic piece is not globally bounds-preserving.'
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

  end subroutine new_remap_2_sbr

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Extrapolate the total mass function.

  function extrap_totmass(ncell, interpcnt, totmass, grid1, xcoord, K) &
       & result(res)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells.
    integer(int32), intent(in) :: interpcnt ! Number of interpolation points.
    real(real64), intent(in)   :: totmass(0:ncell) ! Total mass on the
    ! non-ghost interfaces.
    real(real64), intent(in)   :: grid1(-2:ncell+2) ! Original grid interfaces,
    ! including ghost interfaces.
    real(real64), intent(in)   :: xcoord ! Coordinate we are extrapolating to.
    integer(int32), intent(in) :: K(interpcnt) ! Indices to use for
    ! extrapolation.
    real(real64)               :: res
    real(real64)               :: temps_dp(2) ! Working real values
    integer(int32)             :: kk, kidx, ll, lidx

    res = 0.0_real64
    do kk = 1, interpcnt
       kidx = K(kk)
       temps_dp(1) = totmass(kidx)
       temps_dp(2) = 1.0_real64
       do ll = 1, interpcnt
          lidx = K(ll)
          if (lidx .ne. kidx) then
             temps_dp(2) = temps_dp(2) &
                  & * (xcoord - grid1(lidx)) / (grid1(kidx) - grid1(lidx))
          end if
       end do
       res = res + temps_dp(1) * temps_dp(2)
    end do

    ! Limit the value to the maximum/minimum bounds.
    if (res .gt. 1.0_real64) then
       res = 1.0_real64
    else if (res .lt. 0.0_real64) then
       res = 0.0_real64
    end if

  end function extrap_totmass
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Gets the edge value approximation.
  ! In particular, gets the derivative of a interpolated polynomial for
  ! the cumulative mass function.

  function get_edgeval(ncell, interpcnt, K, grid1, totmassfor, &
       & totmassrev, cell) &
       & result(edgeval)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells.
    integer(int32), intent(in) :: interpcnt ! Order of polynomial
    ! interpolation.
    integer(int32), intent(in) :: K(interpcnt) ! Points for polynomial
    ! intrpolation.
    real(real64), intent(in) :: grid1(-2:ncell+2) ! Grid points of
    ! original grid.
    real(real64), intent(in) :: totmassfor(-2:ncell+2) ! Cumulative mass
    ! function (forward).
    real(real64), intent(in) :: totmassrev(-2:ncell+2) ! Cumulative mass
    ! function (reverse).
    integer(int32), intent(in) :: cell ! Cell we are working with.
    real(real64) :: edgeval ! Approximated edge value.
    real(real64) :: forweight ! Weight for the forward approximation.
    integer(int32) :: kk, kidx, ll, lidx ! For iteration.
    real(real64) :: temps_dp(1) ! Working variables.

    temps_dp = 0.0_real64
    edgeval = 0.0_real64
    forweight = real(cell, real64) / real(ncell, real64)
    forweight = grid1(cell) / grid1(ncell)
    do kk = 1, interpcnt
       kidx = K(kk)
       if (kidx .ne. cell) then
          temps_dp(1) = 1.0_real64
          do ll = 1, interpcnt
             lidx = K(ll)
             if ((lidx .ne. cell) &
                  & .and. (lidx .ne. kidx) ) then
                temps_dp(1) = temps_dp(1) * (grid1(cell) - grid1(lidx)) &
                     & / (grid1(kidx) - grid1(lidx))
             end if
          end do
          edgeval = edgeval &
               & + (forweight * totmassfor(kidx) &
               &    - (1.0_real64 - forweight) * totmassrev(kidx)) &
               &   / (grid1(kidx) - grid1(cell)) * temps_dp(1)

       else if (kidx .eq. cell) then
          temps_dp(1) = 0.0_real64
          do ll = 1, interpcnt
             lidx = K(ll)
             if (lidx .ne. cell) then
                temps_dp(1) = temps_dp(1) &
                     & + (1.0_real64 / (grid1(cell) - grid1(lidx)))
             end if
          end do
          edgeval = edgeval &
               & + (forweight * totmassfor(kidx) &
               &    - (1.0_real64 - forweight) * totmassrev(kidx)) &
               &   * temps_dp(1)
       end if
    end do

  end function get_edgeval
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! First phase of the limiter: enforcing monotnicity of interface values.
  function correct_edgeval(ncell, dp1, avgdens, edgeval, C, cell) &
       & result(corr_edgeval)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells in the grid
    real(real64), intent(in)   :: dp1(-1:ncell+2) ! Width of cells on original
    ! grid, including ghost cells.
    real(real64), intent(in)   :: avgdens(-1:ncell+2) ! Average density of
    ! original cells.
    real(real64), intent(in)   :: edgeval ! Original edge value
    real(real64), intent(in)   :: C ! Scale factor for second derivative
    ! approximations.
    integer(int32), intent(in) :: cell ! Current cell we are correcting the
    ! edge value for.
    real(real64)               :: corr_edgeval ! Corrected edge value
    real(real64)               :: temps_dp(4) ! Holds working real variables

    temps_dp = 0.0_real64

    ! Calculate second derivative approximations.
    !! Centered
    temps_dp(4) = (1.0_real64 / 6.0_real64) &
         & * (dp1(cell+1) + dp1(cell))
    temps_dp(1) = 8.0_real64 * (temps_dp(4))**2 * &
         & ( avgdens(cell) / (dp1(cell) * (dp1(cell+1) + dp1(cell))) &
         &   - edgeval / (dp1(cell) * dp1(cell+1)) &
         &   + avgdens(cell+1) / (dp1(cell+1) * (dp1(cell+1) + dp1(cell))) )

    !! Left
    temps_dp(4) = (1.0_real64 / 6.0_real64) &
         & * (dp1(cell+1) + 3.0_real64 * dp1(cell) + dp1(cell-1))
    temps_dp(2) = 8.0_real64 * (temps_dp(4))**2 * &
         & ( avgdens(cell-1) &
         &     / (dp1(cell-1) * (dp1(cell+1) + dp1(cell) + dp1(cell-1))) &
         &   - avgdens(cell) &
         &     / (dp1(cell-1) * (dp1(cell+1) + dp1(cell))) &
         &   + avgdens(cell+1) &
         &     / ((dp1(cell+1) + dp1(cell)) &
         &         * (dp1(cell+1) + dp1(cell) + dp1(cell-1))) )

    !! Right
    temps_dp(4) = (1.0_real64 / 6.0_real64) &
         & * (dp1(cell+2) + 3.0_real64 * dp1(cell+1) + dp1(cell))
    temps_dp(3) = 8.0_real64 * (temps_dp(4))**2 * &
         & ( avgdens(cell) &
         &     / ((dp1(cell+1) + dp1(cell)) &
         &         * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) &
         &   - avgdens(cell+1) &
         &     / ((dp1(cell+1) + dp1(cell)) * (dp1(cell+2) + dp1(cell+1))) &
         &   + avgdens(cell+2) &
         &     / ((dp1(cell+2) + dp1(cell+1)) &
         &         * (dp1(cell+2) + 2.0_real64 * dp1(cell+1) + dp1(cell))) )

    ! Check if approxs have same sign
    if ((temps_dp(1) * temps_dp(2) .gt. 0) &
         & .and. (temps_dp(2) * temps_dp(3) .gt. 0)) then
       temps_dp(1) = sign(min(abs(temps_dp(1)), C * abs(temps_dp(2)), &
            &                 C * abs(temps_dp(3))), temps_dp(1))
    else
       temps_dp(1) = 0.0_real64
    end if

    ! Set temps_dp(2) to grid spacing for correction

    corr_edgeval = (avgdens(cell) + avgdens(cell+1))/2.0_real64 &
         & + (1.0_real64/3.0_real64) * temps_dp(1)


  end function correct_edgeval
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine correct_parabvals(ncell, dp1, avgdens, parabvals, C, cell)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell ! Number of cells in the grid
    real(real64), intent(in) :: dp1(-1:ncell+2) ! Width of cells on original
    ! grid.
    real(real64), intent(in) :: avgdens(-1:ncell+2) ! Average density of
    ! original cells.
    real(real64), intent(inout) :: parabvals(3) ! Original parabola values
    real(real64), intent(in) :: C ! Scale factor for second derivative
    ! approximations.
    integer(int32), intent(in) :: cell ! Current cell we are correcting the
    ! edge value for.
    real(real64) :: corr_edgeval ! Corrected edge value
    real(real64) :: temps_dp(4) ! Holds working real variables
    integer(int32) :: temps_int(1)


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


  end subroutine correct_parabvals

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

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Checks that a parabolic piece is globally bounds-preserving.
  function parab_piece_global_bnd_preserve_check(ncell, parabvals, cell) result(res)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in) :: ncell
    real(real64), intent(in) :: parabvals(3)
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
    lo_bnd = 0.0_real64
    up_bnd = 1.0_real64

    ! Compare max, min values
    if ((min_parab .le. lo_bnd) .or. (max_parab .ge. up_bnd)) then
       res = 0_int32
    else
       res = 1_int32
    end if

  end function parab_piece_global_bnd_preserve_check

end submodule vertremap_redux_2
