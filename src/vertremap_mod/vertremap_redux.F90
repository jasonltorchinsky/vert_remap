!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! New vertical remapping code, should take in generally the same arguments
! as the original, just simplified.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (vertremap_mod) vertremap_redux

contains

  module subroutine new_remap_sbr(qdp, ncell, dp1, dp2, remap_alg)

    use iso_fortran_env, only: int32, real32, real64

    implicit none

    integer(int32), intent(in)  :: ncell ! Number of cells in the grid
    real(real64), intent(inout) :: qdp(ncell) ! Mass of each cell in original grid
    real(real64), intent(in)    :: dp1(ncell), dp2(ncell) ! Width of cells for
    ! original, new grids
    integer(int32), intent(in)  :: remap_alg ! Algorithm flag to use for
    ! remapping
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
       A = 0.0_real64
       x = 0.0_real64
       b = 0.0_real64
       do kk = 1, interpord ! We can either access the entry of K once...
          kidx = K(kk, jj)
          temps_dp(1) = grid1(kidx)
          do mm = 1, interpord ! or go through this array contiguously.
             ! I'm not sure which saves more time.
             A(kk, mm) = temps_dp(1)**(interpord - mm)
          end do
          b(kk) = totmass(kidx)
       end do
       call dsgesv(interpord, 1, A, interpord, ipiv, b, interpord, x, &
            & interpord, work, swork, temps_int(1), temps_int(2))
       do kk = 1, interpord-1
          edgevals(jj) = edgevals(jj) + (interpord-kk)*x(kk)*(grid1(jj)**(interpord-kk-1)) 
       end do 
    end do

    ! Now we figure out which edge value estimates need to be limited
    limreq(0) = 0 ! Don't limit at the boundaries
    do jj = 1, ncell-1
       if ((avgdens1(jj+1) - edgevals(jj))*(edgevals(jj) - avgdens1(jj)) .gt. 0) then
          ! Edge value is between
          limreq(jj) = 0
       else ! Edge value estimate is outside
          limreq(jj) = 1
       end if
    end do
    limreq(ncell) = 0 ! Don't limit at the boundaries
    limreq = 0

    ! /Interpolating face values./
    ! Use a modified version of the limiter from Colella08.
    ! It's fine on the interior and we don't limit on the boundary, but
    ! near the boundary we don't have quite enough information.
    do jj = 0, ncell
       if (limreq(jj) .eq. 1) then ! Limiting is required, do calculations
          temps_dp = 0.0_real64 ! Here we use this array for the estimates
          ! of the second derivative a-la Colella08 Eqn. 18
          ! Center approximation is good for all cases.
          temps_dp(2) = &
               & avgdens1(jj) / ((ccells1(jj) - grid1(jj))*(ccells1(jj) - ccells1(jj+1))) &
               & + edgevals(jj) / ((ccells1(jj) - grid1(jj))*(ccells1(jj+1) - grid1(jj))) &
               & + avgdens1(jj+1) / ((ccells1(jj+1) - grid1(jj))*(ccells1(jj+1) - ccells1(jj)))
          if (jj .ne. 1) then 
             ! Not on the left edge region, so left approximation is okay.
             temps_dp(1) = 2.0_real64 * &
                  & (- avgdens1(jj-1) / ((ccells1(jj) - ccells1(jj-1))*(ccells1(jj-1) - ccells1(jj+1))) &
                  &  + avgdens1(jj) / ((ccells1(jj) - ccells1(jj-1))*(ccells1(jj) - ccells1(jj+1))) &
                  &  + avgdens1(jj+1) / ((ccells1(jj) - ccells1(jj+1))*(ccells1(jj-1) - ccells1(jj+1))))
          end if
          if (jj .ne. ncell-1) then 
             ! Not on the right edge region, so right approximation is okay.
             temps_dp(3) = 2.0_real64 * &
                  & (- avgdens1(jj) / ((ccells1(jj+1) - ccells1(jj))*(ccells1(jj) - ccells1(jj+2))) &
                  &  + avgdens1(jj+1) / ((ccells1(jj+1) - ccells1(jj))*(ccells1(jj+1) - ccells1(jj+2))) &
                  &  + avgdens1(jj+2) / ((ccells1(jj+1) - ccells1(jj+2))*(ccells1(jj) - ccells1(jj+2))))
          end if

          ! Check if all of the approximations have the same sign or not
          ! Use temps_dp(2) to hold the limiting approximation
          if (jj .eq. 1) then
             ! On left edge, only have two approximations, can multiply to see if signs agree.
             if (temps_dp(2) * temps_dp(3) .gt. 0.0_real64) then
                ! We set C = 1.25 like Colella08
                temps_dp(2) = sign(min(abs(temps_dp(2)), 1.25_real64 * abs(temps_dp(3))), temps_dp(2))
             else
                temps_dp(2) = 0.0_real64
             end if
          else if (jj .eq. ncell-1) then
             ! On right edge, only have two approximations, can multiply to see if the signs agree.
             if (temps_dp(1) * temps_dp(2) .gt. 0.0_real64) then
                !We set C = 1.25 like Colella08
                temps_dp(2) = sign(min(1.25_real64 * abs(temps_dp(1)), abs(temps_dp(2))), temps_dp(2))
             else
                temps_dp(2) = 0.0_real64
             end if
          else
             ! On the interior, so we must check two products
             if ((temps_dp(1) * temps_dp(2) .gt. 0.0_real64) &
                  & .and. (temps_dp(2) * temps_dp(3) .gt. 0.0_real64)) then
                !We set C = 1.25 like Colella08
                temps_dp(2) = sign(min(1.25_real64 * abs(temps_dp(1)), &
                     & abs(temps_dp(2)), 1.25_real64 * abs(temps_dp(3))), temps_dp(2))
             else
                temps_dp(2) = 0.0_real64
             end if
          end if

          ! Adjust the edge value, use analogous formula to Colella08
          ! Use temps_dp(1) to hold the coefficient
          !! TO DO: Change this to reflect Eqn. 2.2.14 of the write-up.
          temps_dp(1) = ((ccells1(jj) - grid1(jj))*(ccells1(jj) - ccells1(jj+1))) &
               & + ((ccells1(jj) - grid1(jj))*(ccells1(jj+1) - grid1(jj))) &
               & + ((ccells1(jj+1) - grid1(jj))*(ccells1(jj+1) - ccells1(jj)))
          edgevals(jj) = (avgdens1(jj) + avgdens1(jj+1))/2.0_real64 &
               & - temps_dp(1) * temps_dp(2)
       end if
    end do

    ! /Constructing the parabolic interpolant./
    ! First we fill in the preliminary parabolic interpolant parameters for
    ! each cell. In the same loop, we correct the values appropriately.
    
    do jj = 1, ncell
       ! Reset temps_dp for new use.
       temps_dp = 0.0_real64
       parabvals(1, jj) = edgevals(jj-1) ! aj,-
       parabvals(2, jj) = edgevals(jj)   ! aj,+
       parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
            & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j
       ! The boundary cells must be handled differently than Colella08.
       if (jj .eq. 1) then ! Leftmost cell
          
          if ((parabvals(2, jj) - avgdens1(jj)) * (avgdens1(jj) - parabvals(1, jj)) &
               & .le. 0.0_real64) then ! We are at a local extremum
             ! Get second derivative estimates, can only do the centered
             ! and right-biased. Store in temps_dp(1), temps_dp(4).
             temps_dp(1) = &
                  & parabvals(1, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj-1) - ccells1(jj))) &
                  & + avgdens1(jj) / ((grid1(jj-1) - ccells1(jj))*(grid1(jj) - ccells1(jj))) &
                  & + parabvals(2, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj) - ccells1(jj)))
             temps_dp(4) = &
                  & avgdens1(jj) / ((ccells1(jj) - ccells1(jj+1))*(ccells1(jj) - ccells1(jj+2))) &
                  & - avgdens1(jj+1) / ((ccells1(jj) - ccells1(jj+1))*(ccells1(jj+1) - ccells1(jj+2))) &
                  & + avgdens1(jj+2) / ((ccells1(jj) - ccells1(jj+2))*(ccells1(jj+1) - ccells1(jj+2)))

             if (temps_dp(1) * temps_dp(4) .gt. 0.0_real64) then
                ! If they have the same sign, the limiting value is the smallest one. We again
                ! use C = 1.25 and store the limiting value in temps_dp(2).
                temps_dp(2) = sign(min(abs(temps_dp(1)), 1.25_real64 * abs(temps_dp(4))), &
                     & temps_dp(1))
             else
                temps_dp(2) = 0.0_real64
             end if

             ! Correct the parabola values
             if (temps_dp(1) .ne. 0.0_real64) then
                parabvals(1, jj) = avgdens1(jj) &
                     & + (parabvals(1, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
                parabvals(2, jj) = avgdens1(jj) &
                     & + (parabvals(2, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
             else
                parabvals(1, jj) = avgdens1(jj)
                parabvals(2, jj) = avgdens1(jj)
             end if
             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j

          else ! We are not at local extremum
             ! Use temps_dp(1), temps_dp(2) to hold alpha_j,+, alpha_j,-
             temps_dp(1) = parabvals(1, jj) - avgdens1(jj)
             temps_dp(2) = parabvals(2, jj) - avgdens1(jj)
             ! Since we are at the leftmost cell, we can't calculcate deltaa
             ! For the alpha_j,- case, so we skip it
             if (abs(temps_dp(1)) .ge. 2.0_real64 * abs(temps_dp(2))) then ! alpha_j,+ case
                ! Use temps_dp(3) for deltaIext, temps_dp(4) for deltaa
                temps_dp(3) = -temps_dp(1)**2 &
                     & / (4.0_real64 * (temps_dp(1) + temps_dp(2)))
                temps_dp(4) = avgdens1(jj+1) - avgdens1(jj)
                ! Get the sign of the cell-average difference
                if (avgdens1(jj+1) .ge. avgdens1(jj)) then ! No cell on the left to compare with
                   temps_int(1) = 1_int32
                else
                   temps_int(1) = -1_int32
                end if
                
                if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                   parabvals(1, jj) = avgdens1(jj) &
                        & - (2.0_real64 * temps_dp(4) &
                        &    + 2.0_real64 * temps_int(1) &
                        &      * sqrt(temps_dp(4)**2 - temps_dp(4) * temps_dp(1)))
                end if
             end if

             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j
             
          end if

       else if (jj .eq. ncell) then ! Rightmost cell

          if ((parabvals(2, jj) - avgdens1(jj)) * (avgdens1(jj) - parabvals(1, jj)) &
               & .le. 0.0_real64) then ! We are at a local extremum
             ! Get second derivative estimates, can only do the centered
             ! and left-biased. Store in temps_dp(1), temps_dp(2).
             temps_dp(1) = &
                  & parabvals(1, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj-1) - ccells1(jj))) &
                  & + avgdens1(jj) / ((grid1(jj-1) - ccells1(jj))*(grid1(jj) - ccells1(jj))) &
                  & + parabvals(2, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj) - ccells1(jj)))
             temps_dp(2) = &
                  & avgdens1(jj-2) / ((ccells1(jj-2) - ccells1(jj-1))*(ccells1(jj-2) - ccells1(jj))) &
                  & - avgdens1(jj-1) / ((ccells1(jj-2) - ccells1(jj-1))*(ccells1(jj-1) - ccells1(jj))) &
                  & + avgdens1(jj) / ((ccells1(jj-2) - ccells1(jj))*(ccells1(jj-1) - ccells1(jj)))

             if (temps_dp(1) * temps_dp(2) .gt. 0.0_real64) then
                ! If they have the same sign, the limiting value is the smallest one. We again
                ! use C = 1.25 and store the limiting value in temps_dp(2).
                temps_dp(2) = sign(min(abs(temps_dp(1)), 1.25_real64 * abs(temps_dp(2))), &
                     & temps_dp(1))
             else
                temps_dp(2) = 0.0_real64
             end if

             ! Correct the parabola values
             if (temps_dp(1) .ne. 0.0_real64) then
                parabvals(1, jj) = avgdens1(jj) &
                     & + (parabvals(1, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
                parabvals(2, jj) = avgdens1(jj) &
                     & + (parabvals(2, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
             else
                parabvals(1, jj) = avgdens1(jj)
                parabvals(2, jj) = avgdens1(jj)
             end if
             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j

          else ! We are not at local extremum
             ! Use temps_dp(1), temps_dp(2) to hold alpha_j,+, alpha_j,-
             temps_dp(1) = parabvals(1, jj) - avgdens1(jj)
             temps_dp(2) = parabvals(2, jj) - avgdens1(jj)
             ! Since we are at the rightmost cell, we can't calculcate deltaa
             ! For the alpha_j,+ case, so we skip it
             if (abs(temps_dp(2)) .ge. 2.0_real64 * abs(temps_dp(1))) then ! alpha_j,- case
                ! Use temps_dp(3) for deltaIext, temps_dp(4) for deltaa
                temps_dp(3) = -temps_dp(2)**2 &
                     & / (4.0_real64 * (temps_dp(1) + temps_dp(2)))
                temps_dp(4) = avgdens1(jj-1) - avgdens1(jj)
                ! Get the sign of the cell-average difference
                if (avgdens1(jj) .ge. avgdens1(jj-1)) then ! No cell on the right to compare with.
                   temps_int(1) = 1_int32
                else
                   temps_int(1) = -1_int32
                end if
                
                if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                   parabvals(1, jj) = avgdens1(jj) &
                        & - (2.0_real64 * temps_dp(4) &
                        &    + 2.0_real64 * temps_int(1) &
                        &      * sqrt(temps_dp(4)**2 - temps_dp(4) * temps_dp(2)))
                end if
             end if

             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j
             
          end if
          
       else ! Interior cell, regular Colella08 works

          if (((parabvals(2, jj) - avgdens1(jj)) * (avgdens1(jj) - parabvals(1, jj)) .le. 0.0_real64) &
               & .or. (((avgdens1(jj-1) - avgdens1(jj)) * (avgdens1(jj) - avgdens1(jj+1)) .le. 0.0_real64))) then ! We are at a local extremum
             ! Get second derivative estimates, can get all four (j, j,L, j,C, j,R)
             temps_dp(1) = &
                  & parabvals(1, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj-1) - ccells1(jj))) &
                  & + avgdens1(jj) / ((grid1(jj-1) - ccells1(jj))*(grid1(jj) - ccells1(jj))) &
                  & + parabvals(2, jj) / ((grid1(jj-1) - grid1(jj))*(grid1(jj) - ccells1(jj)))
             temps_dp(2) = &
                  & avgdens1(jj-2) / ((ccells1(jj-2) - ccells1(jj-1))*(ccells1(jj-2) - ccells1(jj))) &
                  & - avgdens1(jj-1) / ((ccells1(jj-2) - ccells1(jj-1))*(ccells1(jj-1) - ccells1(jj))) &
                  & + avgdens1(jj) / ((ccells1(jj-2) - ccells1(jj))*(ccells1(jj-1) - ccells1(jj)))
             temps_dp(3) = &
                  & avgdens1(jj-1) / ((ccells1(jj-1) - ccells1(jj))*(ccells1(jj-1) - ccells1(jj+1))) &
                  & - avgdens1(jj) / ((ccells1(jj-1) - ccells1(jj))*(ccells1(jj) - ccells1(jj+1))) &
                  & + avgdens1(jj+1) / ((ccells1(jj-1) - ccells1(jj+1))*(ccells1(jj) - ccells1(jj+1)))
             temps_dp(4) = &
                  & avgdens1(jj) / ((ccells1(jj) - ccells1(jj+1))*(ccells1(jj) - ccells1(jj+2))) &
                  & - avgdens1(jj+1) / ((ccells1(jj) - ccells1(jj+1))*(ccells1(jj+1) - ccells1(jj+2))) &
                  & + avgdens1(jj+2) / ((ccells1(jj) - ccells1(jj+2))*(ccells1(jj+1) - ccells1(jj+2)))

             if ((temps_dp(1) * temps_dp(2) .gt. 0.0_real64) &
                  & .and. (temps_dp(2) * temps_dp(3) .gt. 0.0_real64) &
                  & .and. (temps_dp(3) * temps_dp(4) .gt. 0.0_real64)) then
                ! If each approximation has the same sign, the limiting value is the smallest one. We again
                ! use C = 1.25 and store the limiting value in temps_dp(2).
                temps_dp(2) = sign(min(abs(temps_dp(1)), 1.25_real64 * abs(temps_dp(2)), &
                     & 1.25_real64 * abs(temps_dp(3)), 1.25_real64 * abs(temps_dp(4))), &
                     & temps_dp(1))
             else
                temps_dp(2) = 0.0_real64
             end if

             ! Correct the parabola values
             if (temps_dp(1) .ne. 0.0_real64) then
                parabvals(1, jj) = avgdens1(jj) &
                     & + (parabvals(1, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
                parabvals(2, jj) = avgdens1(jj) &
                     & + (parabvals(2, jj) - avgdens1(jj)) * (temps_dp(2)/temps_dp(1))
             else
                parabvals(1, jj) = avgdens1(jj)
                parabvals(2, jj) = avgdens1(jj)
             end if
             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j

          else ! We are not at local extremum
             ! Use temps_dp(1), temps_dp(2) to hold alpha_j,+, alpha_j,-
             temps_dp(1) = parabvals(1, jj) - avgdens1(jj)
             temps_dp(2) = parabvals(2, jj) - avgdens1(jj)

             ! Get the sign of the cell-average difference
             if (avgdens1(jj+1) .ge. avgdens1(jj)) then ! No cell on the left to compare with
                temps_int(1) = 1_int32
             else
                temps_int(1) = -1_int32
             end if

             if (abs(temps_dp(1)) .ge. 2.0_real64 * abs(temps_dp(2))) then ! alpha_j,+ case
                ! Use temps_dp(3) for deltaIext, temps_dp(4) for deltaa
                temps_dp(3) = -temps_dp(1)**2 &
                     & / (4.0_real64 * (temps_dp(1) + temps_dp(2)))
                temps_dp(4) = avgdens1(jj+1) - avgdens1(jj)
                
                if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                   parabvals(1, jj) = avgdens1(jj) &
                        & - (2.0_real64 * temps_dp(4) &
                        &    + 2.0_real64 * temps_int(1) &
                        &      * sqrt(temps_dp(4)**2 - temps_dp(4) * temps_dp(1)))
                end if
             end if

             if (abs(temps_dp(2)) .ge. 2.0_real64 * abs(temps_dp(1))) then ! alpha_j,- case
                ! Use temps_dp(3) for deltaIext, temps_dp(4) for deltaa
                temps_dp(3) = -temps_dp(2)**2 &
                     & / (4.0_real64 * (temps_dp(1) + temps_dp(2)))
                temps_dp(4) = avgdens1(jj-1) - avgdens1(jj)
                ! Get the sign of the cell-average difference
                if (avgdens1(jj) .ge. avgdens1(jj-1)) then ! No cell on the right to compare with.
                   temps_int(1) = 1_int32
                else
                   temps_int(1) = -1_int32
                end if
                
                if (temps_int(1) * temps_dp(3) .ge. temps_int(1) * temps_dp(4)) then
                   parabvals(1, jj) = avgdens1(jj) &
                        & - (2.0_real64 * temps_dp(4) &
                        &    + 2.0_real64 * temps_int(1) &
                        &      * sqrt(temps_dp(4)**2 - temps_dp(4) * temps_dp(2)))
                end if
             end if

             parabvals(3, jj) = 6.0_real64 * avgdens1(jj) &
                  & - 3.0_real64 * (parabvals(1, jj) + parabvals(2, jj)) ! a6,j
             
          end if
          
       end if
       
    end do

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

end submodule vertremap_redux
