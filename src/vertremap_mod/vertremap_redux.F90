!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! New vertical remapping code, should take in generally the same arguments
! as the original, just simplified.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

submodule (vertremap_mod) vertremap_redux

contains

  module subroutine new_remap_sbr(qdp, ncell, dp1, dp2, remap_alg)

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32), intent(in)  :: ncell ! Number of cells in the grid
    real(real64), intent(inout) :: qdp(ncell) ! Mass of each cell
    real(real64), intent(in)    :: dp1(ncell), dp2(ncell) ! Width of cells for
    ! original, new grids
    integer(int32), intent(in)  :: remap_alg ! Algorithm flag to use for
    ! remapping
    real(real64)                :: ogrid(0:ncell) ! Reconstruction of the
    ! original grid.
    real(real64)                :: ocells(ncell) ! Reconstruction of the
    ! cell centers.
    real(real64)                :: totmass(0:ncell) ! Cumulative mass array
    integer(int32)              :: K(5, 0:ncell) ! Indices for each
    ! interpolation polynomial.
    real(real64)                :: edgevals(0:ncell) ! Estimated edge values
    integer(int32)              :: limreq(0:ncell) ! Array for holding whether
    ! an edge estimate needs to be limited.
    real(real64)                :: temps_dp(3) ! Temporary real values used in
    ! calculations.
    integer(int32)              :: jj, kk, kidx, mm, midx, ll, lidx
    ! Counters for do loops


    ! Reconstruct the original grid from 0 to end
    ogrid(0) = 0.0_real64
    do jj = 1, ncell
       ogrid(jj) = ogrid(jj-1) + dp1(jj)
    end do
    ! Get the cell centers explicitly
    do jj = 1, ncell
       ocells(jj) = (ogrid(jj) + ogrid(jj-1))/2.0_real64
    end do
    
    ! Get the cumulative mass function
    totmass(0) = 0
    do jj = 1, ncell
       totmass(jj) = totmass(jj-1) + qdp(jj)*dp1(jj)
    end do

    ! Set the arrays for hold the indices of the points to use for
    ! interpolation
    K(:, 0) = [0, 1, 2, 3, 4]
    K(:, 1) = [0, 1, 2, 3, 4]
    do jj = 2, ncell-2
       K(:, jj) = [jj-2, jj-1, jj, jj+1, jj+2]
    end do
    K(:, ncell-1) = [ncell-4, ncell-3, ncell-2, ncell-1, ncell]
    K(:, ncell  ) = [ncell-4, ncell-3, ncell-2, ncell-1, ncell]

    ! Get the primary estimates for edge values. It's a very messy loop
    ! for Eqn. 2.2.6
    edgevals = 0.0_real64
    do jj = 0, ncell
       temps_dp(1) = 0.0_real64
       do kk = 1, 5
          kidx = K(kk, jj)
          temps_dp(2) = 0.0_real64
          do mm = 1, 5
             midx = K(mm, jj)
             temps_dp(3) = 1.0_real64
             do ll = 1, 5
                lidx = K(ll, jj)
                if ((lidx .ne. midx) .and. (lidx .ne. kidx)) then
                   temps_dp(3) = temps_dp(3) &
                        & * (ogrid(jj) - ogrid(lidx))/(ogrid(kidx) - ogrid(lidx))
                end if
             end do
             temps_dp(2) = temps_dp(2) + temps_dp(3)
          end do
          temps_dp(1) = temps_dp(1) + temps_dp(2)*totmass(kidx)
       end do
       edgevals(jj) = temps_dp(1)
    end do

    ! Now we figure out which edge value estimates need to be limited
    limreq(0) = 0 ! Don't limit at the boundaries
    do jj = 1, ncell-1
       if ((qdp(jj+1) - edgevals(jj))*(edgevals(jj) - qdp(jj)) .gt. 0) then
          ! Edge value is between
          limreq(jj) = 0
       else ! Edge value estimate is outside
          limreq(jj) = 1
       end if
    end do
    limreq(ncell) = 0 ! Don't limit at the boundaries

    ! Use a modified version of the limiter from Colella08.
    ! It's fine on the interior and we don't limit on the boundary, but
    ! near the boundary we don't have quite enough information.
    do jj = 0, ncell
       if (limreq(jj) .eq. 1) then ! Limiting is required, do calculations
          temps_dp = 0.0_real64 ! Here we use this array for the estimates
          ! of the second derivative a-la Colella08 Eqn. 18
          ! Center approximation is good for all cases.
          temps_dp(2) = &
               & qdp(jj) / ((ocells(jj) - ogrid(jj))*(ocells(jj) - ocells(jj+1))) &
               & + edgevals(jj) / ((ocells(jj) - ogrid(jj))*(ocells(jj+1) - ogrid(jj))) &
               & + qdp(jj+1) / ((ocells(jj+1) - ogrid(jj))*(ocells(jj+1) - ocells(jj)))
          if (jj .ne. 1) then 
             ! Not on the left edge region, so left approximation is okay.
             temps_dp(1) = &
                  & qdp(jj-1) / ((ocells(jj-1) - ocells(jj))*(ocells(jj-1) - ocells(jj+1))) &
                  & - qdp(jj) / ((ocells(jj-1) - ocells(jj))*(ocells(jj) - ocells(jj+1))) &
                  & + qdp(jj+1) / ((ocells(jj-1) - ocells(jj+1))*(ocells(jj) - ocells(jj+1)))
          end if
          if (jj .ne. ncell-1) then 
             ! Not on the right edge region, so right approximation is okay.
             temps_dp(3) = &
                  & qdp(jj) / ((ocells(jj) - ocells(jj+1))*(ocells(jj) - ocells(jj+2))) &
                  & - qdp(jj+1) / ((ocells(jj) - ocells(jj+1))*(ocells(jj+1) - ocells(jj+2))) &
                  & + qdp(jj+2) / ((ocells(jj) - ocells(jj+2))*(ocells(jj+1) - ocells(jj+2)))
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
          temps_dp(1) = ((ocells(jj) - ogrid(jj))*(ocells(jj) - ocells(jj+1))) &
               & + ((ocells(jj) - ogrid(jj))*(ocells(jj+1) - ogrid(jj))) &
               & + ((ocells(jj+1) - ogrid(jj))*(ocells(jj+1) - ocells(jj)))
          edgevals(jj) = (qdp(jj) + qdp(jj+1))/2.0_real64 - temps_dp(1) * temps_dp(2)
       end if
    end do

    print *, '~~ ', qdp
    print *, '~~ ', edgevals
    
    print *, '~~ Test complete!'
    stop

    
  end subroutine new_remap_sbr

end submodule vertremap_redux
