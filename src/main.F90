!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main driver for the vertical remap code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TO-DO: Add doc string information
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program vert_remap

  use iso_fortran_env, only: int32, real64
  use vertremap_mod

  implicit none

  integer(int32)            :: ncell, nlev ! Number of cells, and levels for 
                                           ! each grid
  real(real64), allocatable :: dp1(:), dp2(:) ! Grid spacings for each grid
  real(real64), allocatable :: grid1(:), grid2(:) ! Grid at cell boundaries
                                                  ! (0, ..., H)
  real(real64), allocatable :: grid1_stg(:), grid2_stg(:) ! Staggered grid at
                                                          ! cell centers 
                                                          ! (dx/2, ..., H-dx/2)
  real(real64), allocatable :: QOrig(:), QNew(:) ! Average density in each cell
  real(real64), allocatable :: QdpOrig(:), QdpNew(:) ! Mass in each cell
  real(real64), allocatable :: QTrue(:) ! True density on the transformed grid
  character(len=32)         :: arg ! Input argument from command line
  integer(int32)            :: ii ! counter for do loops

  ! Set default input values
  ncell = 3
  nlev = 4
  

  ! Read command line arguments
  do ii = 0, command_argument_count()
    call get_command_argument(ii, arg)
    if (len_trim(arg) .eq. 0) exit
  
    select case(arg)
    case('-h', '--help')
      print *, 'vert_remap flags: '
      print *, '  -h, --help: Returns this message.'
      print *, 'vert_remap arguments: '
      print *, '  nlev: Number of grid levels (cell boundaries).'
      print *, '  ncell: Number of cells (regions between levels).'
    case('ncell')
      call get_command_argument(ii + 1, arg)
      read(arg, *), ncell
    case('nlev')
      call get_command_argument(ii + 1, arg)
      read(arg, *) nlev
    end select
  end do

  ! Check for input errors
  if (ncell .eq. 3) ncell = nlev - 1 ! Assume ncell wasn't set first
  if (nlev .eq. 4) nlev = ncell + 1 ! The assume nlev wasn't set
  if (nlev - ncel .neq 1) then ! Too many/few cells for levels.
    print *, '  ISSUE: Inappropriate number of cells and levels...'
    print *, '  Assuming that the number of levels is correct.'
    ncell = nlev - 1
  end if

  ! Set up grids, these are the boundaries of the cells.
  allocate(grid1(nlev))
  allocate(grid2(nlev))
  
  do ii = 1, nlev ! Fill in grids
    grid2(ii) = real(ii - 1, real64)/real(nlev - 1, real64) ! Target, uniform
    grid1(ii) = grid1_func(grid2(ii)) ! Source, non-uniform
  end do ! NOTE: it would be more efficient to split this loop, but since this
         ! is a small 1-D problem, we leave it like this for readability.
         ! We do this for a few loops here.

  ! Get the cell widths.
  allocate(dp1(ncell))
  allocate(dp2(ncell))

  do ii = 1, ncell
    dp1(ii) = grid1(ii+1) - grid1(ii)
    dp2(ii) = grid2(ii+1) - grid2(ii)
  end do

  ! Calculate cell centers.
  allocate(grid1_stg(ncell))
  allocate(grid2_stg(ncell))

  do ii = 1, ncell ! Cell centers are halfway between their boundaries.
    grid1_stg(ii) = (grid1(ii) + grid1(ii+1))/2.0_real64
    grid2_stg(ii) = (grid2(ii) + grid2(ii+1))/2.0_real64
  end do

  ! Calulate density, mass of each cell.
  allocate(QOrig(ncell))
  allocate(QdpOrig(ncell))
  
  do ii = 1, ncell
    QOrig(ii) = Q_func(grid1_stg(ii))
    QdpOrig(ii) = QOrig(ii) * dp1(ii)
  end do

  ! Remap Qdp to grid 1
  allocate(QdpNew(ncell)) ! Remap is in-place, so we make a new array.
  QdpNew = QdpOrig
  call remap1(QdpNew, 1, ncell, 1, dp1, dp2, 10)

  ! Get density from mass
  allocate(QNew(ncell))
  do ii = 1, nlev
    QNew(ii) = QdpNew(ii) / dp2(ii)
  end do

  ! Get true density on the transformed grid
  allocate(QTrue(ncell))
  do ii = 1, ncell
    QTrue(ii) = Q_func(grid2_stg(ii))
  end do

  ! Output information  
100 format (A, f16.12)

  ! Density
  open(unit=1000, file='Q.txt')
  write(1000, 100, advance='no') '{', QNew(1)
  do ii = 2, ncell
    write(1000, 100, advance='no') ',', QNew(ii)
  end do
  write(1000, '(A)') '}'

  write(1000, 100, advance='no') '{', QTrue(1)
  do ii = 2, ncell
    write(1000, 100, advance='no') ',', QTrue(ii)
  end do
  write(1000, '(A)') '}'
  close(1000)

  ! Cell widths
  open(unit=1001, file='dp.txt')
  write(1001, 100, advance='no') '{', dp2(1)
  do ii = 2, ncell
    write(1001, 100, advance='no') ',', dp2(ii)
  end do
  write(1001, '(A)') '}'
  close(1001)

  ! Cell centers/staggered grid
  open(unit=1002, file='grid_stg.txt')
  write(1002, 100, advance='no') '{', grid2_stg(1)
  do ii = 2, ncell
    write(1002, 100, advance='no') ',', grid2_stg(ii)
  end do
  write(1002, '(A)') '}'

  close(1002)

  contains

    ! Function to get the grid points for grid 1 (source) from 
    ! grid 2 (target, uniform)
    function grid1_func(x) result(y)
    
      implicit none
      
      real(real64), intent(in) :: x ! In [0, 1]
      real(real64)             :: y ! Should be in [0, 1]

      y = x**2

    end function grid1_func

    ! Function to get the value of Qdp at a given grid point
    function Q_func(x) result(y)

      implicit none
      
      real(real64), intent(in) :: x
      real(real64)             :: y

      y = exp(x)

    end function Q_func

end program vert_remap
