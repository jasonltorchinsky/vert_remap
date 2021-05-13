!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main driver for the vertical remap code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TO-DO: Add doc string information
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program vert_remap

  use iso_fortran_env, only: int32, real64
  use vertremap_mod

  implicit none

  integer(int32)            :: nlev ! Number of levels for each grid
  real(real64), allocatable :: dp1(:), dp2(:) ! Grid spacings for each grid
  real(real64), allocatable :: grid1(:), grid2(:) ! Actual grid points
  real(real64), allocatable :: grid1_midpts(:), grid2_midpts(:) ! Midpoints for
                                                          ! each grid cell
  real(real64), allocatable :: QdpOrig(:), QdpNew(:) ! Mass at each grid  
                                    ! mid-point (remapped from grid 2 to grid 1)
  real(real64), allocatable :: QdpDiff(:) ! Difference between remapped Qdp
                                          ! and true Qdp
  character(len=32)         :: arg ! Input argument from command line
  integer(int32)            :: ii ! counter for do loops

  ! Read command line arguments
  do ii = 0, command_argument_count()
    call get_command_argument(ii, arg)
    if (len_trim(arg) .eq. 0) exit
  
    select case(arg)
    case('nlev')
      call get_command_argument(ii+1, arg)
      read(arg,*) nlev
    end select
  end do

  ! Set up grids
  allocate(grid1(nlev+1))
  allocate(grid2(nlev+1))
  
  do ii = 1, nlev+1 ! Fill in grids
    grid2(ii) = real(ii - 1, real64)/real(nlev, real64) ! Target, uniform
    grid1(ii) = grid1_func(grid2(ii)) ! Source, non-uniform
  end do

  allocate(grid1_midpts(nlev))
  allocate(grid2_midpts(nlev))

  do ii = 1, nlev ! Fill in grid midpoints
    grid1_midpts(ii) = (grid1(ii) + grid1(ii+1))/2.0_real64
    grid2_midpts(ii) = (grid2(ii) + grid2(ii+1))/2.0_real64
  end do

  ! Get grid spacings
  allocate(dp1(nlev))
  allocate(dp2(nlev))

  do ii = 1, nlev
    dp1(ii) = grid1(ii+1) - grid1(ii)
    dp2(ii) = grid2(ii+1) - grid2(ii)
  end do

  ! Set Qdp values on grid 1 midpoints (source)
  allocate(QdpOrig(nlev))
  allocate(QdpNew(nlev))
  
  do ii = 1, nlev
    QdpOrig(ii) = Qdp_func(grid1_midpts(ii))
  end do
  QdpNew = QdpOrig

  ! Remap Qdp to grid 1
  call remap_Q_ppm(QdpNew, 1, nlev, 1, dp1, dp2, 10)

  ! Get difference
  allocate(QdpDiff(nlev))
  do ii = 1, nlev
    QdpDiff(ii) = Qdp_func(grid2_midpts(ii)) - QdpNew(ii)
  end do

  ! Output information to .out file for testing
  print *, 'grid 1:'
  print *, grid1

  print *, 'grid 2:'
  print *, grid2

  print *, 'dp1:'
  print *, dp1

  print *, 'dp2:'
  print *, dp2

  print *, 'QdpOrig:'
  print *, QdpOrig

  print *, 'QdpNew:'
  print *, QdpNew

  print *, 'QdpDiff:'
  print *, QdpDiff

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
    function Qdp_func(x) result(y)

      implicit none
      
      real(real64), intent(in) :: x
      real(real64)             :: y

      y = exp(x)

    end function Qdp_func

end program vert_remap
