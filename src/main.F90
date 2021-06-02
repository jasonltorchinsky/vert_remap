!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main driver for the vertical remap code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TO-DO: Add doc string information
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program vert_remap

  use iso_fortran_env, only: int32, real64
  use netcdf
  use output_mod
  use utils_mod
  use vertremap_mod

  implicit none

  integer(int32)            :: ncell, nlev ! Number of cells, and levels for 
  ! each grid
  character(len=8)          :: ogrid, tfunc ! Form of original grid and test
  ! function
  integer(int32)            :: alg ! Turn limiter on (10) or off (11)
  character(len=3)          :: alg_str ! String corresponding to limiter on
  ! or off
  integer(int32)            :: curr_time(8) ! Current time
  integer(int32)            :: seed ! Seed for RNG
  real(real64), allocatable :: dp1(:), dp2(:) ! Grid spacings for each grid
  real(real64), allocatable :: grid1(:), grid2(:) ! Grid at cell boundaries
  ! (0, ..., H)
  real(real64), allocatable :: grid1_stg(:), grid2_stg(:) ! Staggered grid at
  ! cell centers 
  ! (dx/2, ..., H-dx/2)
  real(real64), allocatable :: QOrig(:), QNew(:) ! Average density in each cell
  real(real64), allocatable :: QdpOrig(:), QdpNew(:) ! Mass in each cell
  real(real64), allocatable :: QTrue(:) ! True density on the transformed grid
  type(outfile)             :: dataFile ! File we write the output to
  character(len=50)         :: outdirName, outfileName ! Directory, file for
  ! output
  character(len=32)         :: arg ! Input argument from command line
  integer(int32)            :: ii ! counter for do loops

  ! Set default input values
  ncell = 0
  nlev = 0
  ogrid = 'sqr'
  tfunc = 'exp'
  alg = 10
  alg_str = 'on'
  call date_and_time(VALUES=curr_time)
  seed = 0_int32
  do ii = 4, 8
     seed = seed + curr_time(ii)
  end do
  
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
        print *, '  ogrid: Form of original grid. Options: sqr (x^2),'
        print *, '         cub (x^3), sig (sigmoid-ish), sin (sine),'
        print *, '         uni (uniform).'
        print *, '  tfunc: Test density function. Options: exp (e^x),'
        print *, '         stp (step function), sig (sigmoid-ish),'
        print *, '         wdg (wedge).'
        print *, '  alg: Which algorithm to use. Options: "on" for'
        print *, '       original algorithm with limiter on on the'
        print *, '       boundaries, "off" for original algorithm'
        print *, '       with limit off on the boundaries, "new"'
        print *, '       for new algorithm.'
        print *, '  seed: Seed for RNG.'
        stop
     case('ncell')
        call get_command_argument(ii + 1, arg)
        read(arg, *) ncell
     case('nlev')
        call get_command_argument(ii + 1, arg)
        read(arg, *) nlev
     case('ogrid')
        call get_command_argument(ii + 1, ogrid)
     case('tfunc')
        call get_command_argument(ii + 1, tfunc)
     case('alg')
        call get_command_argument(ii + 1, alg_str)
        if (alg_str .eq. 'on') then
           alg = 10
        else if (alg_str .eq. 'off') then
           alg = 11
        else if (alg_str .eq. 'new') then
           alg = 20
        end if
     case('seed')
        call get_command_argument(ii + 1, arg)
        read(arg, *) seed
     end select
  end do

  ! Check for input errors
  if (ncell .eq. 0) ncell = nlev - 1 ! Assume ncell wasn't set first
  if (nlev .eq. 0) nlev = ncell + 1 ! The assume nlev wasn't set
  if ((nlev - ncell) .ne. 1) then ! Too many/few cells for levels.
     print *, '  ISSUE: Inappropriate number of cells and levels...'
     print *, '  Assuming that the number of levels is correct.'
     ncell = nlev - 1
  end if

  call rand_init(seed)

  write(*,*) '~~ Runnning vert_remap with ncell = ', ncell

  ! Set up grids, these are the boundaries of the cells.
  allocate(grid1(nlev))
  allocate(grid2(nlev))

  do ii = 1, nlev ! Fill in grids
     grid2(ii) = real(ii - 1, real64)/real(nlev - 1, real64) ! Target, uniform
     grid1(ii) = grid1_func(grid2(ii), ogrid, nlev) ! Source, non-uniform
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
     QOrig(ii) = Q_func(grid1_stg(ii), tfunc)
     QdpOrig(ii) = QOrig(ii) * dp1(ii)
  end do

  ! Remap Qdp to grid 1
  allocate(QdpNew(ncell)) ! Remap is in-place, so we make a new array.
  QdpNew = QdpOrig
  call remap1(QdpNew, 1, ncell, 1, dp1, dp2, alg)

  ! Get density from mass
  allocate(QNew(ncell))
  do ii = 1, ncell
     QNew(ii) = QdpNew(ii) / dp2(ii)
  end do

  ! Get true density on the transformed grid
  allocate(QTrue(ncell))
  do ii = 1, ncell
     QTrue(ii) = Q_func(grid2_stg(ii), tfunc)
  end do

  ! Output information
  outdirName = '../output' // repeat(' ', 41)
100 format (A, A, A, A, I0.8, A, A, A)
  write(outfileName, 100) trim(adjustl(ogrid)), '_', trim(adjustl(tfunc)), &
       & '_', ncell, '_', trim(adjustl(alg_str)), '.nc'
  call mkdir(outdirName)
  call netcdf_init_outfile(2, &
       & ['ncell' // repeat(' ', 10), 'nlev' // repeat(' ', 11)], &
       & [ncell, nlev], &
       & 11, &
       & ['grid1' // repeat(' ', 10), 'dp1' // repeat(' ', 12), &
       &  'grid1_stg' // repeat(' ', 6), 'QOrig' // repeat(' ', 10), &
       &  'QdpOrig' // repeat(' ', 8), 'grid2' // repeat(' ', 10), &
       &  'dp2' // repeat(' ', 12), 'grid2_stg' // repeat(' ', 6), &
       &  'QNew' // repeat(' ', 11), 'QdpNew' // repeat(' ', 9), &
       &  'QTrue' // repeat(' ', 10)], &
       & [(nf90_double, ii = 1, 11)], &
       & [2, (1, ii = 2, 5), 2, (1, ii = 7, 11)], &
       & outdirName, outfileName, dataFile)
  call netcdf_add_metadata(dataFile)
  call netcdf_write_output(nlev, ncell, grid1, dp1, grid1_stg, QOrig, &
       & QdpOrig, grid2, dp2, grid2_stg, QNew, QdpNew, QTrue, dataFile)
  call netcdf_close_outfile(dataFile)

contains

  ! Function to get the grid points for grid 1 (source) from 
  ! grid 2 (target, uniform)
  function grid1_func(x, ogrid, nlev) result(y)

    implicit none

    real(real64), intent(in)   :: x ! In [0, 1]
    character(len=8)           :: ogrid
    integer(int32), intent(in) :: nlev
    real(real64)               :: y ! Should be in [0, 1]
    real(real64)               :: pi_dp
    real(real64)               :: uni_spc
    real(real64)               :: rand_dp

    select case(trim(adjustl(ogrid)))
    case('sqr')
       y = x**2
    case('cub')
       y = x**3
    case('sig')
       if ((x .gt. 1.0_real64 / real(nlev - 1, real64)) &
            & .or. (x .lt. 1.0_real64 - 1.0_real64 / real(nlev - 1, real64))) then
          ! Odd looking condition which relaxes the real comparison
          y = 1.0_real64/(1.0_real64 + exp(-30.0_real64*(x-0.5_real64)))
       else
          y = x
       end if
    case('sin')
       pi_dp = 4.0_real64*atan(1.0_real64)
       y = 0.5_real64 * (1 - cos(pi_dp * x))
    case('rng')
       if ((x .lt. 1.0_real64 / real(nlev - 1, real64)) &
            & .or. (x .gt. 1.0_real64 - 1.0_real64 / real(nlev - 1, real64))) then
          y = x
       else
          uni_spc = 1.0_real64 / (nlev - 1.0_real64)
          call random_number(rand_dp)
          y = x + 0.25_real64 * uni_spc * rand_dp
       end if
    case('uni')
       y = x
    end select


  end function grid1_func

  ! Function to get the value of Qdp at a given grid point
  function Q_func(x, tfunc) result(y)

    implicit none

    real(real64), intent(in) :: x
    character(len=8)         :: tfunc
    real(real64)             :: y

    select case(trim(adjustl(tfunc)))
    case('exp')
       y = exp(x)
    case('stp')
       if (x .le. 0.5) then
          y = 0.0_real64
       elseif (x .gt. 0.5) then
          y = 1.0_real64
       end if
    case('sig')
       y = 1.0_real64/(1.0_real64 + exp(-30.0_real64*(x-0.5_real64)))
    case('wdg')
       if (x .le. 0.5) then
          y = 0.5_real64*x
       else if (x .gt. 0.5) then
          y = 1.5_real64*x - 0.5_real64
       end if
    end select

  end function Q_func

  subroutine rand_init(seed)

    use iso_fortran_env, only: int32

    implicit none

    integer(int32), intent(in) :: seed !< the seed to use with the rng
    integer(int32) :: i !< counter for do loops
    integer(int32), dimension(33) :: seedarr !< the array generated for putting
    !! into random_seed

    seedarr = 0_int32
    seedarr(1) = 104729_int32 + 5_int32 * seed ! this is pretty aribtrary

    do i = 2, 33
       seedarr(i) = mod(420_int32 * seedarr(i-1) + 69_int32, 4294967_int32)
       ! this is also pretty arbitrary
    end do

    call random_seed(put = seedarr)

  end subroutine rand_init

end program vert_remap
