!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main driver for the vertical remap code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TO-DO: Add doc string information
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program vert_remap

  use iso_fortran_env, only: int32, real64
  use netcdf
  use conv_comb_mod
  use get_insecs_mod
  use mass_borrow_mod
  use output_mod
  use utils_mod
  use vertremap_mod

  implicit none

  integer(int32)            :: ncell, nlev ! Number of cells, and levels for 
  ! each grid
  character(len=8)          :: ogrid, tfunc ! Form of original grid and test
  ! function
  integer(int32)            :: alg ! Set algorithm. 10 - limiter on. 11 - limiter off.
  ! 20 - new algorithm.
  character(len=3)          :: alg_str ! String corresponding to limiter on
  ! or off
  integer(int32)            :: verbose ! Flag top print debugging messages.
  integer(int32)            :: curr_time(8) ! Current time
  integer(int32)            :: seed ! Seed for RNG
  real(real64), allocatable :: dp1(:), dp2(:) ! Grid spacings for each grid
  real(real64), allocatable :: grid1(:), grid2(:) ! Grid at cell boundaries
  ! (0, ..., H)
  real(real64), allocatable :: grid1_stg(:), grid2_stg(:) ! Staggered grid at
  ! cell centers 
  ! (dx/2, ..., H-dx/2)
  real(real64), allocatable :: Q1(:), Q2(:) ! Average density in each cell
  real(real64), allocatable :: Qdp1(:), Qdp2(:) ! Mass in each cell
  real(real64), allocatable :: Qdp3(:) ! Extra arrays to hold mass in each cell
  ! for extra purposes
  integer(int32), allocatable :: cell_insecs(:,:)
  real(real64), allocatable :: ubnds(:), lbnds(:) ! Local bounds
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
  verbose = 0_int32
  
  ! Read command line arguments
  do ii = 0, command_argument_count()
     call get_command_argument(ii, arg)
     if (len_trim(arg) .eq. 0) exit

     select case(arg)
     case('-h', '--help')
        print *, 'vert_remap flags: '
        print *, '  -h, --help: Returns this message.'
        print *, '  -v, --verbose: Prints debugging messages.'
        print *, 'vert_remap arguments: '
        print *, '  nlev: Number of grid levels (cell boundaries).'
        print *, '  ncell: Number of cells (regions between levels).'
        print *, '  ogrid: Form of original grid. Options: cub (x^3),'
        print *, '         rng (+-h/4), sig (sigmoid-ish), sin (sine),'
        print *, '         sqr (x^2), uni (uniform).'
        print *, '  tfunc: Test density function. Options: ,'
        print *, '         asr (asymmetric parabola), exp (e^x),'
        print *, '         gau (e^(-x^2)), nxp (e^(1-x)), osc (sin(8 pi x)),'
        print *, '         sig (sigmoid-ish), stp (step function),'
        print *, '         sqr (centered parabola), wdg (wedge).'
        print *, '  alg: Which algorithm to use. Options: "on" for'
        print *, '       original algorithm with limiter on everywhere,'
        print *, '       "off" for original algorithm with the limiter'
        print *, '       off on the boundaries, "new" for new algorithm'
        print *, '       "ngh" for new algorithm with ghost cells.'
        print *, '  seed: Seed for RNG.'
        stop
     case('-v', '--verbose')
        verbose = 1_int32
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
        else if (alg_str .eq. 'fff') then
           alg = 12
        else if (alg_str .eq. 'new') then
           alg = 20
        else if (alg_str .eq. 'ngh') then
           alg = 21
        else if (alg_str .eq. 'ngf') then
           alg = 22
        else if (alg_str .eq. 'nmb') then
           alg = 23
        else if (alg_str .eq. 'cs8') then
           alg = 24
        else if (alg_str .eq. 'lmb') then
           alg = 31
        else if (alg_str .eq. 'lco') then
           alg = 32
        else if (alg_str .eq. 'llc') then
           alg = 33
        else if (alg_str .eq. 'gmb') then
           alg = 34
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
  allocate(Q1(ncell))
  allocate(Qdp1(ncell))

  do ii = 1, ncell
     Q1(ii) = Q_func(grid1_stg(ii), tfunc)
     Qdp1(ii) = Q1(ii) * dp1(ii)
  end do

  ! Remap Qdp to grid 1
  allocate(Qdp2(ncell)) ! Remap is in-place, so we make a new array.
  Qdp2 = Qdp1
  if (alg .lt. 30) then
     call remap1(Qdp2, 1, ncell, 1, dp1, dp2, alg, verbose)
  else if (alg .eq. 31) then ! q_alg = 11 +  mass-borrowing
     call remap1(Qdp2, 1, ncell, 1, dp1, dp2, 11, verbose)
     ! Get bounds
     allocate(cell_insecs(2,ncell))
     call get_insecs(ncell, grid1, grid2, cell_insecs)
     do ii = 1, ncell ! Domain of dependence is two greater than intersection
        if (cell_insecs(1,ii) .ge. 3) then
           cell_insecs(1,ii) = cell_insecs(1,ii) - 2
        else if (cell_insecs(1,ii) .eq. 2) then
           cell_insecs(1,ii) = cell_insecs(1,ii) - 1
        end if

        if (cell_insecs(2,ii) .le. ncell-2) then
           cell_insecs(2,ii) = cell_insecs(2,ii) + 2
        else if (cell_insecs(2,ii) .eq. ncell-1) then
           cell_insecs(2,ii) = cell_insecs(2,ii) + 1
        end if
     end do
     
     allocate(ubnds(ncell))
     allocate(lbnds(ncell))
     do ii = 1, ncell
        ubnds(ii) = maxval(Qdp1(cell_insecs(1,ii):cell_insecs(2,ii)) &
             & / dp1(cell_insecs(1,ii):cell_insecs(2,ii)))
        lbnds(ii) = minval(Qdp1(cell_insecs(1,ii):cell_insecs(2,ii)) &
             & / dp1(cell_insecs(1,ii):cell_insecs(2,ii)))
     end do
     call borrow_mass(ncell, Qdp2, dp2, ubnds, lbnds)
  else if (alg .eq. 32) then ! q_alg = 11 + linear combination
     allocate(Qdp3(ncell))
     Qdp3 = Qdp1
     call remap1(Qdp2, 1, ncell, 1, dp1, dp2, 11, verbose)
     call remap1(Qdp3, 1, ncell, 1, dp1, dp2, 10, verbose)
     call conv_comb(ncell, Qdp3, Qdp2, dp2, maxval(Qdp1/dp1), &
          & minval(Qdp1/dp1))
  else if (alg .eq. 33) then
     allocate(Qdp3(ncell))
     Qdp3 = Qdp1
     call remap1(Qdp2, 1, ncell, 1, dp1, dp2, 11, verbose)
     call remap1(Qdp3, 1, ncell, 1, dp1, dp2, 10, verbose)
     call conv_comb(ncell/2, Qdp3(1:ncell/2), Qdp2(1:ncell/2), dp2(1:ncell/2), &
          & maxval(Qdp1(1:ncell/2)/dp1(1:ncell/2)), &
          & minval(Qdp1(1:ncell/2)/dp1(1:ncell/2)))
     call conv_comb(ncell/2, Qdp3(ncell/2+1:ncell), Qdp2(ncell/2+1:ncell), &
          & dp2(ncell/2+1:ncell), &
          & maxval(Qdp1(ncell/2+1:ncell)/dp1(ncell/2+1:ncell)), &
          & minval(Qdp1(ncell/2+1:ncell)/dp1(ncell/2+1:ncell)))
  else if (alg .eq. 34) then
     call remap1(Qdp2, 1, ncell, 1, dp1, dp2, 11, verbose)
     ! Get bounds
     allocate(ubnds(ncell))
     allocate(lbnds(ncell))
     ubnds = maxval(Qdp1/dp1)
     lbnds = minval(Qdp1/dp1)
     call borrow_mass(ncell, Qdp2, dp2, ubnds, lbnds)
  end if
  

  ! Get density from mass
  allocate(Q2(ncell))
  do ii = 1, ncell
     Q2(ii) = Qdp2(ii) / dp2(ii)
  end do

  ! Get true density on the transformed grid
  allocate(QTrue(ncell))
  do ii = 1, ncell
     QTrue(ii) = Q_func(grid2_stg(ii), tfunc)
  end do

  if (verbose .eq. 1) then
     call check_global_bounded(ncell, Q2, maxval(Q1), minval(Q1), ii)
     if (ii .eq. 0) then ! Global bounds violated
        print *, ' ~~ Global bounds have been violated:'
        print *, '    Max extrema difference: ', &
             & max(maxval(Q2) - maxval(Q1), minval(Q1) - minval(Q2))
     else if (ii .eq. 1) then ! Global bounds passed
        continue
     end if
  end if

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
       &  'grid1_stg' // repeat(' ', 6), 'Q1' // repeat(' ', 13), &
       &  'Qdp1' // repeat(' ', 11), 'grid2' // repeat(' ', 10), &
       &  'dp2' // repeat(' ', 12), 'grid2_stg' // repeat(' ', 6), &
       &  'Q2' // repeat(' ', 13), 'Qdp2' // repeat(' ', 11), &
       &  'QTrue' // repeat(' ', 10)], &
       & [(nf90_double, ii = 1, 11)], &
       & [2, (1, ii = 2, 5), 2, (1, ii = 7, 11)], &
       & outdirName, outfileName, dataFile)
  call netcdf_add_metadata(dataFile)
  call netcdf_write_output(nlev, ncell, grid1, dp1, grid1_stg, Q1, &
       & Qdp1, grid2, dp2, grid2_stg, Q2, Qdp2, QTrue, dataFile)
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
    case('sng')
       if ((x .lt. 1.0_real64 / real(nlev - 1, real64)) &
            & .or. (x .gt. 1.0_real64 - 1.0_real64 / real(nlev - 1, real64))) then
          y = x
       else
          uni_spc = 1.0_real64 / (nlev - 1.0_real64)
          call random_number(rand_dp)
          y = x + 3.125e-2_real64 * uni_spc * rand_dp
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
    real(real64)             :: pi_dp
    
    select case(trim(adjustl(tfunc)))
    case('exp')
       y = exp(x) / exp(1.0_real64)
    case('stp')
       if (x .le. 0.5) then
          y = 0.0_real64
       elseif (x .gt. 0.5) then
          y = 1.0_real64
       end if
    case('sig')
       y = 1.0_real64 / (1.0_real64 + exp(-10.0_real64 * x))
    case('wdg')
       if (x .le. 0.5) then
          y = 0.5_real64*x
       else if (x .gt. 0.5) then
          y = 1.5_real64*x - 0.5_real64
       end if
    case('nxp')
       y = exp(1.0_real64 - x) / exp(1.0_real64)
    case('sqr')
       y = 4.0_real64*(x-0.5_real64)**2
    case('osc')
       pi_dp = 4.0_real64 * atan(1.0_real64)
       y = 0.5_real64 * (1.0_real64 + sin(8.0_real64 * pi_dp * x))
    case('gau')
       y = exp(-0.5_real64 * ((x - 0.5_real64)/0.05_real64)**2_int32)
    case('asr')
       y = (16.0_real64 / 9.0_real64) * (x - 0.25_real64)**2_int32
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
