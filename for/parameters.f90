! Define parameters used throughout diablo
! Contains subroutines to read input files and set parameters
! Option to use toml-f to read in the inputs as a toml file

module parameters
#ifdef TOML_INPUT 
  use tomlf, only: toml_table, toml_parse, toml_error, get_value  
#endif
  implicit none
  save

  ! current version - update if code is edited to invalidate old input files
  real :: current_version = 3.4

  ! Specify data-types
  integer, parameter :: single_kind = kind(0.0)
  integer, parameter :: double_kind = kind(0.d0)
  integer, parameter :: rkind = double_kind

  ! Grid
  ! (We hardwire these into the code so that the compiler may perform
  !  optimizations based on the grid size at compile time).
  integer :: Nx, Ny, Nz, N_th
  include 'grid_def'

  ! Parameters set in inputs
  real :: version
  ! scheme
  character(len=35) :: flavor
  logical :: use_mpi ! must be true
  logical :: use_LES
  integer :: num_per_dir ! must be 2
  integer :: time_ad_meth ! must be 1
  integer :: les_model_type
  real(rkind) :: beta
  real(rkind) :: nu_v_scale
  ! physical
  real(rkind) :: Lx, Ly, Lz
  real(rkind) :: Re, nu 
  real(rkind) :: Ro, Ro_inv 
  real(rkind) :: delta
  real(rkind) :: grav_x, grav_y, grav_z
  ! timestepping
  real(rkind) :: wall_time_limit, time_limit
  real(rkind) :: delta_t
  logical :: variable_dt
  real(rkind) :: CFL
  integer :: update_dt
  ! output
  integer :: verbosity
  real(rkind) :: save_flow_dt, save_stats_dt, save_movie_dt
  real(rkind) :: XcMovie, YcMovie, ZcMovie
  ! initial conditions 
  logical :: create_new_flow, reset_time
  integer :: IC_type
  real(rkind) :: kick 
  logical :: physical_noise
  ! forcing
  integer :: f_type
  real(rkind) :: ubulk0, px0, omega0, amp_omega0, force_start
  ! velocity bcs
  integer :: u_BC_Xmin, v_BC_Xmin, w_BC_Xmin
  integer :: u_BC_Xmax, v_BC_Xmax, w_BC_Xmax
  integer :: u_BC_Ymin, v_BC_Ymin, w_BC_Ymin
  integer :: u_BC_Ymax, v_BC_Ymax, w_BC_Ymax
  integer :: u_BC_Zmin, v_BC_Zmin, w_BC_Zmin
  integer :: u_BC_Zmax, v_BC_Zmax, w_BC_Zmax
  real(rkind) :: u_BC_Xmin_c1, v_BC_Xmin_c1, w_BC_Xmin_c1
  real(rkind) :: u_BC_Ymin_c1, v_BC_Ymin_c1, w_BC_Ymin_c1
  real(rkind) :: u_BC_Zmin_c1, v_BC_Zmin_c1, w_BC_Zmin_c1
  real(rkind) :: u_BC_Xmax_c1, v_BC_Xmax_c1, w_BC_Xmax_c1
  real(rkind) :: u_BC_Ymax_c1, v_BC_Ymax_c1, w_BC_Ymax_c1
  real(rkind) :: u_BC_Zmax_c1, v_BC_Zmax_c1, w_BC_Zmax_c1
  ! scalars 
  logical :: create_new_th(1:N_th)
  logical :: filter_th(1:N_th)
  integer :: filter_int(1:N_th)
  real(rkind) :: Ri(1:N_th), Pr(1:N_th)
  integer :: th_BC_Xmin(1:N_th), th_BC_Ymin(1:N_th), th_BC_Zmin(1:N_th)
  integer :: th_BC_Xmax(1:N_th), th_BC_Zmax(1:N_th), th_BC_Ymax(1:N_th)
  real(rkind) :: th_BC_Xmin_c1(1:N_th), th_BC_Ymin_c1(1:N_th), th_BC_Zmin_c1(1:N_th)
  real(rkind) :: th_BC_Xmax_c1(1:N_th), th_BC_Ymax_c1(1:N_th), th_BC_Zmax_c1(1:N_th)

  ! parameters defined in set_parameters
  logical :: homogeneousX
  real(rkind) :: w_BC_Ymax_c1_transient
  real(rkind) :: dTHdX(1:N_th), dTHdZ(1:N_th)

contains


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_inputs
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! read in input parameters
#ifdef TOML_INPUT
    ! read parameters from input.toml
    call read_input_toml
#else
    ! read in parameters from input.dat
    call read_input_dat
    call read_input_chan
#endif

    nu = 1.d0 / Re
    Ro_inv = 1.d0 / Ro

  end


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine log_input_parameters
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! write out the input parameters
    ! should only be called by the rank 0 process

    integer n  
    
    write (*, '("Flavor: ", A35)') flavor
    write (*, '("Nx = ", I10)') Nx
    write (*, '("Ny = ", I10)') Ny
    write (*, '("Nz = ", I10)') Nz
    do n = 1, N_th
      write (*, '("Scalar Number: ", I2)') n
      write (*, '("  Richardson number = ", ES12.5)') Ri(n)
      write (*, '("  Prandtl number    = ", ES12.5)') Pr(n)
    end do
    write (*, '("Use LES: " L1)') use_LES
    write (*, '("Nu   = ", ES12.5)') nu
    write (*, '("Beta = ", ES12.5)') beta

  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_dat
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    integer n
    
    open (11, file='input.dat', form='formatted', status='old')

    ! Read input file.
    !   (Note - if you change the following section of code, update the
    !    current_version number to make obsolete previous input files !)

    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) flavor, version
    if (version /= current_version) stop 'Wrong input data format.'
    read (11, *)
    read (11, *) use_mpi, use_LES
    if (use_mpi .eqv. .false.) stop 'Serial processing has been deprecated in diablo3.'
    read (11, *)
    read (11, *) Re, beta, Lx, Lz, Ly
    read (11, *)
    read (11, *) nu_v_scale
    read (11, *)
    read (11, *) num_per_dir, create_new_flow
    if (num_per_dir /= 2) stop 'DIABLO only supports channel geometry (num_per_dim = 2)'
    read (11, *)
    read (11, *) wall_time_limit, time_limit, delta_t, reset_time, &
      variable_dt, CFL, update_dt
    read (11, *)
    read (11, *) verbosity, save_flow_dt, save_stats_dt, save_movie_dt, XcMovie, ZcMovie, YcMovie
    read (11, *)
    ! Read in the parameters for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) create_new_th(n)
      read (11, *)
      read (11, *) filter_th(n), filter_int(n)
      read (11, *)
      read (11, *) Ri(n), Pr(n)
    end do
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_chan
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    
    integer n
    

    ! Read in input parameters specific for channel flow case
    open (11, file='input_chan.dat', form='formatted', status='old')
    ! Read input file.

    current_version = 3.4
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *)
    read (11, *) version
    if (version /= current_version) &
      stop 'Wrong input data format input_chan'
    read (11, *)
    read (11, *) time_ad_meth
    read (11, *)
    read (11, *) les_model_type
    read (11, *)
    read (11, *) IC_type, kick, physical_noise
    read (11, *)
    read (11, *) Ro
    read (11, *)
    read (11, *) delta
    read (11, *)
    read (11, *) grav_x, grav_z, grav_y
    read (11, *)
    read (11, *) f_type, ubulk0, px0, omega0, amp_omega0, force_start
    read (11, *)
    read (11, *)
    read (11, *) u_BC_Ymin, u_BC_Ymin_c1
    read (11, *)
    read (11, *) w_BC_Ymin, w_BC_Ymin_c1
    read (11, *)
    read (11, *) v_BC_Ymin, v_BC_Ymin_c1
    read (11, *)
    read (11, *) u_BC_Ymax, u_BC_Ymax_c1
    read (11, *)
    read (11, *) w_BC_Ymax, w_BC_Ymax_c1
    read (11, *)
    read (11, *) v_BC_Ymax, v_BC_Ymax_c1
    read (11, *)
    ! Read in boundary conditions and background gradients for the N_th scalars
    do n = 1, N_th
      read (11, *)
      read (11, *) th_BC_Ymin(n), th_BC_Ymin_c1(n)
      read (11, *)
      read (11, *) th_BC_Ymax(n), th_BC_Ymax_c1(n)
    end do

    


    return
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine set_parameters
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|

    ! Compensate no-slip BC in the GS flow direction due to dTHdx
    !   AND also define dTHdx & dTHdz
    if (IC_Type == 4 .or. IC_Type == 5) then ! Infinite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 1.d0
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 1.d0
      end if
      dTHdX(1) = Ro_inv / delta
      dTHdZ(1) = 0.d0

    else if (IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then ! Finite Front
      if (w_BC_Ymin == 1) then
        w_BC_Ymin_c1 = w_BC_Ymin_c1 - 2.d0 * delta / Lx
      end if
      if (w_BC_Ymax == 1) then
        w_BC_Ymax_c1 = w_BC_Ymax_c1 - 2.d0 * delta / Lx
      end if
      dTHdX(1) = 2.d0 / Lx * Ro_inv
      dTHdZ(1) = 0.d0

    end if

    w_BC_Ymax_c1_transient = w_BC_Ymax_c1 ! Mean of surface forcing (compensating for GS flow)


    ! Set the valid averaging directions depending on the IC
    if (IC_Type == 5 .or. IC_Type == 6 .or. IC_Type == 7 .or. IC_Type == 8) then
      homogeneousX = .false.
    else ! Infinite, homogeneous front (or other IC...)
      homogeneousX = .true. ! Assume the x-direction is a valid averaging dimension
    endif
  end

  ! only define toml input option if using the tomlf library
#ifdef TOML_INPUT 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_toml
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  end

#endif 
! TOML_INPUT


end module parameters
