! Define parameters used throughout diablo
! Contains subroutines to read input files and set parameters

module parameters
  ! tomlf library optional
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

  ! Details of the Computational Domain
  ! (We hardwire these into the code so that the compiler may perform
  !  optimizations based on the grid size at compile time).
  integer :: Nx, Ny, Nz, N_th
  include 'grid_def'

  real(rkind) :: time
  integer :: time_step, rk_step
  real(rkind) :: save_flow_time, save_stats_time, save_movie_time

  integer :: previous_time_step

  

  ! Parameters defined in input.dat
  real :: version
  character(len=35)   flavor
  logical             use_mpi
  logical             use_LES
  real(rkind)         nu
  real(rkind)         nu_v_scale
  real(rkind)         beta
  real(rkind)         Lx, Ly, Lz
  real(rkind)         delta_t, dt, delta_t_next_event, kick, ubulk0, px0

  logical             create_new_flow
  real(rkind)         wall_time_limit, time_limit
  real(rkind)         start_wall_time, previous_wall_time, end_wall_time

  logical             variable_dt, first_time
  logical             reset_time
  real(rkind)         CFL
  integer             update_dt

  integer             verbosity
  real(rkind)         save_flow_dt, save_stats_dt
  real(rkind)         save_movie_dt
  real(rkind)         XcMovie, YcMovie, ZcMovie

  logical             create_new_th(1:N_th)
  real(rkind)         Ri(1:N_th), Pr(1:N_th)
  logical             filter_th(1:N_th)
  integer             filter_int(1:N_th)

  integer     num_read_th
  integer     read_th_index(1:N_th)
  real(rkind) dTHdX(1:N_th), dTHdZ(1:N_th)
  real(rkind) dWdX ! Background vorticity


  integer     IC_type, f_type
  logical     physical_noise
  logical     homogeneousX


  ! Periodic
  real(rkind) ek0, ek, epsilon_target
  logical     background_grad(1:N_th)

  ! Rotating Flows
  real(rkind) Ro_inv, grav_x, grav_y, grav_z, delta

  ! Parameters for oscillatory forcing
  real(rkind) omega0, amp_omega0, force_start
  real(rkind) w_BC_Ymax_c1_transient


  ! Numerical parameters
  integer num_per_dir
  integer  time_ad_meth
  integer les_model_type

  ! BCs & Values
  integer :: u_BC_Xmin, v_BC_Xmin, w_BC_Xmin, th_BC_Xmin(1:N_th)
  integer :: u_BC_Xmax, v_BC_Xmax, w_BC_Xmax, th_BC_Xmax(1:N_th)
  integer :: u_BC_Ymin, v_BC_Ymin, w_BC_Ymin, th_BC_Ymin(1:N_th)
  integer :: u_BC_Ymax, v_BC_Ymax, w_BC_Ymax, th_BC_Ymax(1:N_th)
  integer :: u_BC_Zmin, v_BC_Zmin, w_BC_Zmin, th_BC_Zmin(1:N_th)
  integer :: u_BC_Zmax, v_BC_Zmax, w_BC_Zmax, th_BC_Zmax(1:N_th)

  real(rkind) :: u_BC_Xmin_c1, v_BC_Xmin_c1, w_BC_Xmin_c1
  real(rkind) :: u_BC_Ymin_c1, v_BC_Ymin_c1, w_BC_Ymin_c1
  real(rkind) :: u_BC_Zmin_c1, v_BC_Zmin_c1, w_BC_Zmin_c1
  real(rkind) :: th_BC_Xmin_c1(1:N_th), th_BC_Ymin_c1(1:N_th), th_BC_Zmin_c1(1:N_th)
  real(rkind) :: u_BC_Xmax_c1, v_BC_Xmax_c1, w_BC_Xmax_c1
  real(rkind) :: u_BC_Ymax_c1, v_BC_Ymax_c1, w_BC_Ymax_c1
  real(rkind) :: u_BC_Zmax_c1, v_BC_Zmax_c1, w_BC_Zmax_c1
  real(rkind) :: th_BC_Xmax_c1(1:N_th), th_BC_Ymax_c1(1:N_th), th_BC_Zmax_c1(1:N_th)


contains


  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
    ! read in input parameters
#ifdef TOML_INPUT
    ! read parameters from input.toml
    call read_input_toml
#else
    ! read in parameters from input.dat
    call read_input_dat
#endif

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
    real(rkind) Re
    
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
    nu = 1.d0 / Re
    read (11, *)
    read (11, *) nu_v_scale
    read (11, *)
    read (11, *) num_per_dir, create_new_flow
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
    real(rkind) ro

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
    read (11, *) ro
    Ro_inv = 1.d0 / ro
    read (11, *)
    read (11, *) delta
    read (11, *)
    !read (11, *) dWdX ! Background vorticity
    !read (11, *)
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

    if (rank == 0) write (*, '("Ro Inverse = " ES26.18)') Ro_inv


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


    return
  end

  ! only define toml input options if using the tomlf library
#ifdef TOML_INPUT 

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_toml
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  end

  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  subroutine read_input_chan_toml
    !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  end

#endif ! TOML_INPUT


end module parameters
