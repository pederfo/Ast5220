module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point
  real(dp),    allocatable, dimension(:) :: a_eta              ! a values for eta 
  real(dp),    allocatable, dimension(:) :: z_eta              ! z values for eta
  real(dp),    allocatable, dimension(:) :: Omega_rx           ! Omega values for x
  real(dp),    allocatable, dimension(:) :: Omega_mx
  real(dp),    allocatable, dimension(:) :: Omega_bx
  real(dp),    allocatable, dimension(:) :: Omega_lambdax

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2, n3
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, dx2, dx1, da, da2, dx_eta, h1, hmin, eps, dydx, rho_crit, rho_b, rho_m, rho_r, rho_lambda, rho_b0, rho_m0, rho_r0, rho_lambda0, dz, xstart, xstop

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n3          = 100
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation
    h1          = 10.d0**(-8)               ! Stepsize for odesolver
    hmin        = 0.d0                         ! Minimum stepsize
    eps         = 10.d0**(-10)              ! Allowed error


    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))

    ! x grid
    dx = (x_end_rec-x_start_rec) / (n1 -1)
    !x_t(1:200) = [(x_start_rec + ((i-1)*dx), i=1,n1)]
    dx2 = (x_0-x_end_rec) / (n2)
    !x_t(201:500) = [(x_end_rec + ((i-1)*dx2), i=2,n2)]

    x_t(1) = x_start_rec
    do i=2,200
       x_t(i) = x_t(i-1) + dx
    end do
    do i=1, n2
       !x_t(i) = x_t(i-1) + dx2
       x_t(n1+i) = x_end_rec + i*dx2
    end do
    
    !write(*,*) x_t(199), x_t(200), x_t(201), x_t(202)

    !dx1 = (x_start_rec - x_eta1) / (n3 -1)
    !x_t(1:100) = [(x_eta1 + ((i-1)*dx1), i=1,n3)]
    !dx = (x_end_rec-x_start_rec) / (n1 -1)
    !x_t(101:300) = [(x_start_rec + ((i-1)*dx), i=1,n1)]
    !dx2 = (x_0-x_end_rec) / (n2-1)
    !x_t(301:600) = [(x_end_rec + ((i-1)*dx2), i=1,n2)]
    

    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1

    !dx = (xstop - xstart) / real(n_t-1,dp)
    !x_t(1:n_t) = [(xstart + ((i-1)*dx), i=1,n_t)]

    !x_t(500) = 1.d0
    !write(*,*) x_t

    ! a grid
    da = ((1/(1+z_start_rec)) - (1/(1+z_end_rec))) / (n1 -1)
    a_t(1:200) = [(1/(1+z_start_rec) + ((i-1)*da), i=1,n1)]
    da2 = ((1/(1+z_end_rec)) - 1) / (n2 -1)
    a_t(201:500) = [(1/(1+z_end_rec) + ((i-1)*da2), i=1,n2)]
    
        


    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90

    ! allocate arrays for eta
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))
    allocate(a_eta(n_eta))
    allocate(z_eta(n_eta))
    allocate(Omega_bx(n_eta))
    allocate(Omega_mx(n_eta))
    allocate(Omega_rx(n_eta))
    allocate(Omega_lambdax(n_eta))
    
    ! x_eta grid
    dx_eta = (x_eta2 - x_eta1) / (n_eta -1)
    x_eta(1:1000)= [(x_eta1 + ((i-1)*dx_eta), i=1,n_eta)] 

    ! initial eta value
    eta(1)= exp(x_eta(1)) /(H_0*sqrt(Omega_r))

    ! initial a_eta value 
    a_eta(1) = exp(x_eta(1))

    ! initial z_eta value
    !z_eta(1) = z_start_rec

    ! initial rho values
    rho_b0 = Omega_b*rho_c
    rho_m0 = Omega_m*rho_c
    rho_r0 = Omega_r*rho_c
    rho_lambda0 = Omega_lambda*rho_c

    !Calculating eta
    do i=2, n_eta
       eta(i) = eta(i-1)
       a_eta(i) = exp(x_eta(i))
       !z_eta(i) = 1/a_eta(i) - 1
       
       call odeint(eta(i:i), a_eta(i-1), a_eta(i), eps, h1, hmin, derivs_eta, bsstep, output_eta)
    end do
    
    !spline eta
    do i=1, n_eta
       a_eta(i) = exp(x_eta(i))
       call spline(a_eta(i:i), eta(i:i), 1d30, 1d30, eta2(i:i))
    end do

    !Calculating densities
    do i=1, n_eta
       ! critical density
       rho_crit = 3*get_H(x_eta(i))**2/(8*pi*G_grav)

       rho_b = rho_b0*a_eta(i)**(-3)
       rho_m = rho_m0*a_eta(i)**(-3)
       rho_r = rho_r0*a_eta(i)**(-4)
       rho_lambda = rho_lambda0

       Omega_bx(i) = rho_b/rho_crit
       Omega_mx(i) = rho_m/rho_crit
       Omega_rx(i) = rho_r/rho_crit
       Omega_lambdax(i) = rho_lambda/rho_crit


    end do
    
    do i=1, n_eta
       z_eta(i) = 1/a_eta(i) - 1
    end do

  end subroutine initialize_time_mod




  subroutine derivs_eta(x_eta, eta, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_eta
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: a
    a = exp(x_eta)
    
    dydx = c / (a*get_H_p(x_eta))
  end subroutine derivs_eta
  


  subroutine output_eta(a_eta, eta)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: a_eta
    real(dp), dimension(:), intent(in)  :: eta
  end subroutine output_eta


  
  ! Task: Write a function t,hat computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp)             :: a
    a = exp(x)
    get_H = H_0*sqrt((Omega_m + Omega_b)/a**3 + Omega_r/a**4 + Omega_lambda)
  end function get_H



  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp)             :: a
    a = exp(x)
    get_H_p = a*get_H(x)
  end function get_H_p



  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    get_dH_p = 0.5d0*H_0*(-(Omega_m + Omega_b)*(exp(-x))-2.d0*Omega_r*(exp(-2.d0*x)) + 2.d0*Omega_lambda*exp(2.d0*x)) / sqrt((Omega_m+Omega_b)*exp(-x) + Omega_r*exp(-2.d0*x) + Omega_lambda*exp(2.d0*x))

  end function get_dH_p



  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
        
    get_eta = splint(a_eta, eta, eta2, x_in)
    
  end function get_eta

end module time_mod
