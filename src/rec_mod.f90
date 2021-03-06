module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b)                                 :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:)          :: x_rec             ! Grid
  real(dp), allocatable, dimension(:)          :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:)          :: log_tau, log_tau2  
  real(dp), allocatable, dimension(:)          :: n_e, n_e2, logn_e, logn_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:)          :: g, g2, g22        ! Splined visibility function
  real(dp), allocatable, dimension(:)          :: z_rec             ! Redshift
  real(dp), allocatable, dimension(:)          :: X_e ! Fractional electron density, n_e / n_H


contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k, n1, n2
    real(dp)     :: h1, h2, eps, hmin, saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, dx2, f, n_e0, X_e0, xstart, xstop, z_start_rec, z_end_rec, x_start_rec, x_end_rec, a, a_end_rec, a_start, b1, b2, b3, lamba_2s_to_1s, lambda_alpha, C_r, alpha_2, beta, beta_2, phi_2, n_1s
    logical(lgt) :: use_saha
    !real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-8)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo
    !h1         = 1.d-8         ! Stepsize for odesolver
    hmin       = 0.d0         ! Minimum stepsize
    eps        = 1.d-10        ! Allowed error
  


    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(log_tau(n))
    allocate(log_tau2(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(logn_e(n))
    allocate(logn_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    allocate(z_rec(n))

    ! Task: Fill in x (rec) grid

    dx = (xstop - xstart) / real(n-1,dp)
    x_rec(1:n) = [(xstart + ((i-1)*dx), i=1,n)]
    h1 = 0.01d0 *dx
    

    ! Fill in z - redshift
    do i = 1, n
       a = exp(x_rec(i))
       z_rec(i) = 1.d0/a -1
    end do


    ! initial value for tau

    tau(n) = 0.d0

    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n

       a = exp(x_rec(i))
       T_b = T_0/a
       n_b = Omega_b*rho_c/(m_H*a**3)
       

       !Saha's constants
       b1 = 1.d0
       b2 = 1.d0/n_b*(m_e*T_b*k_b/(2.d0*pi))**(1.5d0)*exp(-epsilon_0/(k_b*T_b))*(1.d0/hbar**3)
       b3 = -b2

       if (use_saha) then
          ! Use the Saha equation
          X_e(i) = (-b2 + sqrt(b2**2 + 4.d0*b2))/2.d0



          if (X_e(i) < saha_limit) use_saha = .false.
          !write(*,*) X_e(i)
          !stop
       else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, h1, hmin, derivs_Xe, bsstep, output_Xe) 
          !write(*,*) T_b, X_e(i)
          !stop

       end if
       n_e(i) = X_e(i)*n_b

    end do

    !write(*,*) 'Xe =',X_e(700)
    !stop
    !write(*,*) 'X_rec=', x_rec(1)
    

    ! Task: Compute splined (log of) electron density function
    logn_e = log(n_e)
    call spline(x_rec, logn_e, 1.d30, 1.d30, logn_e2)


    ! Task: Compute optical depth at all grid points
    do i = n-1,1,-1
       tau(i) = tau (i+1)
       call odeint(tau(i:i), x_rec(i+1), x_rec(i), eps, h1, hmin, derivs_tau, bsstep, output_tau)

    end do

 
    ! Task: Compute splined (log of) optical depth
    tau = log(tau)
    call spline(x_rec(1:n-1), tau(1:n-1), 1.d30, 1.d30, tau2(1:n-1))


    
    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec(1:n-1), tau2(1:n-1), 1.d30, 1.d30, tau22(1:n-1))



    ! Task: Compute splined visibility function
    do i = 1, n
       g(i) = -get_dtau(x_rec(i))*exp(-get_tau(x_rec(i)))
       
    end do
    
    call spline(x_rec, g, 1.d30, 1.d30, g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec, g2, 1.d30, 1.d30, g22)


  end subroutine initialize_rec_mod


  subroutine derivs_Xe(x_rec, X_e, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_rec
    real(dp), dimension(:), intent(in)  :: X_e
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: phi_2, alpha_2, beta, beta_2
    real(dp)                            :: lambda_2s_to_1s, T_b, a, n_b, C_r
    real(dp)                            :: lambda_alpha, n_1s, Xe, H

    Xe = X_e(1)
    H = get_H(x_rec)
    a = exp(x_rec)
    T_b = T_0/a
    n_b = Omega_b*rho_c/(m_H*a**3)
    phi_2 = 0.448d0*log(epsilon_0/(k_b*T_b))
    alpha_2 = (64.d0*pi)/sqrt(27.d0*pi)*alpha**2/m_e**2*sqrt(epsilon_0/(k_b*T_b))*phi_2*hbar**2/c
    beta = alpha_2*(m_e*T_b*k_b/(2.d0*pi*hbar**2))**(3.d0/2.d0)*exp(-epsilon_0/(k_b*T_b))
    


    
    if (T_b <= 170) then
       beta_2 = 0.d0
    else
       beta_2 = beta*exp(3.d0*epsilon_0/(4.d0*k_b*T_b))
    end if

    n_1s = (1.d0 - Xe)*n_b
    lambda_2s_to_1s = 8.227d0
    lambda_alpha = H*(3.d0*epsilon_0)**3/((8.d0*pi)**2*n_1s)/(hbar*c)**3
    C_r = (lambda_2s_to_1s + lambda_alpha) / (lambda_2s_to_1s + lambda_alpha + beta_2)
    dydx = C_r/H*(beta*(1.d0-Xe) - n_b*alpha_2*Xe**2)

    !write(*,*) dydx, Xe, C_r !lambda_2s_to_1s, lambda_alpha, beta_2, beta, Xe, alpha_2
    !write(*,*) beta*(1.d0-Xe), n_b*alpha_2*Xe**2

    !stop


  end subroutine derivs_Xe


  subroutine output_Xe(x_rec, X_e)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_rec
    real(dp), dimension(:), intent(in)  :: X_e
  end subroutine output_Xe

  subroutine derivs_tau(x_rec, tau, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_rec
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: Hp, a, ne
    Hp = get_H_p(x_rec)
    a = exp(x_rec)
    ne = get_n_e(x_rec)
    dydx = - c*ne*sigma_T*a/Hp
  end subroutine derivs_tau
  
  subroutine output_tau(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output_tau

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    get_n_e = splint(x_rec,logn_e,logn_e2, x)
    get_n_e = exp(get_n_e)

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau
    integer(i4b)         :: n
    n = size(x_rec)
    if (x>x_rec(n-1)) then 
       get_tau = 0.d0
    else
       get_tau = splint(x_rec(1:n-1), tau(1:n-1), tau2(1:n-1), x)
       !write(*,*) x, get_tau
       get_tau = exp(get_tau)
    end if

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    integer(i4b)         :: n
    n = size(x_rec)
    if (x>x_rec(n-1)) then 
       get_dtau = 0.d0
       !WRITE(*,*) x, x_rec(n-1)
    else
       get_dtau = splint_deriv(x_rec(1:n-1), tau(1:n-1), tau2(1:n-1), x)*get_tau(x)
       
    end if

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    integer(i4b)         :: n
    n = size(x_rec)
    if (x>x_rec(n-1)) then 
       get_ddtau = 0.d0
    else
       get_ddtau = splint(x_rec(1:n-1), tau2(1:n-1), tau22(1:n-1), x)*get_tau(x) - get_dtau(x)**2/get_tau(x)
    end if

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x)

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg


end module rec_mod
