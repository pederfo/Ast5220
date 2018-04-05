module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter          :: a_init   = 1.d-8
  real(dp),     parameter          :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter          :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter          :: lmax_int = 6
  integer(i4b), parameter          :: n_br     = 100
  integer(i4b), parameter          :: n_tot    = 600
  integer(i4b), parameter          :: ns = 5000

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta
  

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

  real(dp), allocatable, dimension(:)     :: x_br  !x grid before recombination 
  real(dp), allocatable, dimension(:)     :: x_tot

  real(dp), allocatable, dimension(:,:)     :: S_lores
  real(dp), allocatable, dimension(:,:,:,:) :: S_coeff


contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), allocatable, dimension(:),   intent(out) :: k, x
    real(dp), allocatable, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j, t
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddH_p, Pi, dPi, ddPi,r,rp,u,up,dx,x0,y0
    !real(dp), allocatable, dimension(:,:) :: S_lores
    !real(dp), allocatable, dimension(:,:,:,:) :: S_coeff

    !ns = 5000

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !

    allocate(x(ns))
    allocate(k(ns))
    allocate(S(ns,ns))
    allocate(S_lores(n_tot, n_k))
    allocate(S_coeff(4,4,n_tot,n_k))

    !write(*,*) x

    do i=1, ns
       !k(i) = k_min + (k_max-k_min)*((i-1)/5000.d0)**2
       k(i) = minval(ks) + (maxval(ks)-minval(ks))*((i-1)/5000.d0)**2
       
    end do

    dx      = (0.d0 - log(a_init))/ (ns -1.d0)
    x(1:ns) = [(log(a_init) + ((i-1)*dx), i=1, ns)] 

    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    do j=1, n_k
       do i=1, n_tot
          !x = x_tot(i)
          tau = get_tau(x_tot(i))
          dt = get_dtau(x_tot(i))
          ddt = get_ddtau(x_tot(i))
          H_p = get_H_p(x_tot(i))
          dH_p = get_dH_p(x_tot(i))

          g = get_g(x_tot(i))!-dt*exp(-tau)
          
          

          !if (i == 16) then
             !write(*,*) g, -dt*exp(-tau), tau, dt
             !stop
          !end if
          dg = get_dg(x_tot(i))
          ddg = get_ddg(x_tot(i))
          Pi = Theta(i,2,j)
          dPi = dTheta(i,2,j)
          ddPi = 2.d0*c*ks(j)/5.d0*(Theta(i,1,j)*dH_p - dTheta(i,1,j)*H_p)/H_p**2 - 3.d0*c*ks(j)/5.d0*(Theta(i,3,j)*dH_p - dTheta(i,3,j)*H_p)/H_p**2 + 9.d0/10.d0*(ddt*Theta(i,2,j) + dt*dTheta(i,2,j))

          !write(*,*) ddPi

          S_lores(i,j) = g*(Theta(i,0,j)+Psi(i,j)+0.25d0*Pi) &
               + exp(-tau)*(dPsi(i,j)-dPhi(i,j)) &
               - 1.d0/(c*ks(j))*(g*H_p*dv_b(i,j)+v_b(i,j)*g*dH_p + v_b(i,j)*H_p*dg) &
               + 0.75d0*(c*ks(j))**2*(0.5d0*H_0**2*((Omega_m + Omega_b)*exp(-x_tot(i)) + 4.d0*Omega_r*exp(-2.d0*x_tot(i)) + 4.d0*Omega_lambda*exp(2.d0*x_tot(i)))*g*Pi + 3.d0*H_p*dH_p*(dg*Pi + g*dPi) + H_p**2*(ddg*Pi + 2*dg*dPi + g*ddPi)) 

          !!(dH_p*(H_p*g*dPi +Pi*g*dH_p + Pi*H_p*dg) + H_p*(H_p*g*ddPi + Pi*g*ddH_p + Pi*H_p*ddg + 2.d0*dPi*dH_p*g + 2.d0*dPi*dg*H_p + 2.d0*dH_p*dg*Pi))

       end do
       
       S_lores(n_tot-50:n_tot,j)=S_lores(n_tot-50,j)
       !stop
    end do
    



    !   2) Then spline this function with a 2D spline

    
    !write(*,*) minval(ks), maxval(ks), k(ns-1)
    !stop

    call splie2_full_precomp(x_tot, ks, S_lores, S_coeff)
    !write(*,*) S_coeff
    !stop

    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

    do j = 1, ns
       do i = 1,ns 
           S(i,j) = splin2_full_precomp(x_tot, ks, S_coeff, x(i), k(j))

           

       end do
       !write(*,*) S(:,j)
       !stop
    end do

    !write(*,*) S(:,8)
    !stop


  end subroutine get_hires_source_function



 
  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    real(dp)     :: dk, dydx, dx, xstop, xstart, h1, dx_br, H_p, dtau
    integer(i4b) :: l, i, n

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))

    do i=1, n_k
       ks(i) = k_min + (k_max-k_min)*((i-1)/99.d0)**2
       
    end do
    
    
    
    

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_tot, 0:lmax_int, n_k))
    allocate(delta(0:n_tot, n_k))
    allocate(delta_b(0:n_tot, n_k))
    allocate(v(0:n_tot, n_k))
    allocate(v_b(0:n_tot, n_k))
    allocate(Phi(0:n_tot, n_k))
    allocate(Psi(0:n_tot, n_k))
    allocate(dPhi(0:n_tot, n_k))
    allocate(dPsi(0:n_tot, n_k))
    allocate(dv_b(0:n_tot, n_k))
    allocate(dTheta(0:n_tot, 0:lmax_int, n_k))
    
    allocate(x_br(n_br))
    allocate(x_tot(n_tot))
    
    dx_br        = (x_t(1)-0.1d0 - log(a_init))/ (n_br -1)
    x_br(1:n_br) = [(log(a_init) + ((i-1)*dx_br), i=1, n_br)] 

    x_tot(1:n_br) = x_br
    x_tot(n_br+1:n_tot) = x_t


    H_p = get_H_p(x_br(1))
    dtau = get_dtau(x_br(1))

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    delta(0,:)   = 1.5d0*Phi(0,:)
    delta_b(0,:) = 1.5d0*Phi(0,:)
       

    do i = 1, n_k
       v(0,i)       = c*ks(i)/(2.d0*H_p)*Phi(0,i)
       v_b(0,i)     = c*ks(i)/(2.d0*H_p)*Phi(0,i)
       Theta(0,0,i) = 0.5d0*Phi(0,i)
       Theta(0,1,i) = -c*ks(i)/(6.d0*H_p)*Phi(0,i)
       Theta(0,2,i) = -20.d0*c*ks(i)/(45.d0*H_p*dtau)*Theta(0,1,i)
       do l = 3, lmax_int
          Theta(0,l,i) = -l/(2.d0*l + 1.d0)*c*ks(i)/(H_p*dtau)*Theta(0,l-1,i)
       end do
       Psi(0,i) = -1.d0
    end do
    
    
  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l, i_tc
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, H_p, dt, t1, t2, x_tc, a_pe, D

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k
       write(*,*) k
       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5
       

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)

       
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       
  


       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       do i=2, n_tot
          if (x_tot(i) < x_tc) then
          
             call odeint(y_tight_coupling, x_tot(i-1), x_tot(i), eps, h1, hmin, derivs_tc, bsstep, output_tc)
             delta(i,k)   = y_tight_coupling(1)
             delta_b(i,k) = y_tight_coupling(2)
             v(i,k)       = y_tight_coupling(3)
             v_b(i,k)     = y_tight_coupling(4)
             Phi(i,k)     = y_tight_coupling(5)
             Theta(i,0,k) = y_tight_coupling(6)
             Theta(i,1,k) = y_tight_coupling(7)
             Theta(i,2,k) = -20.d0/(45.d0*get_dtau(x_tot(i)))*c*k_current/get_H_p(x_tot(i))*Theta(i,1,k)
             do l = 3, lmax_int
                Theta(i,l,k) = -l/(2.d0*l + 1.d0)*c*k_current/(get_H_p(x_tot(i))*get_dtau(x_tot(i)))*Theta(i,l-1,k)
             end do
             
             a_pe = exp(x_tot(i))
             Psi(i,k)     =  -y_tight_coupling(5)-12.d0*H_0**2/(c**2*k_current**2*a_pe**2)*(Omega_r*Theta(i,2,k))
          

          else 
             i_tc = i
             exit
          end if
       end do



       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(i_tc-1,2,k)
       do l = 3, lmax_int
          y(6+l) = Theta(i_tc-1,l,k) 
       end do


       do i = i_tc, n_tot

          
          ! Task: Integrate equations from tight coupling to today

          call odeint(y, x_tot(i-1), x_tot(i), eps, h1, hmin, derivs_tot, bsstep, output_tot)


          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do
          a_pe = exp(x_tot(i))
          D = c*k_current/get_H_p(x_tot(i))
          Psi(i,k)     =  -y(5)-12.d0*H_0**2/(c**2*k_current**2*a_pe**2)*(Omega_r*y(8))    
          
          ! Task: Store derivatives that are required for C_l estimation
          dPhi(i,k)     = Psi(i,k) - D**2/3.d0*y(5) + 0.5d0*H_0**2/get_H_p(x_tot(i))**2*(Omega_m/a_pe*y(1) + Omega_b/a_pe*y(2) +4.d0*Omega_r/a_pe**2*y(6)) 
          dv_b(i,k)     = -y(4) - D*Psi(i,k) + get_dtau(x_tot(i))*(4.d0*Omega_r/(3.d0*Omega_b*a_pe))*(3.d0*y(7)+y(4))
          dTheta(i,0,k) = -D*y(7) - dPhi(i,k)
          dTheta(i,1,k) = D/3.d0*y(6) - 2.d0*D/3.d0*y(8) + D/3.d0*Psi(i,k) + get_dtau(x_tot(i))*(y(7)+1.d0/3.d0*y(4))
          dTheta(i,2,k) = 2.d0/(4.d0 + 1.d0)*D*y(7) - (2.d0+1.d0)/(4.d0 + 1.d0)*D*y(9) + get_dtau(x_tot(i))*(y(8)-1.d0/10.d0*y(8))
          do l=3,lmax_int-1
             dTheta(i,l,k) = l/(2.d0*l + 1.d0)*D*Theta(i,l-1,k) - (l+1.d0)/(2.d0*l + 1.d0)*D*Theta(i,l+1,k) + get_dtau(x_tot(i))*Theta(i,l,k)
          end do
          dTheta(i,6,k) = D*y(11) - D/k_current*(12.d0+1.d0)/get_eta(x_tot(i))*y(12) + get_dtau(x_tot(i))*y(12)
          dPsi(i,k)     = - dPhi(i,k) - dTheta(i,2,k)*12.d0*H_0**2/(c**2*k_current**2*a_pe**2)*Omega_r + 24.d0*H_0**2/(c**2*k_current**2*a_pe**2)*Omega_r*Theta(i,2,l)
          
       end do
       
       

    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time, x, dx
    integer(i4b)          :: i

    
    dx = 1.d-3
    x = log(a_init)

    do while (abs(get_dtau(x)) > 10.d0 .and. x < x_t(1) .and. abs(c*k/(get_H_p(x)*get_dtau(x))) < 0.1d0)
       x = x + dx
    end do
    get_tight_coupling_time = min(x, x_t(1))

  end function get_tight_coupling_time

  subroutine derivs_tc(x_tot, y_tight_coupling, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_tot
    real(dp), dimension(:), intent(in)  :: y_tight_coupling
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: delta, ddelta, delta_b, ddelta_b, v, dv, v_b, dv_b, Phi, dPhi, Theta_0, dTheta_0, Theta_1, dTheta_1, Theta_2, dTheta_2,Theta_l, dTheta_l, Psi, q, R, a, k, x, D, B 
    x = x_tot
    k = k_current
    a = exp(x)

    D = c*k_current/get_H_p(x)
    B = get_dH_p(x)/get_H_p(x)
    !write(*,*) x
    !write(*,*) y_tight_coupling
    !stop

    delta     = y_tight_coupling(1)
    delta_b   = y_tight_coupling(2)
    v         = y_tight_coupling(3)
    v_b       = y_tight_coupling(4)
    Phi       = y_tight_coupling(5)
    Theta_0   = y_tight_coupling(6)
    Theta_1   = y_tight_coupling(7)
    Theta_2   = -20.d0/(45.d0*get_dtau(x))*D*Theta_1
    
    
    R = 4.d0*Omega_r/(3.d0*Omega_b*a)

    Psi = -Phi-12.d0*H_0**2/(c*k_current*a)**2*(Omega_r*Theta_2)    

    dydx(5) = Psi - D**2/3.d0*Phi + 0.5d0*H_0**2/get_H_p(x)**2*(Omega_m/a*delta + Omega_b/a*delta_b +4.d0*Omega_r/a**2*Theta_0)
    !write(*,*) dydx(5)

    dydx(3) = -v-D*Psi
    !write(*,*) dydx(3)

    dydx(2) = D*v_b - 3.d0*dydx(5)

    dydx(1) = D*v - 3.d0*dydx(5)

    dydx(6) = -D*Theta_1 - dydx(5)

    q = -(((1.d0-2.d0*R)*get_dtau(x) + (1.d0+R)*get_ddtau(x))*(3.d0*Theta_1+v_b) - (D*Psi) + (1.d0-B)*D*(-Theta_0 + 2.d0*Theta_2) - D*dydx(6))/((1.d0+R)*get_dtau(x) + B -1.d0)

    dydx(4) = 1.d0/(1.d0+R)*(-v_b-D*Psi + R*(q+D*(-Theta_0 + 2.d0*Theta_2) - D*Psi))

    dydx(7) = 1.d0/3.d0*(q - dydx(4))

    !write(*,*) dydx(5)
    !write(*,*) dydx(3)
    !write(*,*) dydx(2)
    !write(*,*) dydx(1)
    !write(*,*) dydx(6)
    !write(*,*) dydx(4)
    !write(*,*) dydx(7)
    !write(*,*) 'k=',k_current
    !write(*,*) Psi
    !write(*,*) H_0, Omega_r, Theta_2
    !write(*,*) 'dtau=',get_dtau(x)
    !stop    
    
  end subroutine derivs_tc
  
  
  subroutine output_tc(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output_tc

  subroutine derivs_tot(x_tot, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x_tot
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp)                            :: a, R, Psi, k, D, B, x
    integer(i4b)                        :: l

    x = x_tot
    k = k_current
    a = exp(x)
    R = 4.d0*Omega_r/(3.d0*Omega_b*a)
    
    D = c*k_current/get_H_p(x)
    
    Psi = -y(5)-12.d0*H_0**2/(c*k_current*a)**2*(Omega_r*y(8))    

    dydx(5) = Psi - D**2/3.d0*y(5) + 0.5d0*H_0**2/get_H_p(x)**2*(Omega_m/a*y(1) + Omega_b/a*y(2) +4.d0*Omega_r/a**2*y(6))

    dydx(1) = D*y(3) - 3.d0*dydx(5)

    dydx(2) = D*y(4) - 3.d0*dydx(5)

    dydx(3) = -y(3)-D*Psi       

    dydx(4) = -y(4) - D*Psi + get_dtau(x)*R*(3.d0*y(7)+y(4))

    dydx(6) = -D*y(7) - dydx(5)

    dydx(7) = D/3.d0*y(6) - 2.d0*D/3.d0*y(8) + D/3.d0*Psi + get_dtau(x)*(y(7)+1.d0/3.d0*y(4))  

    dydx(8) = 2.d0/(4.d0 + 1.d0)*D*y(7) - (2.d0+1.d0)/(4.d0 + 1.d0)*D*y(9) + get_dtau(x)*(y(8)-1.d0/10.d0*y(8))  

    !dydx(9) = 3.d0/(2.d0*3.d0 + 1.d0)*D*y(8) - (3.d0+1.d0)/(2.d0*3.d0 + 1.d0)*D*y(10) + get_dtau(x)*y(9) 

    !dydx(10) = 4.d0/(2.d0*4.d0 + 1.d0)*D*y(9) - (4.d0+1.d0)/(2.d0*4.d0 + 1.d0)*D*y(11) + get_dtau(x)*y(10)

    !dydx(11) = 5.d0/(2.d0*5.d0 + 1.d0)*D*y(10) - (5.d0+1.d0)/(2.d0*5.d0 + 1.d0)*D*y(12) + get_dtau(x)*y(11)


    do l=3,lmax_int-1
       dydx(l+6) = l/(2.d0*l + 1.d0)*D*y(l+5) - (l+1.d0)/(2.d0*l + 1.d0)*D*y(l+7) + get_dtau(x)*y(l+6)
    end do


    dydx(12) = D*y(11) - D/k_current*(12.d0+1.d0)/get_eta(x)*y(12) + get_dtau(x)*y(12) !dTheta_6


  end subroutine derivs_tot
  
  subroutine output_tot(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output_tot

end module evolution_mod
