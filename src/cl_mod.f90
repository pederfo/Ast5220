module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  integer(i4b) :: i, j, l, l_num, x_num, n_spline

  integer(i4b), allocatable, dimension(:)       :: ls
  real(dp),     allocatable, dimension(:)       :: integrand1, integrand2, ls2
  real(dp),     allocatable, dimension(:,:)     :: j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp, C_l, cl_spline
  real(dp),     allocatable, dimension(:)       :: k, x
  !real(dp),     allocatable, dimension(:,:,:,:) :: S_coeff
  real(dp),     allocatable, dimension(:,:)     :: S, S2
  real(dp),     allocatable, dimension(:,:)     :: Theta_l
  real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires
  real(dp)                                      :: eta_0

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e, dz, a1, a2, h, np, h2, b2, b3
    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    allocate(ls_dp(l_num))
    allocate(ls2(1200))
    allocate(cl_spline(1200))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    do i = 1, l_num
       ls_dp(i) = ls(i)
    end do

    do i = 1, maxval(ls)
       ls2(i) = i
    end do

    ! Task: Get source function from evolution_mod

    call get_hires_source_function(k, x, S)

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    

    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    allocate(integrand1(ns))
    allocate(integrand2(ns))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(Theta_l(l_num, ns))
    allocate(C_l(l_num))

    
    dz      = (3500 - 0.d0)/ (n_spline -1.d0)
    z_spline(1:n_spline) = [(0.d0 + ((i-1)*dz), i=1, n_spline)]
 
    do i = 1, n_spline
       
       do l = 1, l_num
          if (z_spline(i) > 2.d0) then
          
             call sphbes( ls(l), z_spline(i), j_l(i,l))
             
          end if
          
       end do
    end do
    
    

    
    np = 0.96d0
    eta_0 = get_eta(0.d0)


    ! Overall task: Compute the C_l's for each given l
    do l = 1, l_num

       call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
       

       ! Task: Compute the transfer function, Theta_l(k)

       do j = 1, ns
          a1 = x(1)
          a2 = x(ns)
         h = (a2-a1)/(ns-1.d0)

          !write(*,*) a1, a2, h
          !stop

          do i = 1, ns
             integrand1(i) =  S(i,j)*j_lfunc(l, k(j), x(i))
             !write(*,*) j_lfunc(l, k(j), x(i)), S(i,j)
             !stop
          end do


         Theta_l(l,j) = 0.5d0*(integrand1(1) + integrand1(ns))
          
          do i = 2, ns-1
             Theta_l(l,j) = Theta_l(l,j) + integrand1(i)
          end do

          Theta_l(l,j) = Theta_l(l,j)*h
          
          !write(*,*) Theta_l(l,j)
          !stop

       end do


       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's

       h2 = (k(ns) - k(1)) / (ns-1.d0)

       do j=1, ns
          integrand2(j) = (c*k(j)/H_0)**(np-1.d0)*Theta_l(ns, j)**2/k(j)
       end do

       C_l(l) = 0.5d0*(integrand2(1) + integrand2(ns))

       do j = 2, ns-1
          C_l(l) = C_l(l) + integrand2(j)
       end do
       !write(*,*) C_l
       !C_l(l) = C_l(l)*h2

       ! Task: Store C_l in an array. Optionally output to file

       cls(l) = C_l(l)*l*(l+1)/(2.d0*pi)*h2

       

       !write(*,*) clk
       !stop

    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    call spline(ls_dp, cls, 1.d30, 1.d30, cls2)

    do i = 1, maxval(ls)
       cl_spline(i) = splint(ls_dp, cls, cls2, ls2(i))
    end do
    


  end subroutine compute_cls
  
  function j_lfunc(l,k,x)
    implicit none
    integer(i4b), intent(in) :: l
    real(dp),     intent(in) :: k
    real(dp),     intent(in) :: x
    real(dp)             :: j_lfunc

    j_lfunc = splint(z_spline, j_l(:,l), j_l2(:,l), k*(get_eta(0.d0) - get_eta(x)))
  end function j_lfunc


end module cl_mod
