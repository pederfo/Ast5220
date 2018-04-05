program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none
  !integer i,k
  integer, parameter :: out_unit=20
  ! Initialize time grids
  call initialize_time_mod

  call initialize_rec_mod

  call initialize_perturbation_eqns
  call integrate_perturbation_eqns
  call compute_cls

  !call get_hires_source_function

  ! Output to file desired quantities here

  

  open (unit=out_unit,file="results_j.txt",action="write",status="replace")
  do i=1, ns! maxval(ls)
        
     write (out_unit,*) j_lfunc(ls(2), k(4000),x(i)), j_lfunc(ls(17), k(4000),x(i)), j_lfunc(ls(36), k(4000),x(i))!ls2(i), cl_spline(i)!x(i), S(i, 4000), j_lfunc(ls(17), k(4000),x(i)) !Theta_l(30,i), Theta_l(36,i), Theta_l(44,i) !

  end do

  close (out_unit)

  !write (out_unit,'(*(2X, ES14.6))')
  

end program cmbspec
