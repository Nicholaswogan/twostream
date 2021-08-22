program main
  use twostream, only: two_stream
  implicit none
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: nz = 200
  real(real_kind) :: tau(nz)
  real(real_kind) :: w0(nz), u0, Rsfc, amean(nz+1), surface_radiance
  integer :: ierr, i, nt
  real :: end, start

  u0 = 0.6427876096865394d0
  Rsfc = 0.25d0
  
  open(2,file='tau_and_w0.txt',status='old')
  do i = 1,nz
    read(2,*) tau(i), w0(i)
  enddo
  close(2)
  
  call two_stream(nz, tau, w0, u0, Rsfc, amean, surface_radiance,ierr)
  open(2,file='results/amean_fortran.txt',status='replace')
  do i = 1,nz+1
    write(2,*) amean(i)
  enddo
  close(2)
    
  nt = 100000
  call cpu_time(start)
  do i=1,nt
    call two_stream(nz, tau, w0, u0, Rsfc, amean, surface_radiance,ierr)
  enddo
  call cpu_time(end)
  
  print*,'One run in Fortran takes',(end-start)/nt,'seconds'

end program