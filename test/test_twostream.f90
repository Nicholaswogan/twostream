program test_twostream
  use twostream, only: dp
  use twostream, only: two_stream_solar, TwoStreamSolarWrk
  use twostream, only: two_stream_ir, TwoStreamIRWrk
  use test_utils, only: TwostreamTestData
  implicit none

  write(*,"(A)") ""
  write(*,"(A)") "~~~~~~~~~~~~~~~~~~~~"
  write(*,"(A)") "~~~ Memory Tests ~~~"
  write(*,"(A)") "~~~~~~~~~~~~~~~~~~~~"
  write(*,"(A)") ""
  call test_twostream_solar_noalloc()
  call test_twostream_solar_alloc()
  call test_twostream_ir_noalloc()
  call test_twostream_ir_alloc()

  write(*,"(A)") ""
  write(*,"(A)") "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,"(A)") "~~~ Reproducing Toon et al. (1989) ~~~"
  write(*,"(A)") "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  call test_twostream_solar_table_6()

contains

  !~~~~~~~~~~~~~~~~~~!
  !~~ Memory Tests ~~!
  !~~~~~~~~~~~~~~~~~~!

  subroutine test_twostream_solar_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamSolarWrk) :: wrk

    call d%standard_case()
    wrk = TwoStreamSolarWrk(d%nz)
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn, wrk)
    print*,'test_twostream_solar_noalloc:'
    print*,'fup(1) = ',d%fup(1)
  end subroutine

  subroutine test_twostream_solar_alloc()
    type(TwostreamTestData) :: d

    call d%standard_case()
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn)
    print*,'test_twostream_solar_alloc:'
    print*,'fup(1) = ',d%fup(1)
  end subroutine

  subroutine test_twostream_ir_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamIRWrk) :: wrk

    call d%standard_case()
    wrk = TwoStreamIRWrk(d%nz)
    call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn, wrk)
    print*,'test_twostream_ir_noalloc:'
    print*,'fup(1) = ',d%fup(1)
  end subroutine
  
  subroutine test_twostream_ir_alloc()
    type(TwostreamTestData) :: d

    call d%standard_case()
    call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn)
    print*,'test_twostream_ir_alloc:'
    print*,'fup(1) = ',d%fup(1)
  end subroutine
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !~~ Reproducing Toon et al. (1989) ~~!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine test_twostream_solar_table_6()
    type(TwostreamTestData) :: d

    real(dp), parameter :: u0(*) = [0.1_dp, 0.1_dp, 0.1_dp, 0.4_dp, 0.4_dp, 0.4_dp, &
                                    0.92_dp, 0.92_dp, 0.92_dp]
    real(dp), parameter :: tau(*) = [0.02_dp, 0.25_dp, 1.0_dp, 0.02_dp, 0.25_dp, 1.0_dp, &
                                     0.02_dp, 0.25_dp, 1.0_dp]
    real(dp), parameter :: Rsfc(*) = [0.0_dp, 0.25_dp, 0.8_dp]

    character(15) :: tmp_str
    character(:), allocatable :: line, labels, under_labels
    real(dp) :: amean1(size(tau),size(Rsfc)), amean2(size(tau),size(Rsfc))
    integer :: i,j

    call d%allocate(1)

    do i = 1,size(Rsfc)
      do j = 1,size(tau)
        d%tau = tau(j)
        d%w0 = 0.9999999_dp
        d%gt = 0.0_dp
        d%u0 = u0(j)
        d%Rsfc = Rsfc(i)
        call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn)
        amean1(j,i) = d%amean(1)
        amean2(j,i) = d%amean(2)
      enddo
    enddo

    write(*,"(A)") ""
    write(*,"(A)") "**Table 6: Mean Intensities for Conservative Rayleigh Scatter (J(0)/(\pi*F_s))**"

    if (allocated(line)) deallocate(line)
    allocate(character(0)::line)
    tmp_str = "| u0"
    line = line//tmp_str
    tmp_str = "| tau"
    line = line//tmp_str
    tmp_str = "| Rsfc = 0"
    line = line//tmp_str
    tmp_str = "| Rsfc = 0.25"
    line = line//tmp_str
    tmp_str = "| Rsfc = 0.80"
    line = line//tmp_str
    line = line//"|"
    labels = line
    write(*,"(A)") line

    if (allocated(line)) deallocate(line)
    allocate(character(0)::line)
    tmp_str = "| ------------ "
    do i = 1,5
      line = line//tmp_str
    enddo
    line = line//"|"
    under_labels = line
    write(*,"(A)") line
    
    do j = 1,size(tau)
      if (allocated(line)) deallocate(line)
      allocate(character(0)::line)

      call round_format_string(u0(j), tmp_str)
      line = line//tmp_str
      
      call round_format_string(tau(j), tmp_str)
      line = line//tmp_str

      do i = 1,size(Rsfc)
        call round_format_string(amean1(j,i), tmp_str)
        line = line//tmp_str
      enddo
      line = line//"|"
      write(*,"(A)") line
    enddo


    write(*,"(A)") ""
    write(*,"(A)") "**Table 6: Mean Intensities for Conservative Rayleigh Scatter (J(\tau)/(\pi*F_s))**"
    write(*,"(A)") labels
    write(*,"(A)") under_labels
    do j = 1,size(tau)
      if (allocated(line)) deallocate(line)
      allocate(character(0)::line)

      call round_format_string(u0(j), tmp_str)
      line = line//tmp_str
      
      call round_format_string(tau(j), tmp_str)
      line = line//tmp_str

      do i = 1,size(Rsfc)
        call round_format_string(amean2(j,i), tmp_str)
        line = line//tmp_str
      enddo
      line = line//"|"
      write(*,"(A)") line
    enddo

  end subroutine

  subroutine round_format_string(a, b)
    real(dp), intent(in) :: a
    character(*), intent(out) :: b

    integer :: tmp_i
    real(dp) :: tmp_r

    ! rounding
    tmp_i = int(a*1000.0_dp)
    tmp_r = tmp_i/1000.0_dp

    ! formating
    write(b,"(f15.3)") tmp_r
    b = "| "//adjustl(b)

  end subroutine

end program