program benchmark_twostream
  use twostream, only: dp
  use twostream, only: two_stream_solar, TwoStreamSolarWrk
  use twostream, only: two_stream_ir, TwoStreamIRWrk
  use twostream_test_utils, only: TwostreamTestData
  implicit none

  call benchmark_twostream_solar_noalloc()
  call benchmark_twostream_solar_alloc()

  print*,''
  call benchmark_twostream_ir_noalloc()
  call benchmark_twostream_ir_alloc()

contains
  subroutine benchmark_twostream_solar_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamSolarWrk) :: wrk

    integer, parameter :: nt = 100000
    integer :: i
    real(dp) :: t(10)

    call d%standard_case()
    wrk = TwoStreamSolarWrk(d%nz)

    call cpu_time(t(1))
    do i=1,nt
      call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn, wrk)
    enddo
    call cpu_time(t(2))

    print*,'benchmark_twostream_solar_noalloc'
    print*,'fup(1) = ',d%fup(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

  subroutine benchmark_twostream_solar_alloc()
    type(TwostreamTestData) :: d

    integer, parameter :: nt = 100000
    integer :: i
    real(dp) :: t(10)

    call d%standard_case()

    call cpu_time(t(1))
    do i=1,nt
      call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn)
    enddo
    call cpu_time(t(2))

    print*,'benchmark_twostream_solar_alloc'
    print*,'fup(1) = ',d%fup(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

  subroutine benchmark_twostream_ir_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamIRWrk) :: wrk

    integer, parameter :: nt = 100000
    integer :: i
    real(dp) :: t(10)

    call d%standard_case()
    wrk = TwoStreamIRWrk(d%nz)

    call cpu_time(t(1))
    do i=1,nt
      call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn, wrk)
    enddo
    call cpu_time(t(2))

    print*,'benchmark_twostream_ir_noalloc'
    print*,'fup(1) = ',d%fup(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

  subroutine benchmark_twostream_ir_alloc()
    type(TwostreamTestData) :: d

    integer, parameter :: nt = 100000
    integer :: i
    real(dp) :: t(10)

    call d%standard_case()

    call cpu_time(t(1))
    do i=1,nt
      call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn)
    enddo
    call cpu_time(t(2))

    print*,'benchmark_twostream_ir_alloc'
    print*,'fup(1) = ',d%fup(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

end program