program benchmark_twostream
  use twostream, only: two_stream_solar, TwoStreamSolarWrk, dp
  use test_utils, only: TwostreamTestData
  implicit none

  call benchmark_twostream_noalloc()
  call benchmark_twostream_alloc()

contains
  subroutine benchmark_twostream_noalloc()
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

    print*,'benchmark_twostream_noalloc'
    print*,'amean(1) = ',d%amean(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

  subroutine benchmark_twostream_alloc()
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

    print*,'benchmark_twostream_alloc'
    print*,'amean(1) = ',d%amean(1)
    print*,'Time / run = ',(t(2)-t(1))/real(nt,dp),'seconds'
  end subroutine

end program