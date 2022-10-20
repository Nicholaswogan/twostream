program test_twostream
  use twostream, only: two_stream_solar, TwoStreamSolarWrk, dp
  use test_utils, only: TwostreamTestData
  implicit none

  call test_twostream_noalloc()
  call test_twostream_alloc()

contains

  subroutine test_twostream_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamSolarWrk) :: wrk

    call d%standard_case()
    wrk = TwoStreamSolarWrk(d%nz)
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn, wrk)
    print*,'test_twostream_noalloc:'
    print*,'amean(1) = ',d%amean(1)
  end subroutine

  subroutine test_twostream_alloc()
    type(TwostreamTestData) :: d

    call d%standard_case()
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn)
    print*,'test_twostream_alloc:'
    print*,'amean(1) = ',d%amean(1)
  end subroutine

end program