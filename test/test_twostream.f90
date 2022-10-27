program test_twostream
  use twostream, only: dp
  use twostream, only: two_stream_solar, TwoStreamSolarWrk
  use twostream, only: two_stream_ir, TwoStreamIRWrk
  use twostream_test_utils, only: TwostreamTestData, is_close
  implicit none

  call test_twostream_solar_noalloc()
  call test_twostream_solar_alloc()
  call test_twostream_ir_noalloc()
  call test_twostream_ir_alloc()

contains

  !~~~~~~~~~~~~~~~~!
  !~~ Unit Tests ~~!
  !~~~~~~~~~~~~~~~~!

  subroutine test_twostream_solar_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamSolarWrk) :: wrk

    call d%standard_case()
    wrk = TwoStreamSolarWrk(d%nz)
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn, wrk)
    if (.not. is_close(d%fup(1), 3.0752040848372568e-2_dp, tol=1.0e-8_dp)) then
      error stop '"test_twostream_solar_noalloc" failed'
    endif
    print*,'"test_twostream_solar_noalloc" passed'
  end subroutine

  subroutine test_twostream_solar_alloc()
    type(TwostreamTestData) :: d

    call d%standard_case()
    call two_stream_solar(d%nz, d%tau, d%w0, d%gt, d%u0, d%Rsfc, d%amean, d%surface_radiance, d%fup, d%fdn)
    if (.not. is_close(d%fup(1), 3.0752040848372568e-2_dp, tol=1.0e-8_dp)) then
      error stop '"test_twostream_solar_alloc" failed'
    endif
    print*,'"test_twostream_solar_alloc" passed'
  end subroutine

  subroutine test_twostream_ir_noalloc()
    type(TwostreamTestData) :: d
    type(TwoStreamIRWrk) :: wrk

    call d%standard_case()
    wrk = TwoStreamIRWrk(d%nz)
    call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn, wrk)
    if (.not. is_close(d%fup(1), 9.8521506723050778e-9_dp, tol=1.0e-8_dp)) then
      error stop '"test_twostream_ir_noalloc" failed'
    endif
    print*,'"test_twostream_ir_noalloc" passed'
  end subroutine
  
  subroutine test_twostream_ir_alloc()
    type(TwostreamTestData) :: d

    call d%standard_case()
    call two_stream_ir(d%nz, d%tau, d%w0, d%gt, 0.0_dp, d%bplanck, d%fup, d%fdn)
    if (.not. is_close(d%fup(1), 9.8521506723050778e-9_dp, tol=1.0e-8_dp)) then
      error stop '"test_twostream_ir_alloc" failed'
    endif
    print*,'"test_twostream_ir_alloc" passed'
  end subroutine

end program