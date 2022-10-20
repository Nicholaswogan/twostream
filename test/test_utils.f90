module test_utils
  use twostream, only: dp
  implicit none
  private
  public :: TwostreamTestData

  !! This module it to help out the testing and benchmarking functions

  type :: TwostreamTestData
    ! number of layers
    integer :: nz

    ! input
    real(dp) :: u0, Rsfc
    real(dp), allocatable :: tau(:), w0(:), gt(:)

    ! output
    real(dp) :: surface_radiance
    real(dp), allocatable :: amean(:), fup(:), fdn(:)
  contains
    procedure :: allocate => TwostreamTestData_alloc
    procedure :: standard_case
  end type
  interface TwostreamTestData
    procedure :: TwostreamTestData_alloc_1
  end interface

contains

  subroutine TwostreamTestData_alloc(self, nz)
    class(TwostreamTestData), intent(inout) :: self
    integer, intent(in) :: nz

    self%nz = nz
    if (allocated(self%tau)) then
      deallocate(self%tau,self%w0,self%gt)
      deallocate(self%amean,self%fup,self%fdn)
    endif
    allocate(self%tau(nz), self%w0(nz), self%gt(nz))
    allocate(self%amean(nz+1), self%fup(nz+1), self%fdn(nz+1))

  end subroutine

  function TwostreamTestData_alloc_1(nz) result(d)
    integer, intent(in) :: nz
    type(TwostreamTestData) :: d
    call TwostreamTestData_alloc(d, nz)
  end function

  subroutine standard_case(self)
    class(TwostreamTestData), intent(inout) :: self
    integer :: i

    call TwostreamTestData_alloc(self, 200)
    self%u0 = 0.6427876096865394_dp
    self%Rsfc = 0.25_dp
    self%gt(:) = 0.0_dp
    open(2,file='../test/tau_and_w0.txt',status='old')
    do i = 1,self%nz
      read(2,*) self%tau(i), self%w0(i)
    enddo
    close(2)
  end subroutine

end module

