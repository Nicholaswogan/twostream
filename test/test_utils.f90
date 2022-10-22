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
    real(dp), allocatable :: bplanck(:)

    ! output
    real(dp) :: surface_radiance
    real(dp), allocatable :: amean(:), fup(:), fdn(:)
  contains
    procedure :: allocate => TwostreamTestData_alloc
    procedure :: standard_case
  end type

contains

  subroutine TwostreamTestData_alloc(self, nz)
    class(TwostreamTestData), intent(inout) :: self
    integer, intent(in) :: nz

    self%nz = nz
    if (allocated(self%tau)) then
      deallocate(self%tau,self%w0,self%gt,self%bplanck)
      deallocate(self%amean,self%fup,self%fdn)
    endif
    allocate(self%tau(nz), self%w0(nz), self%gt(nz),self%bplanck(nz+1))
    allocate(self%amean(nz+1), self%fup(nz+1), self%fdn(nz+1))

  end subroutine

  subroutine standard_case(self)
    class(TwostreamTestData), intent(inout) :: self
    integer :: i
    real(dp) :: nu, T

    call TwostreamTestData_alloc(self, 200)
    self%u0 = 0.6427876096865394_dp
    self%Rsfc = 0.25_dp
    self%gt(:) = 0.0_dp
    open(2,file='../test/tau_and_w0.txt',status='old')
    do i = 1,self%nz
      read(2,*) self%tau(i), self%w0(i)
    enddo
    close(2)

    ! planck
    nu = 29979245799999.996_dp ! frequency (1/s), corresponding to 10 um
    T = 300.0_dp ! K
    self%bplanck(self%nz+1) = planck_fcn(nu, T)
    do i = 1,self%nz
      self%bplanck(i) = planck_fcn(nu, T)
    enddo

  end subroutine

  function planck_fcn(nu, T) result(B)
    real(dp), intent(in) :: nu ! (1/s) 
    real(dp), intent(in) :: T ! (K)
    real(dp) :: B ! mW sr^−1 m^−2 Hz^-1

    real(dp), parameter :: c_light = 299792458.0_dp ! Speed of light (m / s)
    real(dp), parameter :: k_boltz = 1.380649e-23_dp ! boltzmann's constant si units (J/K)
    real(dp), parameter :: plank = 6.62607004e-34_dp ! planks constant (m2 kg / s)
    
    B = 1.0e3_dp*((2.0_dp*plank*nu**3.0_dp)/(c_light**2.0_dp)) * &
        ((1.0_dp)/(exp((plank*nu)/(k_boltz*T)) - 1.0_dp)) 
    
  end function

end module

