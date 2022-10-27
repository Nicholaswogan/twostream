
module twostream
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: dp
  public :: two_stream_solar, TwoStreamSolarWrk
  public :: two_stream_ir, TwoStreamIRWrk

  type :: TwoStreamSolarWrk
    real(dp), allocatable :: e1(:), e2(:), e3(:), e4(:)
    real(dp), allocatable :: tauc(:), direct(:)
    real(dp), allocatable :: cp0(:), cpb(:), cm0(:), cmb(:)
    real(dp), allocatable :: A(:), B(:), D(:), E(:)
    real(dp), allocatable :: y1(:), y2(:)
  end type
  interface TwoStreamSolarWrk
    procedure :: TwoStreamSolarWrk_create
  end interface

  type :: TwoStreamIRWrk
    real(dp), allocatable :: gam1(:), gam2(:)
    real(dp), allocatable :: lambda(:), cap_gam(:)
    real(dp), allocatable :: e1(:), e2(:), e3(:), e4(:)
    real(dp), allocatable :: cp0(:), cpb(:), cm0(:), cmb(:)
    real(dp), allocatable :: A(:), B(:), D(:), E(:)
    real(dp), allocatable :: y1(:), y2(:)
  end type
  interface TwoStreamIRWrk
    procedure :: TwoStreamIRWrk_create
  end interface

contains
  
  subroutine two_stream_solar(nz, tau, w0, gt, u0, Rsfc, &
                              amean, surface_radiance, fup, fdn, wrk)
    integer, intent(in) :: nz
    real(dp), intent(in) :: tau(nz)
    real(dp), intent(in) :: w0(nz)
    real(dp), intent(in) :: gt(nz)
    real(dp), intent(in) :: u0, Rsfc
    real(dp), intent(out) :: amean(nz+1)
    real(dp), intent(out) :: surface_radiance
    real(dp), intent(out) :: fup(nz+1), fdn(nz+1)
    type(TwoStreamSolarWrk), optional, target, intent(inout) :: wrk
    
    ! local
    real(dp), pointer :: e1(:), e2(:), e3(:), e4(:)
    real(dp), pointer :: tauc(:), direct(:)
    real(dp), pointer :: cp0(:), cpb(:), cm0(:), cmb(:)
    real(dp), pointer :: A(:), B(:), D(:), E(:)
    real(dp), pointer :: y1(:), y2(:)
    ! pointer
    type(TwoStreamSolarWrk), pointer :: w
    logical :: deallocate_w
    
    integer :: i, l, ll
    real(dp) :: gam1, gam2, gam3, gam4
    real(dp) :: lambda, cap_gam
    real(dp) :: wrk_real, facp, facm, et0, etb, denom, Ssfc
    
    real(dp), parameter :: u1 = 1.0_dp/sqrt(3.0_dp) ! (Quadrature)
    real(dp), parameter :: Fs_pi = 1.0_dp

    if (present(wrk)) then
      ! we use memory in wrk
      w => wrk
      deallocate_w = .false.
    else
      ! we allocate memory for the problem
      allocate(w)
      w = TwoStreamSolarWrk(nz)
      deallocate_w = .true.
    endif

    ! Pointers to work memory
    e1 => w%e1
    e2 => w%e2
    e3 => w%e3
    e4 => w%e4
    tauc => w%tauc
    direct => w%direct
    cp0 => w%cp0
    cpb => w%cpb
    cm0 => w%cm0
    cmb => w%cmb
    A => w%A
    B => w%B
    D => w%D
    E => w%E
    y1 => w%y1
    y2 => w%y2

    ! tauc - cumulative optical depth at the top of each layer
    tauc(1) = 0.0_dp
    do i = 2,nz+1
      tauc(i) = tauc(i-1) + tau(i-1)
    enddo

    do i = 1,nz

      ! Quadrature Two-Stream coefficients (Table 1)
      gam1 = sqrt(3.0_dp)*(2.0_dp-w0(i)*(1+gt(i)))/2.0_dp
      gam2 = sqrt(3.0_dp)*w0(i)*(1.0_dp-gt(i))/2.0_dp
      gam3 = (1.0_dp-sqrt(3.0_dp)*gt(i)*u0)/2.0_dp
      gam4 = 1.0_dp - gam3
      ! u1 = 1.0_dp/sqrt(3.0_dp) (parameter)

      ! lambda, and capital Gamma (Equations 21, 22)
      lambda = sqrt(gam1**2.0_dp - gam2**2.0_dp)
      cap_gam = gam2 / (gam1 + lambda)
    
      ! e's (Equation 44)
      wrk_real = exp(-lambda*tau(i))
      e1(i) = 1.0_dp + cap_gam*wrk_real
      e2(i) = 1.0_dp - cap_gam*wrk_real
      e3(i) = cap_gam + wrk_real
      e4(i) = cap_gam - wrk_real
      
      ! C+ and C- (Equation 23, 24)
      ! Fs_pi = 1.0_dp (parameter) 
      ! We take the solar flux at the top of the atmosphere to be = 1
      ! Note: Toon et al. 1989 defines Fs * pi = solar flux. So Fs = solar flux / pi
      facp = w0(i)*Fs_pi*((gam1-1.0_dp/u0)*gam3+gam4*gam2)
      facm = w0(i)*Fs_pi*((gam1+1.0_dp/u0)*gam4+gam2*gam3) 
      et0 = exp(-tauc(i)/u0)
      etb = et0*exp(-tau(i)/u0)
      denom = lambda**2.0_dp - 1.0_dp/(u0**2.0_dp)

      direct(i+1) = u0*Fs_pi*etb
      cp0(i) = et0*facp/denom
      cpb(i) = etb*facp/denom
      cm0(i) = et0*facm/denom
      cmb(i) = etb*facm/denom

    enddo

    direct(1) = u0*Fs_pi
    Ssfc = Rsfc*direct(nz+1)
      
    ! Coefficients of tridiagonal linear system (Equations 39 - 43)
    A(1) = 0.0_dp
    B(1) = e1(1)
    D(1) = -e2(1)
    E(1) = 0.0_dp - cm0(1) ! assumes no downward diffuse flux at top-of-atmosphere
    do i = 1, nz-1
      ! Odd coeficients (Equation 41)
      l = 2*i + 1 ! is this right?
      A(l) = e2(i)*e3(i) - e4(i)*e1(i)
      B(l) = e1(i)*e1(i+1) - e3(i)*e3(i+1)
      D(l) = e3(i)*e4(i+1) - e1(i)*e2(i+1)
      E(l) = e3(i)*(cp0(i+1) - cpb(i)) + e1(i)*(cmb(i) - cm0(i+1))

      ! Even coefficients (Equation 42)
      ll = 2*i ! is this right?
      A(ll) = e2(i+1)*e1(i) - e3(i)*e4(i+1)
      B(ll) = e2(i)*e2(i+1) - e4(i)*e4(i+1)
      D(ll) = e1(i+1)*e4(i+1) - e2(i+1)*e3(i+1)
      E(ll) = e2(i+1)*(cp0(i+1) - cpb(i)) - e4(i+1)*(cm0(i+1) - cmb(i))
    enddo
    l = 2*nz
    A(l) = e1(nz) - Rsfc*e3(nz)
    B(l) = e2(nz) - Rsfc*e4(nz)
    D(l) = 0.0_dp
    E(l) = Ssfc - cpb(nz) + Rsfc*cmb(nz)
    
    ! Solve tridiagonal system. e is solution
    call tridiag(nz*2,a,b,d,e)
    
    ! unpack solution
    do i = 1, nz
      l = 2*i
      y1(i) = e(l-1)
      y2(i) = e(l)
    enddo

    ! amean = integral(J_n d_Omega) = J_n*4*pi
    ! Above is an integration of J_n over a complete sphere, 
    ! and is thus has units of flux density, or irradiance.
    ! amean(i) is for at top of layer i. amean(nz + 1) is the ground.
    
    ! very top edge of atmosphere (Not in paper. Derive from Equation 17 and 31. bit confusing)
    amean(1) = (1.0_dp/u1)*(y1(1)*e3(1)-y2(1)*e4(1)+cp0(1)) + direct(1)/u0
    do i = 1, nz
      ! J_n*4*pi = amean (Equation 49)
      amean(i+1) = (1.0_dp/u1)*(y1(i)*(e1(i)+e3(i))+y2(i)*(e2(i)+e4(i)) &
                   + cpb(i)+cmb(i)) + direct(i+1)/u0
    enddo

    ! Equations 31 and 32. But direct flux added
    fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
    fdn(1) = direct(1)
    do i = 1 , NZ
       fup(i+1) = (y1(i)*e1(i)+y2(i)*e2(i)+cpb(i))
       fdn(i+1) = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i)) + direct(i+1)
    enddo
    
    ! surface radiance (Ranjan and Sasselov 2017). photons actually hitting the ground (i think)
    i = nz
    surface_radiance = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i))/u1 + exp(-tauc(i+1)/u0)

    if (deallocate_w) then
      deallocate(w)
    endif

  end subroutine
  
  subroutine two_stream_ir(nz, tau, w0, gt, Rsfc, bplanck, &
                           fup, fdn, wrk)
    integer, intent(in) :: nz
    real(dp), intent(in) :: tau(nz)
    real(dp), intent(in) :: w0(nz)
    real(dp), intent(in) :: gt(nz)
    real(dp), intent(in) :: Rsfc
    real(dp), intent(in) :: bplanck(nz+1)
    real(dp), intent(out) :: fup(nz+1), fdn(nz+1)
    type(TwoStreamIRWrk), optional, target, intent(inout) :: wrk
    
    ! local
    real(dp), pointer :: gam1(:), gam2(:)
    real(dp), pointer :: lambda(:), cap_gam(:)
    real(dp), pointer :: e1(:), e2(:), e3(:), e4(:)
    real(dp), pointer :: cp0(:), cpb(:), cm0(:), cmb(:)
    real(dp), pointer :: A(:), B(:), D(:), E(:)
    real(dp), pointer :: y1(:), y2(:)
    ! pointer
    type(TwoStreamIRWrk), pointer :: w
    logical :: deallocate_w
    
    integer :: i, l
    real(dp) :: wrk_real, Ssfc
    real(dp) :: b0n, b1n
    
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: u1 = 0.5_dp ! (Hemispheric mean)
    real(dp), parameter :: emissivity = 1.0_dp
    real(dp), parameter :: norm = 2.0_dp*pi*u1

    if (present(wrk)) then
      ! we use memory in wrk
      w => wrk
      deallocate_w = .false.
    else
      ! we allocate memory for the problem
      allocate(w)
      w = TwoStreamIRWrk(nz)
      deallocate_w = .true.
    endif

    ! Pointers to work memory
    gam1 => w%gam1
    gam2 => w%gam2
    lambda => w%lambda
    cap_gam => w%cap_gam
    e1 => w%e1
    e2 => w%e2
    e3 => w%e3
    e4 => w%e4
    cp0 => w%cp0
    cpb => w%cpb
    cm0 => w%cm0
    cmb => w%cmb
    A => w%A
    B => w%B
    D => w%D
    E => w%E
    y1 => w%y1
    y2 => w%y2

    ! Hemispheric mean Two-Stream coefficients (Table 1)
    gam1 = 2.0_dp - w0*(1.0_dp + gt)
    gam2 = w0*(1.0_dp - gt)
    ! u1 = 0.5_dp (parameter)

    ! lambda, and capital Gamma (Equations 21, 22)
    lambda = sqrt(gam1**2.0_dp - gam2**2.0_dp)
    cap_gam = gam2 / (gam1 + lambda)
    ! cap_gam = (gam1-lambda)/gam2 ! this is the same as above
    
    ! e's (Equation 44)
    do i = 1,nz
      wrk_real = exp(-lambda(i)*tau(i))
      e1(i) = 1.0_dp + cap_gam(i)*wrk_real
      e2(i) = 1.0_dp - cap_gam(i)*wrk_real
      e3(i) = cap_gam(i) + wrk_real
      e4(i) = cap_gam(i) - wrk_real
    enddo
    
    ! C+ and C- (Equation 27)
    ! norm = 2.0_dp*pi*u1 (parameter)
    do i=1,nz
      b0n = bplanck(i)
      b1n = (bplanck(i+1) - b0n)/tau(i)
      
      cp0(i) = norm*(b0n + b1n*(1.0_dp/(gam1(i)+gam2(i))))
      cpb(i) = norm*(b0n + b1n*(tau(i) + 1.0_dp/(gam1(i)+gam2(i))))
      cm0(i) = norm*(b0n + b1n*(-1.0_dp/(gam1(i)+gam2(i))))
      cmb(i) = norm*(b0n + b1n*(tau(i) - 1.0_dp/(gam1(i)+gam2(i))))

    enddo

    Ssfc = emissivity*pi*bplanck(nz+1) ! ground
      
    ! Coefficients of tridiagonal linear system (Equations 39 - 43)
    ! Odd coeficients (Equation 41)
    A(1) = 0.0_dp
    B(1) = e1(1)
    D(1) = -e2(1)
    E(1) = 0.0_dp - cm0(1) ! assumes no downward diffuse flux at top-of-atmosphere
    do i = 1, nz-1
      l = 2*i + 1
      A(l) = e2(i)*e3(i) - e4(i)*e1(i)
      B(l) = e1(i)*e1(i+1) - e3(i)*e3(i+1)
      D(l) = e3(i)*e4(i+1) - e1(i)*e2(i+1)
      E(l) = e3(i)*(cp0(i+1) - cpb(i)) + e1(i)*(cmb(i) - cm0(i+1))
    enddo
    
    ! Even coefficients (Equation 42)
    do i = 1, nz-1
      l = 2*i
      A(l) = e2(i+1)*e1(i) - e3(i)*e4(i+1)
      B(l) = e2(i)*e2(i+1) - e4(i)*e4(i+1)
      D(l) = e1(i+1)*e4(i+1) - e2(i+1)*e3(i+1)
      E(l) = e2(i+1)*(cp0(i+1) - cpb(i)) - e4(i+1)*(cm0(i+1) - cmb(i))
    enddo
    l = 2*nz
    A(l) = e1(nz) - Rsfc*e3(nz)
    B(l) = e2(nz) - Rsfc*e4(nz)
    D(l) = 0.0_dp
    E(l) = Ssfc - cpb(nz) + Rsfc*cmb(nz)
    
    ! Solve tridiagonal system. e is solution
    call tridiag(nz*2,a,b,d,e)
    
    ! unpack solution
    do i = 1, nz
      l = 2*i
      y1(i) = e(l-1)
      y2(i) = e(l)
    enddo
    
    ! Equations 31 and 32.
    fup(1) = ((y1(1)*e3(1)-y2(1)*e4(1))+cp0(1))
    fdn(1) = 0.0_dp
    do i = 1,nz
       fup(i+1) = (y1(i)*e1(i)+y2(i)*e2(i)+cpb(i))
       fdn(i+1) = (y1(i)*e3(i)+y2(i)*e4(i)+cmb(i))
    enddo

    if (deallocate_w) then
      deallocate(w)
    endif
    
  end subroutine
  
  subroutine tridiag(n, a, b, c, d)
    integer,intent(in) :: n
    real(dp), intent(inout) :: a(n), b(n), c(n), d(n)

    integer :: i
    
    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)
    
    do i= 2, n-1
      c(i) = c(i)/(b(i) - a(i)*c(i-1))
      d(i) = (d(i) - a(i)*d(i-1)) / (b(i) - a(i)*c(i-1))
    enddo
    
    d(n) = (d(n) - a(n)*d(n-1)) / (b(n) - a(n)*c(n-1))

    do i = n-1,1,-1
      d(i) = d(i) - c(i)*d(i+1)
    enddo
  end subroutine

  function TwoStreamSolarWrk_create(nz) result(wrk)
    integer, intent(in) :: nz
    type(TwoStreamSolarWrk) :: wrk

    allocate(wrk%e1(nz))
    allocate(wrk%e2(nz))
    allocate(wrk%e3(nz))
    allocate(wrk%e4(nz))
    allocate(wrk%tauc(nz+1))
    allocate(wrk%direct(nz+1))
    allocate(wrk%cp0(nz))
    allocate(wrk%cpb(nz))
    allocate(wrk%cm0(nz))
    allocate(wrk%cmb(nz))
    allocate(wrk%A(nz*2))
    allocate(wrk%B(nz*2))
    allocate(wrk%D(nz*2))
    allocate(wrk%E(nz*2))
    allocate(wrk%y1(nz))
    allocate(wrk%y2(nz))

  end function

  function TwoStreamIRWrk_create(nz) result(wrk)
    integer, intent(in) :: nz
    type(TwoStreamIRWrk) :: wrk

    allocate(wrk%gam1(nz))
    allocate(wrk%gam2(nz))
    allocate(wrk%lambda(nz))
    allocate(wrk%cap_gam(nz))
    allocate(wrk%e1(nz))
    allocate(wrk%e2(nz))
    allocate(wrk%e3(nz))
    allocate(wrk%e4(nz))
    allocate(wrk%cp0(nz))
    allocate(wrk%cpb(nz))
    allocate(wrk%cm0(nz))
    allocate(wrk%cmb(nz))
    allocate(wrk%A(nz*2))
    allocate(wrk%B(nz*2))
    allocate(wrk%D(nz*2))
    allocate(wrk%E(nz*2))
    allocate(wrk%y1(nz))
    allocate(wrk%y2(nz))

  end function
  
end module
