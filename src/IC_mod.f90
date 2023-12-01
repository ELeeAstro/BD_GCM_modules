!! This module

module IC_mod
  use, intrinsic :: iso_fortran_env
  use k_Rosseland_mod, only : k_Ross_Freedman, k_Ross_Valencia, gam_Parmentier, Bond_Parmentier
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  public :: IC_profile
  private :: Parmentier_IC, Guillot_IC

contains

  subroutine IC_profile(iIC,corr,nlay,p0,pl,k_IR,Tint,grav,Tl,prc, &
    & kappa)
    implicit none

    !! Input flags
    integer, intent(in) :: iIC
    logical, intent(in) :: corr

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), intent(in) :: p0, Tint, grav, kappa
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), dimension(nlay), intent(in) :: k_IR

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl
    real(dp), intent(out) :: prc

    select case(iIC)

    case(3)
      ! Guillot (2010) analytic RE T-p profile
      call Guillot_IC(nlay,p0,pl,k_IR(1),Tint,grav,Tl)
    case(4)
      ! Parmentier et al. (2014, 2015) IC picket fence IC
      ! check the table_num and metallicity parameters for different options
      call Parmentier_IC(nlay,pl,Tint,0.01_dp,1.0_dp,grav,Tl)
    case default
      print*, 'Invalid IC integer in IC_mod, stopping'
      stop

    end select

    !! Perform adiabatic correction according to Parmentier et al. (2015)
    if (corr .eqv. .True.) then
      call adiabat_correction(nlay,Tl,pl,prc)
    else
      prc = p0
    end if

  end subroutine IC_profile

  subroutine Guillot_IC(nlay,p0,pl,k_IR,Tint,grav,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), intent(in) :: p0, k_IR, Tint, grav
    real(dp), dimension(nlay), intent(in) :: pl

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl

    !! Work variables
    real(dp), dimension(nlay) :: tau_IRl
    real(dp) :: gam, tau0

    tau0 = k_IR/grav * p0

    tau_IRl(:) = (pl(:)/p0 * tau0)

    Tl(:) = ((3.0_dp/4.0_dp) * Tint**4 * (tau_IRl(:) + 2.0_dp/3.0_dp))
    Tl(:) = Tl(:)**(1.0_dp/4.0_dp)

  end subroutine Guillot_IC

  !! This subroutine follows Parmentier & Guillot (2014, 2015) non-grey picket fence scheme
  subroutine Parmentier_IC(nlay,pl,Tint,mu,Tirr,grav,Tl)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), intent(in) :: Tint, mu, Tirr, grav


    real(dp), dimension(nlay), intent(out) :: Tl

    integer, parameter :: table_num = 1
    real(dp), parameter :: met = 0.0_dp

    real(dp) :: Teff0, Teff, Tmu, Bond, Tskin
    real(dp), dimension(3) :: gam_V, Beta_V
    real(dp), dimension(2) :: Beta
    real(dp) :: gam_1, gam_2, gam_P, tau_lim

    integer :: i, j
    real(dp) :: a0, a1, b0, A, B, At1, At2
    real(dp), dimension(3) :: a2, a3, b1, b2, b3, Av1, Av2
    real(dp), dimension(3) :: C, D, E
    real(dp), dimension(nlay+1) :: tau
    real(dp), dimension(nlay) :: kRoss

    !! Effective temperature parameter
    Tmu = (mu * Tirr**4)**(1.0_dp/4.0_dp)

    !! Find Bond albedo of planet - Bond albedo is given by mu = 1/sqrt(3)
    Teff0 = (Tint**4 + (1.0_dp/sqrt(3.0_dp)) * Tirr**4)**(1.0_dp/4.0_dp)
    call Bond_Parmentier(Teff0, grav, Bond)

    Teff = (Tint**4 + (1.0_dp - Bond) * mu * Tirr**4)**(1.0_dp/4.0_dp)

    !! Find the V band gamma, beta and IR gamma and beta ratios for this profile
    ! Passed mu, so make lat = acos(mu) and lon = 0
    call gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim)

    gam_V(:) = gam_V(:) / mu

    !! Hard work starts here - first calculate all the required coefficents
    At1 = gam_1**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_1))
    At2 = gam_2**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_2))
    Av1(:) = gam_1**2*log(1.0_dp + gam_V(:)/gam_1)
    Av2(:) = gam_2**2*log(1.0_dp + gam_V(:)/gam_2)

    a0 = 1.0_dp/gam_1 + 1.0_dp/gam_2

    a1 = -1.0_dp/(3.0_dp * tau_lim**2) * (gam_P/(1.0_dp-gam_P) * (gam_1 + gam_2 - 2.0_dp)/(gam_1 + gam_2) &
    &  + (gam_1 + gam_2)*tau_lim - (At1 + At2)*tau_lim**2)

    a2(:) = tau_lim**2/(gam_P*gam_V(:)**2) &
    &  * ((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(gam_1+gam_2) &
    & - 3.0_dp*gam_V(:)*(6.0_dp*gam_1**2*gam_2**2-gam_V(:)**2*(gam_1**2+gam_2**2))) &
    & / (1.0_dp-gam_V(:)**2 * tau_lim**2)

    a3(:) = -tau_lim**2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(Av2(:)+Av1(:)) &
     &/(gam_P*gam_V(:)**3*(1.0_dp-gam_V(:)**2*tau_lim**2))

    b0 = 1.0_dp/(gam_1*gam_2/(gam_1-gam_2)*(At1-At2)/3.0_dp-(gam_1*gam_2)**2/sqrt(3.0_dp*gam_P)-(gam_1*gam_2)**3 &
    & / ((1.0_dp-gam_1)*(1.0_dp-gam_2)*(gam_1+gam_2)))

    b1(:) = gam_1*gam_2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*tau_lim**2 &
    & / (gam_P*gam_V(:)**2*(gam_V(:)**2*tau_lim**2-1.0_dp))

    b2(:) = 3.0_dp*(gam_1+gam_2)*gam_V(:)**3/((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2))

    b3(:) = (Av2(:)-Av1(:))/(gam_V(:)*(gam_1-gam_2))

    A = 1.0_dp/3.0_dp*(a0+a1*b0)
    B = -1.0_dp/3.0_dp*(gam_1*gam_2)**2/gam_P*b0
    C(:) = -1.0/3.0_dp*(b0*b1(:)*(1.0_dp+b2(:)+b3(:))*a1+a2(:)+a3(:))
    D(:) = 1.0/3.0_dp*(gam_1*gam_2)**2/gam_P*b0*b1(:)*(1.0_dp+b2(:)+b3(:))
    E(:) = (3.0_dp-(gam_V(:)/gam_1)**2)*(3.0_dp-(gam_V(:)/gam_2)**2)/(9.0_dp*gam_V(:)*((gam_V(:)*tau_lim)**2-1.0_dp))

    ! T-p structure calculation - we follow exactly V. Parmentier's method
    ! Estimate the skin temperature by setting tau = 0
    tau(1) = 0.0_dp
    Tskin = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tskin = Tskin**(1.0_dp/4.0_dp)
    ! Estimate the opacity TOA at the skin temperature - assume this is = first layer optacity
    call k_Ross_Freedman(Tskin, pl(1), met, kRoss(1))
    !call k_Ross_Valencia(Tskin, pe(1), met, kRoss(1))

    ! Recalculate the upmost tau with new kappa
    tau(1) = kRoss(1)/grav * pl(1)
    ! More accurate layer T at uppermost layer
    Tl(1) = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tl(1) = Tl(1)**(1.0_dp/4.0_dp)

    ! Now we can loop in optical depth space to find the T-p profile
    do i = 2, nlay
      ! Initial guess for layer
      call k_Ross_Freedman(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      !call k_Ross_Valencia(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
      Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
      & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
      Tl(i) = Tl(i)**(1.0_dp/4.0_dp)
      ! Convergence loop
      do j = 1, 5
        call k_Ross_Freedman(sqrt(Tl(i-1)*Tl(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        !call k_Ross_Valencia(sqrt(Tl(i-1)*Tl(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
        Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
        & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
        Tl(i) = Tl(i)**(1.0_dp/4.0_dp)

      end do

    end do

  end subroutine Parmentier_IC

  !! Subroutine that corrects for adiabatic region following Parmentier & Guillot (2015)
  subroutine adiabat_correction(nlay,Tl,pl,prc)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) ::  pl

    !! Output quantities
    real(dp), dimension(nlay), intent(inout) :: Tl
    real(dp), intent(out) :: prc

    !! Work variables
    integer :: k, iRC, iRC1
    real(dp), dimension(nlay) :: gradrad, gradad

    do k = 1, nlay-1
      gradrad(k) = log(Tl(k)/Tl(k+1))/log(pl(k)/pl(k+1))
      gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
      !print*, k, gradrad(k), gradad(k)
    end do
    gradrad(nlay) = 0.0_dp
    gradad(nlay) = 0.0_dp

    iRC = nlay-1
    iRC1 = nlay-1

    do k = nlay-1, 1, -1
      if (IRC1 <= k+1) then
        if (gradrad(k) > 0.7_dp*gradad(k)) then
          iRC1 = k
        endif
        if (gradrad(k) > 0.98_dp*gradad(k)) then
         iRC = k
         prc = pl(iRC)
        endif
      end if
    end do

    if (iRC < nlay) then
      do k = iRC, nlay-1
        gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
        if (gradad(k) < 0.0_dp) then
          gradad(k) = 0.0_dp
        end if
        !if (pl(k) > 1.0_dp*1e6_dp) then
        Tl(k+1)=Tl(k)*(pl(k+1)/pl(k))**gradad(k)
      !end if
      end do
    end if

  end subroutine adiabat_correction

end module IC_mod