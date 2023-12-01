module BD_clouds_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: rho_c = 3.560_dp
  real(dp), parameter :: t_grow = 150.0_dp
  real(dp), parameter :: t_evap = 150.0_dp

  real(dp), parameter :: reff = 5.0_dp * 1e-4_dp
  real(dp), parameter :: rm = 1.504_dp * 1e-4_dp
  real(dp), parameter :: rV = 8.084_dp * 1e-4_dp
  real(dp), parameter :: sigma = 2.0_dp
  
  real(dp), parameter :: mol_w_sp = 140.6931_dp 

  real(dp), parameter :: bar = 1.0e5_dp ! bar to pa
  real(dp), parameter :: atm = 1.01325e5_dp ! atm to pa
  real(dp), parameter :: dyne = 0.1_dp ! dyne to pa

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R = 8.31446261815324_dp

  real(dp), parameter :: cfl = 0.95_dp

  !! Diameter, LJ potential and molecular weight for background gases ! Do everything in cgs for vf calculation
  real(dp), parameter :: d_OH = 3.06e-8_dp, LJ_OH = 100.0_dp * kb, molg_OH = 17.00734_dp  ! estimate
  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_H2O = 2.641e-8_dp, LJ_H2O = 809.1_dp * kb, molg_H2O = 18.01528_dp
  real(dp), parameter :: d_H = 2.5e-8_dp, LJ_H =  30.0_dp * kb, molg_H = 1.00794_dp
  real(dp), parameter :: d_CO = 3.690e-8_dp, LJ_CO = 91.7_dp * kb, molg_CO = 28.0101_dp
  real(dp), parameter :: d_CO2 = 3.941e-8_dp, LJ_CO2 = 195.2_dp * kb, molg_CO2 = 44.0095_dp
  real(dp), parameter :: d_O = 2.66e-8_dp, LJ_O = 70.0_dp * kb, molg_O = 15.99940_dp
  real(dp), parameter :: d_CH4 = 3.758e-8_dp, LJ_CH4 = 148.6_dp * kb, molg_CH4 = 16.0425_dp
  real(dp), parameter :: d_C2H2 = 4.033e-8_dp, LJ_C2H2 = 231.8_dp * kb, molg_C2H2 = 26.0373_dp
  real(dp), parameter :: d_NH3 = 2.900e-8_dp, LJ_NH3 = 558.3_dp * kb, molg_NH3 = 17.03052_dp
  real(dp), parameter :: d_N2 = 3.798e-8_dp, LJ_N2 = 71.4_dp * kb, molg_N2 = 14.0067_dp
  real(dp), parameter :: d_HCN = 3.630e-8_dp, LJ_HCN = 569.1_dp * kb, molg_HCN = 27.0253_dp
  real(dp), parameter :: d_He = 2.511e-8_dp, LJ_He = 10.22_dp * kb, molg_He = 4.002602_dp


  integer :: ir, iwl
  real(dp),  allocatable, dimension(:) :: wl_lr, rad_lr, nd_lr
  real(dp),  allocatable, dimension(:,:) :: kext_lr, a_lr, g_lr, kext_lrN
  real(dp),  allocatable, dimension(:) :: kext_ln, a_ln, g_ln
  logical :: first_call = .True.

  public :: BD_clouds_chem, BD_clouds_vf, BD_clouds_adv
  private :: p_vap_sp, minmod, first_call

contains 

  subroutine BD_clouds_vf(nlay, nq, Rd_air, grav_in, q_VMR, pl_in, Tl, qc, vf)
    implicit none

    integer, intent(in) :: nlay, nq
    real(dp), intent(in) :: grav_in
    real(dp), dimension(nlay,nq), intent(in) :: q_VMR
    real(dp), dimension(nlay), intent(in) :: pl_in, Tl, Rd_air, qc

    real(dp), dimension(nlay), intent(out) :: vf

    integer :: g, k
    real(dp) :: rho, nu_mix, l_scale, Kn, beta, top, bot, grav
    real(dp), dimension(nq) :: d_g, LJ_g, molg_g, nu_g
    real(dp), dimension(nlay) :: mu, pl

    ! Need to set the arrays manually here for each tracer
    ! representation
    d_g = (/d_OH, d_H2, d_H2O, d_H, d_CO, d_CO2, d_O, d_CH4, d_C2H2, & 
      & d_NH3, d_N2, d_HCN, d_He/)
    LJ_g = (/LJ_OH, LJ_H2, LJ_H2O, LJ_H, LJ_CO, LJ_CO2, LJ_O, LJ_CH4, &
      & LJ_C2H2, LJ_NH3, LJ_N2, LJ_HCN, LJ_He/)
    molg_g = (/molg_OH, molg_H2, molg_H2O, molg_H, molg_CO, molg_CO2, &
      & molg_O, molg_CH4, molg_C2H2, molg_NH3, molg_N2, molg_HCN, molg_He/)

    !d_g = (/d_H2, d_He/)
    !LJ_g = (/LJ_H2, LJ_He/)
    !molg_g = (/molg_H2, molg_He/)

    pl(:) = pl_in(:) * 10.0_dp ! Convert Pa to Dyne

    mu(:) = R/Rd_air(:) * 1000.0_dp ! Atmospheric molecualr weight (g mol-1)

    grav = grav_in * 100.0_dp ! gravity in cm s-2 

    do k = 1, nlay

      if (qc(k) < 1e-20_dp) then
        vf(k) = 0.0_dp
        cycle
      end if

      !! Atmospheric density g cm-3
      rho = (pl(k) * mu(k) * amu)/(kb * Tl(k))
  
      !! Find the dynamical viscosity of each gas
      do g = 1, nq
        ! Dynamical viscosity formula - Rosner (2000/2012) using
        ! Ackerman & Marley (2001) constants
        nu_g(g) = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(g)*amu)*kb*Tl(k))/(pi*d_g(g)**2)) &
        & * (((kb*Tl(k))/LJ_g(g))**(0.16_dp)/1.22_dp)
      end do

      !! Find the mixture dynamical viscosity using the square root
      !mixing law (Rosner 2021)
      !! See Lee (2023) for examination
      top = 0.0_dp
      bot = 0.0_dp
      do g = 1, nq
        top = top + sqrt(molg_g(g))*q_VMR(k,g)*nu_g(g)
        bot = bot + sqrt(molg_g(g))*q_VMR(k,g)
      end do
      nu_mix = top/bot

      !! For consistency, we now use the dynamical viscosity to find the
      !mean free path of the layer
      l_scale = (nu_mix/pl(k)) * sqrt((pi * kb * Tl(k)) / (2.0_dp * mu(k) * amu))

      !! Knudsen number and Cunningham slip factor at mean grain size
      Kn = l_scale / rV
      beta = 1.0_dp + Kn * (1.256_dp + 0.4_dp * exp(-1.1_dp/Kn))

      !! Final v_f is negative (downward)
      vf(k) = -(2.0_dp * beta * rV**2 * grav * (rho_c - rho)) / (9.0_dp * nu_mix)

      !! Convert to pressure coordinates (dyne s-1)
      vf(k) = - rho * grav * vf(k)
      
    end do

    !! Convert to Pa s-1
    vf(:) = vf(:) / 10.0_dp

  end subroutine BD_clouds_vf

  subroutine BD_clouds_chem(nlay, t_step, sp, q_v, q_c, pl, Tl, Rd_air, Kzz, q0, grav, sat)
    implicit none

    integer, intent(in) :: nlay
    character(len=10), intent(in) :: sp
    real(dp), intent(in) :: t_step, q0, grav
    real(dp), dimension(nlay), intent(in) :: pl, Tl, Rd_air, Kzz
    
    real(dp), dimension(nlay), intent(inout) :: q_v, q_c, sat

    integer :: k
    real(dp) :: s, p_vap,  q_s, dqvdt, dqcdt, t_c, tau_deep, Hp
    real(dp) :: k1v, k2v, k3v, k4v,  k1c, k2c, k3c, k4c

    !! Calculate stability coefficent and change in vapour/cloud mass
    !fractions
    do k = 1, nlay

      !! Vapour pressure and saturation vapour pressure
      p_vap = p_vap_sp(sp, Tl(k))
      sat(k) = (q_v(k) * pl(k))/p_vap

      !! Give stability coefficent
      if (sat(k) < 0.9999_dp) then
        ! Evaporate the cloud mass portion
        s = 0.0_dp
        t_c = t_evap
      else if (sat(k) > 1.0001_dp) then
        ! Condense the cloud vapour portion
        s = 1.0_dp
        t_c = t_grow
      else
        ! Don't do anything if at saturation
        cycle
      end if

      !! Equilibrium (sat = 1) vapour pressure fraction
      q_s = max(1.0e-30_dp,p_vap/pl(k))
      q_s = min(1.0_dp, q_s)

      !! Calculate tracer tendencies
      if (k == nlay) then
        !! Include lower boundary replenishment rate at tau_deep relaxtion timescale
        Hp = (Rd_air(k) * Tl(k)) / grav
        tau_deep = Hp**2/Kzz(k)
        dqvdt = (1.0_dp - s)*min(q_s - q_v(k), q_c(k))/t_c  &
          & - s*(q_v(k) - q_s)/t_c - (q_v(k) - q0)/tau_deep
      else
        !! Calculate qv tracer tendency
        dqvdt = (1.0_dp - s)*min(q_s - q_v(k), q_c(k))/t_c  &
          & - s*(q_v(k) - q_s)/t_c
      end if

      !! Calculate qc tracer tendency
      dqcdt = s*(q_v(k) - q_s)/t_c  & 
        & - (1.0_dp - s)*min(q_s - q_v(k), q_c(k))/t_c


      !! Timestep vapour and cloud tracer values
      !! check if integration would go negative for vapour and give
      !limits
      if ((q_v(k) + dqvdt * t_step) <= 0.0_dp) then
        q_v(k) = 1e-30_dp
      else
        q_v(k) = q_v(k) + dqvdt * t_step
      end if
      if ((q_c(k) + dqcdt * t_step) <= 0.0_dp) then
        q_c(k) = 1e-30_dp
      else
        q_c(k) = q_c(k) + dqcdt * t_step
      end if

      !! Add limiter to 1
      q_v(k) = min(1.0_dp,q_v(k))
      q_c(k) = min(1.0_dp,q_c(k))   

    end do

  end subroutine BD_clouds_chem

  real(dp) function p_vap_sp(sp, T)
    implicit none

    character(len=10), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC


    ! Return vapour pressure in pa
    select case(sp)
    case('C')
      p_vap_sp = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
    !case('TiC')
    !case('SiC')
    !case('CaTiO3')
    case('TiO2')
      ! NIST 5 param fit
      p_vap_sp = exp(-7.70443e4_dp/T + 4.03144e1_dp - 2.59140e-3_dp*T &
      & + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! NIST 5 param fit
      p_vap_sp = exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
      & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)
    case('Al2O3')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-73503.0_dp/T + 22.01_dp) * atm
    case('Fe')
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      if (T > 1800.0_dp) then
        p_vap_sp = exp(9.86_dp - 37120.0_dp/T) * bar
      else
        p_vap_sp = exp(15.71_dp - 47664.0_dp/T) * bar
      end if
    case('FeO')
      ! NIST 5 param fit
      p_vap_sp = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
      & + 3.18636e-7_dp*T**2 - 0.0_dp*T**3)
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-62279_dp/T + 20.944_dp) * atm
    case('MgSiO3')
      ! A&M (2001)
      p_vap_sp = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('SiO2')
      ! NIST 5 param fit
      p_vap_sp = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3)
    case('SiO')
      ! Gail et al. (2013)
      p_vap_sp = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! NIST 5 param fit
      p_vap_sp = exp(-4.78455e4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T &
      & - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
    case('MnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap_sp =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
    case('ZnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(12.812_dp - 15873.0_dp/T) * bar
    case('KCl')
      ! NIST 5 param fit
      p_vap_sp = exp(-2.69250e+4_dp/T + 3.39574e1_dp - 2.04903e-3_dp*T &
      & - 2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
    case('NaCl')
      ! NIST 5 param fit
      p_vap_sp = exp(-2.79146e+4_dp/T + 3.46023e1_dp - 3.11287e-3_dp*T &
      & + 5.30965e-07_dp*T**2 - 2.59584e-12_dp*T**3)
    case('NH4Cl')
      p_vap_sp = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      ! Ackerman & Marley et al. (2001) H2O liquid & ice vapour pressure expressions
      TC = T - 273.15_dp
      if (T > 1048.0_dp) then
        p_vap_sp = 6.0e8_dp
      else if (T < 273.16_dp) then
        p_vap_sp = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      else
        p_vap_sp = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp))
      end if
    case('NH3')
      ! Ackerman & Marley et al. (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap_sp = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap_sp = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      !--- Stull (1947) ---
      if (T < 30.0_dp) then ! Limiter for very cold T
        p_vap_sp = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
      end if
      if (T < 212.8_dp) then
        p_vap_sp = 10.0_dp**(4.43681_dp - 829.439_dp/(T-25.412_dp)) * bar
      else
        p_vap_sp = 10.0_dp**(4.52887_dp - 958.587_dp/(T-0.539_dp)) * bar
      end if
    case('S2')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap_sp = exp(27.0_dp - 18500.0_dp/T) * bar
      else
        p_vap_sp = exp(16.1_dp - 14000.0_dp/T) * bar
      end if
    case('S8')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap_sp = exp(20.0_dp - 11800.0_dp/T) * bar
      else
        p_vap_sp = exp(9.6_dp - 7510.0_dp/T) * bar
      end if
    case default
      print*, 'Saturation: dust species not found: ', trim(sp), 'Stopping!'
      stop
    end select

  end function p_vap_sp

  subroutine BD_clouds_adv(nlay, nq, tend, q, vf, delp)
    implicit none

    integer, intent(in) :: nlay, nq
    real(dp), dimension(nlay,nq), intent(inout) :: q
    real(dp), dimension(nlay), intent(in) :: vf, delp
    real(dp), intent(in) :: tend

    real(dp), dimension(nlay) :: qc
    real(dp), dimension(nlay) :: sig, c, vf_c
    real(dp) :: D, tnow, dt
    integer :: i, iit, n

    vf_c(:) = vf(:)

    dt = tend
    do i = 1, nlay
      dt = min(dt,cfl*(delp(i))/abs(vf_c(i)))
    enddo

    tnow = 0.0_dp
    iit = 1

    do while ((tnow < tend) .and. (iit < 10000))

      ! If next time step overshoots - last time step is equal tend
      if (tnow + dt > tend) then
        dt = tend - tnow
      end if

      !! Find the courant number
      c(:) = abs(vf_c(:)) * dt / delp(:)

      do n = 1, nq

        ! Apply boundary conditions
        q(1,n) = 1e-30_dp
        q(nlay,n) = 1e-30_dp

        q(:,n) = max(q(:,n),1e-30_dp)

        !! Find the minmod limiter
        call minmod(nlay,q(:,n),delp,sig)

        !! Perform McCormack step
        qc(:) = q(:,n)
        qc(1:nlay-1) = q(1:nlay-1,n) - sig(1:nlay-1)*c(1:nlay-1)*(q(2:nlay,n) - q(1:nlay-1,n))
        q(2:nlay,n) = 0.5_dp * (q(2:nlay,n) + qc(2:nlay) - c(2:nlay)*(qc(2:nlay) - qc(1:nlay-1)))

        q(:,n) = max(q(:,n),1e-30_dp)

        ! Apply boundary conditions
        q(1,n) = 1e-30_dp
        q(nlay,n) = 1e-30_dp

      end do

      tnow = tnow + dt
      iit = iit + 1

    end do

    ! Apply boundary conditions
    q(1,:) = 1e-30_dp
    q(nlay,:) = 1e-30_dp

  end subroutine BD_clouds_adv


  subroutine minmod(nlay,q,dz,sig)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q, dz

    real(dp), dimension(nlay), intent(out) :: sig

    integer :: i
    real(dp) :: de_minus, de_plus

    sig(1) = 0.0_dp
    do i = 2, nlay-1
      de_minus = (q(i) - q(i-1)) / dz(i)
      de_plus = (q(i+1) - q(i)) / dz(i)
      if ((de_minus > 0.0_dp) .and. (de_plus > 0.0_dp)) then
        sig(i) = min(de_minus, de_plus)
      else if ((de_minus < 0.0_dp) .and. (de_plus < 0.0_dp)) then
        sig(i) = max(de_minus, de_plus)
      else
        sig(i) = 0.0_dp
      end if
    end do
    sig(nlay) = 0.0_dp

  end subroutine minmod

end module BD_clouds_mod

