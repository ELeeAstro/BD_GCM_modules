module BD_clouds_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: rho_c = 3.560_dp
  real(dp), parameter :: t_grow = 10.0_dp
  real(dp), parameter :: t_evap = 10.0_dp

  real(dp), parameter :: reff = 0.15_dp * 1e-4_dp
  real(dp), parameter :: rm = 0.1_dp * 1e-4_dp
  real(dp), parameter :: rV = 10_dp * 1e-4_dp
  real(dp), parameter :: sigma = 1.5_dp

  real(dp), parameter :: mol_w_sp = 60.08430_dp

  real(dp), parameter :: bar = 1.0e5_dp ! bar to pa
  real(dp), parameter :: atm = 1.01325e5_dp ! atm to pa
  real(dp), parameter :: dyne = 0.1_dp ! dyne to pa

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R_gas = 8.31446261815324_dp

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


  integer :: ir, iwl, iT
  real(dp),  allocatable, dimension(:) :: wl_lr, rad_lr, nd_lr, nd_lr_norm
  real(dp),  allocatable, dimension(:,:) :: kext_lr, a_lr, g_lr, kext_lrN
  real(dp),  allocatable, dimension(:) :: kext_ln, a_ln, g_ln, T_Ross
  real(dp), allocatable, dimension(:,:) :: Qext_Ross, Qsca_Ross, gg_Ross
  real(dp), allocatable, dimension(:) :: iQext, iQsca, igg, ifunc_k, ifunc_a, ifunc_g
  logical :: first_call = .True.

  public :: BD_clouds_chem, BD_clouds_vf, BD_clouds_Ross_opac
  private :: p_vap_sp, first_call, locate, trapz

contains 

  subroutine BD_clouds_Ross_opac(nlay,sp,Rd_air,pl_in,Tl,qc,k_ext_cld,a_cld,g_cld)
    implicit none

    integer, intent(in) :: nlay
    character(len=20), intent(in) :: sp
    real(dp), dimension(nlay), intent(in) :: qc, Rd_air, pl_in, Tl

    real(dp), dimension(nlay), intent(out) :: k_ext_cld, a_cld, g_cld

    integer :: k, r, iT1, iT2
    real(dp), dimension(nlay) :: pl
    real(dp) :: rho, mu, mu_ratio, expterm, N0

    if (first_call .eqv. .True.) then
      call read_Ross_table(sp)
      first_call = .False.
      allocate(iQext(ir),iQsca(ir),igg(ir),ifunc_k(ir), ifunc_a(ir), ifunc_g(ir))
    end if

    pl(:) = pl_in(:) * 10.0_dp ! Convert Pa to Dyne

    do k = 1, nlay

      if (qc(k) < 1.0e-8) then
        k_ext_cld(k) = 1.0e-30_dp
        a_cld(k) = 1.0e-30_dp
        g_cld(k) = 1.0e-30_dp 
        cycle
      end if

      mu = R_gas/Rd_air(k) * 1000.0_dp

      rho = (pl(k) * mu * amu)/ (kb * Tl(k)) ! g cm-3

      mu_ratio = mol_w_sp/mu

      expterm = exp(-9.0_dp/2.0_dp * log(sigma)**2)

      N0 = (3.0_dp*qc(k)*mu_ratio*rho)/(4.0_dp*pi*rho_c*rm**3) * expterm ! cm-3

      ! Find number density for each size and calculate kappa
      do r = 1, ir
        nd_lr(r) = (N0  / (rad_lr(r) * sqrt(twopi) * log(sigma))) * &
          & exp(-(log(rad_lr(r)/rm))**2/(2.0_dp * log(sigma)**2))
        nd_lr_norm(r) = nd_lr(r)/N0
      end do
      
      ! Interpolate tables to find Rosseland mean weighted values at the layer temperature
      if (Tl(k) <= T_Ross(1)) then
        iT1 = 1 
        iT2 = iT1 + 1
      else if (Tl(k) >= T_Ross(iT)) then
        iT1 = iT - 1
        iT2 = iT
      else
        call locate(T_Ross(:), iT, Tl(k), iT1)
        iT2 = iT1 + 1
      end if

      do r = 1, ir
        call linear_interp(Tl(k), T_Ross(iT1), T_Ross(iT2), Qext_Ross(r,iT1), Qext_Ross(r,iT2), iQext(r))
        call linear_interp(Tl(k), T_Ross(iT1), T_Ross(iT2), Qsca_Ross(r,iT1), Qsca_Ross(r,iT2), iQsca(r))
        call linear_interp(Tl(k), T_Ross(iT1), T_Ross(iT2), gg_Ross(r,iT1), gg_Ross(r,iT2), igg(r))
      end do

      !! Calculate weighted opacity values
      do r = 1, ir
        ifunc_k(r) = iQext(r) *  pi * rad_lr(r)**2 * nd_lr(r)
        ifunc_a(r) = iQsca(r)/iQext(r)
        ifunc_g(r) = igg(r)
      end do


      ! Use trapezoid rule - function in [cm-3 cm-1]
      ! Total extinction = integral over all sizes
      k_ext_cld(k) = trapz(rad_lr(:),ifunc_k(:))
      ! effective SSA = integral for k_sca / k_ext =  integral(k_ext * a) / k_ext
      a_cld(k) = trapz(rad_lr(:),ifunc_k(:) * ifunc_a(:)) ! Store intermediate result for cl_out_g
      ! effective g = itergral for k_sca * g / k_sca = integeral(k_sca * g) / k_sca
      g_cld(k) = trapz(rad_lr(:),ifunc_k(:) * ifunc_a(:) * ifunc_g(:))/a_cld(k)

      a_cld(k) = a_cld(k)/ k_ext_cld(k) ! Albedo is scattering/extinction

      k_ext_cld(k) = k_ext_cld(k)/rho/10.0_dp
 
    end do


  end subroutine BD_clouds_Ross_opac

  subroutine read_Ross_table(sp)
    implicit none
    
    character(len=20), intent(in) :: sp

    integer :: u, r

    open(newunit=u,file='cloud_data/RTtable.txt',action='read')
    read(u,*) ir, iT
    allocate(rad_lr(ir),nd_lr(ir),nd_lr_norm(ir),T_Ross(iT))
    read(u,*) rad_lr(:) 
    rad_lr(:) = rad_lr(:) * 1.0e-4_dp
    read(u,*) T_Ross(:)
    close(u)

    open(newunit=u,file='cloud_data/'//trim(sp)//'_rosselandMean_qext.txt',action='read')
    read(u,*)
    allocate(Qext_Ross(ir,iT))
    do r = 1, ir
      read(u,*) Qext_Ross(r,:)
    end do
    close(u)

    open(newunit=u,file='cloud_data/'//trim(sp)//'_rosselandMean_qscat.txt',action='read')
    read(u,*)
    allocate(Qsca_Ross(ir,iT))
    do r = 1, ir
      read(u,*) Qsca_Ross(r,:)
    end do
    close(u)

   
    open(newunit=u,file='cloud_data/'//trim(sp)//'_rosselandMean_gg.txt',action='read')
    read(u,*)
    allocate(gg_Ross(ir,iT))
    do r = 1, ir
      read(u,*) gg_Ross(r,:)
    end do
    close(u)

  end subroutine read_Ross_table

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

    mu(:) = R_gas/Rd_air(:) * 1000.0_dp ! Atmospheric molecualr weight (g mol-1)

    grav = grav_in * 100.0_dp ! gravity in cm s-2 

    do k = 1, nlay

      if (qc(k) < 1e-30_dp) then
        vf(k) = -1.0e-10_dp
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
      !vf(k) = - rho * grav * vf(k)
      
    end do

    !! Convert to Pa s-1
    !vf(:) = vf(:) / 10.0_dp

  end subroutine BD_clouds_vf

  subroutine BD_clouds_chem(nlay, t_end, sp, q_v, q_c, pl, Tl, Rd_air, Kzz, q0, grav)
    implicit none

    integer, intent(in) :: nlay
    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: t_end, q0, grav
    real(dp), dimension(nlay), intent(in) :: pl, Tl, Rd_air, Kzz
    
    real(dp), dimension(nlay), intent(inout) :: q_v, q_c

    integer :: k
    real(dp) :: p_vap, q_s, tau_deep, Hp

    real(dp) :: t_now

    integer :: itol, iout, idid
    integer, dimension(2) :: ipar
    real(dp) :: rtol, atol
    real(dp), dimension(5) :: rpar

    integer, parameter :: n = 2
    integer, parameter :: lwork = 37 !8*2
    integer, parameter :: liwork = 21
    real(dp), dimension(lwork) :: work
    integer, dimension(liwork) :: iwork

    real(dp), dimension(n) :: y

    itol = 0
    rtol = 1e-4_dp
    atol = 1e-30_dp

    iout = 0

    ! Calculate 
    Hp = (Rd_air(nlay) * Tl(nlay)) / grav
    tau_deep = Hp**2/Kzz(nlay)

    do k = 1, nlay

      !! Vapour pressure and saturation vapour pressure
      p_vap = p_vap_sp(sp, Tl(k))

      !! Equilibrium (sat = 1) vapour pressure fraction
      q_s = min(1.0_dp, p_vap/pl(k))


      t_now = 0.0_dp

      y(1) = q_v(k)
      y(2) = q_c(k)

      ipar(1) = k
      ipar(2) = nlay

      !! Set rpar to send through variables into dopri5 to calculate dqdt
      rpar(1) = pl(k)
      rpar(2) = p_vap
      rpar(3) = q_s
      rpar(4) = q0
      rpar(5) = tau_deep


      work(:) = 0.0_dp
      iwork(:) = 0

      call dopri5(n,dqdt,t_now,y,t_end, &
        &         rtol,atol,itol, &
        &         solout,iout, &
        &         work,lwork,iwork,liwork,rpar,ipar,idid)

      q_v(k) = y(1)
      q_c(k) = y(2)

      if (idid /= 1) then
        print*, k, idid
      end if

    end do

  end subroutine BD_clouds_chem

  subroutine dqdt(n,x,y,f,rpar,ipar)
    implicit none

    integer, intent(in) :: n
    integer, dimension(*), intent(in) :: ipar

    real(dp), intent(in) :: x 
    real(dp), dimension(n), intent(inout) :: y
    real(dp), dimension(*), intent(in) :: rpar

    real(dp), dimension(n), intent(out) :: f

    real(dp) :: sat, t_c

    !y(1) = q_v, y(2) = q_c
    !rpar(1) = pl(k) , rpar(2) = p_vap, rpar(3) = q_s , rpar(4) = q0 ,  rpar(5) = tau_deep

    y(1) = max(y(1),1e-30_dp)
    y(2) = max(y(2),1e-30_dp)

    !! Calculate dqdt for vapour and condensate
    sat = (y(1) * rpar(1))/rpar(2)

    !! Calculate dqdt given the supersaturation ratio
    if (sat < 0.99_dp) then
      ! Evaporate the cloud mass portion
      t_c = t_evap

      ! Evaporate from q_c
      f(1) = min(rpar(3) - y(1), y(2))/t_c

    else if (sat > 1.01_dp) then

      ! Condense the cloud vapour portion
      t_c = t_grow

      ! Condense q_v toward the saturation ratio
      f(1) = -(y(1) - rpar(3))/t_c

    else
      f(1) = 0.0_dp
    end if

    f(2) = -f(1)

    ! Add replenishment to lower boundary at the tau_deep rate
    if (ipar(1) == ipar(2)) then
      f(1) = f(1) - (y(1) - rpar(4))/rpar(5)
    end if

  end subroutine dqdt

  subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
    implicit none
    integer, intent(in) :: nr, n, nd, irtrn, ipar(*)
    real(dp), intent(in) ::  Y(*),CON(*),ICOMP(*),RPAR(*), xold, x
  end subroutine solout


  real(dp) function p_vap_sp(sp, T)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC


    ! Return vapour pressure in pa
    select case(sp)
    case('C')
      p_vap_sp = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp)) * dyne
    !case('TiC')
    !case('SiC')
    !case('CaTiO3')
    case('TiO2')
      ! Woitke & Helling (2004)
      p_vap_sp = exp(35.8027_dp - 74734.7_dp/T) * dyne
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
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-62279_dp/T + 20.944_dp) * atm
    case('MgSiO3','MgSiO3_amorph')
      ! Ackerman & Marley (2001)
      p_vap_sp = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('SiO2','SiO2_amorph')
      ! NIST 5 param fit
      p_vap_sp = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) * bar
    case('SiO')
      ! Gail et al. (2013)
      p_vap_sp = exp(-49520.0_dp/T + 32.52_dp) * dyne
    case('Cr')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(7.490_dp - 20592_dp/T) * bar
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
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(7.611_dp - 11382_dp/T) * bar
    case('NaCl')
      ! Stull (1947)
      if (T < 83.0_dp) then ! Limiter for very cold T
        p_vap_sp = 10.0_dp**(5.07184_dp - 8388.497_dp / (83.0_dp - 82.638_dp)) * bar
      else
        p_vap_sp = 10.0_dp**(5.07184_dp - 8388.497_dp / (T - 82.638_dp)) * bar
      end if
    case('NH4Cl')
      p_vap_sp = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
      TC = T - 273.15_dp
      if (T > 1048.0_dp) then
        p_vap_sp = 6.0e8_dp * dyne
      else if (T < 273.16_dp) then
        p_vap_sp = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp)) * dyne
      else
        p_vap_sp = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) * dyne
      end if
    case('NH3')
      ! Ackerman & Marley (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap_sp = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      if (T < 0.5_dp) then ! Limiter for very very very cold T
        p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(0.5_dp-0.49_dp)) * bar       
      else
        p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
      end if
    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap_sp = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      !--- Stull (1947) ---
      if (T < 30.0_dp) then ! Limiter for very cold T
        p_vap_sp = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
      else if (T < 212.8_dp) then
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

  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp), intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  pure function trapz(x, y) result(r)
      !! Calculates the integral of an array y with respect to x using the trapezoid
      !! approximation. Note that the mesh spacing of x does not have to be uniform.
      real(kind=dp), intent(in)  :: x(:)         !! Variable x
      real(kind=dp), intent(in)  :: y(size(x))   !! Function y(x)
      real(kind=dp)              :: r            !! Integral ∫y(x)·dx

      ! Integrate using the trapezoidal rule
      associate(n => size(x))
        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2.0_dp
      end associate
    end function trapz

end module BD_clouds_mod

