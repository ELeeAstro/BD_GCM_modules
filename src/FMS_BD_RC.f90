!!!

!!!

program FMS_BD_RC
  use, intrinsic :: iso_fortran_env
  use k_Rosseland_mod, only : k_Ross_Freedman
  use IC_mod, only : IC_profile
  use CE_mod, only : CE_interpolate_Bezier
  use mini_ch_i_seulex, only: mini_ch_seulex
  use mini_ch_read_reac_list, only : read_react_list
  use BD_kappa_poly_mod, only : calc_kappa
  use BD_MLT_mod, only : BD_MLT
  use BD_clouds_mod, only : BD_clouds_chem, BD_clouds_vf, BD_clouds_adv
  use BD_vert_diff_mod, only: BD_vert_diff
  use ts_Toon_scatter_mod, only : ts_Toon_scatter
  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp

  integer :: n, i, k, u, j, inan
  integer :: nstep, nlay, nlev
  real(dp) :: t_step, t_tot
  real(dp) :: Tint, Fint, pref, pu, met_R
  real(dp) :: tau_IRref
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe
  real(dp), allocatable, dimension(:) :: k_IRl
  real(dp), allocatable, dimension(:) :: tau_IRe
  real(dp), allocatable, dimension(:) :: dT_rad, net_F

  real(dp), allocatable, dimension(:) :: lw_a, lw_g

  real(dp), allocatable, dimension(:) :: vf, sat

  real(dp) :: cp_air, grav, k_IR, kappa_air, Rd_air
  real(dp) :: lw_ac, lw_gc

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  real(dp) :: olr

  logical :: Bezier

  !! kapp_prime calculation
  real(dp), allocatable, dimension(:) :: mu, Rd_bar, cp_bar, k_prime


  !! MLT variables
  real(dp), allocatable, dimension(:) :: dt_conv, Kzz, w_osh

  !! Tracer arrays
  integer :: nq
  real(dp), allocatable, dimension(:,:) :: q

  integer :: ua, ub
  character(len=50) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a, b

  real(dp) :: start, finish

  character(len=50) :: ts_scheme, opac_scheme, adj_scheme, diff_scheme, cloud_chem_scheme, chem_scheme
  character(len=50) :: kappa_prime_calc, cloud_opac_scheme

  !! Mini-chem variables
  integer :: nsp
  character(len=200) :: VMR_tab_sh
  character(len=10), allocatable, dimension(:) :: sp_list
  character(len=200) :: data_file, sp_file, network,  net_dir, met
  real(dp), allocatable, dimension(:) :: mw, q0
  real(dp) :: p_CE

  !! Cloud variables
  integer :: ncld
  character(len=10), allocatable, dimension(:) :: cld_sp
  real(dp), allocatable, dimension(:) :: q0_cld



  integer :: u_nml

  !! Shortwave arrays and variabiles, not needed yet
  real(dp), allocatable, dimension(:) :: mu_z, tau_Ve, sw_a, sw_g
  real(dp) :: F0, sw_a_surf, asr, AB, lw_a_surf

  namelist /FMS_BD_RC_nml/ ts_scheme, opac_scheme, adj_scheme, diff_scheme, cloud_chem_scheme, chem_scheme, &
    & kappa_prime_calc, cloud_opac_scheme, &
    & nlay, a_sh, b_sh, pref, &
    & t_step, nstep, Rd_air, cp_air, grav, Tint, k_IR, met_R, &
    & iIC, corr, lw_ac, lw_gc, Bezier, lw_a_surf, &
    & &
    & nsp, ncld

  namelist /FMS_BD_RC_mc_nml/ VMR_tab_sh, sp_list, data_file, sp_file, network, net_dir, met, &
    & mw, p_CE

  namelist /FMS_BD_RC_cld_nml/ q0_cld, cld_sp

  !! Read input variables from main namelist
  open(newunit=u_nml, file='FMS_BD_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_BD_RC_nml)
  close(u_nml)

  !! Allocate and read variabiles for mini-chem namelist
  allocate(sp_list(nsp),mw(nsp))
  open(newunit=u_nml, file='FMS_BD_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_BD_RC_mc_nml)
  close(u_nml)

  !! Allocate and read variabiles for cloud scheme namelist
  allocate(cld_sp(ncld),q0_cld(ncld))
  open(newunit=u_nml, file='FMS_BD_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_BD_RC_cld_nml)
  close(u_nml)


  !! Number of layer edges (levels)
  nlev = nlay + 1

  !! Read in hybrid sigma grid values
  open(newunit=ua,file=trim(a_sh), action='read', status='old')
  open(newunit=ub,file=trim(b_sh), action='read', status='old')
  allocate(a(nlev),b(nlev))
  do k = 1, nlev
    read(ua,*) a(k)
    read(ub,*) b(k)
  end do
  close(ua); close(ub)

  !! Contruct pressure array [pa] at the levels using the hybrid sigma formula
  ! Reference surface pressure [pa] is pref
  allocate(pe(nlev))
  do k = 1, nlev
    pe(k) = a(k) + b(k)*pref
  end do
  pu = pe(1)

  !! Pressure at the layers
  allocate(pl(nlay),dpe(nlay))
  do k = 1, nlay
    dpe(k) = pe(k+1) - pe(k)
    pl(k) = dpe(k) / log(pe(k+1)/pe(k))
  end do

  !! Allocate other arrays we need
  allocate(Tl(nlay), dT_rad(nlay), dT_conv(nlay), net_F(nlev))
  allocate(tau_IRe(nlev), k_IRl(nlay))
  allocate(Kzz(nlay), mu(nlay), w_osh(nlay), Rd_bar(nlay), cp_bar(nlay), k_prime(nlay))
  Kzz(:) = 1e1_dp

  !! shortwave arrays (not needed here yet)
  allocate(mu_z(nlev),tau_Ve(nlev),sw_a(nlay),sw_g(nlay))

  !! Allocate cloud properties - assume constant for all layers in this test code
  !! Adjust these parts for your own situations
  allocate(lw_a(nlay), lw_g(nlay))
  lw_a(:) = lw_ac
  lw_g(:) = lw_gc

  allocate(vf(nlay), sat(nlay))
  vf(:) = 1e-30_dp

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp

  Fint = sb * Tint**4 ! Internal flux

  print*, 'Tint ', 'pref ', 'pu ', 'grav '
  print*, Tint, pref/1e5_dp, pu/1e5_dp, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_IRl(:) = k_IR

  select case(opac_scheme)
  case default
    tau_IRref = (k_IR * pref) / grav
  end select

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_IRl,Tint,grav,Tl,prc, &
  & kappa_air)

  !! Print initial T-p profile
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i)
  end do
  print*, '--------------'

  !! Prepare tracer scheme
  nq = nsp+ncld*2
  allocate(q(nlay,nq),q0(nq))
  q(:,:) = 0.0_dp
  q0(:) = 1e-30_dp

  !! Make CE initial conditions for chemistry
  do k = 1, nlay
    call CE_interpolate_Bezier(pl(k)/1e5_dp,Tl(k),sp_list(:),nsp,VMR_tab_sh,q(k,1:nsp),mu(k))
    q(k,1:nsp) = q(k,1:nsp)/sum(q(k,1:nsp))
    !print*, k, q(k,:), mu(k)
  end do

  !! Initialise mini-chem if needed
  select case(chem_scheme)
  case('mini-chem')
    call read_react_list(data_file, sp_file, net_dir, met)
  case default
  end select  


  !! Recaclulate mu for each layer
  do k = 1, nlay
    mu(k) = 0.0_dp
    do i = 1, nsp
      mu(k) = mu(k) + q(k,i) * mw(i)
    end do
  end do

  !! Write out initial conditions
  open(newunit=u,file='FMS_RC_ic.out',action='readwrite')
  write(u,*) nlay, nsp, ncld 
  do k = 1, nlay
    write(u,*) k, pl(k), Tl(k), mu(k), q(k,:)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep


    select case(chem_scheme)
      case('mini-chem')
        do k = 1, nlay
          if (pl(k) > p_CE) then
            call CE_interpolate_Bezier(pl(k)/1e5_dp,Tl(k),sp_list(:),nsp,VMR_tab_sh,q(k,1:nsp),mu(k))
          else
            q(k,1:nsp) = q(k,1:nsp)/sum(q(k,1:nsp))
            call mini_ch_seulex(Tl(k), pl(k), t_step, q(k,1:nsp), network)
          end if
          q(k,1:nsp) = q(k,1:nsp)/sum(q(k,1:nsp))
        end do
      case('CE') 
        do k = 1, nlay
          call CE_interpolate_Bezier(pl(k)/1e5_dp,Tl(k),sp_list(:),nsp,VMR_tab_sh,q(k,1:nsp),mu(k))
          q(k,1:nsp) = q(k,1:nsp)/sum(q(k,1:nsp))
        end do
      case default
    end select

    !! Recaclulate mu for each layer
    do k = 1, nlay
      mu(k) = 0.0_dp
      do i = 1, nsp
        mu(k) = mu(k) + q(k,i) * mw(i)
      end do
    end do

    select case(kappa_prime_calc)
    case('NASA9')
      do k = 1, nlay
        call calc_kappa(nsp, sp_list, Tl(k), q(k,1:nsp), mu(k), Rd_bar(k), cp_bar(k), k_prime(k))
        !print*, k, Rd_bar(k), cp_bar(k), k_prime(k)
      end do
    case('Constant')
      Rd_bar(:) = Rd_air
      cp_bar(:) = cp_air
      k_prime(:) = kappa_air
    end select

    select case(cloud_chem_scheme)
    case('Tracer')
      call BD_clouds_chem(nlay, t_step, cld_sp(1), q(:,nsp+1), q(:,nsp+2), pl, Tl, Rd_bar, Kzz, q0_cld(1), grav, sat)
      call BD_clouds_vf(nlay, nsp, Rd_bar, grav, q(:,1:nsp), pl, Tl, q(:,nsp+2), vf)
      call BD_clouds_adv(nlay, 1, t_step, q(:,nsp+2), vf, dpe)
    case('None')
    case default
      print*, 'Invalid cloud_scheme: ', trim(cloud_chem_scheme)
      stop
    end select

    !! Convective adjustment scheme
    select case(adj_scheme)
    case('MLT')
      ! Mixing length theory
      call BD_MLT(nlay, nlev, t_step, Tl, pl, pe, Rd_bar, cp_bar, k_prime, &
         & grav, dT_conv, Kzz, w_osh)
    case('None')
      dT_conv(:) = 0.0_dp
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    !! Tracer diffusion scheme
    select case(diff_scheme)
    case('Explicit')
      ! 1st order explicit scheme
      q0(:) = q(nlay,:)
      call BD_vert_diff(nlay, nlev, t_step, Rd_bar, grav, Tl, pl, pe, Kzz, nq-1, q(:,1:nq-1), q0(1:nq-1))
    case('None')
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    select case(opac_scheme)
    case('Freedman')
      ! Calculate optical depth structure for Freedman et al. (2014) Rosseland mean fitting function scheme
      call k_Ross_Freedman(Tl(1), pe(1), met_R, k_IRl(1))
      k_IRl(1) = max(k_IRl(1),1.0e-5_dp)
      tau_IRe(1) = (k_IRl(1) * pe(1)) / grav
      do k = 1, nlay
        call k_Ross_Freedman(Tl(k), pl(k), met_R, k_IRl(k))
        k_IRl(k) = max(k_IRl(k),1.0e-5_dp)
        tau_IRe(k+1) = tau_IRe(k) + (k_IRl(k) * dpe(k)) / grav
      end do
    case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select

    select case(cloud_opac_scheme)
    case('Rosseland')

    case('None')

    case default
      print*,'Invalid cloud opacity scheme: ', cloud_opac_scheme
    end select 


    !! Two stream radiative transfer step
    select case(ts_scheme)
    case('Toon_scatter')
      tau_Ve(:) = 0.0_dp
      mu_z = 0.0_dp
      F0 = 0.0_dp
      AB = 0.0_dp
      sw_a(:) = 0.0_dp
      sw_g(:) = 0.0_dp
      sw_a_surf = 0.0_dp
      call ts_Toon_scatter(Bezier, nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
        & sw_a, sw_g, lw_a, lw_g, sw_a_surf, lw_a_surf, net_F, olr, asr)
    case('None')
      net_F(:) = 0.0_dp
    case default
      print*, 'Invalid ts_scheme: ', trim(ts_scheme)
      stop
    end select

    !! Calculate the temperature tendency due to radiation
    do k = 1, nlay
      dT_rad(k) = (grav/cp_bar(k)) * (net_F(k+1)-net_F(k))/(dpe(k))
      !print*, n, i, dT_rad(i)
    end do

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Check for NaN's in the temperature and exit the simulation if detected
    do i = 1, nlay
      if (ieee_is_nan(Tl(i)) .eqv. .True.) then
        do j = 1, nlay
          print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
        end do
        print*, nlev, net_F(nlev), olr
        inan = 1
        exit
      end if
    end do
    if (inan == 1) then
      exit
    end if

    !! Increase the total time simulated
    t_tot = t_tot + t_step

  end do

  !! cpu timer end
  call cpu_time(finish)

  !! Output the results
  print*, 'sec: ', 'hours: ', 'days: '
  print*, t_tot, t_tot/60.0_dp/60.0_dp,t_tot/60.0_dp/60.0_dp/24.0_dp

  print*, 'For profile properties: '
  print*, 'Tint: ', 'Pref: '
  print*, Tint, pref

  print*, 'OLR [W m-2], T_b [K]:'
  print*, olr, (olr/sb)**(0.25_dp)

  print*, 'Outputting results: '
  open(newunit=u,file='FMS_RC_pp_1.out', action='readwrite')
  write(u,*) nlay, nsp, ncld 
  do k = 1, nlay
    write(u,*) k, pl(k), Tl(k), mu(k), q(k,:)
  end do
  close(u)

  open(newunit=u,file='FMS_RC_pp_2.out', action='readwrite')
  write(u,*) nlay, nsp, ncld 
  do k = 1, nlay
    write(u,*) k, pl(k), Tl(k), Rd_bar(k), cp_bar(k), k_prime(k), Kzz(k), w_osh(k)
  end do
  close(u)

  print*, n, 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start

end program FMS_BD_RC
