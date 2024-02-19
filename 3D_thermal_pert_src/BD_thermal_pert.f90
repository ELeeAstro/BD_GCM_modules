! E.K.H Lee - Feb 2023
! Exo-FMS version of Showman and Tan's boundary thermal perturbation
!

module thermal_perturb
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64
  integer, parameter :: i64 = INT64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

  integer, parameter :: mforce = 10
  integer, parameter :: nforce = 29
  integer, parameter :: delta_n = 2

  integer, parameter :: mmax = 40
  integer, parameter :: nmax = 40

  integer, parameter :: Nrforctop = 40
  real(dp), parameter :: t_storm = 1e5_dp
  real(dp), parameter :: t_amp = 1.5e-4_dp
  real(dp), parameter :: p_rcb = 1e6_dp

  real(dp), allocatable, dimension(:,:) :: bforce, bturb
  real(dp), allocatable, dimension(:,:,:) :: PM

  integer :: Nspecmodes

  integer, parameter :: iseed = -10
  integer(i64), dimension(4) :: u4
  integer, parameter :: nburn = 100

  logical :: first_call = .True.

  public :: thermal_pert
  private :: LPMN_NORM, set_seed, ran2, rotl

contains

  subroutine thermal_pert(nlat, nlon, npz, lat, lon, pl, t_step, thermpert)
    implicit none

    integer, intent(in) :: nlat, nlon, npz

    real(dp), intent(in) :: t_step
    real(dp), dimension(nlat), intent(in) :: lat
    real(dp), dimension(nlon), intent(in) :: lon
    real(dp), dimension(nlat, nlon, npz), intent(in) :: pl

    real(dp), dimension(nlat, nlon, npz), intent(out) :: thermpert

    integer :: i, j, k, n, m
    real(dp) :: r, slat
    real(dp) :: rphase, thet_over_t
    real(dp), dimension(npz) :: forc_vert_profile

    r = 1.0_dp - t_step/t_storm

    if (first_call .eqv. .True.) then
      allocate(bforce(nlat,nlon))
      allocate(bturb(nlat,nlon))

      bforce(:,:) = 0.0_dp

      call set_seed()

      allocate(PM(nlat,0:mmax,0:nmax))
      do i = 1, nlat
        slat = sin(lat(i))
        call LPMN_NORM(mmax, mmax, nmax, slat, PM(i,:,:))
      end do

      Nspecmodes = 0
      do n = nforce,nforce+delta_n-1
        Nspecmodes = Nspecmodes + n
      end do

      first_call = .False.

      ! We need to 'burn in' the baseline thermal perturbation
      ! to make the map truly random at the start of the simulation
      do i = 1, nburn

        bturb(:,:) = 0.0_dp

        do n = nforce, nforce+delta_n-1
          do m = 1, n
            rphase = ran2() * 2.0_dp * pi
            do j = 1, nlon
              bturb(:,j) = bturb(:,j) + PM(:,m,n) * cos(real(m,dp) * (lon(j) + rphase))
            end do
          end do
        end do

        bforce(:,:) = r*bforce(:,:) + sqrt(1.0_dp - r**2) * (t_amp*bturb(:,:))

      end do

    end if

    bturb(:,:) = 0.0_dp

    do n = nforce, nforce+delta_n-1
      do m = 1, n
        rphase = ran2() * 2.0_dp * pi
        do j = 1, nlon
          bturb(:,j) = bturb(:,j) + PM(:,m,n) * cos(real(m,dp) * (lon(j) + rphase))
        end do
      end do
    end do

    bforce(:,:) = r*bforce(:,:) + sqrt(1.0_dp - r**2) * (t_amp*bturb(:,:))

    do i = 1, nlat
      do j = 1, nlon
        do k = 1, npz
          thet_over_t = 1.0_dp!(p0/pl(i,j,k))**kappa
          if (pl(i,j,k) <= p_rcb) then
            forc_vert_profile(k) = (pl(i,j,k)/p_rcb)**2
          else
            forc_vert_profile(k) = (p_rcb/pl(i,j,k))**2
          endif
          if (pl(i,j,k) <= p_rcb*7.389_dp .or. pl(i,j,k) >= p_rcb/7.389_dp) then
            thermpert(i,j,k) = thet_over_t * forc_vert_profile(k) * bforce(i,j)
          else
            thermpert(i,j,k) = 0.0_dp
          end if
        end do
      end do
    end do

  end subroutine thermal_pert


  subroutine LPMN_NORM(mm,m,n,x,PM)
    implicit none
! C       =====================================================
! C       Purpose: Compute the normalized associated Legendre functions 
! C                Pmn(x) over the full grid of m and n at a 
! C                specific value of argument x.  The 2D array
! C                of Pmn (at the given x) is returned in PM.
! C
! C       
! C       Input :  x  --- Argument of Pmn(x)
! C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
! C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
! C                mm --- Physical dimension of PM and PD
! C       Output:  PM(m,n) --- Gives associated Legendre 
! C                        polynomial for degree n and order m,
! C                        EVALUATED AT X, i.e. Pmn(x).
! C
! C       Built by Adam Showman from the routine mlpmn.f that calculates
! C       the non-normalized associated Legendre functions... obtained 
! C       from the website of Jianming Jin (Illinois) who provides
! C       it for general use; original routine copyright by him.
! C       From his book "Computation of Special Functions".
! C       (Original routine also calculated the derivative of 
! C       Pmn but I don't need that so I removed it.)
! C
! C       Note: the original version of this code had a bug:
! C       the array bounds on the "DO 35" line was originally
! C       I=0,M.  But since that loop sets PM(I,I+1) it means
! C       that PM can be accessed beyond where it is defined
! C       in the second element.  One solution would be to
! C       change array bounds to I=0,M-1.  This would work fine
! C       if N=M (triangular truncation) but not if N>M
! C       (rhomboidal truncation) which would then leave some
! C       gaps where certain Pmn were not evaluated.  To 
! C       fix the bug while allowing for N>M, 
! C       I kept the original array bounds but simply added a
! C       conditional to only assign PM if I+1 <= N.
! C       =====================================================

    integer, intent(in) :: mm, m, n
    real(dp), intent(in) :: x

    real(dp), dimension(0:mm,0:n), intent(out) :: PM

    integer :: i, j
    real(dp) :: ls, xq, xs

    do i = 0,n
      do j = 0,m
        PM(j,i) = 0.0_dp
      end do
    end do
    PM(0,0) = 0.5_dp*sqrt(2.0_dp)

    if (abs(x) == 1.0_dp) then
      do i = 1,n
        PM(0,i) = sqrt(real(i,dp)+0.5_dp) * x**i
      end do
      return
    end if

    ls = 1.0_dp
    if (abs(x) > 1.0_dp) then
      ls = -1.0_dp
    end if
    xq = sqrt(ls*(1.0_dp - x**2))
    xs = ls*(1.0_dp-x**2)
    do i = 1,m
      PM(i,i) = ls * &
        & sqrt((2.0_dp*real(i,dp)+1.0_dp)/(2.0_dp*real(i,dp))) * &
        & xq*PM(i-1,i-1)
    end do
    do i = 0,m
! c          Showman added the following if conditional to
! c          prevent PM from being accessed outside of array
! c          bounds when N does not exceed M.
      if (i+1 <= n) then
! c     PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        PM(i,i+1) = sqrt(2.0_dp*real(i,dp)+3.0_dp) * x * PM(i,i)
      endif
    end do
    do i = 0,m
      do j = i+2, n
! c      PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
! c        & (I+J-1.0D0)*PM(I,J-2))/(J-I)
        PM(i,j) = sqrt((2.0_dp*real(j,dp)+1.0_dp)*(2.0_dp*real(j,dp)-1.0_dp) / &
          & real((j-i)*(j+i),dp))*x*PM(i,j-1) - &
          & sqrt((2.0_dp*real(j,dp)+1.0_dp)*(real(j,dp)-real(i,dp)-1.0_dp)*(real(j,dp)+real(i,dp)-1.0_dp) / &
          & ((2.0_dp*real(j,dp)-3.0_dp)*(real(j,dp)-real(i,dp))*(real(j,dp)+real(i,dp))))*PM(i,j-2)
      end do
    end do

  end subroutine LPMN_NORM

  subroutine set_seed()
    ! Note: this should be called with a different seed in each thread or
    ! process. But the following is usually good enough.
    ! We use the xoshiro256+ random number generator, which again, is usually good enough for this application.
    implicit none

    integer(i64) :: iseed_i64

    iseed_i64 = iseed
    u4(1) = 101 * iseed_i64
    u4(2) = 303 * iseed_i64
    u4(3) = 505 * iseed_i64
    u4(4) = 707 * iseed_i64

  end subroutine set_seed
  
  function ran2()
    implicit none

    real(dp) :: ran2

    integer(i64) :: res, t
    real(dp) :: rl

    integer(i64), dimension(4) :: s4

    s4(:) = u4(:)

    t = shiftl(s4(2),17)

    s4(3) = ieor(s4(3),s4(1))
    s4(4) = ieor(s4(4),s4(2))
    s4(2) = ieor(s4(2),s4(3))
    s4(1) = ieor(s4(1),s4(4))

    s4(3) = ieor(s4(3),t)

    s4(4) = rotl(s4(4), 45)

    res = s4(1) + s4(4)
    !res = rotl(s4(1) + s4(4),23) + s4(1)

    res = ior(shiftl(1023_i64, 52), shiftr(res, 12))

    ran2 = transfer(res,rl) - 1.0_dp

    u4(:) = s4(:)

  end function ran2

  function rotl(x, k)
    implicit none

    integer(i64), intent(in) :: x
    integer, intent(in) :: k

    integer(i64) :: rotl

    rotl =  ior(shiftl(x,k), shiftr(x,(64 - k)))

  end function rotl

end module thermal_perturb

program test_pert
  use thermal_perturb
  implicit none

  integer, parameter :: nlat = 96
  integer, parameter :: nlon = 96
  integer, parameter :: npz = 1

  double precision ::  t_step
  double precision, dimension(nlat) :: lat 
  double precision, dimension(nlon) :: lon
  double precision, dimension(nlat, nlon, npz) :: pl

  double precision, dimension(nlat, nlon, npz) :: thermpert

  integer :: i, j, n, n_pert, u
  double precision :: dlat, dlon

  dlat = pi/real(nlat+1,dp)
  dlon = (2.0d0 * pi)/real(nlon+1,dp)

  lat(1) = -pi/2.0d0 + dlat
  do j = 2, nlat
    lat(j) = lat(j-1) + dlat
  end do

  lon(1) = dlon
  do i = 2, nlon
    lon(i) = lon(i-1) + dlon
  end do

  pl(:,:,:) = 1e5

  t_step = 30d0

  n_pert = 100

  do n = 1, n_pert
    call thermal_pert(nlat, nlon, npz, lat, lon, pl, t_step, thermpert)
  end do

  open(newunit=u,file='bforce.txt',action='readwrite')
  write(u,*) nlat, lon(:)
  do j = 1, nlat
    write(u,*) lat(j), bforce(j,:) 
  end do
    
  close(u)

end program test_pert