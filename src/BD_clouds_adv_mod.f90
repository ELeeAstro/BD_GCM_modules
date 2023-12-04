module BD_clouds_adv_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: cfl = 0.95_dp
  real(dp), parameter :: kb = 1.380649e-23_dp

  public :: BD_clouds_adv
  private :: minmod

  contains
 
  subroutine BD_clouds_adv(nlay, nlev, Rd_air, grav, Tl, pl, pe, nq, tend, q, vf)
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), dimension(nlay,nq), intent(inout) :: q
    real(dp), dimension(nlay), intent(in) :: vf, Rd_air, Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), intent(in) :: tend, grav

    real(dp), dimension(nlay,nq) :: qc
    real(dp), dimension(nlay) :: sig, c, vf_c, delz, delz_mid
    real(dp), dimension(nlev) :: alte, vf_e
    real(dp) :: tnow, dt
    integer :: k, iit, n

    ! Convert to m s-1
    vf_c(:) = vf(:)/100.0_dp 

    ! Find velocity at levels
    vf_e(1) = vf_c(1)
    do k = 2, nlev
      vf_e(k) = (vf_c(k-1) + vf_c(k))/2.0_dp
    end do

    !! First calculate the vertical height (m) assuming hydrostatic equilibrium
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (Rd_air(k)*Tl(k))/grav * log(pe(k+1)/pe(k))
      delz(k) = alte(k) - alte(k+1)
    end do

    !! Find differences between layers directly
    do k = 1, nlay-1
      delz_mid(k) = (alte(k) + alte(k+1))/2.0_dp - (alte(k+1) + alte(k+2))/2.0_dp
    end do
    delz_mid(nlay) = delz(nlay)

    !! Find minimum timestep that allows the CFL condition
    dt = tend
    do k = 1, nlay
      dt = min(dt,cfl*(delz_mid(k))/abs(vf_e(k+1)))
    enddo

    do while ((tnow < tend) .and. (iit < 1000000))

      ! If next time step overshoots - last time step is equal tend
      if (tnow + dt > tend) then
        dt = tend - tnow
      end if

      !! Find the courant number
      c(:) = abs(vf_e(2:nlev)) * dt / delz_mid(:)

      do n = 1, nq

        !! Find the minmod limiter
        call minmod(nlay,q(:,n),delz_mid,sig)

        !! Perform McCormack step
        qc(:,n) = q(:,n)
        qc(:nlay-1,n) = q(:nlay-1,n) - sig(:nlay-1)*c(:nlay-1)*(q(2:nlay,n) - q(:nlay-1,n))
        q(2:nlay,n) = 0.5_dp * (q(2:nlay,n) + qc(2:nlay,n) - c(2:nlay)*(qc(2:nlay,n) - qc(:nlay-1,n)))

        q(:,n) = max(q(:,n),1e-30_dp)

      end do

      tnow = tnow + dt
      iit = iit + 1

      q(nlay,:) = 1e-30_dp

    end do

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

end module BD_clouds_adv_mod