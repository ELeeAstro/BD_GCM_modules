module BD_vert_diff_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: CFL = 0.40_dp
  real(dp), parameter :: kb = 1.380649e-23_dp

  public :: BD_vert_diff
  private :: linear_interp

contains

  subroutine BD_vert_diff(nlay, nlev, t_end, Rd_air, grav, Tl, pl, pe, Kzz, nq, q, q0)
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav
    real(dp), dimension(nlay), intent(in) :: Tl, pl, Kzz, Rd_air
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q

    integer :: k
    real(dp) :: h1, h2, d1Kzz
    real(dp), dimension(nlev) :: alte, alte_r, lpe, Kzz_e, Kzz_er, Te, nde, nde_r
    real(dp), dimension(nlay) :: dh

    real(dp), dimension(nlay) :: lKzz, Kzz_r, lpl, lTl, nd, nd_r
    real(dp), dimension(nlay,nq) :: q_r, flx
    real(dp), dimension(nq) :: phi_jp, phi_jm

    integer :: n_it
    real(dp) :: dt, t_now

    !! First interpolate Kzz at the layer to Kzz at the edges
    !! Log pressure grid
    lpe(:) = log10(pe(:))
    lpl(:) = log10(pl(:))
    lTl(:) = log10(Tl(:))
    lKzz(:) = log10(Kzz(:))
    !! Perform interpolation to level edges using log-linear
    !interpolation
    do k = 2, nlay
      call linear_interp(lpe(k), lpl(k-1), lpl(k), lKzz(k-1), lKzz(k), Kzz_e(k))
      Kzz_e(k) = 10.0_dp**Kzz_e(k)
    end do
    
    Kzz_e(1) = Kzz(1)
    Kzz_e(nlev) = Kzz(nlay)

    do k = 2, nlay
      call linear_interp(lpe(k), lpl(k-1), lpl(k), lTl(k-1), lTl(k), Te(k))
    end do
    ! Edges are linearly interpolated
    Te(1) = (log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/10.0_dp**Te(2)))
    Te(nlev) = (log10(Tl(nlay)) &
    &  + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/10.0_dp**Te(nlay)))
    ! De-log the temperature levels (edges)
    Te(:) = 10.0_dp**Te(:)

    ! Number density at layers and levels
    nd(:) = pl(:)/(kb * Tl(:))
    nde(:) = pe(:)/(kb * Te(:)) 

    !! First calculate the vertical height (m) assuming hydrostatic equilibrium
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (Rd_air(k)*Tl(k))/grav * log(pe(k+1)/pe(k))
    end do

    !! Now reverse alte, Kzz, Kzz_e and q so indexing starts at 1->nlay from 0 altitude
    do k = 1, nlev
      alte_r(k) = alte(nlev-k+1)
      Kzz_er(k) = Kzz_e(nlev-k+1)
      nde_r(k) = nde(nlev-k+1)
    end do
    do k = 1, nlay
      Kzz_r(k) = Kzz(nlay-k+1)
      q_r(k,:) = q(nlay-k+1,:)
      nd_r(k) = nd(nlay-k+1)
    end do

    !! We now follow the 1st order VULCAN scheme

    !! Altitude differences
    do k = 1, nlay
      dh(k) = alte_r(k+1) - alte_r(k)
    enddo

    !! Prepare timestepping routine
    !! Find minimum timestep that satifies the CFL condition
    dt = t_end
    do k = 1, nlay
      dt = min(dt,CFL*(alte_r(k+1)-alte_r(k))**2/Kzz_r(k))
    end do

    !! Begin timestepping routine
    t_now = 0.0_dp
    n_it = 1

    do while ((t_now < t_end) .and. (n_it < 100000))

      !! If next time step overshoots - last time step is equal tend
      if ((t_now + dt) > t_end) then
        dt = t_end - t_now
      end if

      !! Apply tracer lower boundary conditions
      q_r(1,:) = q0(:)
      q_r(nlay,:) = q_r(nlay-1,:)


      !! Find flux between layers
      do k = 2, nlay-1
        phi_jp(:) = -Kzz_er(k+1) * nde_r(k+1) * (q_r(k+1,:) - q_r(k,:))/dh(k)
        phi_jm(:) = -Kzz_er(k) * nde_r(k) * (q_r(k,:) - q_r(k-1,:))/dh(k)
        flx(k,:) = -(phi_jp(:) - phi_jm(:))/dh(k)
        !print*, k,   phi_jp(:),  phi_jm(:), flx(k,:), dh(k)
      end do

      !! Perform tracer timestepping
      do k = 2, nlay-1
        q_r(k,:) = q_r(k,:) + flx(k,:)/nd_r(k)*dt
      end do

      !! Apply tracer lower boundary conditions
      q_r(1,:) = q0(:)
      q_r(nlay,:) = q_r(nlay-1,:)

      !! Increment time and iteration number
      t_now = t_now + dt
      n_it = n_it + 1

    end do

    !! Reverse q_r and pass back to q
    do k = 1, nlay
      q(k,:) = q_r(nlay-k+1,:)
    end do

  end subroutine BD_vert_diff

  ! Perform linear interpolation
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

end module BD_vert_diff_mod

