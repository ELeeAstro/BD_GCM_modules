module BD_MLT_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: alp = 1.0_dp, beta = 2.2_dp
  real(dp), parameter :: dt_max = 0.5_dp
  real(dp), parameter :: Kzz_min = 1e1_dp ! Minimum Kzz (typically 1e5 cm2 s-1)
  real(dp), parameter :: Kzz_max = 1e8_dp ! Maximum Kzz

  real(dp), parameter :: f1 = 1.0_dp/8.0_dp, f2 = 1.0_dp/2.0_dp

  public :: BD_MLT
  private :: linear_interp

contains

  subroutine BD_MLT(nlay, nlev, t_end, Tl, pl, pe, Rd_air, cp_air, kappa_air, &
    & grav, dT_conv, Kzz, w_osh)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Rd_air, cp_air, kappa_air
    real(dp) ::  grav, t_end
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), dimension(nlev), intent(in) :: pe

    real(dp), dimension(nlay), intent(in) :: Tl
    real(dp), dimension(nlay), intent(out) :: dT_conv, Kzz, w_osh

    integer :: k
    real(dp) :: Hp, L, rho, t_now, dt, w, w_over, p_rcb
    real(dp), dimension(nlay) :: gradrad, Fc, lpl, Tl_c, lTl_c
    real(dp), dimension(nlev) :: lpe, Fc_e, Te

    !! Save the input temperature into copy
    Tl_c(:) = Tl(:)

    !! Loop over time until we reach t_end
    t_now = 0.0_dp
    do while (t_now < t_end)

      !! Calculate the timestep
      if ((t_now + dt_max) >= t_end) then
        dt = t_end - t_now
      else
        dt = dt_max
      end if

      do k = 2, nlay
        call linear_interp(pe(k), pl(k-1), pl(k), Tl_c(k-1), Tl_c(k), Te(k))
      end do
      ! Edges are linearly interpolated
      Te(1) = 10.0_dp**(log10(Tl_c(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl_c(1)/Te(2)))
      Te(nlev) = 10.0_dp**(log10(Tl_c(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl_c(nlay)/Te(nlay)))


      !! Calculate the regions that are convectivly unstable
      do k = nlay, 1, -1 
        gradrad(k) = log(Te(k)/Te(k+1))/log(pe(k)/pe(k+1))
      end do

      !! Follow the Joyce & Tayar (2023), Marley & Robinson (2015) and
      !Robinson & Marley (2014) formalisims
      do k = 1, nlay
        if (gradrad(k) > kappa_air(k)) then
          ! Convective region, calculate convective flux
          ! Pressure scale height and mixing length
          Hp = (Rd_air(k) * Tl_c(k)) / grav
          L = alp * Hp
          ! Buoyancy characteristic vertical velocity Robinson & Marley (2014)
          w = (f1 * L**2 * grav * (gradrad(k) - kappa_air(k))) / Hp
          w = sqrt(w)
          ! Vertical convective thermal flux (Joyce & Tayar 2023)
          rho = pl(k)/(Rd_air(k) * Tl_c(k))
          Fc(k) = f2 * (rho * cp_air(k) * w * Tl_c(k) * L * (gradrad(k) - kappa_air(k)))/Hp
        else
          ! Convective flux is not adjusted
          Fc(k) = 0.0_dp
        end if
      end do

      !! Perform interpolation to level edges using log-linear
      !interpolation
      do k = 2, nlay
        call linear_interp(pe(k), pl(k-1), pl(k), Fc(k-1), Fc(k), Fc_e(k))
      end do

      !! Edges are linearly extrapolted (or equal 0)
      Fc_e(1) = 0.0_dp !(Fc(1) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * (Fc(1) - Fc_e(2)))
      Fc_e(nlev) = 0.0_dp !(Fc(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * (Fc(nlay) - Fc_e(nlay))

      !! Calculate the temperature tendency due to radiation
      do k = 1, nlay
        dT_conv(k) = (grav/cp_air(k)) * (Fc_e(k+1)-Fc_e(k))/(pe(k+1) - pe(k))
        !print*, n, i, dT_conv(k)
      end do

      !! Forward march the temperature change copy   
      Tl_c(:) = Tl_c(:) + dt * dT_conv(:)

      !! Add dt to the time
      t_now = t_now + dt

    end do

    do k = 2, nlay
      call linear_interp(pe(k), pl(k-1), pl(k), Tl_c(k-1), Tl_c(k), Te(k))
    end do
    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl_c(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl_c(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl_c(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl_c(nlay)/Te(nlay)))

    !! Perform Kzz calculation after last iteration is compelte 
    !! so we don't calculate it inside the loop
    do k = nlay, 1, -1 
      gradrad(k) = log(Te(k)/Te(k+1))/log(pe(k)/pe(k+1))
    end do

    !! Find rcb
    p_rcb = pl(nlay)
    do k = nlay, 1, -1
      if (gradrad(k) > kappa_air(k)) then
        p_rcb = pl(k)
      else
        p_rcb = pl(k)
        exit
      end if
    end do

    !! Follow the Joyce & Tayar (2023), Marley & Robinson (2015) and Robinson & Marley (2014) formalisims
    w_osh(:) = 0.0_dp
    Kzz(:) = Kzz_min
    do k = nlay, 1, -1

      if (gradrad(k) > kappa_air(k)) then
        ! Convective region, calculate Kzz
        ! Pressure scale height and mixing length
        Hp = (Rd_air(k) * Tl_c(k)) / grav
        L = alp * Hp
        ! Buoyancy characteristic vertical velocity Robinson & Marley (2014)
        w = (f1 * L**2 * grav * (gradrad(k) - kappa_air(k))) / Hp
        w = sqrt(w)
        ! Kzz calculation
        Kzz(k) = max(w * L, Kzz_min)
      end if

      if (pl(k) <= p_rcb) then
        if (Kzz(k) == Kzz_min) then
          ! Kzz is given by convective overshoot following Woitke & Helling (2004), Woitke et al. (2020)
          w_over = exp(log(w) - beta * max(0.0_dp, log(p_rcb/pl(k))))
          Hp = (Rd_air(k) * Tl_c(k)) / grav
          Kzz(k) = max(Hp * w_over, Kzz_min)
          w_osh(k) = w_over
        end if
      end if

      !! TODO: Calculate molecular diffusion and add to Kzz

    end do

    !! Make sure Kzz is smaller than the maximum value
    Kzz(:) = min(Kzz(:),Kzz_max)

    !! Pass back the net temperature tendency
    dT_conv(:) = (Tl_c(:) - Tl(:))/t_end

  end subroutine BD_MLT

  ! Perform linear interpolation
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

end module BD_MLT_mod
