module BD_kappa_poly_mod
  use, intrinsic :: iso_fortran_env
  implicit none
 
  integer, parameter :: dp = REAL64

  type species
    integer :: id, n_a
    character(len=10) :: c
    real(dp) :: mw, nd
    real(dp), allocatable, dimension(:) :: a_l, a_h
  end type species

  
  type(species), allocatable, dimension(:) :: g_sp
  real(dp), parameter :: R = 8.31446261815324_dp
  logical :: first_call = .True.

  integer :: n_sp
  integer, allocatable, dimension(:) :: idx


  public :: calc_kappa
  private :: read_poly, first_call

  contains

  subroutine calc_kappa(n_mol, m_list, T_in, VMR, mmw, R_bar, cp_bar, k_prime)
    implicit none

    integer, intent(in) :: n_mol
    character(len=10), dimension(n_mol), intent(in) :: m_list
    real(dp), intent(in) :: T_in, mmw
    real(dp), dimension(n_mol), intent(in) :: VMR

    real(dp), intent(out) :: R_bar, cp_bar, k_prime

    integer :: n, i
    real(dp) :: mmr, cp_val, T, T2, T3, T4

    if (first_call .eqv. .True.) then
      call read_poly()
      first_call = .False.

      ! Find index of each input species
      allocate(idx(n_mol))
      idx(:) = 0 

      do n = 1, n_mol

        do i = 1, n_sp
          if (trim(m_list(n)) == trim(g_sp(i)%c)) then
            idx(n) = i
          end if
        end do

        if (idx(n) == 0) then
          print*, 'idx = 0, species not found in polynomial list: ', trim(m_list(n))
          stop
        end if

      end do
  
    end if

    T = T_in
    T2 = T**2
    T3 = T**3
    T4 = T**4

    !! Main loop calculations
    R_bar = 0.0_dp
    cp_bar = 0.0_dp
    do n = 1, n_mol

      i = idx(n)

      ! Mass mixing ratio = VMR * molecular weight / mean molecular weight
      mmr = VMR(n) * g_sp(i)%mw/mmw

      ! Contribution to R_bar = mmr * specific gas constant
      R_bar = R_bar + mmr * R/g_sp(i)%mw

      if (T <= 1000.0_dp) then
        cp_val = g_sp(i)%a_l(1)/T2 + g_sp(i)%a_l(2)/T + g_sp(i)%a_l(3) + g_sp(i)%a_l(4)*T + &
          & g_sp(i)%a_l(5)*T2 + g_sp(i)%a_l(6)*T3 + g_sp(i)%a_l(7)*T4
      else
        cp_val = g_sp(i)%a_h(1)/T2 + g_sp(i)%a_h(2)/T + g_sp(i)%a_h(3) + g_sp(i)%a_h(4)*T + &
          & g_sp(i)%a_h(5)*T2 + g_sp(i)%a_h(6)*T3 + g_sp(i)%a_h(7)*T4
      end if

      cp_val = cp_val*R

      ! Contribution to cp_bar = mmr * cp_val / molecular weight
      cp_bar = cp_bar + mmr * cp_val/g_sp(i)%mw

    end do

    ! Convert R_bar [J g-1 K-1] to SI units [J kg-1 K-1]
    R_bar = R_bar * 1000.0_dp

    ! Convert cp_bar [J g-1 K-1] to SI units [J kg-1 K-1]
    cp_bar = cp_bar * 1000.0_dp

    ! kappa_prime evaluation
    k_prime = R_bar / cp_bar


  end subroutine calc_kappa


  subroutine read_poly()
    implicit none

    integer :: u, i

    open(newunit=u,file='data/NASA_9_poly.txt',action='read',status='old',form='formatted')

    do i = 1, 6
      read(u,*)
    end do
    read(u,*) n_sp
    read(u,*)

    allocate(g_sp(n_sp))

    do i = 1, n_sp
      g_sp(i)%id = i
      read(u,*) g_sp(i)%c, g_sp(i)%mw, g_sp(i)%n_a
      allocate(g_sp(i)%a_l(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_l(:)
      allocate(g_sp(i)%a_h(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_h(:)
      
      ! print*, g_sp(i)%id, g_sp(i)%c, g_sp(i)%mw,  g_sp(i)%n_a
      ! print*, g_sp(i)%a_l(:)
      ! print*, g_sp(i)%a_h(:)
    end do

    close(u)

  end subroutine read_poly

end module BD_kappa_poly_mod


! program test_kappa_poly
!   use BD_kappa_poly_mod!, only : clac_kappa
!   implicit none

!   integer, parameter :: n_mol = 3
!   character(len=4), dimension(n_mol) :: m_list
!   double precision :: T_in, mmw, R_bar, cp_bar, k_prime
!   double precision, dimension(n_mol) :: VMR
!   logical :: restart

!   ! Test for Earth air - gets very close to the wikipedia quoted values for cp air and R air
!   ! m_list = (/'N2  ','O2  ','Ar  ', 'CO2 ', 'H2O '/)
!   ! VMR = (/0.78d0,0.21d0,0.0093d0,365d-6,1d-4/)
!   ! T_in = 298.15d0
!   ! mmw = 28.9645d0

!   m_list = (/'H2','He','H '/)
!   VMR = (/0.854d0,0.145d0,9.70d-6/)
!   T_in = 1000.0d0
!   mmw = 2.33d0

!   call calc_kappa(n_mol, m_list, T_in, VMR, mmw, R_bar, cp_bar, k_prime)

!   print*, m_list(:)
!   print*, T_in, mmw, VMR(:)
!   print*, R_bar, cp_bar, k_prime

! end program test_kappa_poly
