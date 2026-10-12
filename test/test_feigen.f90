! SPDX-License-Identifier: MIT
program test_feigen
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures, only: cubature, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs, only: shapefunc, SHP_LAGRANGE, SHP_SERENDIPITY
  use feigen, only: material, stiffness, modes, name_modes, plot_modes
  implicit none

  real(dp), parameter :: eps = epsilon(1._dp)
  real(dp), parameter :: tol = 2**19*eps   ! eigenvalues and eigenvectors of K
  integer :: ncheck = 0, nfail = 0
  real(dp), allocatable :: lambda(:)

  ! Material ------------------------------------------------------------
  call check(all(material(1, 2._dp, 0.3_dp) == 2._dp), 'material 1D: D = E')
  call check(close([material(2, 1._dp, 0._dp)], &
                 & [1._dp, 0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 0.5_dp]), &
           & 'material 2D, nu = 0: diag(E, E, E/2)')

  ! Zero modes: rigid body motions, plus hourglass modes under reduced integration
  call zeros(CUB_LIN, 1, 2, SHP_LAGRANGE, 1, 'lin2')
  call zeros(CUB_QUA, 1, 2, SHP_LAGRANGE, 3, 'qua4')
  call zeros(CUB_QUA, 1, 1, SHP_LAGRANGE, 5, 'qua4 one point')
  call zeros(CUB_QUA, 2, 4, SHP_SERENDIPITY, 3, 'qua8')
  call zeros(CUB_TRI, 2, 4, SHP_LAGRANGE, 3, 'tri6')
  call zeros(CUB_HEX, 1, 2, SHP_LAGRANGE, 6, 'hex8')
  call zeros(CUB_HEX, 1, 1, SHP_LAGRANGE, 18, 'hex8 one point')
  call zeros(CUB_TET, 1, 2, SHP_LAGRANGE, 6, 'tet4')
  call zeros(CUB_WED, 2, 4, SHP_LAGRANGE, 6, 'wed18')

  ! Eigenvalues with known values ---------------------------------------
  ! quadratic bar, L = 2: K = E/(3L) [7 1 -8; 1 7 -8; -8 -8 16]
  lambda = eig(CUB_LIN, 2, 4, SHP_LAGRANGE, 1._dp, 0.3_dp, 'bar3')
  call check(close(lambda, [0._dp, 1._dp, 4._dp]), 'bar3 eigenvalues 0, 1, 4')
  ! 4-node quad, one point: matches the Finite_Elements/brick2_1_g1 prototype
  lambda = eig(CUB_QUA, 1, 1, SHP_LAGRANGE, 800._dp, 0.2_dp, 'quad4 reduced')
  call check(close(lambda(6:), [2000/3._dp, 2000/3._dp, 10000/9._dp]), 'quad4 reduced eigenvalues')
  ! 8-node hex: dilatation mode x/|x| has energy V*3E/(1-2nu)/24
  lambda = eig(CUB_HEX, 1, 2, SHP_LAGRANGE, 1._dp, 0.3_dp, 'hex8')
  call check(close(lambda(24:), [8*3/0.4_dp/24]), 'hex8 dilatation eigenvalue')

  ! Named modes: every mode of linear line/quad/hex is a pure catalog field --
  call named(CUB_LIN, 2, 'tension x  u=(x)')
  call named(CUB_QUA, 2, 'expansion  u=(x,y)')
  call named(CUB_QUA, 1, 'expansion  u=(x,y)')
  call named(CUB_HEX, 2, 'expansion  u=(x,y,z)')
  call named(CUB_HEX, 1, 'expansion  u=(x,y,z)')

  ! Plot ----------------------------------------------------------------
  call plotted()
  ! gfortran does not free main-program allocatables, and LeakSanitizer would report it
  deallocate (lambda)

  write(*, '(i0,a,i0,a)') ncheck - nfail, ' of ', ncheck, ' checks passed'
  if (nfail > 0) error stop 1

contains

subroutine check(cond, name)
  logical, intent(in) :: cond
  character(*), intent(in) :: name
  ncheck = ncheck + 1
  if (.not. cond) then
    nfail = nfail + 1
    print '(a)', 'FAIL: ' // name
  end if
end subroutine check

pure logical function close(a, b)
  !! a and b agree to tol relative to the larger of 1 and max|b|
  real(dp), intent(in) :: a(:), b(:)
  close = all(abs(a - b) <= tol*max(1._dp, maxval(abs(b))))
end function close

function eig(elm, order, degree, family, E, nu, name) result(lambda)
  !! Eigenvalues of the reference element; checks that K is symmetric
  integer, intent(in) :: elm, order, degree, family
  real(dp), intent(in) :: E, nu
  character(*), intent(in) :: name
  real(dp), allocatable :: lambda(:), K(:,:), phi(:,:)
  type(cubature) :: q
  type(shapefunc) :: s

  q = cubature(elm, degree)
  s = shapefunc(q, order, family=family)
  K = stiffness(s, q, s%coords, material(s%dim, E, nu))
  call check(maxval(abs(K - transpose(K))) <= 2**12*eps*maxval(abs(K)), name // ': K symmetric')
  allocate (lambda(size(K,1)), phi(size(K,1), size(K,1)))
  call modes(K, lambda, phi)
end function eig

subroutine zeros(elm, order, degree, family, expected, name)
  integer, intent(in) :: elm, order, degree, family, expected
  character(*), intent(in) :: name
  real(dp), allocatable :: lambda(:)

  lambda = eig(elm, order, degree, family, 1._dp, 0.3_dp, name)
  call check(count(lambda == 0) == expected, name // ': number of zero modes')
end subroutine zeros

subroutine named(elm, degree, last)
  integer, intent(in) :: elm, degree
  character(*), intent(in) :: last
  real(dp), allocatable :: K(:,:), lambda(:), phi(:,:), match(:)
  character(40), allocatable :: names(:)
  character(40) :: name
  type(cubature) :: q
  type(shapefunc) :: s
  integer :: n

  q = cubature(elm, degree)
  s = shapefunc(q, 1)
  K = stiffness(s, q, s%coords, material(s%dim, 1._dp, 0.3_dp))
  n = size(K,1)
  allocate (lambda(n), phi(n,n), names(n), match(n))
  call modes(K, lambda, phi)
  call name_modes(s%coords, lambda, phi, names, match)
  write(name, '(a,i0,a,i0)') 'named elm ', elm, ' degree ', degree
  call check(.not. any(names == '-') .and. all(match >= 1 - tol), &
           & trim(name) // ': every mode named, 100% match')
  call check(names(n) == last, trim(name) // ': stiffest mode is ' // last)
  call check(maxval(abs(matmul(K, phi) - phi*spread(lambda, 1, n))) <= tol*maxval(lambda), &
           & trim(name) // ': named modes are eigenvectors')
end subroutine named

subroutine plotted()
  character(*), parameter :: file = 'test_feigen.svg'
  real(dp), allocatable :: K(:,:), lambda(:), phi(:,:)
  type(cubature) :: q
  type(shapefunc) :: s
  logical :: exists
  integer :: u

  q = cubature(CUB_QUA, 2)
  s = shapefunc(q, 1)
  K = stiffness(s, q, s%coords, material(2, 1._dp, 0.3_dp))
  allocate (lambda(8), phi(8,8))
  call modes(K, lambda, phi)
  call plot_modes(s, s%coords, lambda, phi, spread('-', 1, 8), 0.2_dp, file)
  inquire (file=file, exist=exists)
  call check(exists, 'plot_modes writes ' // file)
  if (exists) then
    open (newunit=u, file=file)
    close (u, status='delete')
  end if
end subroutine plotted

end program test_feigen
