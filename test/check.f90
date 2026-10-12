program check
  use cubatures, only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs, only: shapefunc, SHP_LAGRANGE, SHP_SERENDIPITY
  use feigen, only: material, stiffness, modes, name_modes
  implicit none

  real(rk), allocatable :: lambda(:)

  ! zero modes: rigid body motions, plus hourglass modes under reduced integration
  call zeros(CUB_LIN, 1, 2, SHP_LAGRANGE, 1)
  call zeros(CUB_QUA, 1, 2, SHP_LAGRANGE, 3)
  call zeros(CUB_QUA, 1, 1, SHP_LAGRANGE, 5)
  call zeros(CUB_QUA, 2, 4, SHP_SERENDIPITY, 3)
  call zeros(CUB_TRI, 2, 4, SHP_LAGRANGE, 3)
  call zeros(CUB_HEX, 1, 2, SHP_LAGRANGE, 6)
  call zeros(CUB_HEX, 1, 1, SHP_LAGRANGE, 18)
  call zeros(CUB_TET, 1, 2, SHP_LAGRANGE, 6)
  call zeros(CUB_WED, 2, 4, SHP_LAGRANGE, 6)

  ! quadratic bar, L = 2: K = E/(3L) [7 1 -8; 1 7 -8; -8 -8 16]
  lambda = eig(CUB_LIN, 2, 4, SHP_LAGRANGE, 1.0_rk, 0.3_rk)
  call close(lambda, [0.0_rk, 1.0_rk, 4.0_rk], 'bar3')

  ! 4-node quad, one point: matches Finite_Elements/brick2_1_g1 prototype
  lambda = eig(CUB_QUA, 1, 1, SHP_LAGRANGE, 800.0_rk, 0.2_rk)
  call close(lambda(6:), [2000/3.0_rk, 2000/3.0_rk, 10000/9.0_rk], 'quad4 reduced')

  ! 8-node hex: dilatation mode x/|x| has energy V*3E/(1-2nu)/24
  lambda = eig(CUB_HEX, 1, 2, SHP_LAGRANGE, 1.0_rk, 0.3_rk)
  call close(lambda(24:), [8*3/0.4_rk/24], 'hex8 dilatation')

  ! every mode of linear line/quad/hex is a pure catalog field, and stays an eigenvector
  call named(CUB_LIN, 2, 'tension x  u=(x)')
  call named(CUB_QUA, 2, 'expansion  u=(x,y)')
  call named(CUB_QUA, 1, 'expansion  u=(x,y)')
  call named(CUB_HEX, 2, 'expansion  u=(x,y,z)')
  call named(CUB_HEX, 1, 'expansion  u=(x,y,z)')

  print '(a)', 'all checks passed'

contains

  function eig(elm, order, degree, family, E, nu) result(lambda)
    integer, intent(in) :: elm, order, degree, family
    real(rk), intent(in) :: E, nu
    real(rk), allocatable :: lambda(:), K(:,:), phi(:,:)
    type(cubature) :: q
    type(shapefunc) :: s

    q = cubature(elm, degree)
    s = shapefunc(q, order, family=family)
    K = stiffness(s, q, s%coords, material(s%dim, E, nu))
    if (maxval(abs(K - transpose(K))) > 1e-12_rk*maxval(abs(K))) error stop 'K not symmetric'
    allocate (lambda(size(K,1)), phi(size(K,1), size(K,1)))
    call modes(K, lambda, phi)
  end function eig

  subroutine named(elm, degree, last)
    integer, intent(in) :: elm, degree
    character(*), intent(in) :: last
    real(rk), allocatable :: K(:,:), lambda(:), phi(:,:), match(:)
    character(40), allocatable :: names(:)
    type(cubature) :: q
    type(shapefunc) :: s
    integer :: n

    q = cubature(elm, degree)
    s = shapefunc(q, 1)
    K = stiffness(s, q, s%coords, material(s%dim, 1.0_rk, 0.3_rk))
    n = size(K,1)
    allocate (lambda(n), phi(n,n), names(n), match(n))
    call modes(K, lambda, phi)
    call name_modes(s%coords, lambda, phi, names, match)
    if (any(names == '-') .or. any(match < 1 - 1e-10_rk)) error stop 'unnamed or impure mode'
    if (names(n) /= last) error stop 'wrong name for the stiffest mode'
    if (maxval(abs(matmul(K, phi) - phi*spread(lambda, 1, n))) > 1e-10_rk*maxval(lambda)) &
      error stop 'named mode is not an eigenvector'
  end subroutine named

  subroutine zeros(elm, order, degree, family, expected)
    integer, intent(in) :: elm, order, degree, family, expected
    real(rk), allocatable :: lambda(:)
    lambda = eig(elm, order, degree, family, 1.0_rk, 0.3_rk)
    if (count(lambda == 0) /= expected) then
      print '(a,i0,a,i0,a,i0)', 'element ', elm, ': ', count(lambda == 0), ' zero modes, expected ', expected
      error stop
    end if
  end subroutine zeros

  subroutine close(a, b, name)
    real(rk), intent(in) :: a(:), b(:)
    character(*), intent(in) :: name
    if (any(abs(a - b) > 1e-10_rk*max(1.0_rk, maxval(abs(b))))) then
      print *, name, ': ', a, ' expected ', b
      error stop
    end if
  end subroutine close

end program check
