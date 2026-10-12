program demo_hourglass
  !! Reduced integration of a 4-node quadrilateral: one cubature point adds two
  !! zero-energy hourglass modes to the three rigid-body modes
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use cubatures
  use shapefuncs
  use feigen
  implicit none
  type(cubature) :: q
  type(shapefunc) :: s
  real(dp) :: K(8,8), lambda(8), phi(8,8), match(8)
  character(40) :: names(8)
  integer :: degree, nzero, i

  do degree = 2, 1, -1                          ! 2x2 points, then 1 point
    q = cubature(CUB_QUA, degree)
    s = shapefunc(q, 1)
    K = stiffness(s, q, s%coords, material(s%dim, 1._dp, 0.3_dp))
    call modes(K, lambda, phi)
    call name_modes(s%coords, lambda, phi, names, match)
    nzero = count(lambda == 0)                  ! 3, then 5
    write(*, '(a,i0,a,i0)') 'points: ', q%npoints, ', zero modes: ', nzero
    write(*, '(4x,a)') (trim(names(i)), i = 1, nzero)
  end do
  ! translation x, translation y, rotation z;
  ! with one point also bending x u=(xy,0) and bending y u=(0,xy): the hourglass modes
end program demo_hourglass
