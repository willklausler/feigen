program main
  use flap, only: command_line_interface
  use cubatures, only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs, only: shapefunc, SHP_LAGRANGE, SHP_SERENDIPITY
  use feigen, only: material, stiffness, modes, name_modes, plot_modes
  implicit none

  type(command_line_interface) :: cli
  type(cubature) :: q
  type(shapefunc) :: s
  character(99) :: element, family, coords, output
  integer :: order, degree, elm, fam, error, u, i
  real(rk) :: E, nu, scale
  real(rk), allocatable :: X(:,:), K(:,:), lambda(:), phi(:,:), match(:)
  character(40), allocatable :: names(:)

  call cli%init(progname='feigen', description='Stiffness eigenmodes of a single finite element')
  call cli%add(switch='--element', switch_ab='-e', help='Element type', def='qua', &
               choices='lin,tri,qua,tet,hex,wed', act='store', error=error)
  call cli%add(switch='--order', switch_ab='-p', help='Shape function order', def='1', act='store', error=error)
  call cli%add(switch='--family', switch_ab='-f', help='Shape function family', def='lagrange', &
               choices='lagrange,serendipity', act='store', error=error)
  call cli%add(switch='--degree', switch_ab='-d', help='Cubature degree of exactness, default 2*order', &
               def='-1', act='store', error=error)
  call cli%add(switch='--E', help="Young's modulus", def='1.0', act='store', error=error)
  call cli%add(switch='--nu', help="Poisson's ratio", def='0.3', act='store', error=error)
  call cli%add(switch='--coords', switch_ab='-x', &
               help='Nodal coordinates file, one node per line, default reference element', &
               def='', act='store', error=error)
  call cli%add(switch='--scale', switch_ab='-s', help='Mode amplitude relative to element size', &
               def='0.2', act='store', error=error)
  call cli%add(switch='--output', switch_ab='-o', help='Plot file (.png, .svg, .pdf)', &
               def='modes.png', act='store', error=error)
  call cli%parse(error=error)
  if (error /= 0) stop 1, quiet=.true.
  call cli%get(switch='--element', val=element, error=error)
  call cli%get(switch='--order', val=order, error=error)
  call cli%get(switch='--family', val=family, error=error)
  call cli%get(switch='--degree', val=degree, error=error)
  call cli%get(switch='--E', val=E, error=error)
  call cli%get(switch='--nu', val=nu, error=error)
  call cli%get(switch='--coords', val=coords, error=error)
  call cli%get(switch='--scale', val=scale, error=error)
  call cli%get(switch='--output', val=output, error=error)

  select case (element)
  case ('lin'); elm = CUB_LIN
  case ('tri'); elm = CUB_TRI
  case ('qua'); elm = CUB_QUA
  case ('tet'); elm = CUB_TET
  case ('hex'); elm = CUB_HEX
  case ('wed'); elm = CUB_WED
  end select
  fam = merge(SHP_SERENDIPITY, SHP_LAGRANGE, family == 'serendipity')
  if (order < 1) error stop 'feigen: order must be at least 1'
  if (degree < 0) degree = 2*order

  q = cubature(elm, degree)
  s = shapefunc(q, order, family=fam)
  X = s%coords
  if (coords /= '') then
    open (newunit=u, file=coords, status='old', action='read')
    read (u, *) X
    close (u)
  end if

  K = stiffness(s, q, X, material(s%dim, E, nu))
  allocate (lambda(size(K,1)), phi(size(K,1), size(K,1)))
  call modes(K, lambda, phi)
  allocate (names(size(lambda)), match(size(lambda)))
  call name_modes(X, lambda, phi, names, match)

  call s%summary()
  print '(a,i0,a,i0)', 'cubature points: ', q%npoints, ', zero modes: ', count(lambda == 0)
  print '(a)', ' mode    eigenvalue  match  name'
  print '(i5,es14.5,f6.0,"%  ",a)', (i, lambda(i), 100*match(i), trim(names(i)), i = 1, size(lambda))

  call plot_modes(s, X, lambda, phi, names, scale, trim(output))
  print '(a)', 'wrote '//trim(output)
end program main
