!> Stiffness eigenmodes of a single finite element.
module feigen
  use cubatures, only: cubature, rk, CUB_LIN, CUB_TRI, CUB_QUA, CUB_TET, CUB_HEX, CUB_WED
  use shapefuncs, only: shapefunc
  use stdlib_linalg, only: det, inv, eigh
  implicit none
  private

  public :: material, stiffness, modes, name_modes, plot_modes

contains

  !> Isotropic elasticity matrix in Voigt notation, engineering shear strains.
  !> 1D: bar, 2D: plane strain, 3D: order xx, yy, zz, xy, xz, yz.
  pure function material(dim, E, nu) result(D)
    integer, intent(in) :: dim
    real(rk), intent(in) :: E, nu
    real(rk) :: D(dim*(dim+1)/2, dim*(dim+1)/2)
    integer :: i

    D = 0
    if (dim == 1) then
      D(1,1) = E
      return
    end if
    D(1:dim, 1:dim) = nu
    do i = 1, dim
      D(i,i) = 1 - nu
    end do
    do i = dim+1, size(D,1)
      D(i,i) = (1 - 2*nu)/2
    end do
    D = E/((1 + nu)*(1 - 2*nu))*D
  end function material

  !> Stiffness matrix K = sum_g B^T D B det(J) w_g, dofs ordered node by node.
  function stiffness(s, q, X, D) result(K)
    type(shapefunc), intent(in) :: s
    type(cubature), intent(in) :: q
    real(rk), intent(in) :: X(:,:)   !< nodal coordinates (dim, nnodes)
    real(rk), intent(in) :: D(:,:)
    real(rk) :: K(s%dim*s%nnodes, s%dim*s%nnodes)
    real(rk) :: jac(s%dim, s%dim), dN(s%dim, s%nnodes), B(size(D,1), s%dim*s%nnodes)
    integer :: g, n, i

    K = 0
    do g = 1, q%npoints
      jac = matmul(X, transpose(s%derv(:,:,g)))        ! dx/dxi
      dN = matmul(transpose(inv(jac)), s%derv(:,:,g))  ! dN/dx
      B = 0
      do n = 1, s%nnodes
        i = s%dim*(n - 1)
        select case (s%dim)
        case (1)
          B(1, i+1) = dN(1,n)
        case (2)
          B(1, i+1) = dN(1,n)
          B(2, i+2) = dN(2,n)
          B(3, i+1) = dN(2,n); B(3, i+2) = dN(1,n)
        case (3)
          B(1, i+1) = dN(1,n)
          B(2, i+2) = dN(2,n)
          B(3, i+3) = dN(3,n)
          B(4, i+1) = dN(2,n); B(4, i+2) = dN(1,n)
          B(5, i+1) = dN(3,n); B(5, i+3) = dN(1,n)
          B(6, i+2) = dN(3,n); B(6, i+3) = dN(2,n)
        end select
      end do
      K = K + matmul(transpose(B), matmul(D, B))*det(jac)*q%weights(g)
    end do
  end function stiffness

  !> Eigenvalues (ascending) and eigenvectors of symmetric K.
  !> Eigenvalues below tol*max|lambda| are set to exactly zero.
  subroutine modes(K, lambda, phi, tol)
    real(rk), intent(in) :: K(:,:)
    real(rk), intent(out) :: lambda(:), phi(:,:)
    real(rk), intent(in), optional :: tol
    real(rk) :: A(size(K,1), size(K,2)), t

    t = 1e-10_rk
    A = K   ! eigh overwrites its input
    if (present(tol)) t = tol
    call eigh(A, lambda, vectors=phi)
    where (abs(lambda) < t*maxval(abs(lambda))) lambda = 0
  end subroutine modes

  !> Replace the eigenvectors of each group of equal eigenvalues by the projections
  !> of the best-matching catalog fields onto that eigenspace, signed to agree with
  !> the field. These are exact eigenvectors, but not orthogonal within a group.
  !> match is the cosine between field and eigenspace; unmatched modes keep an
  !> orthonormal remainder of the eigenspace, named "-".
  subroutine name_modes(X, lambda, phi, names, match, tol)
    real(rk), intent(in) :: X(:,:), lambda(:)
    real(rk), intent(inout) :: phi(:,:)
    character(len=*), intent(out) :: names(:)
    real(rk), intent(out) :: match(:)
    real(rk), intent(in), optional :: tol
    real(rk), parameter :: min_match = 0.5_rk, min_independent = 0.1_rk
    real(rk), allocatable :: F(:,:), V(:,:), P(:,:), Q(:,:), frac(:), c(:)
    character(len=40), allocatable :: tnames(:)
    logical, allocatable :: avail(:)
    real(rk) :: t
    integer :: i0, i1, k, n, j

    t = 1e-8_rk
    if (present(tol)) t = tol
    call catalog(X, F, tnames)
    names = '-'
    match = 0
    i0 = 1
    do while (i0 <= size(lambda))
      i1 = i0
      do while (i1 < size(lambda))
        if (abs(lambda(i1+1) - lambda(i0)) > t*maxval(abs(lambda))) exit
        i1 = i1 + 1
      end do
      k = i1 - i0 + 1
      V = phi(:, i0:i1)
      P = matmul(transpose(V), F)   ! fields in eigenspace coordinates
      frac = norm2(P, 1)
      avail = frac > min_match
      allocate (Q(k, 0))            ! orthonormal basis of the accepted fields
      n = 0
      ! greedy: best match first, catalog order breaks ties
      do while (n < k .and. any(avail))
        j = findloc(avail .and. frac >= maxval(frac, mask=avail) - 1e-8_rk, .true., dim=1)
        avail(j) = .false.
        c = P(:,j)/frac(j)
        c = c - matmul(Q, matmul(transpose(Q), c))
        if (norm2(c) < min_independent) cycle   ! dependent on accepted fields
        Q = reshape([Q, c/norm2(c)], [k, n+1])
        n = n + 1
        phi(:, i0+n-1) = matmul(V, P(:,j))/frac(j)
        names(i0+n-1) = tnames(j)
        match(i0+n-1) = frac(j)
      end do
      ! fill the rest with the orthonormal complement of the accepted fields
      do j = 1, k
        if (n == k) exit
        c = spread(0.0_rk, 1, k)
        c(j) = 1
        c = c - matmul(Q, matmul(transpose(Q), c))
        if (norm2(c) < min_independent) cycle
        Q = reshape([Q, c/norm2(c)], [k, n+1])
        n = n + 1
        phi(:, i0+n-1) = matmul(V, Q(:,n))
      end do
      deallocate (Q)
      i0 = i1 + 1
    end do
  end subroutine name_modes

  !> Named displacement fields of linear line, quad and hex elements, in coordinates
  !> centred on the element, as unit columns of F (dofs node by node).
  !> Bending is named by the fibres that bend: (xy,0) bends x-fibres through y.
  subroutine catalog(X, F, names)
    real(rk), intent(in) :: X(:,:)
    real(rk), allocatable, intent(out) :: F(:,:)
    character(len=40), allocatable, intent(out) :: names(:)
    real(rk), dimension(size(X,2)) :: cx, cy, cz, o

    cx = X(1,:) - sum(X(1,:))/size(X,2)
    o = 1
    allocate (F(size(X), 0), names(0))
    select case (size(X,1))
    case (1)
      call add('translation x  u=(1)', o)
      call add('tension x  u=(x)', cx)
    case (2)
      cy = X(2,:) - sum(X(2,:))/size(X,2)
      call add('translation x  u=(1,0)', o, 0*o)
      call add('translation y  u=(0,1)', 0*o, o)
      call add('rotation z  u=(-y,x)', -cy, cx)
      call add('expansion  u=(x,y)', cx, cy)
      call add('tension x/compression y  u=(x,-y)', cx, -cy)
      call add('shear xy  u=(y,x)', cy, cx)
      call add('bending x  u=(xy,0)', cx*cy, 0*o)
      call add('bending y  u=(0,xy)', 0*o, cx*cy)
    case (3)
      cy = X(2,:) - sum(X(2,:))/size(X,2)
      cz = X(3,:) - sum(X(3,:))/size(X,2)
      call add('translation x  u=(1,0,0)', o, 0*o, 0*o)
      call add('translation y  u=(0,1,0)', 0*o, o, 0*o)
      call add('translation z  u=(0,0,1)', 0*o, 0*o, o)
      call add('rotation x  u=(0,-z,y)', 0*o, -cz, cy)
      call add('rotation y  u=(z,0,-x)', cz, 0*o, -cx)
      call add('rotation z  u=(-y,x,0)', -cy, cx, 0*o)
      call add('expansion  u=(x,y,z)', cx, cy, cz)
      call add('tension x/compression y  u=(x,-y,0)', cx, -cy, 0*o)
      call add('tension y/compression z  u=(0,y,-z)', 0*o, cy, -cz)
      call add('tension z/compression x  u=(-x,0,z)', -cx, 0*o, cz)
      call add('shear xy  u=(y,x,0)', cy, cx, 0*o)
      call add('shear xz  u=(z,0,x)', cz, 0*o, cx)
      call add('shear yz  u=(0,z,y)', 0*o, cz, cy)
      call add('bending yz, synclastic  u=x(0,y,z)', 0*o, cx*cy, cx*cz)
      call add('bending xz, synclastic  u=y(x,0,z)', cx*cy, 0*o, cy*cz)
      call add('bending xy, synclastic  u=z(x,y,0)', cx*cz, cy*cz, 0*o)
      call add('bending yz, anticlastic  u=x(0,y,-z)', 0*o, cx*cy, -cx*cz)
      call add('bending xz, anticlastic  u=y(-x,0,z)', -cx*cy, 0*o, cy*cz)
      call add('bending xy, anticlastic  u=z(x,-y,0)', cx*cz, -cy*cz, 0*o)
      call add('torsion x  u=x(0,-z,y)', 0*o, -cx*cz, cx*cy)
      call add('torsion y  u=y(z,0,-x)', cy*cz, 0*o, -cx*cy)
      call add('torsion z  u=z(-y,x,0)', -cy*cz, cx*cz, 0*o)
      call add('warping  u=(yz,xz,xy)', cy*cz, cx*cz, cx*cy)
      call add('hourglass x  u=(xyz,0,0)', cx*cy*cz, 0*o, 0*o)
      call add('hourglass y  u=(0,xyz,0)', 0*o, cx*cy*cz, 0*o)
      call add('hourglass z  u=(0,0,xyz)', 0*o, 0*o, cx*cy*cz)
    end select

  contains

    subroutine add(name, ux, uy, uz)
      character(len=*), intent(in) :: name
      real(rk), intent(in) :: ux(:)
      real(rk), intent(in), optional :: uy(:), uz(:)
      real(rk) :: u(size(X,1), size(X,2))
      u(1,:) = ux
      if (present(uy)) u(2,:) = uy
      if (present(uz)) u(3,:) = uz
      F = reshape([F, reshape(u, [size(u)])/norm2(u)], [size(u), size(F,2)+1])
      names = [character(len=40) :: names, name]
    end subroutine add

  end subroutine catalog

  !> Vertex pairs of the element edges, classical numbering.
  pure function edge_vertices(elm) result(ev)
    integer, intent(in) :: elm
    integer, allocatable :: ev(:,:)

    select case (elm)
    case (CUB_LIN)
      ev = reshape([1,2], [2,1])
    case (CUB_TRI)
      ev = reshape([1,2, 2,3, 3,1], [2,3])
    case (CUB_QUA)
      ev = reshape([1,2, 2,3, 3,4, 4,1], [2,4])
    case (CUB_TET)
      ev = reshape([1,2, 2,3, 3,1, 1,4, 2,4, 3,4], [2,6])
    case (CUB_HEX)
      ev = reshape([1,2, 2,3, 3,4, 4,1, 5,6, 6,7, 7,8, 8,5, 1,5, 2,6, 3,7, 4,8], [2,12])
    case (CUB_WED)
      ev = reshape([1,2, 2,3, 3,1, 4,5, 5,6, 6,4, 1,4, 2,5, 3,6], [2,9])
    end select
  end function edge_vertices

  !> One PNG/SVG with a subplot per mode: undeformed (black), deformed (green),
  !> nodal displacements (blue). Displacements are scaled to scale*element size.
  subroutine plot_modes(s, X, lambda, phi, names, scale, filename)
    use fplotlib, only: figure, subplot, plot, plot3d, scatter, scatter3d, title, &
                        set_aspect, axis, savefig, dp
    type(shapefunc), intent(in) :: s
    real(rk), intent(in) :: X(:,:), lambda(:), phi(:,:), scale
    character(len=*), intent(in) :: names(:), filename
    real(rk) :: U(size(X,1), size(X,2)), h
    integer, allocatable :: ev(:,:), path(:)
    integer :: ndof, rows, cols, m, e, i, nv, ne
    character(len=120) :: label

    ndof = size(lambda)
    cols = ceiling(sqrt(real(ndof)))
    rows = (ndof + cols - 1)/cols
    h = maxval(maxval(X, 2) - minval(X, 2))
    ev = edge_vertices(s%elm)
    nv = maxval(ev)
    ne = max(s%order - 1, 0)   ! interior nodes per edge
    call figure(figsize=[4.5_dp*cols, 3.5_dp*rows])

    do m = 1, ndof
      call subplot(rows, cols, m)
      U = X + scale*h*reshape(phi(:,m), shape(X))/maxval(abs(phi(:,m)))
      write (label, '("$\lambda_{",i0,"}$=",es9.2,2x,a)') m, lambda(m), trim(names(m))
      call title(trim(label), fontsize=8.0_dp)

      ! 1D: undeformed at y = 0, deformed at y = h/2
      if (s%dim == 1) then
        call draw2d(reshape([X, 0*X], [2, size(X,2)], order=[2,1]), &
                    reshape([U, 0*U + h/2], [2, size(U,2)], order=[2,1]))
      else if (s%dim == 2) then
        call draw2d(X, U)
      else
        do e = 1, size(ev, 2)
          path = edge_path(e)
          call plot3d(X(1,path), X(2,path), X(3,path), "k-")
          call plot3d(U(1,path), U(2,path), U(3,path), "g-")
        end do
        do i = 1, size(X, 2)
          call plot3d([X(1,i), U(1,i)], [X(2,i), U(2,i)], [X(3,i), U(3,i)], "b-", lw=0.5_dp)
        end do
        call scatter3d(X(1,:), X(2,:), X(3,:), c="k", s=10.0_dp)
        call scatter3d(U(1,:), U(2,:), U(3,:), c="k", s=10.0_dp)
      end if
    end do
    do m = ndof + 1, rows*cols
      call subplot(rows, cols, m)
      call axis("off")
    end do
    call savefig(filename)

  contains

    !> Nodes along edge e: first vertex, edge nodes, second vertex.
    function edge_path(e) result(p)
      integer, intent(in) :: e
      integer :: p(ne + 2), k
      p(1) = ev(1,e)
      p(ne+2) = ev(2,e)
      p(2:ne+1) = [(nv + (e - 1)*ne + k, k = 1, ne)]
    end function edge_path

    subroutine draw2d(A, C)
      real(rk), intent(in) :: A(:,:), C(:,:)
      integer :: e, i
      do e = 1, size(ev, 2)
        path = edge_path(e)
        call plot(A(1,path), A(2,path), "k-")
        call plot(C(1,path), C(2,path), "g-")
      end do
      do i = 1, size(A, 2)
        call plot([A(1,i), C(1,i)], [A(2,i), C(2,i)], "b-", lw=0.5_dp)
      end do
      call scatter(A(1,:), A(2,:), c="k", s=10.0_dp)
      call scatter(C(1,:), C(2,:), c="k", s=10.0_dp)
      call set_aspect(1.0_dp)
    end subroutine draw2d

  end subroutine plot_modes

end module feigen
