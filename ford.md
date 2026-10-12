---
project: feigen
summary: Stiffness eigenmodes of a single finite element, named and plotted.
author: Will Klausler
github: https://github.com/willklausler
project_github: https://github.com/willklausler/feigen
license: mit
src_dir: ./src
output_dir: ./build/doc
page_dir: ./build/pages
preprocess: false
display: public
source: true
search: true
graph: false
---

`use feigen` gives `material`, `stiffness`, `modes`, `name_modes` and
`plot_modes`. The element comes from fortran_cubatures (`cubature`) and
fortran_shapefuncs (`shapefunc`). The `feigen` program wraps them in a
command line. Start with the [user manual](page/index.html), which is also
available as plain text in [llms.txt](llms.txt).

```fortran
q = cubature(CUB_QUA, 1)                    ! one point: reduced integration
s = shapefunc(q, 1)                         ! 4-node quadrilateral
K = stiffness(s, q, s%coords, material(2, 1._dp, 0.3_dp))
call modes(K, lambda, phi)
print *, count(lambda == 0)                 ! 5: 3 rigid-body + 2 hourglass
```
