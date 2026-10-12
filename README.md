# feigen

Stiffness eigenmodes of a single finite element: solves K·φ = λφ and plots every mode.
Rigid-body and hourglass modes show up as zero eigenvalues.

```sh
fpm run -- --element hex --order 1 --degree 1 --E 200 --nu 0.3 -o hex8.png
fpm run -- -e qua -p 2 -f serendipity -x coords.txt    # distorted 8-node quad
fpm run -- --help
fpm test
```

- Elements: `lin tri qua tet hex wed`, any Lagrange order, serendipity quad/hex orders 1–3
  ([fortran_shapefuncs](https://github.com/willklausler/fortran_shapefuncs)).
- Integration: `--degree` is the degree of exactness, default `2*order`
  ([fortran_cubatures](https://github.com/willklausler/fortran_cubatures)).
- Material: isotropic linear elastic; 2D is plane strain.
- Geometry: reference element by default. `--coords` reads one node per line in the
  library's [node numbering](https://github.com/willklausler/fortran_shapefuncs/blob/main/manual/elements.md).
- Mode names: within each group of equal eigenvalues, modes are replaced by the projections
  of the best-matching named fields (translation, rotation, expansion, tension/compression,
  shear, bending, torsion, warping, hourglass), with a match %. These are exact eigenvectors
  but not orthogonal within a group. The catalog covers linear line, quad and hex elements;
  other modes are named `-`. 3D bending is named by the plane of the fibres that bend,
  synclastic `z(x,y,0)` or anticlastic `z(x,-y,0)`.
- Output: eigenvalues and names on stdout, one PNG/SVG/PDF with a subplot per mode
  ([fplotlib](https://github.com/certik/fplotlib)).
