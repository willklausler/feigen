# Changelog

All notable changes are listed here. The project follows [Semantic Versioning](https://semver.org).

## [0.1.0] - 2026-10-12

First public release.

- Module `feigen`: `material` (isotropic elasticity, plane strain in 2D), `stiffness`, `modes` (eigenvalues with exact zeros), `name_modes` (named rigid-body, deformation and hourglass fields) and `plot_modes` (one PNG, SVG or PDF with a subplot per mode).
- App `feigen`: command-line front end for every element of fortran_shapefuncs; the default plot file is named after the element, `<element>_p<order>_<fam>.png`.
- Depends on fortran_cubatures v0.2.0, fortran_shapefuncs v0.2.0, stdlib, fplotlib and FLAP.
- FORD documentation, man page, `llms.txt` API reference for coding agents.
