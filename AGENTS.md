# AGENTS.md

Instructions for coding agents working on this repository. The public API and
its conventions are in `llms.txt`; read it before changing behaviour.

## Layout

| Path | Role |
|---|---|
| `src/feigen.f90` | **Source of truth**: module `feigen` |
| `app/main.f90` | The `feigen` command-line program (`fpm run`) |
| `test/test_feigen.f90` | The test program |
| `example/*.f90` | Demo programs; each is copied verbatim into `llms.txt` |
| `llms.txt` | API reference for agents, and the online user manual |
| `man/feigen.3`, `README.md`, `ford.md` | Human documentation |
| `images/qua4_full.png`, `images/hex8_full.png` | README images: `fpm run -- -e qua` and `-e hex`, renamed so later runs cannot overwrite them |
| `images/logo.png` | README logo and FORD favicon: photo by Tijana Drndarski, Unsplash License (credited in README) |
| `tools/check_docs.py` | CI check: every public name and example is in `llms.txt` |
| `tools/doc.sh` | FORD build, used by the Docs workflow |

## Commands

```sh
pip install ford                     # once; the sources need no preprocessor
fpm test --profile debug             # fast check, gfortran
fpm test --profile release
fpm run --example --all
fpm run -- -e hex -p 2               # the app; arguments after --
python3 tools/check_docs.py
./ci.sh                              # gfortran and ifx, if installed
tools/doc.sh                         # FORD HTML in build/doc
```

Done means: the tests pass for every compiler you have in both profiles, and
`tools/check_docs.py` passes.

## Decided conventions (do not change without the maintainer)

- fpm only, no Makefile: stdlib, fplotlib and FLAP are multi-file fpm packages that a
  Makefile could not bundle in a few lines. `tools/doc.sh` replaces `make doc`.
- Kinds are not exported: everything is `real64` (`dp`) from `iso_fortran_env`.
- Degrees of freedom node by node; Voigt order xx, yy, zz, xy, xz, yz with engineering
  shear strains; 2D is plane strain.
- `modes` sets eigenvalues below `tol*max|lambda|` to exactly zero, so zero modes are
  counted with `== 0`.
- `name_modes` replaces eigenvectors within a group of equal eigenvalues by projections
  of catalog fields: exact eigenvectors, not orthogonal within the group.
- Errors: `feigen` does not check its arguments. Invalid elements and orders stop in
  fortran_cubatures (`cubature%set:`) and fortran_shapefuncs (`shapefunc%set:`), a singular
  Jacobian stops in stdlib's `inv`, and the app stops with messages that start with `feigen:`.
  Errors cannot be caught.
- The default plot file is `<element>_p<order>_<fam>_d<degree>.png`, `fam` the first three
  letters of the family; plot files in the root are gitignored.
- The README images are committed generated output, an agreed exception. `images/logo.png` is a
  256x256 crop of an Unsplash photo; keep the credit in the README License section.

## Adding or changing a procedure

1. Implement it in `src/feigen.f90`, with a `!!` line after its declaration.
2. Add it to the `public ::` statement.
3. Add checks to `test/test_feigen.f90`.
4. Document it in `llms.txt` (the CI check requires the name), the man page and `CHANGELOG.md`.
5. Run the tests in both profiles and `python3 tools/check_docs.py`.

## Test rules

- Debug profiles trap floating-point exceptions (ifx `-fpe0`, nvfortran `-Ktrap=fp`):
  tests must not create NaN or Inf through arithmetic.
- Tolerances are relative to `epsilon`: `2**19*eps` for eigenvalues and eigenvectors,
  which come from an iterative solver, and `2**12*eps` for the symmetry of K.
- The plot check writes `test_feigen.svg` in the working directory and deletes it.

## Compiler quirks that shaped the code

- LLVM flang 21 is not in CI: it cannot resolve stdlib's kind parameter `ilp` in an
  `equivalence` statement of the submodule `stdlib_lapack_eigv_gen2`.
- ifx on Windows is not in CI: setup-fortran installs no `icx` there, which fpm needs for
  stdlib's C sources.
- nvfortran is not in CI: 26.3 rejects stdlib's `implicit none(type, external)` and links
  fplotlib's module `fplotlib_contour` and procedure `contour` to the same symbol.
- fpm compiles stdlib's C sources with `icx` for ifx; `ci.sh` sets `FPM_CC=gcc` when `icx`
  is missing.
- ifx 2026.1.1 miscompiles an internal function with an automatic result sized by a host
  variable when a sibling internal procedure calls it: `edge_path` in `plot_modes` has an
  allocatable result instead.
- gfortran does not deallocate main-program allocatables at the end, and the debug
  profile's LeakSanitizer reports them: the app and the test deallocate them explicitly.

## Commits

Commit as the maintainer, Will Klausler, with the GitHub no-reply address
`226649647+willklausler@users.noreply.github.com` (set in the repository's git config).
