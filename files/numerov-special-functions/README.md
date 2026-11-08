# Numerov method & special functions of mathematical physics

Course notes, codes and exercises (Python 3 and Fortran) for solving the
Schrödinger equation with the **Numerov method**, and for computing the
**special functions** that appear when it is separated in three dimensions.

| Volume | Contents | Document |
|---|---|---|
| **1. Numerov method** (`vol1_numerov/`) | free particle, square barrier, harmonic oscillator; normalisation, error analysis, convergence; 10 exercises with worked solutions | `Numerov_Schrodinger_Course.docx`, `Numerov_Exercise_Solutions.docx` |
| **2. Special functions** (`vol2_special_functions/`) | Laguerre → hydrogen atom, Legendre → angular equation, Bessel $J_m$ → circular well, spherical Bessel → spherical well; stable recurrences, Miller's algorithm, Gauss quadrature, the Frobenius start | `Special_Functions_Numerov_Course.docx` |

The built Word documents are attached to each
[GitHub Release](../../releases); the Markdown sources are in `*/docs/`.

## Quick start

```bash
git clone https://github.com/<your-user>/numerov-special-functions.git
cd numerov-special-functions
python3 -m pip install -r requirements.txt      # numpy scipy matplotlib mpmath sympy pytest
make test                                        # ~1 minute of regression tests
make fortran                                     # needs gfortran (Fortran 2008)
make python                                      # regenerate every figure and table (slow)
make docs                                        # rebuild the .docx files into dist/ (needs pandoc)
```

Run a single program:

```bash
cd vol2_special_functions/python && python3 ch1_laguerre_hydrogen.py
cd vol2_special_functions/fortran && gfortran -O2 -o ch1 ch1_laguerre_hydrogen.f90 && ./ch1
```

## Repository layout

```
vol1_numerov/
  python/      p1_free_particle.py  p2_square_barrier.py  p3_harmonic_oscillator.py (+ outputs, figures)
  fortran/     p1 … p3 *.f90, plot_fortran_results.py (+ .dat files)
  exercises/   ex01_05.py  ex06_09.py  ex10_precision.F90  ex10_plot.py
  figures/     illustrations and logic charts (Graphviz .dot sources)
  docs/        course/ and solutions/: Markdown sources + build scripts
vol2_special_functions/
  python/      ch1_laguerre_hydrogen.py … ch4_spherical_bessel_well.py
  fortran/     ch1 … ch4 *.f90, plot_fortran_results.py
  figures/     illustrations and logic charts
  docs/        part0.md … part5.md, refs.py (bibliography), cite.py (numbered citations), build.py, style.py
tests/         fast pytest checks against exact results
```

## What each program verifies

Every program prints its own verification tables: comparison with exact
eigenvalues and eigenfunctions, normalisation, orthogonality, expectation values
(virial theorem, ⟨r⟩, ⟨r²⟩, ⟨cos²θ⟩), independent identities (Jacobi–Anger, cross
products, quadrature sums), and convergence studies with the observed order of
accuracy, Richardson extrapolation and the round-off floor.

## Citing

See [`CITATION.cff`](CITATION.cff) — GitHub shows a "Cite this repository" button.
The literature used is listed with numbered citations in each document; the
bibliographies live in `vol1_numerov/docs/course/refs_numerov.py` and
`vol2_special_functions/docs/refs.py`.

## Licence

Code: MIT (see `LICENSE`). Course texts: CC BY 4.0.
