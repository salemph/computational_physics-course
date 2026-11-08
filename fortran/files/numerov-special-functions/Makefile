# Numerov method & special functions — build everything from source.
#   make test      quick pytest checks (≈1 min)
#   make fortran   compile and run all Fortran programs (writes .dat / output files)
#   make python    run all Python programs (regenerates figures and outputs; slow)
#   make figures   regenerate illustrations and logic charts (needs graphviz `dot`)
#   make docs      rebuild the three Word documents into dist/ (needs pandoc)
#   make clean
PY  ?= python3
FC  ?= gfortran
FFLAGS ?= -O2 -Wall
V1 = vol1_numerov
V2 = vol2_special_functions

.PHONY: all test fortran python figures docs clean
all: fortran python docs

test:
	$(PY) -m pytest -q tests

fortran:
	cd $(V1)/fortran && for p in p1_free_particle p2_square_barrier p3_harmonic_oscillator; do \
	  $(FC) $(FFLAGS) -o $${p%%_*} $$p.f90 && ./$${p%%_*} > $${p%%_*}_output.txt; done && $(PY) plot_fortran_results.py
	cd $(V1)/exercises && $(FC) -O2 -cpp -DWP=4 -o ex10_sp ex10_precision.F90 && ./ex10_sp > ex10_sp.txt \
	  && $(FC) -O2 -cpp -DWP=8 -o ex10_dp ex10_precision.F90 && ./ex10_dp > ex10_dp.txt && $(PY) ex10_plot.py
	cd $(V2)/fortran && for c in 1 2 3 4; do f=$$(ls ch$$c*.f90); \
	  $(FC) $(FFLAGS) -o ch$$c $$f && ./ch$$c > ch$${c}_output.txt; done && $(PY) plot_fortran_results.py

python:
	cd $(V1)/python && for f in p1_free_particle p2_square_barrier p3_harmonic_oscillator; do $(PY) $$f.py > $${f%%_*}_output.txt; done
	cd $(V1)/exercises && $(PY) ex01_05.py > ex01_05_output.txt && $(PY) ex06_09.py > ex06_09_output.txt
	cd $(V2)/python && for f in ch1_laguerre_hydrogen ch2_legendre_rotor ch3_bessel_circular ch4_spherical_bessel_well; do \
	  $(PY) $$f.py > $${f%%_*}_output.txt; done

figures:
	cd $(V1)/figures && $(PY) make_illustrations.py && for d in *.dot; do dot -Tpng -Gdpi=170 $$d -o $${d%.dot}.png; done
	cd $(V2)/figures && $(PY) make_illustrations.py && for d in *.dot; do dot -Tpng -Gdpi=170 $$d -o $${d%.dot}.png; done

docs:
	mkdir -p dist
	cd $(V1)/docs/course && $(PY) build.py && $(PY) style.py && cp Numerov_Schrodinger_Course.docx ../../../dist/
	cd $(V1)/docs/solutions && $(PY) build.py && $(PY) style.py && cp Numerov_Exercise_Solutions.docx ../../../dist/
	cd $(V2)/docs && cat part0.md part1.md part2.md part3.md part4.md part5.md > course.md && $(PY) build.py && $(PY) style.py \
	  && cp Special_Functions_Numerov_Course.docx ../../dist/

clean:
	rm -f $(V1)/fortran/p1 $(V1)/fortran/p2 $(V1)/fortran/p3 $(V2)/fortran/ch[1234] $(V1)/exercises/ex10_sp $(V1)/exercises/ex10_dp
	find . -name '*.mod' -delete
	find . -type d -name '__pycache__' -prune -exec rm -rf {} +
	rm -rf dist
