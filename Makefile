
OVERVIEWTEXFILES=$(wildcard overview*.tex)

OVERVIEWFILES=$(OVERVIEWTEXFILES:.tex=.pdf)

LECTURETEXFILES=$(wildcard lecture*.tex)

# Keeping manual list of lectures
#LECTUREFILES=$(LECTURETEXFILES:.tex=.pdf)

LECTUREFILES=lecture_00_fem_introduction.pdf \
	     lecture_00_python_programming.pdf \
             lecture_01_fenics_installation.pdf \
             lecture_01_fenics_introduction.pdf \
             lecture_02_static_linear_pdes.pdf \
             lecture_02_static_linear_pdes_alt.pdf \
             lecture_03_static_nonlinear_pdes.pdf \
             lecture_04_time_dependent_pdes.pdf \
             lecture_04_time_dependent_pdes_alt.pdf \
             lecture_05_happy_hacking.pdf \
             lecture_06_static_hyperelasticity.pdf \
             lecture_07_dynamic_hyperelasticity.pdf \
             lecture_08_stokes.pdf \
             lecture_08_stokes_alt.pdf \
             lecture_09_navier_stokes.pdf \
             lecture_09_navier_stokes_alt.pdf \
             lecture_10_discontinuous_galerkin.pdf \
             lecture_11_error_control.pdf \
             lecture_12_sensitivities.pdf \
             lecture_13_introduction_to_dolfin_adjoint.pdf \
             lecture_14_from_sensitivities_to_optimisation.pdf \
             lecture_15_one_shot_optimisation.pdf \
             lecture_16_optimisation_challenge.pdf \
             lecture_17_cpp_programming.pdf \
             lecture_18_cpp_fenics.pdf \
             lecture_19_fenics_implementation.pdf \
	     lecture_21_tools_for_online_collaboration.pdf \
             lecture_22_linear_elasticity.pdf

# List of files for convenient 'make course' while working on a
# specific course, modify this list to the lectures you want to include
CURRENTCOURSEFILES=overview_simula_2016.pdf \
      lecture_02_static_linear_pdes.pdf \
      lecture_03_static_nonlinear_pdes.pdf \
      lecture_04_time_dependent_pdes.pdf \
      lecture_05_happy_hacking.pdf \
      lecture_06_static_hyperelasticity.pdf \
      lecture_08_stokes.pdf \
      lecture_09_navier_stokes.pdf \


.PHONY: all
all: overviews lectures

.PHONY: overviews
overviews: $(OVERVIEWFILES)

.PHONY: lectures
lectures: $(LECTUREFILES)

.PHONY: course
course: $(CURRENTCOURSEFILES)
	mkdir -p course
	rm course/*.pdf
	cp $(CURRENTCOURSEFILES) course

course/course.pdf: course
	(cd course && pdfjoin --outfile course.pdf $(CURRENTCOURSEFILES))

.PHONY: clean
clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *~
	(cd output && rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *~)
	rm -f `find -name \*.pyc` `find -name \*~`
	rm -f `find -name \*.pvd` `find -name \*.vtu`

.PHONY: purge
purge: clean
	rm -rf output course
	rm -f $(OVERVIEWFILES) $(LECTUREFILES)

# Optimal dependencies from %.tex file: (don't know how to make these dependencies of %.pdf)
#DEPS=$(shell grep '\\input{' %.tex | cut -d { -f 2 | cut -d } -f 1)

%.pdf: %.tex $(wildcard slides/*.tex)
	mkdir -p output
	pdflatex --output-directory=output $*.tex
	mv output/$*.pdf .
