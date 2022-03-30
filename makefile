F90 = gfortran
ACMLDIR :=
FLAGS :=   -fbounds-check
LIBACML :=
all:  rns.x

rns.x:  DR_rns.f90 comp_f_p.f90 TOV_guess.f90 iteration.f90 comp.f90 spin_up_freq.f90 kepler.F m0_models.f90 


	$(F90) $(LIBACML)   $(FLAGS) DR_rns.f90 comp.f90 comp_f_p.f90 TOV_guess.f90 iteration.f90 spin_up_freq.f90 kepler.F m0_models.f90 -o rns.x


clean:
	rm *.x
