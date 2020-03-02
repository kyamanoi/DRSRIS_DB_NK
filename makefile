TARGET = maccormack
OBJECTS = maccormack.o Initial_condition.o corrector.o boundary_condition_c.o fwrite_asc.o boundary_condition_p.o fwrite_xyz.o predictor.o check_cfl.o variables.o dry_or_wet.o
MOD_FILES  = Initial_condition.mod corrector.mod boundary_condition_c.mod fwrite_asc.mod boundary_condition_p.mod fwrite_xyz.mod predictor.mod check_cfl.mod variables.mod dry_or_wet.mod
# FC = ifort
# FC = ifort
# FC = mpifrtpx
FC = mpif90

FLAGS =

# for gfortran
ifeq (${FC},gfortran)
	# FLAGS += 
	# FLAGS += -fimplicit-none -fbounds-check
	#-fopenmp 
endif

# for ifort
ifeq (${FC},ifort)
	FLAGS += -fopenmp -parallel
	# FLAGS += -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback
	# FLAGS +=  -check all -fpe0 -ftrapuv -traceback
endif

ifeq (${FC},mpifrtpx)
	# for development
	# FLAGS += -Kopenmp -NRtrap -O1
	# for optimisation
	#FLAGS += -Kopenmp -Kparallel -Koptmsg=2 -Nlst=t -O3 -Kdalign -Kprefetch_conditional -Komitfp
	FLAGS += -Kopenmp -Kparallel -Koptmsg=2 -Nlst=t -Kfast -Ksimd=2 -Kocl -Kpreex -Karray_private
	# for debug
	# FLAGS += mpifrtpx -NRtrap -Nquickdbg -Nsetvalue -Kopenmp -Kparallel -Koptmsg=2 -Nlst=t 
	# for fastest_compile
	# FLAGS += -Kopenmp -O0
endif

maccormack: maccormack.o Initial_condition.o corrector.o boundary_condition_c.o fwrite_asc.o boundary_condition_p.o fwrite_xyz.o predictor.o check_cfl.o variables.o dry_or_wet.o sub_mpi.o
	${FC} ${FLAGS} -o maccormack maccormack.o Initial_condition.o corrector.o boundary_condition_c.o fwrite_asc.o boundary_condition_p.o fwrite_xyz.o predictor.o check_cfl.o variables.o dry_or_wet.o sub_mpi.o
maccormack.o: maccormack.f90 Initial_condition.o corrector.o boundary_condition_c.o fwrite_asc.o boundary_condition_p.o fwrite_xyz.o predictor.o check_cfl.o variables.o dry_or_wet.o sub_mpi.o
	${FC} ${FLAGS} -c maccormack.f90
variables.o: variables.f90
	${FC} ${FLAGS} -c variables.f90
Initial_condition.o: Initial_condition.f90 variables.o dry_or_wet.o
	${FC} ${FLAGS} -c Initial_condition.f90
corrector.o: corrector.f90 variables.o dry_or_wet.o
	${FC} ${FLAGS} -c corrector.f90
boundary_condition_c.o: boundary_condition_c.f90 variables.o
	${FC} ${FLAGS} -c boundary_condition_c.f90
fwrite_asc.o: fwrite_asc.f90 variables.o
	${FC} ${FLAGS} -c fwrite_asc.f90
boundary_condition_p.o: boundary_condition_p.f90 variables.o
	${FC} ${FLAGS} -c boundary_condition_p.f90
fwrite_xyz.o: fwrite_xyz.f90 variables.o
	${FC} ${FLAGS} -c fwrite_xyz.f90
predictor.o: predictor.f90 variables.o dry_or_wet.o
	${FC} ${FLAGS} -c predictor.f90
check_cfl.o: check_cfl.f90  variables.o
	${FC} ${FLAGS} -c check_cfl.f90
dry_or_wet.o: dry_or_wet.f90 
	${FC} ${FLAGS} -c dry_or_wet.f90
sub_mpi.o: sub_mpi.f90 variables.o
	${FC} ${FLAGS} -c sub_mpi.f90 
clean:
	rm -f maccormack maccormack.o Initial_condition.o corrector.o boundary_condition_c.o fwrite_asc.o boundary_condition_p.o fwrite_xyz.o predictor.o check_cfl.o variables.o dry_or_wet.o
