# Makefile for RAMSES utils
#F90=gfortran
#F90 = ifort 
F90 = mpif90
BINDIR=.
#FLAGS=-ffree-line-length-none
FFLAGS = -cpp -fast

MODOBJ=random.o io_ramses.o

# Make all targets
all: amr2prof amr2cylprof ramses2tipsy amr2map part2map part2prof part2birth
amr2prof: $(MODOBJ) amr2prof.o
	$(F90) $(MODOBJ) $(FLAGS) amr2prof.o -o $(BINDIR)/amr2prof
amr2cylprof: $(MODOBJ) amr2cylprof.o
	$(F90) $(MODOBJ) $(FLAGS) amr2cylprof.o -o $(BINDIR)/amr2cylprof
ramses2tipsy: $(MODOBJ) ramses2tipsy.o
	$(F90) $(MODOBJ) $(FLAGS) ramses2tipsy.o -o $(BINDIR)/ramses2tipsy
amr2map: $(MODOBJ) amr2map.o
	$(F90) $(FLAGS) amr2map.o -o $(BINDIR)/amr2map
part2map: $(MODOBJ) part2map.o
	$(F90) $(FLAGS) part2map.o -o $(BINDIR)/part2map
part2prof: $(MODOBJ) part2prof.o
	$(F90) $(FLAGS) part2prof.o -o $(BINDIR)/part2prof
part2birth: $(MODOBJ) part2birth.o
	$(F90) $(FLAGS) part2birth.o -o $(BINDIR)/part2birth
part2sfr: $(MODOBJ) part2sfr.o
	$(F90) $(FLAGS) part2sfr.o -o $(BINDIR)/part2sfr
part2cube: $(MODOBJ) part2cube.o
	$(F90) $(FLAGS) part2cube.o -o $(BINDIR)/part2cube

# Make a specific object file
%.o: %.f90
	$(F90) $(FLAGS) -c $^
clean:
	rm *.o *.mod
