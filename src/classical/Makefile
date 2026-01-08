FC = gfortran
FFLAGS = -O2 -Wall -Wextra -std=f2008

OBJS = constants.o standard_map.o diagnostics.o main_standardmap.o
EXE  = standardmap

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

constants.o: constants.f90
	$(FC) $(FFLAGS) -c $<

standard_map.o: standard_map.f90 constants.o
	$(FC) $(FFLAGS) -c $<

diagnostics.o: diagnostics.f90 constants.o standard_map.o
	$(FC) $(FFLAGS) -c $<

main_standardmap.o: main_standardmap.f90 constants.o standard_map.o diagnostics.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE)

run: $(EXE)
	./$(EXE)
