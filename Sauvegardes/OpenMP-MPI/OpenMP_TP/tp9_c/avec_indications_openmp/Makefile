include ../../arch/make_arch_inc
include ../parametres/make_tp_inc

# Build $OBJ from $SRC
TMPOBJ   ?= $(SRC:.F90=.o)
OBJ      ?= $(TMPOBJ:.c=.o)

OBJ_MONO ?= $(OBJ:.o=_mono.o)
OBJ_OMP  ?= $(OBJ:.o=_omp.o)

LINKER   ?= $(CC)

default: mono para

mono: $(OBJ_MONO) $(DEPLIB)
	$(LINKER) $(OBJ_MONO) -o $(CMD)_mono.exe $(LDFLAGS) $(CLDFLAGS) $(LDFLAGS_MONO)

para: $(OBJ_OMP) $(DEPLIB)
	$(LINKER) $(OBJ_OMP) -o $(CMD)_omp.exe $(LDFLAGS) $(CLDFLAGS) $(LDFLAGS_OMP)

%_mono.o: %.F90
	$(FC) $(FFLAGS) $(FFLAGS_MONO) -c $< -o $@

%_omp.o: %.F90
	$(FC) $(FFLAGS) $(FFLAGS_OMP) -c $< -o $@

%_mono.o: %.c
	$(CC) $(CFLAGS) $(CFLAGS_MONO) -c $< -o $@

%_omp.o: %.c
	$(CC) $(CFLAGS) $(CFLAGS_OMP) -c $< -o $@

clean:
	rm -f *.o core.*

cleanall: clean $(DEPCLEAN)
	rm -f *.exe *.res

visu:
	../../scripts/speedup.sh
