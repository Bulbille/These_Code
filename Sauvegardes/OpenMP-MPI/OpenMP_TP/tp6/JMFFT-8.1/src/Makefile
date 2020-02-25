############################ -*- Mode: Makefile -*- ###########################
## Makefile --- 
## 
## Auteur          : Jean-Marie Teuler (CNRS) <Jean-Marie.Teuler@lcp.u-psud.fr.fr>
## Créé le         : Thu May 13 11:52:48 2004
## Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Dern. mod. le   : Thu May 13 14:49:58 2004
###############################################################################

include ../arch/Make.inc

LIBNAME = jmfft
LIB = lib$(LIBNAME).a

.SUFFIXES:
.SUFFIXES: .f90 .o .a

include Makefile.objs

all: $(OBJS)

.f90.a:
	$(FC) $(FFLAGS) -c $(*).f90
	$(AR) $(ARFLAGS) $(LIB) $(*).o
	$(RM) $(RMFLAGS) $(*).o

clean:
	$(RM) $(RMFLAGS) *.f *.o *.a
