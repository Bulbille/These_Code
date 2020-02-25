############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Construction de la bibliotheque JMFFT
##
## Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Créé le         : Tue Dec  8 18:18:21 1998
## Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
## Dern. mod. le   : Wed Jun 12 11:20:14 2013
###############################################################################

SHELL = /bin/ksh

include arch/Make.inc

all: jmfft

jmfft: lib/libjmfft.a

lib/libjmfft.a:
	@echo
	@echo "*********** Construction de libjmfft.a ***********"
	@echo
	(cd src; make)
	mv src/libjmfft.a lib/

clean:
	(cd src; make clean)
	rm -f lib/libjmfft.a