# $Id: GNUmakefile,v 1.2 2000/10/19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := orb
G4TARGET := $(name)
G4EXLIB := true

#ifndef G4INSTALL
#  G4INSTALL = ../../..
#endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

# Show activity with keyboard LED under X (with xset led 3):
CPPFLAGS += -DACTIVITY_LED

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

