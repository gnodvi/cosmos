
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

include ./Make.mak


AA-56M_DIR = ../aa-56m/R
COMMON_DIR = ../common
ABMD_DIR   = ../iaaras/abmd
NSTARS_DIR = ../nstars/R
DE118I_DIR = ../de118i
EPHOLD_DIR = ../iaaras/eph-old-access
STARCH_DIR = ../starch/R
NAVALO_DIR = ../navalo/R

#-------------------------------------------------------------------------------

all:  ABMD__ NSTARS__ AA-56M__ DE118I__ EPHOLD__ STARCH__ NAVALO__

#-------------------------------------------------------------------------------


NSTARS__ :
	(cd $(NSTARS_DIR); make)

ABMD__  :
	(cd $(ABMD_DIR);   make)

AA-56M__ :
	(cd $(AA-56M_DIR); make)

DE118I__ :
	(cd $(DE118I_DIR); make)

EPHOLD__ :
	(cd $(EPHOLD_DIR); make)

STARCH__ :
	(cd $(STARCH_DIR); make)

NAVALO__ :
	(cd $(NAVALO_DIR); make)

#-------------------------------------------------------------------------------

clean : 
	(cd $(ABMD_DIR);   make clean)
	(cd $(NSTARS_DIR); make clean)
	(cd $(AA-56M_DIR); make clean)
	(cd $(DE118I_DIR); make clean)
	(cd $(EPHOLD_DIR); make clean)
	(cd $(STARCH_DIR); make clean)
	(cd $(NAVALO_DIR); make clean)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
