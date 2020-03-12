
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

include ./Make.mak


AA-56M_DIR = ../aa-56m/R
COMMON_DIR = ../common
DE118I_DIR = ../de118i/R
ABMD_DIR   = ../iaaras/abmd
NSTARS_DIR = ../nstars/R
EPHOLD_DIR = ../iaaras/eph-old-access
STARCH_DIR = ../starch/R
NAVALO_DIR = ../navalo/R

#-------------------------------------------------------------------------------

all:  AA-56M__ DE118I__ ABMD__ EPHOLD__ NAVALO__ NSTARS__ STARCH__ 

#-------------------------------------------------------------------------------


AA-56M__ :
	(cd $(AA-56M_DIR); make)

DE118I__ :
	(cd $(DE118I_DIR); make)

ABMD__  :
	(cd $(ABMD_DIR);   make)

EPHOLD__ :
	(cd $(EPHOLD_DIR); make)

NAVALO__ :
	(cd $(NAVALO_DIR); make)

NSTARS__ :
	(cd $(NSTARS_DIR); make)

STARCH__ :
	(cd $(STARCH_DIR); make)


#-------------------------------------------------------------------------------

clean : 
	(cd $(AA-56M_DIR); make clean)
	(cd $(DE118I_DIR); make clean)
	(cd $(ABMD_DIR);   make clean)
	(cd $(EPHOLD_DIR); make clean)
	(cd $(NSTARS_DIR); make clean)
	(cd $(STARCH_DIR); make clean)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
