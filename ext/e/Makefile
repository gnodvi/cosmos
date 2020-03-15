
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


AA-56M_DIR = ../aa-56m/R
COMMON_DIR = ../common
DE118I_DIR = ../de118i/R
IAARAS_DIR = ../iaaras
NSTARS_DIR = ../nstars/R
STARCH_DIR = ../starch/R
NAVALO_DIR = ../navalo/R

#-------------------------------------------------------------------------------

all:  AA-56M__ DE118I__ IAARAS__ NAVALO__ NSTARS__ STARCH__ 

#-------------------------------------------------------------------------------


AA-56M__ :
	(cd $(AA-56M_DIR); make)

DE118I__ :
	(cd $(DE118I_DIR); make)

IAARAS__  :
	(cd $(IAARAS_DIR); make)

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
	(cd $(IAARAS_DIR); make clean)
	(cd $(NSTARS_DIR); make clean)
	(cd $(STARCH_DIR); make clean)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
