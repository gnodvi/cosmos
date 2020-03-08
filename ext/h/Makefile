
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

MAIN_SRC_DIR = .

include $(MAIN_SRC_DIR)/Make.mak


COMMON_DIR = $(MAIN_SRC_DIR)/../common
ABMD_DIR   = $(MAIN_SRC_DIR)/../iaaras/abmd
NSTARS_DIR = $(MAIN_SRC_DIR)/../nstars/R
AA-56M_DIR = $(MAIN_SRC_DIR)/../aa-56m
DE118I_DIR = $(MAIN_SRC_DIR)/../de118i
EPHOLD_DIR = $(MAIN_SRC_DIR)/../iaaras/eph-old-access
STARCH_DIR = $(MAIN_SRC_DIR)/../starch/R
NAVALO_DIR = $(MAIN_SRC_DIR)/../navalo/R

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
