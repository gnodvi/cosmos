
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


ABMD_DIR   = ./abmd/R
EPHOLD_DIR = ./eph-old-access/R
EPHNEW_DIR = ./eph-new-access/R

#-------------------------------------------------------------------------------

all:  ABMD__ EPHOLD__ EPHNEW__ 

#-------------------------------------------------------------------------------


ABMD__  :
	(cd $(ABMD_DIR);   make)

EPHOLD__ :
	(cd $(EPHOLD_DIR); make)

EPHNEW__ :
	(cd $(EPHNEW_DIR); make)

#-------------------------------------------------------------------------------

clean : 
	(cd $(ABMD_DIR);   make clean)
	(cd $(EPHOLD_DIR); make clean)
	(cd $(EPHNEW_DIR); make clean)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
