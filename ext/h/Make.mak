#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

ARCHPATH := $(shell dpkg-architecture -q DEB_BUILD_MULTIARCH)

TCLTK_CONF_DIR?=/usr/lib/$(ARCHPATH)/tcl8.6

#-------------------------------------------------------------------------------

AR = @ar ru
CP = @-cp
RM = @-rm
MV = mv
EDIT_PROG  = xemacs
MAKE = make

SOFT_DIR = $(HOME)/.MSOFT

#-------------------------------------------------------------------------------

_OS := $(shell uname -s)

ifeq ($(_OS), Linux)
else
endif

_NODENAME := $(shell uname -n)


#------------------------------------------------

include $(TCLTK_CONF_DIR)/tclConfig.sh

#-------------------------------------------------------------------------------

TCL_INCS=$(TCL_INCLUDE_SPEC) 

TCL_LIB_SPEC_clean	= $(subst ', , $(TCL_LIB_SPEC))

TCL_LIBS=$(TCL_LIB_SPEC_clean) 

TCL_STUB_LIB_SPEC_clean = $(subst ', , $(TCL_STUB_LIB_SPEC))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
