#===============================================================================
# SuitOPT/Makefile:
#===============================================================================

default: all

#-------------------------------------------------------
# all - Makes SuiteSparse and SuiteOPT shared libraries
#-------------------------------------------------------

all:
	( cd SuiteOPT ; $(MAKE) py )

#-------------------------------------
# Purge
#-------------------------------------

purge:
	( cd SuiteOPT ; $(MAKE) purge )

#-------------------------------------
# Clean
#-------------------------------------

clean:
	( cd SuiteOPT ; $(MAKE) clean )