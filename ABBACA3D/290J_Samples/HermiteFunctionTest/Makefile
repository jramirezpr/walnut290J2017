#
# This is the "meta" makefile that builds the sample executables
# by invoking their respective makescripts.
#
#
# Usage :
#
# make release                   === builds release versions of all executables
# make debug                     === builds debug versions of all executables
#
# make .NameOfSamplesDir         === builds the release version of the executables in the samples directory specified
# make .NameOfSamplesDir_debug   === builds the debug version of the executables in the samples directory specified
#
# make clean                     === removes object files and executables
# make cleanall                  === removes object files, executables and any supporting libraries 
#
# Jan. 8, 2016
#
SHELL=/bin/sh

# Location of 290J source files and Makescripts files

290J_Dir=../../290J_2015
MAKESCRIPTS_Dir=$(290J_Dir)/MakeScripts

# Use BaseCommonConfig.mk if it exists, otherwise use BaseCommonConfig_Default.mk 

ifneq ("$(wildcard $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk)","")
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk
else
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig_Default.mk
endif


ifeq ($(MAKECMDGOALS),release)
BUILD_TYPE=release
endif
ifeq ($(MAKECMDGOALS),debug)
BUILD_TYPE=debug
endif

release : .Starting   .Hermite1dTest   .Hermite2dTest .Hermite3dTest  .Finished 
debug   : .Starting   .Hermite1dTest_debug .Hermite2dTest_debug .Hermite3dTest_debug .Finished  

.Starting :
	###########################################################
	#            Build Started
	########################################################### 
	$(QUIET) echo "Compilation Date : " `date` 
ifeq ($(OpenMP),0)
	$(QUIET) echo "Executable Type  :  Single-threaded" 
else
	$(QUIET) echo "Executable Type  :  OpenMP based multi-threaded" 
endif

 

.Hermite1dTest :
	###########################################################
	#               .Hermite1dTest
	###########################################################
	$(QUIET)$(MAKE) -f Hermite1dTest.mk  release 
.Hermite1dTest_debug :
	###########################################################
	#               .Hermite1dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Hermite1dTest.mk  debug

.Hermite2dTest :
	###########################################################
	#               .Hermite2dTest
	###########################################################
	$(QUIET)$(MAKE) -f Hermite2dTest.mk  release 
.Hermite2dTest_debug :
	###########################################################
	#               .Hermite2dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Hermite2dTest.mk  debug


.Hermite3dTest :
	###########################################################
	#               .Hermite3dTest
	###########################################################
	$(QUIET)$(MAKE) -f Hermite3dTest.mk  release 
.Hermite3dTest_debug :
	###########################################################
	#               .Hermite3dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Hermite3dTest.mk  debug


.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f Hermite1dTest.mk                  clean 
	$(QUIET)$(MAKE) -f Hermite2dTest.mk                  clean
	$(QUIET)$(MAKE) -f Hermite3dTest.mk                  clean




