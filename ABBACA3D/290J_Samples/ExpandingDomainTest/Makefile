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
# make clanall                   === removes object files, executables and any supporting libraries 
#
# Nov. 12,2015
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


release : .Starting   .Poisson1dTest   .Poisson2dTest .Poisson3dTest .ScreenedPoisson2dTest .ScreenedPoisson3dTest .Finished 

debug   : .Starting   .Poisson1dTest_debug .Poisson2dTest_debug  .Poisson3dTest_debug  .ScreenedPoisson2dTest_debug .ScreenedPoisson3dTest  .Finished  

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

 

.Poisson1dTest :
	###########################################################
	#               .Poisson1dTest
	###########################################################
	$(QUIET)$(MAKE) -f Poisson1dTest.mk  release 
.Poisson1dTest_debug :
	###########################################################
	#               .Poisson1dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Poisson1dTest.mk  debug

.Poisson2dTest :
	###########################################################
	#               .Poisson2dTest
	###########################################################
	$(QUIET)$(MAKE) -f Poisson2dTest.mk  release 
.Poisson2dTest_debug :
	###########################################################
	#               .Poisson2dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Poisson2dTest.mk  debug


.Poisson3dTest :
	###########################################################
	#               .Poisson3dTest
	###########################################################
	$(QUIET)$(MAKE) -f Poisson3dTest.mk  release 
.Poisson3dTest_debug :
	###########################################################
	#               .Poisson3dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f Poisson3dTest.mk  debug


.ScreenedPoisson2dTest :
	###########################################################
	#               .ScreenedPoisson2dTest
	###########################################################
	$(QUIET)$(MAKE) -f ScreenedPoisson2dTest.mk  release 
.ScreenedPoisson2dTest_debug :
	###########################################################
	#               .ScreenedPoisson2dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f ScreenedPoisson2dTest.mk  debug

.ScreenedPoisson3dTest :
	###########################################################
	#               .ScreenedPoisson3dTest
	###########################################################
	$(QUIET)$(MAKE) -f ScreenedPoisson3dTest.mk  release 
.ScreenedPoisson3dTest_debug :
	###########################################################
	#               .ScreenedPoisson3dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f ScreenedPoisson3dTest.mk  debug



.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f Poisson1dTest.mk       clean 
	$(QUIET)$(MAKE) -f Poisson2dTest.mk       clean
	$(QUIET)$(MAKE) -f Poisson3dTest.mk       clean
	$(QUIET)$(MAKE) -f ScreenedPoisson2dTest.mk       clean
	$(QUIET)$(MAKE) -f ScreenedPoisson3dTest.mk       clean


