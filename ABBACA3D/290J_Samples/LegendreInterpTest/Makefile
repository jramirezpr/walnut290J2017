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
#
# To build FFTW 3d samples that use the multi-threaded capabilities, add the specification
# _FFTW_OPENMP=1 to the make invocation, e.g. make .FFTWsinSolveTest _FFTW_OPENMP=1 
# 
#
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


release : .Starting   .LegendreInterp1dTest   .LegendreInterp2dTest .LegendreInterp3dTest .Finished 

debug   : .Starting   .LegendreInterp1dTest_debug .LegendreInterp2dTest_debug  .LegendreInterp3dTest_debug .Finished  

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

 

.LegendreInterp1dTest :
	###########################################################
	#               .LegendreInterp1dTest
	###########################################################
	$(QUIET)$(MAKE) -f LegendreInterp1dTest.mk  release 
.LegendreInterp1dTest_debug :
	###########################################################
	#               .LegendreInterp1dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f LegendreInterp1dTest.mk  debug

.LegendreInterp2dTest :
	###########################################################
	#               .LegendreInterp2dTest
	###########################################################
	$(QUIET)$(MAKE) -f LegendreInterp2dTest.mk  release 
.LegendreInterp2dTest_debug :
	###########################################################
	#               .LegendreInterp2dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f LegendreInterp2dTest.mk  debug


.LegendreInterp3dTest :
	###########################################################
	#               .LegendreInterp3dTest
	###########################################################
	$(QUIET)$(MAKE) -f LegendreInterp3dTest.mk  release 
.LegendreInterp3dTest_debug :
	###########################################################
	#               .LegendreInterp3dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f LegendreInterp3dTest.mk  debug

.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f LegendreInterp1dTest.mk       clean 
	$(QUIET)$(MAKE) -f LegendreInterp2dTest.mk       clean
	$(QUIET)$(MAKE) -f LegendreInterp3dTest.mk       clean



