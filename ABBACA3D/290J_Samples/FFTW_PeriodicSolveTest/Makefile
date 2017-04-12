#
# This is the "meta" makefile that builds the sample executables
# by invoking their respective makescripts.
#
# Usage :
#
# make release                   === builds release versions of all executables
# make debug                     === builds debug versions of all executables
#
# make clean                     === removes object files and executables
# make cleanall                  === removes object files, executables and temporary directories
#
#
# To build FFTW 3d samples that use the multi-threaded capabilities, add the specification
# _FFTW_OPENMP=1 to the make invocation, e.g. make .FFTWsinSolveTest _FFTW_OPENMP=1 
# 
# Dec. 5, 2015
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


release : .Starting   .fftw_1dSolveTest   .fftw_2dSolveTest .fftw_3dSolveTest .Finished 

debug   : .Starting   .fftw_1dSolveTest_debug .fftw_2dSolveTest_debug  .fftw_3dSolveTest_debug .Finished  

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

 
.fftw_1dSolveTest :
	###########################################################
	#               .fftw_1dSolveTest
	###########################################################
	$(QUIET)$(MAKE) -f fftw_1dSolveTest.mk  release 
.fftw_1dSolveTest_debug :
	###########################################################
	#               .fftw_1dSolveTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_1dSolveTest.mk  debug

.fftw_2dSolveTest :
	###########################################################
	#               .fftw_2dSolveTest
	###########################################################
	$(QUIET)$(MAKE) -f fftw_2dSolveTest.mk  release 
.fftw_2dSolveTest_debug :
	###########################################################
	#               .fftw_2dSolveTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_2dSolveTest.mk  debug


.fftw_3dSolveTest :
	###########################################################
	#               .fftw_3dSolveTest
	###########################################################
	$(QUIET)$(MAKE) -f fftw_3dSolveTest.mk  release 
.fftw_3dSolveTest_debug :
	###########################################################
	#               .fftw_3dSolveTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_3dSolveTest.mk  debug
	

.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f fftw_1dSolveTest.mk       clean 
	$(QUIET)$(MAKE) -f fftw_2dSolveTest.mk       clean
	$(QUIET)$(MAKE) -f fftw_3dSolveTest.mk       clean



