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

release : .Starting  .fftw3_XdTest  .fftw3_sinXdTest  .fftw3_cosXdTest  .fftw_scaleTest .Finished 

debug   : .Starting  .fftw3_XdTest_debug .fftw3_sinXdTest_debug .fftw3_cosXdTest_debug .fftw_scaleTest_debug .Finished  

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

.fftw3_XdTest :
	###########################################################
	#                 .fftw3_XdTest 
	###########################################################
	$(QUIET)$(MAKE) -f fftw3_1dTest.mk  release 
	$(QUIET)$(MAKE) -f fftw3_2dTest.mk  release
	$(QUIET)$(MAKE) -f fftw3_3dTest.mk  release

.fftw3_XdTest_debug :
	###########################################################
	#                 .fftw3_XdTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f fftw3_1dTest.mk  debug
	$(QUIET)$(MAKE) -f fftw3_2dTest.mk  debug
	$(QUIET)$(MAKE) -f fftw3_3dTest.mk  debug

.fftw3_sinXdTest :
	###########################################################
	#                 .fftw3_sinXdTest 
	###########################################################
	$(QUIET)$(MAKE) -f fftw3_sin1dTest.mk  release 
	$(QUIET)$(MAKE) -f fftw3_sin2dTest.mk  release
	$(QUIET)$(MAKE) -f fftw3_sin3dTest.mk  release


.fftw3_sinXdTest_debug :
	###########################################################
	#                 .fftw3_sinXdTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f fftw3_sin1dTest.mk  debug
	$(QUIET)$(MAKE) -f fftw3_sin2dTest.mk  debug
	$(QUIET)$(MAKE) -f fftw3_sin3dTest.mk  debug

.fftw3_cosXdTest :
	###########################################################
	#                 .fftw3_cosXdTest 
	###########################################################
	$(QUIET)$(MAKE) -f fftw3_cos1dTest.mk release


.fftw3_cosXdTest_debug :
	###########################################################
	#                 .fftw3_cosXdTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f fftw3_cos1dTest.mk  debug

.fftw_scaleTest :
	###########################################################
	#                 .fftw_scaleTest
	###########################################################
	$(QUIET)$(MAKE) -f fftw_scalingTest.mk    release 
	$(QUIET)$(MAKE) -f fftw_scaledSinTest.mk  release
	$(QUIET)$(MAKE) -f fftw_scaledCosTest.mk  release 

.fftw_scaleTest_debug :
	###########################################################
	#                 .fftw_scaleTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f fftw_scalingTest.mk    debug
	$(QUIET)$(MAKE) -f fftw_scaledSinTest.mk  debug
	$(QUIET)$(MAKE) -f fftw_scaledCosTest.mk  debug

.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*

clean: 
	$(QUIET)$(MAKE) -f fftw3_1dTest.mk               clean 
	$(QUIET)$(MAKE) -f fftw3_2dTest.mk               clean
	$(QUIET)$(MAKE) -f fftw3_3dTest.mk               clean
	$(QUIET)$(MAKE) -f fftw3_sin1dTest.mk            clean 
	$(QUIET)$(MAKE) -f fftw3_sin2dTest.mk            clean
	$(QUIET)$(MAKE) -f fftw3_sin3dTest.mk            clean
	$(QUIET)$(MAKE) -f fftw3_cos1dTest.mk            clean
	$(QUIET)$(MAKE) -f fftw_scalingTest.mk           clean 
	$(QUIET)$(MAKE) -f fftw_scaledSinTest.mk         clean
	$(QUIET)$(MAKE) -f fftw_scaledCosTest.mk         clean




