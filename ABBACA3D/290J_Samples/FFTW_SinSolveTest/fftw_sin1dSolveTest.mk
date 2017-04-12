#
# fftw_sin1dSolveTest program makefile using MakeScripts
#
# This makefile assumes a directory structure
#
# Course Directory
# --------------------------------------------------------------- 
#       |                |                |            |                
# Project_1 Project_2   ....         290J_Samples   290J_2015
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


CPPfiles   +=  fftw_sin1dSolveTest.cpp
INCLUDES    = -I./
INCLUDES   += -I$(290J_Dir)

# Specifying FFTW3 libraries and defines when OpenMP version is to be used. 

LIBS     += $(FFTW_LIB)  
LIB_PATH += $(FFTW_PATH)

ifeq ($(_FFTW_OPENMP),1)
CXXDEFINES   += -D_FFTW_OPENMP
endif


RELEASE_EXEC  = fftw_sin1dSolveTest.exe
DEBUG_EXEC    = fftw_sin1dSolveTest_debug.exe

RELEASE_DIR  = ./_releasefftw_sin1dSolveTest
DEBUG_DIR    = ./_debugfftw_sin1dSolveTest

include $(290J_Dir)/MakeScripts/ExecutableMake.mk

