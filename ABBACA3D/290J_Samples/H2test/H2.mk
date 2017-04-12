#
# H2.exe program makefile using MakeScripts
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


CPPfiles   +=  H2.cpp
CPPfiles   +=  ./$(290J_Dir)/MultiParticleHamiltonian/MultiParticleHamiltonianMatrix.cpp
CPPfiles   +=  ./$(290J_Dir)/MultiParticleHamiltonian/SlaterDeterminantReference.cpp
INCLUDES   += -I./
INCLUDES   += -I$(290J_Dir)

LIBS     += $(FFTW_LIB)  
LIB_PATH += $(FFTW_PATH)


RELEASE_EXEC  =  H2.exe
DEBUG_EXEC    =  H2_debug.exe

RELEASE_DIR  = ./_releaseH2
DEBUG_DIR    = ./_debugH2

include $(290J_Dir)/MakeScripts/ExecutableMake.mk

