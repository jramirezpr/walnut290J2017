#
# LegendreInterp3dTest program makefile using MakeScripts
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

CPPfiles     +=   LegendreInterp3dTest.cpp

INCLUDES      += -I./
INCLUDES      += -I$(290J_Dir)

RELEASE_EXEC  = LegendreInterp3dTest.exe
DEBUG_EXEC    = LegendreInterp3dTest_debug.exe


RELEASE_DIR  = ./_releaseLegendreInterp3dTest
DEBUG_DIR    = ./_debugLegendreInterp3dTest

include $(290J_Dir)/MakeScripts/ExecutableMake.mk

