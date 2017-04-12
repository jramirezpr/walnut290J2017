#
# RitzMethod2d.exe program makefile using MakeScripts
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


CPPfiles   +=  RitzMethod2d.cpp
INCLUDES   += -I./
INCLUDES   += -I../LaplaceOp2dTest
INCLUDES   += -I$(290J_Dir)

# Add flag for compilation with c++11 constructs

CXXDEFINES   += -std=c++11

RELEASE_EXEC  =  RitzMethod2d.exe
DEBUG_EXEC    =  RitzMethod2d_debug.exe

RELEASE_DIR  = ./_releaseRitzMethod2d
DEBUG_DIR    = ./_debugRitzMethod2d

include $(290J_Dir)/MakeScripts/ExecutableMake.mk

