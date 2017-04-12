#
# ParameterListTest Program Makefile
#
# This makefile assumes a directory structure
#
# Course Directory
# --------------------------------------------------------------- 
#       |                |                |            |                
# Project_1 Project_2   ....         290J_Samples   290J_2015
#
#
# The make file depends upon the XML_ParameterListLib being created
# prior to invoking this make file. This library can be 
# manually created by entering the 290J_2015/XML_ParameterList 
# directory and executing 
#
# make -f XML_ParameterListLib.mk release 
#
#
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


CPPfiles  +=  ParameterListTest.cpp
INCLUDES   = -I./
INCLUDES  += -I$(290J_Dir)/

LIB_PATH  += -L$(290J_Dir)/XML_ParameterList/lib

ifeq ($(MAKECMDGOALS),debug)
LIBS      += -lXML_ParameterList_debug
endif
ifeq ($(MAKECMDGOALS),release)
LIBS      += -lXML_ParameterList
endif


RELEASE_EXEC  = ParameterListTest.exe
DEBUG_EXEC    = ParameterListTest_debug.exe

RELEASE_DIR  = ./_releaseParameterListTest
DEBUG_DIR    = ./_debugParameterListTest

include $(290J_Dir)/MakeScripts/ExecutableMake.mk

