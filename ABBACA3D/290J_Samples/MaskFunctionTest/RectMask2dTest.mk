#
# RectMask2dTest program makefile using MakeScripts
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

CPPfiles      +=   RectMask2dTest.cpp

INCLUDES   += -I./
INCLUDES   += -I$(290J_Dir)


RELEASE_EXEC  = RectMask2dTest.exe
DEBUG_EXEC    = RectMask2dTest_debug.exe


RELEASE_DIR  = ./_releaseRectMask2dTest
DEBUG_DIR    = ./_debugRectMask2dTest

include $(MAKESCRIPTS_Dir)/ExecutableMake.mk

