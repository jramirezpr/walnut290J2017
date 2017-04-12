#
# SqTest program makefile using MakeScripts
#
# The name and location of the sqlite3 libraries 
# must be set by specifying 
#
# SQLITE_PATH and SQLITE_LIB in the BaseCommonConfig.mk
#
# If sqlite3 is not installed on your system, you can 
# build the library from source using the makescript
# SIMdb_SQLite/sqlite3lib.mk. After building, you can
# then set the name and paths to library constructed. 
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


CPPfiles +=   SqTest.cpp
 
INCLUDES +=  -I$(290J_Dir)/SimDb_SQLite
INCLUDES +=  -I$(290J_Dir)/SimDb_SQLite/sqlite

LIB_PATH +=  $(SQLITE_PATH)
LIBS     +=  $(SQLITE_LIB)

RELEASE_EXEC  = SqTest.exe
DEBUG_EXEC    = SqTestd.exe


RELEASE_DIR  = ./_releaseSqTest
DEBUG_DIR    = ./_debugSqTest

include $(290J_Dir)/MakeScripts/ExecutableMake.mk


