#
# Sample_2 : A makefile for creating an executable 
#            using the MakeScripts makefile scripts 
#            This sample demonstrates the specification
#            of additional include paths and incorporation
#            of multiple source files. 
#
# To build the executable using this script use the
# command
#
# make -f Sample_2.mk release 
#
# or 
#
# make -f Sample_2.mk debug
#
# To see a verbose output of the build add VERBOSE=1 to the make invocation commands above 
#
#
# To clean up use 
#
# make -f Sample_2.mk clean 
#
# C. Anderson July 30, 2015
#
# Note: If you are using the makefiles within a repository, you'll want to 
# make sure to exclude the executables, the *.o and *.d files created. 
# For this script, this can be accomplished by excluding *.exe, _release*,
# and _debug* (directories that contain the created *.o and *.d files).  
#

SHELL=/bin/sh

# Specifying the location of the MakeScripts *.mk files

SCRIPTS_PATH=../../../MakeScripts

# Use BaseCommonConfig.mk if it exists otherwise use BaseCommonConfig_Default.mk 

ifneq ("$(wildcard $(SCRIPTS_PATH)/BaseCommonConfig.mk)","")
	include $(SCRIPTS_PATH)/BaseCommonConfig.mk
else
	include $(SCRIPTS_PATH)/BaseCommonConfig_Default.mk
endif


# Specifying the cpp files to use. For additional files make sure to 
# use the += operator.

CPPfiles  +=  ./Sample_2.cpp
CPPfiles  +=  ./Subdirectory/Sample_2_Class.cpp

# Specifying the include files to use. Additional directories
# are specified using 
#
# INCLUDES += -I/[directory to include]
#
#
INCLUDES  +=   -I./
INCLUDES  +=   -I./Subdirectory

# Specifing the release and debug directories for the *.o (object files) and 
# *.d (depenency) files for the different builds 

RELEASE_DIR  = ./_releaseSample_2
DEBUG_DIR    = ./_debugSample_2d

# Specifying the name of the executables to be created 

RELEASE_EXEC  = Sample_2.exe
DEBUG_EXEC    = Sample_2_debug.exe


# Includng the script that creates the executable

include $(SCRIPTS_PATH)/ExecutableMake.mk


 
