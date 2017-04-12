#
# This is the high level makefile that first builds the
# XML_ParameterListLib library (archive) and then 
# builds the ParameterListTest executable. 
#
# Note: The indentation in the specification of commands
# following a target are tabs --- not spaces! Tabs are 
# significant in make files (arrgh). 
#
SHELL=/bin/sh

# Location of 290J source files 

290J_Dir=../../290J_2015

include $(290J_Dir)/MakeScripts/BaseCommonConfig.mk


ifeq ($(MAKECMDGOALS),release)
BUILD_TYPE=release
endif
ifeq ($(MAKECMDGOALS),debug)
BUILD_TYPE=debug
endif

# Define targets for release and debug 

release : .Starting  .XML_ParameterListLib .ParameterListTest .Finished

debug   : .Starting .XML_ParameterListLib_debug .ParameterListTest_debug .Finished  

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


.XML_ParameterListLib :
	###########################################################
	#       Creating XML_ParameterListLib 
	###########################################################
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=$(290J_Dir)/XML_ParameterList release 

.XML_ParameterListLib_debug :
	###########################################################
	#       Creating XML_ParameterListLib_debug
	###########################################################
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=$(290J_Dir)/XML_ParameterList debug

.ParameterListTest :
	###########################################################
	#            Creating ParameterListTest.exe
	###########################################################
	$(QUIET)$(MAKE) -f ParameterListTest.mk --directory=./ release 
.ParameterListTest_debug :
	###########################################################
	#            Creating ParameterListTest_debug.exe
	########################################################### 
	$(QUIET)$(MAKE) -f ParameterListTest.mk --directory=./ debug
	
.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  :
	rm -rf ./_*
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk   --directory=$(290J_Dir)/XML_ParameterList cleanall 
	$(QUIET)$(MAKE) -f ParameterListTest.mk     --directory=./ clean
clean: 
	$(QUIET)$(MAKE) -f ParameterListTest.mk   --directory=./ clean 

