#
# This is the "meta" makefile that builds the sample executables
# by invoking their respective makescripts.
#
#
# Usage :
#
# make release                   === builds release versions of all executables
# make debug                     === builds debug versions of all executables
#
# make .NameOfSamplesDir         === builds the release version of the executables in the samples directory specified
# make .NameOfSamplesDir_debug   === builds the debug version of the executables in the samples directory specified
#
# make clean                     === removes object files and executables
# make clanall                   === removes object files, executables and any supporting libraries 
#
#
# To build FFTW 3d samples that use the multi-threaded capabilities, add the specification
# _FFTW_OPENMP=1 to the make invocation, e.g. make .FFTWsinSolveTest _FFTW_OPENMP=1 
# 
#
#
# Nov. 12,2015
#
SHELL=/bin/sh

# Location of 290J source files and Makescripts files

290J_Dir=../290J_2015
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

# Define targets with dependencies by checking for specification of supporting libraries.  

# FFTW Targets built if FFTW_LIB is specified 

ifneq ("$(FFTW_LIB)","")
FFTW_TARGETS= .FFTW_PeriodicSolveTest .FFTW_SinSolveTest .FFTW_InterfaceTest .FFTW_DiscreteSinSolveTest 
FFTW_TARGETS_DEBUG= .FFTW_PeriodicSolveTest_debug .FFTW_SinSolveTest_debug .FFTW_InterfaceTest_debug .FFTW_DiscreteSinSolveTest_debug 
else
FFTW_TARGETS=
endif


ifneq ("$(SQLITE_LIB)","")
SQLITE_TARGETS= .SimDb_SQLiteTest  
SQLITE_TARGETS_DEBUG= .SimDb_SQLiteTest_debug 
else
SQLITE_TARGETS=
endif



release : .Starting  .XML_ParameterListLib .ParameterListTest .GridFunction3dTest .LaplaceOp1dTest .LaplaceOp2dTest .RayleighChebyshevTest \
.MollifiedNuclearPotTest .RitzMethodTest .LegendreInterpTest $(FFTW_TARGETS)  $(SQLITE_TARGETS)  .ExpandingDomainTest  .HermiteFunctionTest  .Finished 

debug   : .Starting .XML_ParameterListLib_debug .ParameterListTest_debug .GridFunction3dTest_debug \
.LaplaceOp1dTest_debug .LaplaceOp2dTest_debug .RayleighChebyshevTest_debug  .RitzMethodTest_debug .LegendreInterpTest_debug .MollifiedNuclearPotTest_debug   $(FFTW_TARGETS_DEBUG)\
$(SQLITE_TARGETS) .ExpandingDomainTest_debug  .HermiteFunctionTest_debug .Finished  

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
	#              .XML_ParameterListLib 
	###########################################################
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=$(290J_Dir)/XML_ParameterList release 

.XML_ParameterListLib_debug :
	###########################################################
	#             .XML_ParameterListLib_debug
	###########################################################
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=$(290J_Dir)/XML_ParameterList debug

.ParameterListTest :
	###########################################################
	#               .ParameterListTest
	###########################################################
	$(QUIET)$(MAKE) -f ParameterListTest.mk --directory=XML_ParameterListTest release 
.ParameterListTest_debug :
	###########################################################
	#               .ParameterListTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f ParameterListTest.mk --directory=XML_ParameterListTest debug


.GridFunction3dTest :
	###########################################################
	#               .GridFunction3dTest
	###########################################################
	$(QUIET)$(MAKE) -f GridFunction3dTest.mk --directory=GridFunction3dTest release 
.GridFunction3dTest_debug :
	###########################################################
	#               .GridFunction3dTest_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f GridFunction3dTest.mk --directory=GridFunction3dTest debug

.LaplaceOp1dTest :
	###########################################################
	#                 .LaplaceOp1dTest 
	###########################################################
	$(QUIET)$(MAKE) -f LaplaceOp1dTest.mk --directory=LaplaceOp1dTest release 
.LaplaceOp1dTest_debug :
	###########################################################
	#                .LaplaceOp1dTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f LaplaceOp1dTest.mk --directory=LaplaceOp1dTest debug


.LaplaceOp2dTest :
	###########################################################
	#                 .LaplaceOp2dTest 
	###########################################################
	$(QUIET)$(MAKE) -f LaplaceOp2dTest.mk --directory=LaplaceOp2dTest release 
.LaplaceOp2dTest_debug :
	###########################################################
	#                .LaplaceOp2dTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f LaplaceOp2dTest.mk --directory=LaplaceOp2dTest debug


.RayleighChebyshevTest :
	###########################################################
	#                 .RayleighChebyshevTest 
	###########################################################
	$(QUIET)$(MAKE) -f LaplaceOp2dEigTest.mk --directory=RayleighChebyshevTest release 
.RayleighChebyshevTest_debug :
	###########################################################
	#                 .RayleighChebyshevTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f LaplaceOp2dEigTest.mk --directory=RayleighChebyshevTest debug
	
.MollifiedNuclearPotTest :
	###########################################################
	#                 .MollifiedNuclearPotTest 
	###########################################################
	$(QUIET)$(MAKE) -f MollifiedNuclearPotTest.mk --directory=MollifiedNuclearPotTest release 
.MollifiedNuclearPotTest_debug :
	###########################################################
	#                 .MollifiedNuclearPotTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f MollifiedNuclearPotTest.mk --directory=MollifiedNuclearPotTest debug

.FFTW_SinSolveTest :
	###########################################################
	#                 .FFTW_SinSolveTest
	###########################################################
	$(QUIET)$(MAKE) -f Makefile  --directory=FFTW_SinSolveTest release 
.FFTW_SinSolveTest_debug :
	###########################################################
	#                 .FFTW_SinSolveTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_SinSolveTest debug

.FFTW_PeriodicSolveTest :
	###########################################################
	#                 .FFTW_PeriodicSolveTest
	###########################################################
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_PeriodicSolveTest release 

.FFTW_PeriodicSolveTest_debug :
	###########################################################
	#                 .FFTW_PeriodicSolveTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_PeriodicSolveTest debug


.FFTW_InterfaceTest :
	###########################################################
	#                 .FFTW_InterfaceTest 
	###########################################################
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_InterfaceTest release 

.FFTW_InterfaceTest_debug :
	###########################################################
	#                 .FFTW_InterfaceTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_InterfaceTest debug

	
.FFTW_DiscreteSinSolveTest :
	###########################################################
	#                 .FFTW_DiscreteSinSolveTest   
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_DiscreteSinSolveTest release

.FFTW_DiscreteSinSolveTest_debug : 
	###########################################################
	#                .FFTW_DiscreteSinSolveTest_debug  
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=FFTW_DiscreteSinSolveTest debug
	
.SimDb_SQLiteTest :
	###########################################################
	#                 .SimDb_SQLiteTest
	########################################################### 
	$(QUIET)$(MAKE) -f SqTest.mk --directory=SimDb_SQLiteTest release
	
.SimDb_SQLiteTest_debug : 
	###########################################################
	#                 .SimDb_SQLiteTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f SqTest.mk --directory=SimDb_SQLiteTest debug	


.RitzMethodTest :
	###########################################################
	#                 .RitzMethodTest
	########################################################### 
	$(QUIET)$(MAKE) -f RitzMethod2d.mk --directory=RitzMethodTest release
	
.RitzMethodTest_debug : 
	###########################################################
	#                 .RitzMethodTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f RitzMethod2d.mk --directory=RitzMethodTest debug	


.LegendreInterpTest :
	###########################################################
	#                 .LegendreInterpTest
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile  --directory=LegendreInterpTest release
	
.LegendreInterpTest_debug : 
	###########################################################
	#                 .LegendreInterpTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=LegendreInterpTest debug	


.ExpandingDomainTest :
	###########################################################
	#                 .ExpandingDomainTest
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile  --directory=ExpandingDomainTest release
	
.ExpandingDomainTest_debug : 
	###########################################################
	#                 .ExpandingDomainTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=ExpandingDomainTest debug	

	
.HermiteFunctionTest :
	###########################################################
	#                 .HermiteFunctionTest
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile  --directory=HermiteFunctionTest release
	
.HermiteFunctionTest_debug : 
	###########################################################
	#                 .HermiteFunctionTest_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=HermiteFunctionTest debug	

.H2test :
	###########################################################
	#                 .H2test
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile  --directory=HermiteFunctionTest release
	
.HermiteFunctionTest_debug : 
	###########################################################
	#                 .H2test_debug
	########################################################### 
	$(QUIET)$(MAKE) -f Makefile --directory=H2test debug	
	
.Finished :
	###########################################################
	#            Build Completed
	########################################################### 
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk       --directory=$(290J_Dir)/XML_ParameterList cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=LegendreInterpTest            cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=ExpandingDomainTest           cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=HermiteFunctionTest           cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_DiscreteSinSolveTest     cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_InterfaceTest            cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_SinSolveTest             cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_PeriodicSolveTest        cleanall
	$(QUIET)$(MAKE) -f Makefile                      --directory=H2test                        cleanall
clean: 
	$(QUIET)$(MAKE) -f ParameterListTest.mk          --directory=XML_ParameterListTest clean 
	$(QUIET)$(MAKE) -f GridFunction3dTest.mk         --directory=GridFunction3dTest clean
	$(QUIET)$(MAKE) -f LaplaceOp1dTest.mk            --directory=LaplaceOp1dTest clean
	$(QUIET)$(MAKE) -f LaplaceOp2dTest.mk            --directory=LaplaceOp2dTest clean
	$(QUIET)$(MAKE) -f LaplaceOp2dEigTest.mk         --directory=RayleighChebyshevTest clean
	$(QUIET)$(MAKE) -f MollifiedNuclearPotTest.mk    --directory=MollifiedNuclearPotTest clean
	$(QUIET)$(MAKE) -f SqTest.mk                     --directory=SimDb_SQLiteTest clean
	$(QUIET)$(MAKE) -f RitzMethod2d.mk               --directory=RitzMethodTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=LegendreInterpTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=ExpandingDomainTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=HermiteFunctionTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_DiscreteSinSolveTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_SinSolveTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_PeriodicSolveTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=FFTW_InterfaceTest clean
	$(QUIET)$(MAKE) -f Makefile                      --directory=H2test             clean


