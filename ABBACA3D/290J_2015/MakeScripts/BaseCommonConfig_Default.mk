# 
# A "basic" configuration file used by the MakeScript make files
# 
# This file should be copied to a file named BaseCommonConfig.mk 
# and then customized as needed for a particular operating system.
#
# 
# Fri 07 Aug 2015 02:20:29 PM PDT 
# 
# To obtain verbose compilation output, define VERBOSE on the makefile invocation
# 
#
# For OpenMP : The default is now to build the executable with OpenMP. 
# To disable OpenMP specify OpenMP=0 on the command line invocation of the makefile
#
# Echo use of configuration default configuration file 

ifeq ("$(wildcard ./BaseCommonConfig.mk)","")
$(info # Using BaseCommonConfig_Default.mk)
endif

LIB_DIR      =  ./lib
BIN_DIR      =  ./bin

# Extracts the host name to be used for machine dependent make statements

HOSTNAME := $(shell hostname)

#
# This code sets the make variable _USE_EXTERN_LAPACK_ when lapack=1 is specified
# on the command line. The variable _USE_EXTERN_LAPACK_ defined to be 1 is then
# available to all makefiles that include the common configuration 
#

ifdef lapack
_USE_EXTERN_LAPACK_ = 1
endif

ifndef VERBOSE
QUIET     :=@
MAKEFLAGS += --no-print-directory
MAKEFLAGS := $(sort $(MAKEFLAGS))
endif

# Compiler specifications

CC   := gcc
CXX  := g++
AR   := ar

ifeq ($(COMPILER),MinGW)
AR  :=C:/MinGW/bin/ar
CC  :=C:/MinGW/bin/gcc
CXX :=C:/MinGW/bin/g++
endif

#
# Compiler options: Specify C and C++ flags and defines here (except openMP) flags
# 

ifeq ($(MAKECMDGOALS),release)
CFLAGS         +=-O2  -fno-gcse -fno-optimize-sibling-calls -Wno-write-strings
CXXFLAGS       +=-O2  -fno-gcse -fno-optimize-sibling-calls -Wno-write-strings -std=c++11
CDEFINES       += 
CXXDEFINES     +=  
endif


ifeq ($(MAKECMDGOALS),debug)
CFLAGS      +=-g -Wall -fno-inline
CDEFINES    +=-D_DEBUG
CXXFLAGS    +=-g -Wall -fno-inline -std=c++11
CXXDEFINES  +=-D_DEBUG 
endif


# For OpenMP : The default is now to build the executable with OpenMP. 
# To disable OpenMP specify OpenMP=0 on the command line invocation of the makefile

ifeq ($(OpenMP),0)
else
CXXFLAGS       += -fopenmp
LINK_ADDITIONS += -fopenmp  
endif

# ====== fftw3 library spcification ======== 

# Modify appropriately to specify the name and 
# location of the FFTW3 libraries. 
#
# The current default is the location for a Linux machine
# without OpenMP. 

FFTW_LIB= -lfftw3  
FFTW_PATH= 

# For Linux machines with OpenMP you may need to use 
# FFTW_LIB= -lfftw3_omp -lfftw3 -lm

# ====== sqlite3 library spcification ======== 

# The default specification here refers to a Linux installation
# with sqlite3 and libsqlite3-dev packages installed. 
#
# The source for sqlite3 is contained in a subdirectory of
# SIMdb_SQLite, so one can build the library using 
# the SIMdb_SQLite/sqlite3lib.mk makescript, or by downloading
# the source and build tools from https://www.sqlite.org/ 
# and then set the path and libraries to link appropriately.
#
# The UCLA Math department linux systems do not have the sqlite3 
# header files installed, so you'll have to create the library
# yourself if you are using Joshua. 
 

SQLITE_PATH=  -L/usr/include
SQLITE_LIB=   -lsqlite3  

#On some Linux machines you may need to add -ldl, e.g. 
#SQLITE_LIB=  -lsqlite3 -ldl
#
#On the UCLA Math department machines you need to add -ldl 

# For static builds 

ifdef STATIC 
LINK_ADDITIONS += -static 
endif 

# Revisions:
# (11/25.2015):
# Added specification of SQLite libraries 
#
# (11/24/2015): 
# Modified specification of FFTW libraries, so they are not
# included in every build that uses BaseCommonConfig
#
# (11/18/2015): 
# Added -std=c++11 to CXXFLAGS   
# Added specification of fftw3 library: Modify appropriately
# for different fftw3 installations, or comment out to avoid
# their inclusion. 

