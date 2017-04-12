
# This makefile builds the sqlite3lib support library. 
#

SHELL=/bin/sh

# Parameters for library construction script 

CXX        := g++
LIB_DIR    = ./lib
Cfiles     = ./sqlite/sqlite3.c
INCLUDES  += -I./

# Create library directory if it doesn't exist 

create_lib_dir := \
$(shell if [ ! -e $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi) 

# specify make parameters

ifndef VERBOSE
QUIET     :=@
MAKEFLAGS += --no-print-directory
MAKEFLAGS := $(sort $(MAKEFLAGS))
endif

ifeq ($(MAKECMDGOALS),release)

CFLAGS      :=-O2  
CXXFLAGS    :=-O2 
CDEFINES    += 
CXXDEFINES  += 
endif

ifeq ($(MAKECMDGOALS),debug)
CFLAGS      =-g -Wall 
CDEFINES    =-D_DEBUG
CXXFLAGS    =-g -Wall 
CXXDEFINES  +=-D_DEBUG 
endif

RELEASE_DIR  = ./_releaseSQLite3
DEBUG_DIR    = ./_debugSQLite3

DEBUG_LIBRARY    = $(LIB_DIR)/libsqlite3_debug.a
RELEASE_LIBRARY  = $(LIB_DIR)/libsqlite3.a

include ./ComponentMake.mk

