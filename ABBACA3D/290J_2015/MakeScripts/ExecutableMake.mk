#
#############################################################################
#                    ExecutableMake.mk 
#!!!!!!!!!!!!!!!!!!!!!!! Do Not Edit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#            All variables required by this script are specfied 
#               before this script is included. 
#############################################################################
#
# A common source file for creating executables. 
#
# Variables that must be specified prior to including this
# file are
#
# SRC_DIR         : Directories containing source for the 
#                   library. All .c and .cpp files contained
#                   in these directories will be included in 
#                   the executable. May be left unspecified
#                   in which case only those files specified by 
#                   CPPfiles will be included. 
#
# CPPfiles        : Specification of individual C++ source files
#                   to be included in the executable. These files are
#                   added to those specified by SRC_DIR. 
#
# Cfiles          : Specification of individual C source files
#                   to be included in the executable. These files are
#                   added to those specified by SRC_DIR. 
#
# INCLUDES : list of include directories (in format -Idir1 -Idir2 ...)
#
# LIB_PATH : directories containing libraries (in format -Ldir1 -Ldir2 ...)
# LIBS     : list of libraries (in format -llib1 -llib2 ...)
#
# RELEASE_EXEC : Name of the release executable
# DEBUG_EXEC   : Name of the debug executable
#
# RELEASE_DIR  : location of release object (*.o) and dependency (*.d) files 
# DEBUG_DIR    : location of debug object (*.o) and dependency (*.d) files  
#                (Optional, default = ./_release and ./_debug)
#
# STATIC       : If defined then this will link in static versions of the
#                libraries. 
#
# Additional compiler parameters may be appended, e.g. for CXDEFINES,
# CXXDEFINES  += -D__NoGateGeo__ -D__NOGRAPHICS__
#
#############################################################################
#
# Copyright 2007-2015 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
# Verson: Wed Nov 25 19:56:12 2015 -0800
# Collection of files to be used in the executable
#
Cfiles   +=    $(foreach i,$(SRC_DIR),$(wildcard $(i)/*.c))
CPPfiles +=    $(foreach i,$(SRC_DIR),$(wildcard $(i)/*.cpp))

VPATH += $(dir $(CPPfiles))
VPATH += $(dir $(Cfiles))

LIB_PATH_NAMES = $(patsubst -L%,%,$(LIB_PATH))
VPATH += $(LIB_PATH_NAMES)

VPATH := $(sort $(VPATH))

ifndef RELEASE_DIR
RELEASE_DIR = ./_release
endif

ifndef DEBUG_DIR
DEBUG_DIR = ./_debug
endif

ifdef STATIC
LINK_ADDITIONS += -static -static-libgcc 
endif

ifndef VERBOSE
QUIET     :=@
MAKEFLAGS += --no-print-directory
MAKEFLAGS := $(sort $(MAKEFLAGS))
endif

ifeq ($(MAKECMDGOALS),release)
BUILD_EXEC          =$(RELEASE_EXEC)
OBJ_OUTPUT          =$(RELEASE_DIR)
EXISTING_OBJ_NAMES  =$(patsubst %.o,%,$(notdir $(wildcard $(RELEASE_DIR)/*.o)))
endif
ifeq ($(MAKECMDGOALS),debug)
BUILD_EXEC          =$(DEBUG_EXEC)
OBJ_OUTPUT          =$(DEBUG_DIR)
EXISTING_OBJ_NAMES  =$(patsubst %.o,%,$(notdir $(wildcard $(DEBUG_DIR)/*.o)))
endif

#
# This code creates the output directory for object and
# dependency files if it doesn't already exist
#
create_output_dir := \
$(shell if [ ! -e $(OBJ_OUTPUT) ]; then mkdir $(OBJ_OUTPUT); fi) 
#
# This following code removes .o and .d files from the 
# $(DEBUG_DIR) and $(RELEASE_DIR) directories if the associated 
# source files are not found. 
#
# If a file is removed then the executlable is  removed to
# force a recreation of the executable.  
#
# This code insures that if a source file is deleted from the 
# source list then it's corresponding object and dependency file are 
# deleted from the system, and hence the object file will
# not be included into the recreated executable. 
#
#
BUILD_NAMES    = $(notdir $(patsubst %.cpp,%,$(CPPfiles)))
BUILD_NAMES   += $(notdir $(patsubst %.c,%,$(Cfiles)))

ifneq ($(words $(EXISTING_OBJ_NAMES)),0)
del_orphaned_objects := \
$(shell for k in $(EXISTING_OBJ_NAMES);              \
    do                                               \
      result=0;                                      \
      for j in $(BUILD_NAMES);                       \
        do if [ $$j = $$k ]; then result=1; fi       \
      done;                                          \
      if [ $$result -eq 0 ]; then                    \
         rm --force "$(OBJ_OUTPUT)/$$k".o;           \
         rm --force "$(OBJ_OUTPUT)/$$k".d;           \
         rm --force  $(BUILD_EXEC); fi               \
    done)
endif
#
# This code creates the list of object files with a prefix of
# either $(RELEASE_DIR)/ or $(DEBUG_DIR)/. 
#
# The object files are be placed in subdirectories so that
# we can simultaneously work with debug and release
# versions. 
#
OBJECT_NAMES  = $(notdir $(patsubst %.c,%.o,$(Cfiles)))
OBJECT_NAMES += $(notdir $(patsubst %.cpp,%.o,$(CPPfiles)))

ReleaseObjects  = $(patsubst %.o,$(RELEASE_DIR)/%.o,$(OBJECT_NAMES))
DebugObjects    = $(patsubst %.o,$(DEBUG_DIR)/%.o,$(OBJECT_NAMES))

ifeq ($(MAKECMDGOALS),release)
Objects    = $(ReleaseObjects)
DFILES     = $(patsubst %.d,$(RELEASE_DIR)/%.d,$(notdir $(patsubst %.c,%.d,$(Cfiles))))
DFILES    += $(patsubst %.d,$(RELEASE_DIR)/%.d,$(notdir $(patsubst %.cpp,%.d,$(CPPfiles))))
endif
ifeq ($(MAKECMDGOALS),debug)
Objects     = $(DebugObjects)
DFILES      = $(patsubst %.d,$(DEBUG_DIR)/%.d,$(notdir $(patsubst %.c,%.d,$(Cfiles))))
DFILES     += $(patsubst %.d,$(DEBUG_DIR)/%.d,$(notdir $(patsubst %.cpp,%.d,$(CPPfiles))))
endif

release :  $(BUILD_EXEC)
debug   :  $(BUILD_EXEC)

EXEC_NAME = $(notdir $(BUILD_EXEC))
	
$(BUILD_EXEC): $(Objects) $(LIBS)
	$(QUIET)################################################
	# Linking executable $(EXEC_NAME)
	$(QUIET)################################################
	$(QUIET)$(CXX) $(LIB_PATH) $(Objects) $(LIBS) $(LINK_ADDITIONS) -o $(BUILD_EXEC) 
	$(QUIET)################################################
	$(QUIET)# Done building $(EXEC_NAME)
	$(QUIET)################################################
	
$(OBJ_OUTPUT)/%.o : %.cpp
	$(QUIET)################################################
	# Compiling $< 
	$(QUIET)################################################
	$(QUIET)$(call make-depend-cpp,$<,$@,$(subst .o,.d,$@))
	$(QUIET)$(CXX) -c $(CXXFLAGS) $(CXXDEFINES) \
$(INCLUDES) \
$< -o $@

$(OBJ_OUTPUT)/%.o : %.c
	$(QUIET)################################################
	# Compiling $< 
	$(QUIET)################################################
	$(QUIET)$(call make-depend-c,$<,$@,$(subst .o,.d,$@))
	$(QUIET)$(CC) -c $(CFLAGS) $(CDEFINES) \
$(INCLUDES) \
$< -o $@	
	
.PHONY : clean
clean  :
	rm -f $(RELEASE_EXEC)
	rm -f $(DEBUG_EXEC)
	rm -f $(RELEASE_DIR)/*.o
	rm -f $(DEBUG_DIR)/*.o
	rm -f $(RELEASE_DIR)/*.d
	rm -f $(DEBUG_DIR)/*.d
#
# This code creates the dependency files. This is a slight
# modification of the code given in the GNU make book, in
# that I've created a dependency list when the object files
# are contained in either the $(DEBUG_DIR) or $(RELEASE_DIR) 
# subdirectories. 
#
ifeq ($(MAKECMDGOALS),debug)
-include $(DFILES)
endif
ifeq ($(MAKECMDGOALS),release)
-include $(DFILES)
endif

#$(call make-depend-cpp,source-file,object-file,depend-file)
define make-depend-cpp
	$(CXX) -MM -MF $3 -MP -MT $2 $(CXXFLAGS) $(CXXDEFINES) $(INCLUDES)  $1
endef

#$(call make-depend-c,source-file,object-file,depend-file)
define make-depend-c
	$(CC) -MM -MF $3 -MP -MT $2 $(CFLAGS) $(CDEFINES) $(INCLUDES)  $1
endef




