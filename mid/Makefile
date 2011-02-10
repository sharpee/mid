#
# Copyright (C) CNL-T, The Salk Institute for Biological Studies
# Author: Anirvan Nandy
#
# THE PROGRAMS ARE PROVIDED "AS IS" WITHOUT WARRANTY OF MERCANTABILITY
# OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
# IMPLIED. IN NO EVENT SHALL THE SALK INSTITUTE OR THE AUTHOR(S)
# BE LIABLE FOR ANY DIRECT OR CONSEQUENTIAL DAMAGES RESULTING
# FROM USE OF THE PROGRAMS.  THE USER BEARS THE ENTIRE RISK FOR USE OF
# THE PROGRAMS.
#

#
# M.I.D. Version
#
VERSION			=	3.0


#
# MXML definitions
#
MXML_DIR		=	./mxml-2.6
MXML_LIB_DIR	=	-L$(MXML_DIR)/lib
MXML_LIBS		=	-lmxml
MXML_INC_DIR	=	-I$(MXML_DIR)/include

#
# Compiler tool definitions
#
CC				=	g++ 
CFLAGS			=	-c -D_THREAD_SAFE -D_REENTRANT
LDFLAGS			=	-fstrict -O3 -D_THREAD_SAFE -D_REENTRANT
LIBS			=	$(MXML_LIBS) -lpthread
LIB_DIR			=	$(MXML_LIB_DIR)
INC_DIR			=	$(MXML_INC_DIR)

#
# Targets
#
TARG_1D			=	mid1d
TARG_ND			=	midnd
TARGETS			=	$(TARG_1D) $(TARG_ND)

#
# Source files
#
SRC_1D			=	$(TARG_1D).cpp
SRC_ND			=	$(TARG_ND).cpp
SRC_COMMON_CXX	=	ConfigFile.cpp
SRC_COMMON_C	=	nrutil.c ran2.c mnbrakN.c dbrentsaN.c \
					dlinminsaN.c f1dimN.c Proj1d.c myutil.c neginform.c \
					df1dimN.c readin.c derivative.c readinspikes.c \
					twoneginform.c choose.c trysa.c
SRCS			=	$(SRC_1D) $(SRC_ND) $(SRC_COMMON_CXX) $(SRC_COMMON_C)
HDRS            =   ConfigFile.h myutil.h nrutil.h patch.h prototypes.h

#
# Other files that go into the tar distribution
#
MAKE_FILE_NAME  =   Makefile
MISC_FILES		=	README.txt

#
# Distribution tar file
#
TAR_FILE_PREFIX =   MID-v$(VERSION)-$(shell date +%F)
TAR_FILE		=	$(TAR_FILE_PREFIX).tar
TAR_FILE_LIST	=	$(SRCS) $(HDRS) $(MISC_FILES) $(MAKE_FILE_NAME)

#
# Objects
#
OBJ_1D 			= 	$(TARG_1D).o
OBJ_ND 			= 	$(TARG_ND).o
OBJ_COMMON		=	$(SRC_COMMON_CXX:.cpp=.o) $(SRC_COMMON_C:.c=.o)
OBJS			=	$(OBJ_1D) $(OBJ_ND) $(OBJ_COMMON)

#
# Object dependency files
#
DEPS            =   $(OBJS:.o=.d)

#
# Make everything
#
all:            $(TARGETS)

#
# Target dependencies
#
$(TARG_1D):	$(OBJ_1D) $(OBJ_COMMON) $(MAKE_FILE_NAME)
	$(CC) $(LDFLAGS) -o $@ $(OBJ_1D) $(OBJ_COMMON) $(LIB_DIR) $(LIBS)

$(TARG_ND):	$(OBJ_ND) $(OBJ_COMMON) $(MAKE_FILE_NAME)
	$(CC) $(LDFLAGS) -o $@ $(OBJ_ND) $(OBJ_COMMON) $(LIB_DIR) $(LIBS)

#
# Automatically generate pre-requisites for source code
# (i.e. .c & .cpp files). This will generate a bunch of mini
# dependency files (*.d) corresponding to each object file.
# This saves one the trouble of explicitly listing the header
# files as dependencies for individual object files.
# Borrowed from the GNU Make Manual
# (http://www.gnu.org/software/make/manual/make.html)
#
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# read in the dependency make files
include $(DEPS)

#
# Implicit rules
#
.c.o: 
	$(CC) $(CFLAGS) $(INC_DIR) $< -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INC_DIR) $< -o $@

#
# Tar distribution
#
tar: 
	mkdir $(TAR_FILE_PREFIX)
	cp $(TAR_FILE_LIST) ./$(TAR_FILE_PREFIX)
	tar cvf $(TAR_FILE) ./$(TAR_FILE_PREFIX)
	rm -rf $(TAR_FILE_PREFIX)

targz: tar
	gzip -9 $(TAR_FILE)

#
# Cleanup
# 
clean:
	rm -f $(DEPS) $(OBJS)

purge:
	rm -f $(DEPS) $(OBJS) $(TARGETS) $(TAR_FILE)*

