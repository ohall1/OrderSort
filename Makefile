# **********************************
# Makefile for AIDAsort-v2cg program
#      C. Griffin 04/10/16
#   Adapted to compile on ubuntu
#       O. Hall 27/03/27
# **********************************

.KEEP_STATE:

OS=Linux64

# Compiler.
CC          	= g++
# Flags for compiler.
CFLAGS		=  -c -O3 -Wall -Wextra -std=c++0x -DUNIX -DLINUX -DLINUX64 -DPOSIX #-DDEBUG_UNP -DDEBUG_CAL -DDEBUG_ANA
#If getting seg faults and can't work out why, use the below CFLAGS for extra debug options and help
#CFLAGS		= -g -Q -v -c -O3 -Wall -Wextra -DUNIX -DLINUX -DLINUX64 -DPOSIX -std=c++0x
# Flags for linker.
LDFLAGS		=-lrt

ROOTCONFIG   := $(shell which root-config)
ROOTCINT     := $(shell which rootcint)

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --glibs)
LIBS         := $(ROOTLIBS)
HOST      = $(shell hostname)

# Add root flags.
CFLAGS 		+= $(ROOTCFLAGS)
# Add root.
LDFLAGS 	+= $(ROOTLDFLAGS) $(ROOTLIBS)

#INCLUDES=  -I/MIDAS/DataPackage/DataXferLib/V4_TCP -I/MIDAS/DataPackage/DataSpyLib #!!! on aidas1 !!!
#LDLIBS=   -L/usr/ucblib -lrt -lpthread -L/MIDAS/Linux/lib64 -lxfer -ldataspy #!!! on aidas1 !!!

ifeq '$(HOST)' 'vorbis.ph.ed.ac.uk'
	INCLUDES=  -I/Disk/ds-sopa-group/np/RIKEN/AIDAsort/DataXferLib/V4_TCP -I/Disk/ds-sopa-group/np/RIKEN/AIDAsort/DataSpyLib
	LDLIBS=   -L/usr/ucblib -lrt -lpthread -L/home/s1668713/Documents/Code/AIDALib/MIDAS/Linux/lib64 -lxfer -ldataspy
	COMP= $($(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^)
else
	INCLUDES=  -I/mnt/c/Users/oscar/Linux/Code/AIDA/AIDALib/DataXferLib/V4_TCP -I/mnt/c/Users/oscar/Linux/Code/AIDA/AIDALib/DataSpyLib
	LDLIBS=   -L/usr/ucblib -lrt -lpthread -L/mnt/c/Users/oscar/Linux/Code/AIDA/AIDALib/MIDAS/Linux/lib64 -lxfer -ldataspy
	COMP= $($(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS))
endif


# The object files.
OBJECTS =   DataSource.o Unpacker.o Calibrator.o Analysis.o

# Specify the file to search for .cc/.h files.
vpath %.cpp ./src
vpath %.h ./include

AIDAsort-v2cg.exe : main.o $(OBJECTS)
	$(COMP)
main.o : main.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -Iinclude/ $^

DataSource.o : DataSource.h DataSource.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -Iinclude/ $^

Unpacker.o : Unpacker.h Unpacker.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -Iinclude/ $^

Calibrator.o : Calibrator.h Calibrator.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -Iinclude/ $^

Analysis.o : Analysis.h Analysis.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -Iinclude/ $^

clean:
	rm -vf *.exe *.o *~ include/*.gch
