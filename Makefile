#!/bin/make

CC = g++
#CC = icpc 
OPTS =  -O2 -ffast-math 
#OPTS = -g
#OPTS = -O3
SRC = Main.cpp EGSPhant.cpp DicomReader.cpp ./DICOMParser/DICOMParser.cxx ./DICOMParser/DICOMFile.cxx ./DICOMParser/DICOMAppHelper.cxx
OBJS = $(SRC:.cpp=.o)
INC = 
LIBS = 
PROFILE = 

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(OPTS) -c $< -o $@ $(PROFILE)

CTCombine: $(OBJS)
	$(CC) $(OPTS) $(OBJS) -lm -o CTCombine $(LIBS) $(PROFILE)
	@echo Binary created!!

clean:
	set nonomatch; rm -f CTCombine *.o

