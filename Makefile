#!/bin/make
#CC = icpc
CC = g++ 
OPTS =  -O2 -ffast-math 
#OPTS = -g
#OPTS = -O3
SRC = Main.cpp DicomReader.cpp EGSPhant.cpp ./DICOMParser/DICOMParser.cxx ./DICOMParser/DICOMFile.cxx ./DICOMParser/DICOMAppHelper.cxx
OBJS = $(SRC:.cpp=.o)
BIN = CTCombine
INC = 
LIBS = -lm
PROFILE = 

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(OPTS) -c $< -o $@ $(PROFILE)

$(BIN): $(OBJS)
	$(CC) $(OPTS) $(OBJS) -o $(BIN) $(LIBS) $(PROFILE)
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(BIN) *.o

