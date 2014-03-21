#!/bin/make
# My first Makefile

CC = /usr/bin/g++
#CC = icpc
OPTS = -O2 -ffast-math 
SRC = Main.cpp 
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

