# Boost base path
BOOST=/usr/local/include/boost
# Include flags
#INC=-I ${BOOST}/include -I ../include
INC=-I ../include
# Library path flags
LIBS=
#${BOOST}/lib/libboost_timer.a ${BOOST}/lib/libboost_system.a ${BOOST}/lib/libboost_chrono.a
# Remove file command
RM=rm -f
# Compiler
#CC=/Library/Developer/CommandLineTools/usr/bin/c++
CC=g++
# Compiler flags
CCFLAGS=-O3 -std=c++11 $(INC)
#-Wno-deprecated -Wall -W -Wextra -fpic -O3 -static $(INC)
# Source files
SRCS=utils.cpp variable.cpp domain.cpp factor.cpp trieset.cpp setfactor.cpp valuation.cpp operations.cpp io.cpp graph.cpp inference.cpp limid.cpp 
# Object files
OBJS=$(subst .cpp,.o,$(SRCS))
# Binary files
BIN=../bin/
BSRC=solve_limid.cpp run_ve.cpp
BINS=$(subst .cpp,,$(BSRC))
BOBJ=$(subst .cpp,.o,$(BSRC))

all: ${BINS}

%.o: %.cpp
	$(CC) $(CCFLAGS) -c -o $@ $<

run_ve: run_ve.o $(OBJS)
	$(CC) $(CCFLAGS) -o ${BIN}$@ $^ $(LIBS)

solve_limid: solve_limid.o $(OBJS)
	$(CC) $(CCFLAGS) -o ${BIN}solve_limid $^ $(LIBS)


# CLEAN

.PHONY: clean

clean:
	$(RM) $(OBJS) ${BOBJ}

dist-clean: clean
	$(RM) ${BIN}${BINS} ${OBJS} ${BOBJ}
