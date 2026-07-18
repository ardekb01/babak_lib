
CFLAGS = -O3  -Wall -Wextra -fstrict-overflow -fstack-protector-all

CC = g++
INC = -I../include
LIBS =  -L../lib -lbabak -lniftiio -lznz -lm -lz -lc -llapack -lblas

all: getDirectoryName_test

# compilation of executables
getDirectoryName_test: getDirectoryName_test.cxx
	$(CC) $(CFLAGS) -o getDirectoryName_test getDirectoryName_test.cxx $(INC) $(LIBS)
