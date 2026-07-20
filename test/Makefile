
CFLAGS = -O3  -Wall -Wextra -fstrict-overflow -fstack-protector-all

CC = g++
INC = -I../include
LIBS =  -L../lib -lbabak -lniftiio -lznz -lm -lz -lc -llapack -lblas

all: get_directory_name_test get_nifti_filename_test

# compilation of executables
get_directory_name_test: get_directory_name_test.cxx
	$(CC) $(CFLAGS) -o get_directory_name_test get_directory_name_test.cxx $(INC) $(LIBS)

get_nifti_filename_test: get_nifti_filename_test.cxx
	$(CC) $(CFLAGS) -o get_nifti_filename_test get_nifti_filename_test.cxx $(INC) $(LIBS)
