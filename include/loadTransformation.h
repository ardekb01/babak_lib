#ifndef LOADTRANSFORMATION_H
#define LOADTRANSFORMATION_H

#define LOADTRANSFORM_OK 0
#define LOADTRANSFORM_NULL_POINTER 1
#define LOADTRANSFORM_FILE_ERROR 2
#define LOADTRANSFORM_PARSE_ERROR 3

int loadTransformation(const char *filename, float *T);

#endif
