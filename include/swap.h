#ifndef SWAP_H
#define SWAP_H

bool bigEndian();
bool swapByteOrder(char *in, size_t N);
bool swapN(char *in, size_t N);
bool swap_float_array(float *x, size_t n);
bool swap_double_array(double *x, size_t n);
bool swap_int_array(int *x, size_t n); 

#endif
