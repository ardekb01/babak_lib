#ifndef MSPTRANSFORATION_H 
#define MSPTRANSFORATION_H 

bool MSPtransformation(const char *filename, const char *orient, float4 *Tmsp, DIM &dim);
bool MSPtransformation(const char *filename, const char *orient, const char *lmfile, float4 *Tmsp, DIM &dim);

#endif
