#ifndef SETLOWHIGH_H
#define SETLOWHIGH_H

void setLowHigh(const short *image, int nv, int &low, int &high, float percent);

// image is modified in-place
void setMX(short *image, short *msk, int nv, int &high, float percent);

#endif
