#ifndef GETNIFTIIMAGEORIENTATOIN_H
#define GETNIFTIIMAGEORIENTATOIN_H

bool getNiftiImageOrientation(const char *filename, char *orientation);

bool getNiftiImageOrientation(nifti_1_header hdr,
                              char *orientation);

#endif
