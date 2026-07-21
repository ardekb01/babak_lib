#ifndef GETNIFTIIMAGEORIENTATOIN_H
#define GETNIFTIIMAGEORIENTATOIN_H

// orientation must point to a buffer of at least 4 bytes.

bool getNiftiImageOrientation(const char *filename,
                              char *orientation);

bool getNiftiImageOrientation(nifti_1_header hdr,
                              char *orientation);

#endif
