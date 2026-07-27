#ifndef READ_NIFTI_IMAGE_H
#define READ_NIFTI_IMAGE_H

char *read_nifti_image(const char *filename, nifti_1_header *hdr);
bool read_nifti_hdr(const char *filename, nifti_1_header *hdr);

#endif
