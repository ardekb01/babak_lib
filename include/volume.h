///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Copyright (C) 2024 Babak A. Ardekani, PhD - All Rights Reserved.  //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// This is the HEADER FILE volume.h. This is The INTERFACE for the class
// VOLUME, which is an ADT for 3D medical images.

#ifndef VOLUME_H 
#define VOLUME_H

class VOLUME
{

public:
	VOLUME();	// constructor
	~VOLUME();	// destructor

	// reads a Siemens MOSAIC image
	// nz is the number of slices in the MOSAIC image      
	void readmosaic(char *pathname);

	// Input - full path of an image in Siemens Vision format
	// Behaviour - Examines all the files in the same directory as the 
	// input image.  
	// Selects the readable files that have the same "patient id",
	// "patient name", "study number", TE and size as the input image.
	// The images in the files that pass all tests are read into the 
	// data member of the VOLUME class along with matrix and voxels 
	// sizes and other information. 
	// The image slices are arranged according to their "image number".
	void readsiemensvision(char *pathname);

	int readsiemensvision(int N, char **filename);

	int read_analyze(const char *pathname);

	void read_ge(char *pathname);
	void read_gelx(char *pathname);
	void read_dicom(char *pathname);

	void read_image(char *pathname);
	void read_smis_header(char *pathname);
	void read_smis_image(char *pathname);
	void read_vision_header(char *pathname);
	void read_vision_image(char *pathname);

	// saves the image volume in CABI format
	void savenki(char *pathname);

	int nx;	// number of image columns
	int ny; // number of image rows
	int nz; // number of slices
	int np; // nx*ny
	int nv; // nx*ny*nz
	int nt; // number of frames
	int image_number;	// image number in multi-slice data

	double dx; // voxel size row (distance between columns) in mm
	double dy; // voxel size column (distance between rows) in mm
	double dz; // voxel size slice (distance between slices) in mm
	double dt; // time between frames (seconds)

	short min;
	short max;

	// position vector pointing to the center of the first slice
	double centervec[3]; 

	double normalvec[3]; // unit normal vector

	// unit vector pointing in the row (left to right) direction
	double rowvec[3]; 

	// unit vector pointing in the column (top to bottom) direction
	double columnvec[3];

	short *data;
	unsigned char *uc_data;

private:
	int check_F_R_permission(char *file);
	void read_analyze_hdr(struct dsr *hdr, char *filename);
	void get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img);
	void get_subimage(short *im1, int nx1, int ny1, int x0, int y0, short *im2);
	void minmax();
	void scale_image();
	int isSMIS(char *pathname);
	int isVision(char *pathname);
};

#endif // VOLUME_H
