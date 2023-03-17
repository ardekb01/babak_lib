#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <spm_analyze.h>
#include <babak_lib.h>
#include <minmax.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int extension_is_hdr(const char *filename)
{
   int L; // number of characters in the filename

   L=strlen(filename);

   if(L<=4)
   {
      return(0);
   }

   if( filename[L-1]=='r' && filename[L-2]=='d' && filename[L-3]=='h' && filename[L-4]=='.')
   {
      return(1);
   }

   return(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function was last updated on Sep. 20, 2002.
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img)
{
	int n;

	/* find the number of character in the filename */
	n=(int)strlen(filename);

	/* make out the file names */
	if(n<4) 
	{
		/* filename is less than 4 characters. Therefore, it
		cannot have .hdr or .img extensions. It must be the
		base name alone. */   
      
		sprintf(basename_hdr,"%s.hdr",filename);
		sprintf(basename_img,"%s.img",filename);

		return;
	}

	if( filename[n-1]=='r' && filename[n-2]=='d' && filename[n-3]=='h' && filename[n-4]=='.') 
	{
		/* filename is the .hdr file */

		sprintf(basename_hdr,"%s",filename);
		sprintf(basename_img,"%s",filename);
		basename_img[n-3]='i'; basename_img[n-2]='m'; basename_img[n-1]='g';

		return;
	}

	if( filename[n-1]=='g' && filename[n-2]=='m' && filename[n-3]=='i' && filename[n-4]=='.') 
	{
		/* filename is the .img file */

		sprintf(basename_img,"%s",filename);
		sprintf(basename_hdr,"%s",filename);
		basename_hdr[n-3]='h'; basename_hdr[n-2]='d'; basename_hdr[n-1]='r';

		return;
     }

	// filename is the basename
	sprintf(basename_hdr,"%s.hdr",filename);
	sprintf(basename_img,"%s.img",filename);

	return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_analyze_hdr(struct dsr *hdr, char *filename)
{
   FILE *fp;
   int n;

   fp=fopen(filename,"r");
   if(fp==NULL) file_open_error(filename);

   n=(int)fread(hdr,sizeof(struct dsr),1,fp);

   if(n!=1) 
   {
      printf("error: reading from file %s !\n",filename);
      fclose(fp);
      exit(0);
   }

   fclose(fp);
}

void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType)
{
	if ( hdr.hk.sizeof_hdr != 348 )
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );

	*dataType = hdr.dime.datatype;

	printf("\n\tData type = %d",*dataType);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);

		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	*ny=hdr.dime.dim[2];
	*nx=hdr.dime.dim[1];
	*nz=hdr.dime.dim[3];

	printf("\n\tMatrix size = %d x %d x %d (voxels)", *nx, *ny, *nz);

	*dx=hdr.dime.pixdim[1];
	*dy=hdr.dime.pixdim[2];
	*dz=hdr.dime.pixdim[3];

	printf("\n\tVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", *dx,*dy,*dz);
}

void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz, short *dataType)
{
	if ( hdr.hk.sizeof_hdr != 348 )
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );

	*dataType = hdr.dime.datatype;

	printf("\n\tData type = %d",*dataType);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);

		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	*ny=hdr.dime.dim[2];
	*nx=hdr.dime.dim[1];
	*nz=hdr.dime.dim[3];

	printf("\n\tMatrix size = %d x %d x %d (voxels)", *nx, *ny, *nz);

	*dx=hdr.dime.pixdim[1];
	*dy=hdr.dime.pixdim[2];
	*dz=hdr.dime.pixdim[3];

	printf("\n\tVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", *dx,*dy,*dz);
}

void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType, int v)
{
	if ( hdr.hk.sizeof_hdr != 348 )
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );

	*dataType = hdr.dime.datatype;

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);

		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	*ny=hdr.dime.dim[2];
	*nx=hdr.dime.dim[1];
	*nz=hdr.dime.dim[3];

	*dx=hdr.dime.pixdim[1];
	*dy=hdr.dime.pixdim[2];
	*dz=hdr.dime.pixdim[3];

	if(v) 
	{
		printf("\n\tMatrix size = %d x %d x %d (voxels)", *nx, *ny, *nz);
		printf("\n\tVoxel size = %8.6f x %8.6f x %8.6f (mm3)", *dx,*dy,*dz);
		printf("\n\tData type = %d\n",*dataType);
	}

	return;
}

void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, short *dataType, int v)
{
	if ( hdr.hk.sizeof_hdr != 348 )
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );

	*dataType = hdr.dime.datatype;

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);

		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	*ny=hdr.dime.dim[2];
	*nx=hdr.dime.dim[1];
	*nz=hdr.dime.dim[3];
	*nt=hdr.dime.dim[4];

	if(v) printf("\n\tMatrix size = %d x %d x %d x %d (voxels)", *nx, *ny, *nz, *nt);

	*dx=hdr.dime.pixdim[1];
	*dy=hdr.dime.pixdim[2];
	*dz=hdr.dime.pixdim[3];

	if(v) printf("\n\tVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", *dx,*dy,*dz);
}

char *read_image(char *file, int n)
{
	FILE *fp;
	char *im;
	int nread;

	im =(char *)calloc(n,1);

	if(im==NULL) return(NULL);
	
	fp = fopen(file,"r");
    if(fp==NULL) file_open_error(file);
	nread=(int)fread(im,1,n,fp);
	fclose(fp);

	if(nread!=n) { free(im); return(NULL);}

	return(im);
}

void read_image(char *file, char *im, int n)
{
   FILE *fp;

   if(im==NULL) return;
	
   fp = fopen(file,"r");
   if(fp==NULL) file_open_error(file);
   if (fread(im,1,n,fp)==(size_t)n ) {};
   fclose(fp);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
void read_analyze_image(const char *filename, short *im)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	int nx,ny,nz;
	float dx,dy,dz;
	short dataType;

	struct dsr analyzehdr;

	if(im==NULL) return;

	get_analyze_file_names(filename, hdrfile, imgfile);

    if( !check_F_R_permission(hdrfile) )
    {
        printf("\nError: cannot read %s\n",hdrfile);
		return;
    }

    if( !check_F_R_permission(imgfile) )
    {
        printf("\nError: cannot read %s\n",imgfile);
		return;
    }

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, &nx, &ny, &nz, &dx, &dy, &dz, &dataType, 0);

	nv = nx*ny*nz;

	if(dataType!=4) return;

	read_image(imgfile, (char *)im, nv*2);

	if ( analyzehdr.hk.sizeof_hdr != 348 )
		swapN( (char *)im, 2*nv);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	short dataType;
	char *im;

	struct dsr analyzehdr;

	get_analyze_file_names(filename, hdrfile, imgfile);

    if( !check_F_R_permission(hdrfile) )
    {
        printf("\nError: cannot read %s\n",hdrfile);
		return(NULL);
    }

    if( !check_F_R_permission(imgfile) )
    {
        printf("\nError: cannot read %s\n",imgfile);
		return(NULL);
    }

	printf("\nReading:");
	printf("\n\t%s",hdrfile);
	printf("\n\t%s",imgfile);

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, nx, ny, nz, dx, dy, dz, &dataType);

	nv = (*nx)*(*ny)*(*nz);

	if(dataType==4)
	{
		im= read_image(imgfile,nv*2);

		if(im==NULL) return(NULL);

		if ( analyzehdr.hk.sizeof_hdr != 348 )
			swapN( (char *)im, 2*nv);

		return(im);
	}

	if(dataType==2)
	{
		im = read_image(imgfile,nv);
		return(im);
	}

	printf("\n\nSorry, but this program cannot handle this data type! Aborting ...\n\n");
	return(NULL);
}

void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int datatype, float dx, float dy, float dz)
{
	/* required elements of the header_key substructure */
	hdr->hk.sizeof_hdr=348;
	hdr->hk.extents=16384;
  	hdr->hk.regular='r';
	hdr->hk.hkey_un0='0';

	hdr->dime.dim[0]=4;
	hdr->dime.dim[1]=nx;
	hdr->dime.dim[2]=ny;
	hdr->dime.dim[3]=nz;
  	hdr->dime.dim[4]=1;
  	hdr->dime.dim[5]=0;
  	hdr->dime.dim[6]=0;
  	hdr->dime.dim[7]=0;

	strcpy(hdr->dime.vox_units,"mm");
	strcpy(hdr->dime.cal_units,""); 
	hdr->dime.unused1=0; 

  	hdr->dime.datatype=datatype;

	switch(hdr->dime.datatype) {
		case 2: 
		{
			hdr->dime.bitpix=sizeof(unsigned char)*8;
			break;
		}
		case 4: 
		{
			hdr->dime.bitpix=sizeof(short)*8;
			break;
		}
		case 8: 
		{
			hdr->dime.bitpix=sizeof(int)*8;
			break;
		}
		case 16: 
		{
			hdr->dime.bitpix=sizeof(float)*8;
			break;
		}
  	}

	hdr->dime.pixdim[0]=0.0;
	hdr->dime.pixdim[1]=dx;
	hdr->dime.pixdim[2]=dy;
	hdr->dime.pixdim[3]=dz;
	hdr->dime.pixdim[4]=0.0;
	hdr->dime.pixdim[5]=0.0;
	hdr->dime.pixdim[6]=0.0;
	hdr->dime.pixdim[7]=0.0;

	hdr->dime.vox_offset=0.0; // The images could not be displayed in SPM
                                  // without this.

	hdr->dime.roi_scale=1.0;        // SPM SCALE variable
	hdr->dime.funused1=0.0;
	hdr->dime.funused2=0.0;
	
	hdr->dime.cal_max = 0.0;
	hdr->dime.cal_min = 0.0;

	hdr->dime.compressed = 0;
	hdr->dime.verified = 0;

	strcpy(hdr->hist.aux_file,"none");
	hdr->hist.orient='\0';

	for(int i=0; i<10; i++)
		hdr->hist.originator[i]='\0';
	
	strcpy(hdr->hist.generated,"");
	strcpy(hdr->hist.scannum,"");
	strcpy(hdr->hist.patient_id,"");
	strcpy(hdr->hist.exp_date,"");
	strcpy(hdr->hist.exp_time,"");
	strcpy(hdr->hist.hist_un0,"");
	hdr->hist.views=0;
	hdr->hist.vols_added=0;
	hdr->hist.start_field=0;
	hdr->hist.field_skip=0;
	hdr->hist.omax=0;
	hdr->hist.omin=0;
	hdr->hist.smax=0;
	hdr->hist.smin=0;
}

void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int nt, int datatype, float dx, float dy, float dz)
{
	/* required elements of the header_key substructure */
	hdr->hk.sizeof_hdr=348;
	hdr->hk.extents=16384;
  	hdr->hk.regular='r';
	hdr->hk.hkey_un0='0';

	hdr->dime.dim[0]=4;
	hdr->dime.dim[1]=nx;
	hdr->dime.dim[2]=ny;
	hdr->dime.dim[3]=nz;
  	hdr->dime.dim[4]=nt;
  	hdr->dime.dim[5]=0;
  	hdr->dime.dim[6]=0;
  	hdr->dime.dim[7]=0;

	strcpy(hdr->dime.vox_units,"mm");
	strcpy(hdr->dime.cal_units,""); 
	hdr->dime.unused1=0; 

  	hdr->dime.datatype=datatype;

	switch(hdr->dime.datatype) {
		case 2: 
		{
			hdr->dime.bitpix=sizeof(unsigned char)*8;
			break;
		}
		case 4: 
		{
			hdr->dime.bitpix=sizeof(short)*8;
			break;
		}
		case 8: 
		{
			hdr->dime.bitpix=sizeof(int)*8;
			break;
		}
		case 16: 
		{
			hdr->dime.bitpix=sizeof(float)*8;
			break;
		}
  	}

	hdr->dime.pixdim[0]=0.0;
	hdr->dime.pixdim[1]=dx;
	hdr->dime.pixdim[2]=dy;
	hdr->dime.pixdim[3]=dz;
	hdr->dime.pixdim[4]=0.0;
	hdr->dime.pixdim[5]=0.0;
	hdr->dime.pixdim[6]=0.0;
	hdr->dime.pixdim[7]=0.0;

	hdr->dime.vox_offset=0.0; // The images could not be displayed in SPM
                                  // without this.

	hdr->dime.roi_scale=1.0;        // SPM SCALE variable
	hdr->dime.funused1=0.0;
	hdr->dime.funused2=0.0;
	
	hdr->dime.cal_max = 0.0;
	hdr->dime.cal_min = 0.0;

	hdr->dime.compressed = 0;
	hdr->dime.verified = 0;

	strcpy(hdr->hist.aux_file,"none");
	hdr->hist.orient='\0';

	for(int i=0; i<10; i++)
		hdr->hist.originator[i]='\0';
	
	strcpy(hdr->hist.generated,"");
	strcpy(hdr->hist.scannum,"");
	strcpy(hdr->hist.patient_id,"");
	strcpy(hdr->hist.exp_date,"");
	strcpy(hdr->hist.exp_time,"");
	strcpy(hdr->hist.hist_un0,"");
	hdr->hist.views=0;
	hdr->hist.vols_added=0;
	hdr->hist.start_field=0;
	hdr->hist.field_skip=0;
	hdr->hist.omax=0;
	hdr->hist.omin=0;
	hdr->hist.smax=0;
	hdr->hist.smin=0;
}

void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	short min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 4, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=min;
	hdr.dime.glmax=max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(short),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz, int v)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	short min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if( v ) printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 4, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=min;
	hdr.dime.glmax=max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(short),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	float min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 16, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=(int)min;
	hdr.dime.glmax=(int)max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(float),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz, int v)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	float min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if(v) printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 16, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=(int)min;
	hdr.dime.glmax=(int)max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(float),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

void write_analyze_image(const char *filename, unsigned char *im, int nx, int ny, int nz, float dx, float dy, float dz)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	unsigned char min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 2, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=min;
	hdr.dime.glmax=max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(unsigned char),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

void write_analyze_image(const char *filename, unsigned char *im, int nx, int ny, int nz, float dx, float dy, float dz, int v)
{
	FILE *fp;
	char hdrfile[1024];
	char imgfile[1024];
	unsigned char min=0, max=0;
	int nv;

	struct dsr hdr;

	nv = nx*ny*nz;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if(v) printf("\nWriting %s ...\n",imgfile);

	create_analyze_hdr(&hdr, nx, ny, nz, 2, dx, dy, dz);

	minmax(im, nv, min, max);
	hdr.dime.glmin=min;
	hdr.dime.glmax=max;

   fp=fopen(imgfile,"w");
   if(fp==NULL) file_open_error(imgfile);
   fwrite(im,sizeof(unsigned char),nv,fp);
   fclose(fp);

   fp=fopen(hdrfile,"w");
   if(fp==NULL) file_open_error(hdrfile);
   fwrite(&hdr,sizeof(struct dsr),1,fp);
   fclose(fp);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	short dataType;
	char *im;

	struct dsr analyzehdr;

	get_analyze_file_names(filename, hdrfile, imgfile);

    if( !check_F_R_permission(hdrfile) )
    {
        printf("\nError: cannot read %s\n",hdrfile);
		return(NULL);
    }

    if( !check_F_R_permission(imgfile) )
    {
        printf("\nError: cannot read %s\n",imgfile);
		return(NULL);
    }

	printf("\nReading:");
	printf("\n\t%s",hdrfile);
	printf("\n\t%s",imgfile);

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, nx, ny, nz, dx, dy, dz, &dataType);

	nv = (*nx)*(*ny)*(*nz);

	*type = dataType;

	if(dataType==16)
	{
		im = read_image(imgfile,(int)(nv*sizeof(float)));

		if(im==NULL) return(NULL);

		if ( analyzehdr.hk.sizeof_hdr != 348 )
		for(int i=0; i<nv; i+=sizeof(float) )
			swapByteOrder( im+i, sizeof(float));

		return(im);
	}

	if(dataType==4)
	{
		im= read_image(imgfile,nv*2);

		if(im==NULL) return(NULL);

		if ( analyzehdr.hk.sizeof_hdr != 348 ) swapN( (char *)im, 2*nv);

		return(im);
	}

	if(dataType==2)
	{
		im = read_image(imgfile,nv);
		return(im);
	}

	printf("\n\nSorry, but this program cannot handle this data type! Aborting ...\n\n");
	return(NULL);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type, int v)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	short dataType;
	char *im=NULL;

	struct dsr analyzehdr;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if( !check_F_R_permission(hdrfile) )
	{
		printf("\nError: cannot read %s\n",hdrfile);
		return(NULL);
	}

	if( !check_F_R_permission(imgfile) )
	{
		printf("\nError: cannot read %s\n",imgfile);
		return(NULL);
	}

	if(v)
	{
		printf("\nReading:");
		printf("\n\t%s",hdrfile);
		printf("\n\t%s",imgfile);
	}

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, nx, ny, nz, dx, dy, dz, &dataType,v);

	nv = (*nx)*(*ny)*(*nz);

	*type = dataType;

	if(dataType==16)
	{
		im = read_image(imgfile,(int)(nv*sizeof(float)));

		if ( im!=NULL && analyzehdr.hk.sizeof_hdr != 348 )
		for(int i=0; i<nv; i += sizeof(float))
			swapByteOrder( im+i, sizeof(float));
	}
	else if(dataType==4)
	{
		im= read_image(imgfile,nv*sizeof(short));

		if (im!=NULL &&  analyzehdr.hk.sizeof_hdr != 348 )
			swapN( im, 2*nv);
	}
	else if(dataType==2)
	{
		im = read_image(imgfile,nv);
	}
	else
	{
		printf("\n\nSorry, but this program cannot handle this data type! Aborting ...\n\n");
	}

	return(im);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
char *read_analyze_image(const char *filename, DIM *dim, int *type, int v)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	short dataType;
	char *im=NULL;

	struct dsr analyzehdr;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if( !check_F_R_permission(hdrfile) )
	{
		printf("\nError: cannot read %s\n\n",hdrfile);
		return(NULL);
	}

	if( !check_F_R_permission(imgfile) )
	{
		printf("\nError: cannot read %s\n\n",imgfile);
		return(NULL);
	}

	if(v)
	{
		printf("\nReading:");
		printf("\n\t%s",hdrfile);
		printf("\n\t%s",imgfile);
	}

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, &(dim->nx), &(dim->ny), &(dim->nz), &(dim->dx), &(dim->dy), &(dim->dz), &dataType,v);

	nv = (dim->nx)*(dim->ny)*(dim->nz);

	*type = dataType;

	if(dataType==16)
	{
		im = read_image(imgfile,(int)(nv*sizeof(float)));

		if ( im!=NULL && analyzehdr.hk.sizeof_hdr != 348 )
		for(int i=0; i<nv; i += sizeof(float))
			swapByteOrder( im+i, sizeof(float));
	}
	else if(dataType==4)
	{
		im= read_image(imgfile,nv*sizeof(short));

		if (im!=NULL &&  analyzehdr.hk.sizeof_hdr != 348 )
			swapN( im, 2*nv);
	}
	else if(dataType==2)
	{
		im = read_image(imgfile,nv);
	}
	else
	{
		printf("\nSorry, but this program cannot handle this data type! Aborting ...\n\n");
		return(NULL);
	}

	return(im);
}

// assumes that the *.img and *.hdr files exist and we have read permission.
char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, int *type, int v)
{
	char hdrfile[1024];
	char imgfile[1024];
	int nv;
	short dataType;
	char *im;

	struct dsr analyzehdr;

	get_analyze_file_names(filename, hdrfile, imgfile);

	if( !check_F_R_permission(hdrfile) )
	{
		printf("\nError: cannot read %s\n",hdrfile);
		return(NULL);
	}

	if( !check_F_R_permission(imgfile) )
	{
		printf("\nError: cannot read %s\n",imgfile);
		return(NULL);
	}

	if(v)
	{
		printf("\nReading:");
		printf("\n\t%s",hdrfile);
		printf("\n\t%s",imgfile);
	}

	read_analyze_hdr(&analyzehdr, hdrfile);
	setDimensions(analyzehdr, nx, ny, nz, nt, dx, dy, dz, &dataType,v);

	nv = (*nx)*(*ny)*(*nz)*(*nt);

	*type = dataType;

	if(dataType==16)
	{
		im = read_image(imgfile,(int)(nv*sizeof(float)));

		if(im==NULL) return(NULL);

		if ( analyzehdr.hk.sizeof_hdr != 348 )
		for(int i=0; i<nv; i+=sizeof(float) )
			swapByteOrder( im+i, sizeof(float));

		return(im);
	}

	if(dataType==4)
	{
		im= read_image(imgfile,nv*2);

		if(im==NULL) return(NULL);

		if ( analyzehdr.hk.sizeof_hdr != 348 )
			swapN( (char *)im, 2*nv);

		return(im);
	}

	if(dataType==2)
	{
		im = read_image(imgfile,nv);
		return(im);
	}

	printf("\n\nSorry, but this program cannot handle this data type! Aborting ...\n\n");
	return(NULL);
}

int read_nt(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);
	}

	return( hdr.dime.dim[4] );
}

float read_dx(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	return( hdr.dime.pixdim[1] );
}

float read_dy(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	return( hdr.dime.pixdim[2] );
}

float read_dz(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		for(int i=0; i<8; i++)
			swapByteOrder( (char *) &hdr.dime.pixdim[i], sizeof(float) );
	}

	return( hdr.dime.pixdim[3] );
}

int read_nx(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);
	}

	return( hdr.dime.dim[1] );
}

int read_ny(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);
	}

	return( hdr.dime.dim[2] );
}

int read_nz(const char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapN( (char *) hdr.dime.dim , 16);
	}

	return( hdr.dime.dim[3] );
}

int read_datatype(char *filename)
{
	struct dsr hdr;
	char hdrfile[1024];
	char imgfile[1024];

	get_analyze_file_names(filename, hdrfile, imgfile);

	read_analyze_hdr(&hdr, hdrfile);

	if ( hdr.hk.sizeof_hdr != 348 )
	{
		swapByteOrder( (char *) &hdr.dime.datatype, sizeof(short) );
	}

	return( (int)(hdr.dime.datatype) );
}
