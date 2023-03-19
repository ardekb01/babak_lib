#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>   gcc compiler on mac did not like this
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>
#include <ge.h>
#include <volume.h>
#include <babak_lib.h>
#include <analyze.h>

static void get_directory(char *, char *);

// destructor
VOLUME::~VOLUME()
{
	if(data!=NULL) free(data);
	if(uc_data!=NULL) free(uc_data);
}

// constructor
VOLUME::VOLUME()
{
	// make sure that the data is initially NULL
	data=NULL;
	uc_data=NULL;
}

int VOLUME::isVision(char *pathname)
{
	int val=0;
        FILE *fp;
        char key[1025];

        fp=fopen(pathname,"r");

        if(fp==NULL)
        {
                printf("Error: Cannot open %s!\n",pathname);
                data=NULL;
                return(0);
        }

        fseek(fp,5742,0);
        if( fread(key,1,2,fp) != 2 )
        {
		val=0;
        }
        else if( strncmp(key,"TE",2)==0 )
        {
		val=1;
        }

        fclose(fp);

        return(val);
}

int VOLUME::isSMIS(char *pathname)
{
   int val=0;
   FILE *fp;
   char key[1025];
   char dum[1025];


   // Is the first file SMIS?
   fp=fopen(pathname,"r");

	if(fp==NULL)
	{
		printf("Error: Cannot open %s!\n",pathname);
		data=NULL;
		return 0;
	}

        while( fgets(dum,1024,fp) != NULL )
        {
                sscanf(dum,"%s",key);
                if( strcmp(key,":IM_TOTAL_SLICES")==0 )
                {
                        val = 1;
                        break;
                }
        }
        fclose(fp);

	return(val);
}

void VOLUME::scale_image()
{
	double scale;

	uc_data = (unsigned char *)calloc(np,sizeof(unsigned char));
                
	scale=255.0/(max-min);

	for(int i=0;i<np;i++)
		uc_data[i]=(unsigned char)rint((data[i]-min)*scale);
}

void VOLUME::read_image(char *pathname)
{
	if( isSMIS(pathname) )
	{
		read_smis_header(pathname);
		read_smis_image(pathname);
	}
	else if( isVision(pathname) )
	{
		read_vision_header(pathname);
		read_vision_image(pathname);
	}
	else
	{
		printf("Error: Cannot detect the image format!\n");
		data=NULL;
		return;
	}
	minmax();
	scale_image();
}

void VOLUME::read_smis_header(char *pathname)
{
	FILE *fp;
	char key[1025];
	char dum[1025];
	float fov;

	fp = fopen(pathname,"r");

	if(fp==NULL)
	{
		printf("Error: Cannot open %s!\n",pathname);
		data=NULL;
		return;
	}

        if( fread(&nx,sizeof(int),1,fp) != 1 ) {};
        if( fread(&ny,sizeof(int),1,fp) != 1 ) {};
        fclose(fp);

        if(bigEndian())
        {
                swapByteOrder( (char *)&nx, sizeof(int));
                swapByteOrder( (char *)&ny, sizeof(int));
        }

	np = nx*ny;

   fp=fopen(pathname,"r");
   if(fp==NULL) file_open_error(pathname);
   fseek(fp,np*2+512,0);
   while( fgets(dum,1024,fp) != NULL )
   {
                sscanf(dum,"%s",key);
                if( strcmp(key,":FOV")==0 )
                {
                        sscanf(dum,"%s %f",key,&fov);
                        continue;
                }
   }
   fclose(fp);

	dx=fov/nx;
	dy=fov/ny;
}

void VOLUME::read_vision_header(char *pathname)
{
	FILE *fp;
	double FOVrow, FOVcol;
	double slicethickness;
	double slicegapfactor;

	fp = fopen(pathname,"r");

	if(fp==NULL)
	{
		printf("Error: Cannot open %s!\n",pathname);
		data=NULL;
		return;
	}

	fseek(fp, 5000, SEEK_SET);
	if( fread(&dx, sizeof(double), 1, fp) != 1 ) {};
	if( !bigEndian() ) swapByteOrder((char *)&dx, sizeof(double) );

	fseek(fp, 5008, SEEK_SET);
	if( fread(&dy, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dy, sizeof(double) );

	fseek(fp, 3744, SEEK_SET);
	if( fread(&FOVrow, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVrow, sizeof(double) );

	fseek(fp, 3752, SEEK_SET);
	if( fread(&FOVcol, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVcol, sizeof(double) );

	nx = (int)(FOVrow/dx + 0.5);
	ny = (int)(FOVcol/dy + 0.5);
	np = nx*ny;

	fseek(fp, 1544, SEEK_SET);
	if( fread(&slicethickness, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicethickness, sizeof(double) );
 
	fseek(fp, 4136, SEEK_SET);
	if( fread(&slicegapfactor, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicegapfactor, sizeof(double) );
        
	// Ensure that the slice gap factor is a reasonable number
	if(slicegapfactor<0.0 || slicegapfactor>10000.0)
                slicegapfactor=0.0;
        
	dz = slicethickness * (1.0 + slicegapfactor);

	fseek(fp, 3212, SEEK_SET);
	if( fread(&image_number,sizeof(int),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&image_number, sizeof(int) );

	fclose(fp);
}

void VOLUME::read_smis_image(char *pathname)
{
   FILE *fp;

   fp = fopen(pathname,"r");

	if(fp==NULL)
	{
		printf("Error: Cannot open %s!\n",pathname);
		data=NULL;
		return;
	}

	if(data!=NULL) free(data);
	data = (short *)calloc(np,sizeof(short));
	fseek(fp, 512, SEEK_SET);
	if( fread(data,sizeof(short),np,fp) != (size_t)np ) {};
        fclose(fp);

	if(bigEndian()) swapN( (char *)data, np*2);
}

void VOLUME::read_vision_image(char *pathname)
{       
	FILE *fp;
        
	fp = fopen(pathname,"r");

	if(fp==NULL)
	{
		printf("Error: Cannot open %s!\n",pathname);
		data=NULL;
		return;
	}

	if(data!=NULL) free(data);
        data = (short *)calloc(np,sizeof(short));

	fseek(fp, 6144, SEEK_SET);
	if( fread(data, sizeof(short), np, fp) != (size_t)np ) {};
	if( !bigEndian() ) swapN((char *)(data), np*2);

	fclose(fp);
}


// Input - full path of an image in Siemens Vision format
// Behaviour - Examines all the files in the same directory as the input
// image.  Selects the readable files that have the same "patient id", 
// "patient name", "study number", TE and size as the input image.
// The images in the files that pass all tests are read into the data member
// of the VOLUME class along with matrix and voxels sizes and other 
// information. The image slices are arranged according to their "image number".
void VOLUME::readsiemensvision(char *pathname)
{
	FILE *fp;
	DIR *dp;		  // constant - directory pointer
	struct dirent *dir;	  // variable - directory entry
	struct stat buf;	  // variable - file stats returned by stat()
	off_t filesize;		  // constant - size of the input file (bytes)

	char patientID[13];       // constant - patient ID of the input file
	char patientID2[13];      // variable - patient ID's
	char patientname[26];     // constant - patient name in the input file
	char patientname2[26];    // variable - patient names
	char dirname[128];        // directory where files are looked for
	char filename[DEFAULT_STRING_LENGTH];       // variable - full path of files
	char filenames[1024][512]; // names of the files that pass the test

	double slicethickness;
	double slicegapfactor;
	double FOVrow, FOVcol;
	double echotime; 	// constant - TE of the input image
	double echotime2;   	// variable - TE's

	int studynumber;       // constant - study number of the input file
	int studynumber2;      // variable - study numbers 
	int imagenumber;       // variable - image numbers
	int imagenumbermin;    // minimum image number
	int imagenumbermax;    // maximum image number
	int numberoffiles=0;   // number of files that pass the test
	int imagenumbers[512]; // image numbers of the files that pass the test
	int i;

	// check to see if the file exist
	if( access(pathname,F_OK) == -1 )
	{
		printf("Error: %s does not exist\n",pathname);
		return;
	}

	// check read permission
	if( access(pathname,R_OK) == -1 )
	{
		printf("Error: No read permission for %s\n",pathname);
		return;
	}

	// obtain filesize
	stat(pathname,&buf);
	filesize=buf.st_size;

	// determine the patient ID and name of the input image
	// we will look for all files which have the same name and ID
   fp=fopen(pathname,"r");
   if(fp==NULL) file_open_error(pathname);

   nt=1;
   dt=0;

/*
fseek(fp, 2864, SEEK_SET);
if( fread(&nx, sizeof(int), 1, fp);
ny=nx;
np=nx*ny;
*/

	fseek(fp, 5000, SEEK_SET);
	if( fread(&dx, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dx, sizeof(double) );

	fseek(fp, 5008, SEEK_SET);
	if( fread(&dy, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dy, sizeof(double) );

	fseek(fp, 3744, SEEK_SET);
	if( fread(&FOVrow, sizeof(double), 1, fp) != 1 ) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVrow, sizeof(double) );

	fseek(fp, 3752, SEEK_SET);
	if( fread(&FOVcol, sizeof(double), 1, fp) != 1 ) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVcol, sizeof(double) );

	nx = (int)(FOVrow/dx + 0.5);
	ny = (int)(FOVcol/dy + 0.5);
	np = nx*ny;

	fseek(fp, 1544, SEEK_SET);
	if( fread(&slicethickness, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicethickness, sizeof(double) );
        
	fseek(fp, 4136, SEEK_SET);
	if( fread(&slicegapfactor, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicegapfactor, sizeof(double) );

	// Ensure that the slice gap factor is a reasonable number
	if(slicegapfactor<0.0 || slicegapfactor>10000.0)
		slicegapfactor=0.0;
        
	dz = slicethickness * (1.0 + slicegapfactor);

	fseek(fp, 795, SEEK_SET);
	if( fread(patientID,sizeof(char),12,fp) != 12) {};
	patientID[12]='\0';
	printf("Patient ID: %s\n",patientID);

	fseek(fp, 768, SEEK_SET);
	if( fread(patientname,sizeof(char),25,fp) != 25) {};
	patientname[25]='\0';
	printf("Patient Name: %s\n",patientname);

	fseek(fp, 3200, SEEK_SET);
	if( fread(&studynumber,sizeof(int),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&studynumber, sizeof(int) );
	printf("Study Number: %d\n",studynumber);

	fseek(fp, 3212, SEEK_SET);
	if( fread(&imagenumber,sizeof(int),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&imagenumber, sizeof(int) );
	imagenumbermin=imagenumber;
	imagenumbermax=imagenumber;

	fseek(fp, 1568, SEEK_SET);
	if( fread(&echotime, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&echotime, sizeof(double) );

	fseek(fp, 3792, SEEK_SET);
	if( fread(normalvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&normalvec[0], sizeof(double) );
		swapByteOrder((char *)&normalvec[1], sizeof(double) );
		swapByteOrder((char *)&normalvec[2], sizeof(double) );
	}

	fseek(fp, 3832, SEEK_SET);
	if( fread(rowvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&rowvec[0], sizeof(double) );
		swapByteOrder((char *)&rowvec[1], sizeof(double) );
		swapByteOrder((char *)&rowvec[2], sizeof(double) );
	}

	fseek(fp, 3856, SEEK_SET);
	if( fread(columnvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&columnvec[0], sizeof(double) );
		swapByteOrder((char *)&columnvec[1], sizeof(double) );
		swapByteOrder((char *)&columnvec[2], sizeof(double) );
	}

	fclose(fp);

	// determine the image directory
	get_directory(pathname,dirname);

	// open the directory and associate the directory stream dp with it
	dp=opendir(dirname);

	// If opendir() returns a NULL pointer, then dirname cannot be accessed
	// or if it cannot malloc() enough memory to hold the whole thing
	if(dp==NULL)
	{
		printf("Error: Cannot access %s\n", dirname);
		return;
	}

	// Get a pointer to the next directory entry using readdir().  
	// If readdir() returns NULL, the end of the directory has been 
	// reached.  
	while( (dir=readdir(dp)) != NULL )
	{
		// create the full path of the file to be examined
		snprintf(filename,sizeof(filename),"%s/%s",dirname,dir->d_name);

		// ensure that the file is readable
		if( access(filename,R_OK) == -1 ) continue;

		// obtain information about the file pointed to by filename
		stat(filename,&buf);

		// make sure that the file is a regular file
		if(!S_ISREG(buf.st_mode)) continue;

		// only process files of size filesize
		if( buf.st_size != filesize ) continue;

		// determine the patient ID and name of the file
		// we will look for all files which have the same name and ID
		fp=fopen(filename,"r");
                if(fp==NULL) file_open_error(filename);
		fseek(fp, 795, SEEK_SET);
		if( fread(patientID2,sizeof(char),12,fp) != 12) {};
		patientID2[12]='\0';
		        
		fseek(fp, 768, SEEK_SET);
		if( fread(patientname2,sizeof(char),25,fp) != 25) {};
		patientname2[25]='\0';
		        
		fseek(fp, 3200, SEEK_SET);
		if( fread(&studynumber2,sizeof(int),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&studynumber2, sizeof(int) );

		fseek(fp, 3212, SEEK_SET);
		if( fread(&imagenumber,sizeof(int),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&imagenumber, sizeof(int) );

		fseek(fp, 1568, SEEK_SET);
		if( fread(&echotime2, sizeof(double), 1, fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&echotime2, sizeof(double) );

		fclose(fp);

		if( strcmp(patientID, patientID2) != 0 ||
		strcmp(patientname, patientname2) != 0 ||
		studynumber != studynumber2 || 
		echotime != echotime2 )  continue;

		if(imagenumber<imagenumbermin) imagenumbermin=imagenumber;
		if(imagenumber>imagenumbermax) imagenumbermax=imagenumber;

		imagenumbers[numberoffiles]=imagenumber;
		snprintf(filenames[numberoffiles++],512,"%s",filename);

	};

	// close the directory stream and free the structure associated with dp
	closedir(dp);

	nz=imagenumbermax-imagenumbermin+1;

	if(data!=NULL) free(data);
	data=(short *)calloc(nx*ny*nz,sizeof(short));

	printf("\nReading image data ");
	for(i=0;i<numberoffiles;i++)
	{
		imagenumber=imagenumbers[i]-imagenumbermin;

		fp=fopen(filenames[i],"r");
        if(fp==NULL) file_open_error(filenames[i]);

		if(imagenumber==0)
		{
			fseek(fp, 3768, SEEK_SET);
			if( fread(centervec, sizeof(double), 3, fp) != 3) {};
			if( !bigEndian() ) 
			{
				swapByteOrder((char *)&centervec[0], sizeof(double) );
				swapByteOrder((char *)&centervec[1], sizeof(double) );
				swapByteOrder((char *)&centervec[2], sizeof(double) );
			}
		}

		fseek(fp, 6144, SEEK_SET);
		if( fread(data+imagenumber*np, sizeof(short), np, fp) != (size_t)np) {};
		if( !bigEndian() ) swapN((char *)(data+imagenumber*np), np*2);

		fclose(fp);

		printf(". "); fflush(NULL);
	}
	printf("\n\n");

	return;
}

int VOLUME::readsiemensvision(int N, char **filename)
{
	off_t filesize;		  // constant - size of the input file (bytes)

	FILE *fp;
	struct stat buf;	  // variable - file stats returned by stat()

	char patientID[13];       // constant - patient ID of the input file
	char patientID2[13];      // variable - patient ID's
	char patientname[26];     // constant - patient name in the input file
	char patientname2[26];    // variable - patient names

	double slicethickness;
	double slicegapfactor;
	double FOVrow, FOVcol;
	double echotime; 	// constant - TE of the input image
	double echotime2;   	// variable - TE's

	int studynumber;       // constant - study number of the input file
	int studynumber2;      // variable - study numbers 
	int imagenumber;       // variable - image numbers
	int imagenumbermin;    // minimum image number
	int imagenumbermax;    // maximum image number

	nt=1;
	dt=0.0;

	// ensure all input image files exist
	for(int i=0; i<N; i++)
	if( access(filename[i],F_OK) == -1 )
	{
		printf("Error: %s does not exist.\n",filename[i]);
		return(0);
	}

	// ensure read permission for all input image files
	for(int i=0; i<N; i++)
	if( access(filename[i],R_OK) == -1 )
	{
		printf("Error: No read permission for %s.\n",filename[i]);
		return(0);
	}

	// 1. ensure all input image files have the same size
	// 2. ensure that all files are regular files
	stat(filename[0],&buf);
	filesize=buf.st_size;
	for(int i=0; i<N; i++)
	{
		stat(filename[i],&buf);
		if( buf.st_size != filesize )
		{
			printf("Error: All input files must have the same size.\n");
			return(0);
		}

		if(!S_ISREG(buf.st_mode))
		{
			printf("Error: %s is not a regular file.\n",filename[i]);
			return(0);
		}
	}

	// determine the patient ID and name of the input image
	// we will look for all files which have the same name and ID
	fp=fopen(filename[0],"r");
    if(fp==NULL) file_open_error(filename[0]);

	fseek(fp, 5000, SEEK_SET);
	if( fread(&dx, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dx, sizeof(double) );

	fseek(fp, 5008, SEEK_SET);
	if( fread(&dy, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dy, sizeof(double) );

	fseek(fp, 3744, SEEK_SET);
	if( fread(&FOVrow, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVrow, sizeof(double) );

	fseek(fp, 3752, SEEK_SET);
	if( fread(&FOVcol, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVcol, sizeof(double) );

	nx = (int)(FOVrow/dx + 0.5);
	ny = (int)(FOVcol/dy + 0.5);
	np = nx*ny;

	fseek(fp, 1544, SEEK_SET);
	if( fread(&slicethickness, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicethickness, sizeof(double) );
        
	fseek(fp, 4136, SEEK_SET);
	if( fread(&slicegapfactor, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicegapfactor, sizeof(double) );

	// Ensure that the slice gap factor is a reasonable number
	if(slicegapfactor<0.0 || slicegapfactor>10000.0) slicegapfactor=0.0;
        
	dz = slicethickness * (1.0 + slicegapfactor);

	fseek(fp, 795, SEEK_SET);
	if(fread(patientID,sizeof(char),12,fp)!=12){};
	patientID[12]='\0';
	printf("Patient ID: %s\n",patientID);

	fseek(fp, 768, SEEK_SET);
	if(fread(patientname,sizeof(char),25,fp)!=25){};
	patientname[25]='\0';
	printf("Patient Name: %s\n",patientname);

	fseek(fp, 3200, SEEK_SET);
	if( fread(&studynumber,sizeof(int),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&studynumber, sizeof(int) );
	printf("Study Number: %d\n",studynumber);

	fseek(fp, 3212, SEEK_SET);
	if( fread(&imagenumber,sizeof(int),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&imagenumber, sizeof(int) );
	imagenumbermin=imagenumber;
	imagenumbermax=imagenumber;

	fseek(fp, 1568, SEEK_SET);
	if( fread(&echotime, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&echotime, sizeof(double) );

	fseek(fp, 3792, SEEK_SET);
	if( fread(normalvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&normalvec[0], sizeof(double) );
		swapByteOrder((char *)&normalvec[1], sizeof(double) );
		swapByteOrder((char *)&normalvec[2], sizeof(double) );
	}

	fseek(fp, 3832, SEEK_SET);
	if( fread(rowvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&rowvec[0], sizeof(double) );
		swapByteOrder((char *)&rowvec[1], sizeof(double) );
		swapByteOrder((char *)&rowvec[2], sizeof(double) );
	}

	fseek(fp, 3856, SEEK_SET);
	if( fread(columnvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&columnvec[0], sizeof(double) );
		swapByteOrder((char *)&columnvec[1], sizeof(double) );
		swapByteOrder((char *)&columnvec[2], sizeof(double) );
	}

	fclose(fp);

	for(int i=0; i<N; i++)
	{
		// determine the patient ID and name of the file
		// ensure all files which have the same name and ID
		fp=fopen(filename[i],"r");
        if(fp==NULL) file_open_error(filename[i]);

		fseek(fp, 795, SEEK_SET);
		if(fread(patientID2,sizeof(char),12,fp)!=12){};
		patientID2[12]='\0';
		        
		fseek(fp, 768, SEEK_SET);
		if(fread(patientname2,sizeof(char),25,fp)!=25){};
		patientname2[25]='\0';
		        
		fseek(fp, 3200, SEEK_SET);
		if( fread(&studynumber2,sizeof(int),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&studynumber2, sizeof(int) );

		fseek(fp, 3212, SEEK_SET);
		if( fread(&imagenumber,sizeof(int),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&imagenumber, sizeof(int) );

		fseek(fp, 1568, SEEK_SET);
		if( fread(&echotime2, sizeof(double), 1, fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&echotime2, sizeof(double) );

		fclose(fp);

		if( strcmp(patientID, patientID2) != 0 )
		{
			printf("Error: All images must have the same Patient ID.\n");
			return(0);
		}

		if( strcmp(patientname, patientname2) != 0 )
		{
			printf("Error: All images must have the same Patient Name.\n");
			return(0);
		}

		if( studynumber != studynumber2 )
		{
			printf("Error: All images must have the same Study Number.\n");
			return(0);
		}

		//if( echotime != echotime2 )
		//{
		//	printf("Error: All images must have the same Echo Time.\n");
		//	return(0);
		//}

		if(imagenumber<imagenumbermin) imagenumbermin=imagenumber;
		if(imagenumber>imagenumbermax) imagenumbermax=imagenumber;
	};

	nz=imagenumbermax-imagenumbermin+1;

	if(data!=NULL) free(data);
	data=(short *)calloc(nx*ny*nz,sizeof(short));

	if(data==NULL)
	{
		printf("Error: Cound not allocate enough memory.\n");
		return(0);
	}

	// printf("\nReading image data: ");
	for(int i=0;i<N;i++)
	{
		fp=fopen(filename[i],"r");
        if(fp==NULL) file_open_error(filename[i]);

		fseek(fp, 3212, SEEK_SET);
		if( fread(&imagenumber,sizeof(int),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&imagenumber, sizeof(int) );

		imagenumber -=  imagenumbermin;

		if(imagenumber==0)
		{
	        fseek(fp, 3768, SEEK_SET);
        	if( fread(centervec, sizeof(double), 3, fp) != 3) {};
			if( !bigEndian() ) 
			{
				swapByteOrder((char *)&centervec[0], sizeof(double) );
				swapByteOrder((char *)&centervec[1], sizeof(double) );
				swapByteOrder((char *)&centervec[2], sizeof(double) );
			}
		}

		fseek(fp, 6144, SEEK_SET);
		if( fread(data+imagenumber*np, sizeof(short), np, fp) != (size_t)np) {};
		if( !bigEndian() ) swapN((char *)(data+imagenumber*np), np*2);

		fclose(fp);

		printf("."); fflush(NULL);
	}
	printf("\n\n");

	return(1);
}

void VOLUME::readmosaic(char *pathname)
{
	FILE *fp;
	double FOVrow;
	double FOVcol;
	double slicethickness;
	double slicegapfactor;
	short *mosaic;
	int nix,niy;
	int i=0; // tile row index 0,...,512/nx-1
	int j=0; // tile column index 0,...,512/ny-1
	int m; // image index 0,1,...,nz-1
	int n; // nz/2 or (nz+1)/2
	int k; // index from 0 to nz-1
	int np;
	int W;

	nt=1;
	dt=0;

	fp=fopen(pathname,"r");
    if(fp==NULL) file_open_error(pathname);

	fseek(fp, 5000, SEEK_SET);
	if( fread(&dx, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&dx, sizeof(double) );

	fseek(fp, 3744, SEEK_SET);
	if( fread(&FOVrow, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVrow, sizeof(double) );

	W = (int)(FOVrow/dx + 0.5);

	fseek(fp, 2864, SEEK_SET);
	if( fread(&nx, sizeof(int), 1, fp) != 1) {} 
	if( !bigEndian() ) swapByteOrder((char *)&nx, sizeof(int) );
	ny=nx;
	np=nx*ny;
	nix=W/nx;
	niy=W/ny;

	fseek(fp, 3752, SEEK_SET);
	if( fread(&FOVcol, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&FOVcol, sizeof(double) );

	// dx and dy may be read directly from locations 5000 and 5008 
	// However, for MOSAIC format they contain FOVrow/512 and FOVcol/512
	// which not the real pixel size.
	dx = FOVrow/nx;
	dy = FOVcol/ny;

	fseek(fp, 1544, SEEK_SET);
	if( fread(&slicethickness, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicethickness, sizeof(double) );

	fseek(fp, 4136, SEEK_SET);
	if( fread(&slicegapfactor, sizeof(double), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&slicegapfactor, sizeof(double) );

	fseek(fp, 3768, SEEK_SET);
	if( fread(centervec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&centervec[0], sizeof(double) );
		swapByteOrder((char *)&centervec[1], sizeof(double) );
		swapByteOrder((char *)&centervec[2], sizeof(double) );
	}

	fseek(fp, 3792, SEEK_SET);
	if( fread(normalvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&normalvec[0], sizeof(double) );
		swapByteOrder((char *)&normalvec[1], sizeof(double) );
		swapByteOrder((char *)&normalvec[2], sizeof(double) );
	}

	fseek(fp, 3832, SEEK_SET);
	if( fread(rowvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&rowvec[0], sizeof(double) );
		swapByteOrder((char *)&rowvec[1], sizeof(double) );
		swapByteOrder((char *)&rowvec[2], sizeof(double) );
	}

	fseek(fp, 3856, SEEK_SET);
	if( fread(columnvec, sizeof(double), 3, fp) != 3) {};
	if( !bigEndian() ) 
	{
		swapByteOrder((char *)&columnvec[0], sizeof(double) );
		swapByteOrder((char *)&columnvec[1], sizeof(double) );
		swapByteOrder((char *)&columnvec[2], sizeof(double) );
	}

	// Ensure that the slice gap factor is a reasonable number 
	if(slicegapfactor<0.0 || slicegapfactor>10000.0)
		slicegapfactor=0.0;

	dz = slicethickness * (1.0 + slicegapfactor);

	mosaic=(short *)calloc(W*W,sizeof(short));
	fseek(fp, 6144, SEEK_SET);
	if( fread(mosaic, sizeof(short), W*W, fp) != (size_t)W*W) {};
	if( !bigEndian() ) swapN((char *)mosaic, 512*512*2);

	fclose(fp);

	// this part automatically determines nz
	k=0;
	j=W;
	while(k==0 && j>0)
	{
		j--;
 		for(i=W-1;i>=0;i--)	
		{
			k=mosaic[W*j+i];
			if( k != 0) break;
		}
	};
	nz = i/nx + 1 + (j/nx)*nix ;

	if(data!=NULL) free(data);
	data=(short *)calloc(nx*ny*nz,sizeof(short));

	if(nz%2)
		n=(nz+1)/2;
	else
		n=nz/2;

	k=0;
	for(j=0;j<niy;j++)
	for(i=0;i<nix;i++)
	if(k<nz)
	{
		if(k<n) m = 2*k;
		else m = 2*(k-n)+1;

		get_subimage(mosaic,W,W,i*nx,j*ny,data+m*np);
		k++;
	}
	else
	{
		break;
	}

	free(mosaic);
}

void VOLUME::savenki(char *pathname)
{
	FILE *fp;
	short id=1;
	short hdrsize=1024;
	char U=0;
	char C=0;
	char DT=1;

	fp=fopen(pathname,"w");
    if(fp==NULL) file_open_error(pathname);
	fwrite("NKI", sizeof(char), 3, fp);
	fwrite(&id, sizeof(short), 1, fp);
	fwrite(&hdrsize, sizeof(short), 1, fp);
	fwrite(&nx, sizeof(int), 1, fp);
	fwrite(&ny, sizeof(int), 1, fp);
	fwrite(&nz, sizeof(int), 1, fp);
	fwrite(&nt, sizeof(int), 1, fp);
	fwrite(&dx, sizeof(double), 1, fp);
	fwrite(&dy, sizeof(double), 1, fp);
	fwrite(&dz, sizeof(double), 1, fp);
	fwrite(&dt, sizeof(double), 1, fp);
	fwrite(centervec, sizeof(double), 3, fp);
	fwrite(normalvec, sizeof(double), 3, fp);
	fwrite(rowvec, sizeof(double), 3, fp);
	fwrite(columnvec, sizeof(double), 3, fp);
	fwrite(&DT, sizeof(char), 1, fp);
	fwrite(&U, sizeof(char), 1, fp);
	fwrite(&C, sizeof(char), 1, fp);
	fseek(fp,hdrsize,SEEK_SET);
	fwrite(data, sizeof(short), nx*ny*nz, fp);
	fclose(fp);
}

void VOLUME::get_subimage(short *im1, int nx1, int ny1, int x0, int y0, short *im2)
{
	int i,j;
	int offset1;
	int offset2;

	if(x0<0 || (x0+nx)>nx1) return;
	if(y0<0 || (y0+ny)>ny1) return;

	for(j=0;j<ny;j++)
	{
		offset1 = nx1*(j+y0) + x0;
		offset2 = nx*j;

		for(i=0;i<nx;i++)
			im2[offset2+i]=im1[offset1+i];
	}
}

static void get_directory(char *pathname, char *dirname)
{
	int n;
	int i;

	n = (int)strlen(pathname);

	for(i=n-1;i>=0;i--)
	if( pathname[i] == '/') break;

	if(i==-1)
	{
		dirname[0]='.';
		dirname[1]='\0';
	}
	else
	{
		strncpy(dirname, pathname, i);
		dirname[i]='\0';
	}

	return;
}

void VOLUME::minmax()
{
	min=max=data[0];
	for(int i=0;i<np;i++) 
	{
		if(data[i]<min) 
			min=data[i];
		else if(data[i]>max) 
			max=data[i];
	}
	return;
}

// Input - full path of an image in GE format
// Behaviour - Examines all the files in the same directory as the input
// image.  Selects the readable files that have the same "patient id", 
// "patient name", "study number", TE and size as the input image.
// The images in the files that pass all tests are read into the data member
// of the VOLUME class along with matrix and voxels sizes and other 
// information. The image slices are arranged according to their "image number".
void VOLUME::read_ge(char *pathname)
{
	FILE *fp;
	DIR *dp;		  // constant - directory pointer
	struct dirent *dir;	  // variable - directory entry
	struct stat buf;	  // variable - file stats returned by stat()
	off_t filesize;		  // constant - size of the input file (bytes)

	char patientID[13];       // constant - patient ID of the input file
	char patientID2[13];      // variable - patient ID's
	char patientname[26];     // constant - patient name in the input file
	char patientname2[26];    // variable - patient names
	char dirname[128];        // directory where files are looked for
	char filename[DEFAULT_STRING_LENGTH];       // variable - full path of files
	char filenames[1024][512]; // names of the files that pass the test

	double ddum;

	int echotime; 	// constant - TE of the input image
	int echotime2;   	// variable - TE's
	int im_hdr_offset;
	int series_hdr_offset;
	int exam_hdr_offset;
	int numberoffiles=0;   // number of files that pass the test
	int i;
	int hdr_size;

	float TLHC[3];
	float TRHC[3];
	float BRHC[3];
	float cntr[3];
	float fdum;

	short sdum;
	short seriesnumber;       
	short seriesnumber2;     

	// check to see if the file exist
	if( access(pathname,F_OK) == -1 )
	{
		printf("Error: %s does not exist\n",pathname);
		data=NULL;
		return;
	}

	// check read permission
	if( access(pathname,R_OK) == -1 )
	{
		printf("Error: No read permission for %s\n",pathname);
		data=NULL;
		return;
	}

	// obtain filesize
	stat(pathname,&buf);
	filesize=buf.st_size;

	// determine the patient ID and name of the input image
	// we will look for all files which have the same name and ID
	fp=fopen(pathname,"r");
    if(fp==NULL) file_open_error(pathname);

	nt=1;
	dt=0;

	fseek(fp,4,0);
	if( fread(&hdr_size,sizeof(int),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&hdr_size, sizeof(int));

	fseek(fp,148,0);
	if( fread(&im_hdr_offset,sizeof(int),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&im_hdr_offset, sizeof(int));
	
	fseek(fp,132,0);
	if( fread(&exam_hdr_offset,sizeof(int),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&exam_hdr_offset, sizeof(int));
	
	fseek(fp,140,0);
	if( fread(&series_hdr_offset,sizeof(int),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&series_hdr_offset, sizeof(int));

	fseek(fp,im_hdr_offset+50,0);
	if( fread(&fdum,sizeof(float),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
	dx = fdum;

	fseek(fp,im_hdr_offset+54,0);
	if( fread(&fdum,sizeof(float),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
	dy = fdum;

	fseek(fp,im_hdr_offset+30,0);
	if( fread(&sdum,sizeof(short),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&sdum, sizeof(short));
	nx = sdum;

	fseek(fp,im_hdr_offset+32,0);
	if( fread(&sdum,sizeof(short),1,fp) != 1) {};
	if(!bigEndian() ) swapByteOrder( (char *)&sdum, sizeof(short));
	ny = sdum;

	np = nx*ny;

	fseek(fp,exam_hdr_offset+84,0);
	if(fread(patientID,sizeof(char),12,fp)!=12) {};
	patientID[12]='\0';
	printf("Patient ID: %s\n",patientID);

	fseek(fp,exam_hdr_offset+97,0);
	if(fread(patientname,sizeof(char),25,fp)!=25) {};
	patientname[25]='\0';
	printf("Patient Name: %s\n",patientname);

	fseek(fp,im_hdr_offset+202,0);
	if( fread(&echotime, sizeof(int), 1, fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&echotime, sizeof(int) );
	printf("Echo time=%d (microseconds)\n",echotime);

	fseek(fp, series_hdr_offset+10, 0);
	if( fread(&seriesnumber,sizeof(short),1,fp) != 1) {};
	if( !bigEndian() ) swapByteOrder((char *)&seriesnumber, sizeof(short) );
	printf("Series number: %d\n",seriesnumber);

	fseek(fp,im_hdr_offset+154,0);
	if( fread(TLHC,sizeof(float),3,fp) != 3) {};
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&TLHC[0], sizeof(float));
		swapByteOrder( (char *)&TLHC[1], sizeof(float));
		swapByteOrder( (char *)&TLHC[2], sizeof(float));
	}
	TLHC[0] *= -1.0; TLHC[2] *= -1.0; // RAS to FAL system converstion

	fseek(fp,im_hdr_offset+166,0);
	if( fread(TRHC,sizeof(float),3,fp) != 3) {};
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&TRHC[0], sizeof(float));
		swapByteOrder( (char *)&TRHC[1], sizeof(float));
		swapByteOrder( (char *)&TRHC[2], sizeof(float));
	}
	TRHC[0] *= -1.0; TRHC[2] *= -1.0; // RAS to FAL system converstion

	fseek(fp,im_hdr_offset+178,0);
	if( fread(BRHC,sizeof(float),3,fp) != 3) {};
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&BRHC[0], sizeof(float));
		swapByteOrder( (char *)&BRHC[1], sizeof(float));
		swapByteOrder( (char *)&BRHC[2], sizeof(float));
	}
	BRHC[0] *= -1.0; BRHC[2] *= -1.0; // RAS to FAL system converstion

	for(int i=0; i<3; i++) rowvec[i] = TRHC[i]-TLHC[i];
	for(int i=0; i<3; i++) columnvec[i] = BRHC[i]-TRHC[i];

	ddum = sqrt( columnvec[0]*columnvec[0] + columnvec[1]*columnvec[1] + columnvec[2]*columnvec[2] );
	for(int i=0; i<3; i++) columnvec[i]/=ddum;

	ddum = sqrt( rowvec[0]*rowvec[0] + rowvec[1]*rowvec[1] + rowvec[2]*rowvec[2] );
	for(int i=0; i<3; i++) rowvec[i]/=ddum;

	fclose(fp);

	// determine the image directory
	get_directory(pathname,dirname);

	// open the directory and associate the directory stream dp with it
	dp=opendir(dirname);

	// If opendir() returns a NULL pointer, then dirname cannot be accessed
	// or if it cannot malloc() enough memory to hold the whole thing
	if(dp==NULL)
	{
		printf("Error: Cannot access %s\n", dirname);
		data=NULL;
		return;
	}

	// Get a pointer to the next directory entry using readdir().  
	// If readdir() returns NULL, the end of the directory has been 
	// reached.  
	while( (dir=readdir(dp)) != NULL )
	{
		// create the full path of the file to be examined
		snprintf(filename,sizeof(filename),"%s/%s",dirname,dir->d_name);

		// ensure that the file is readable
		if( access(filename,R_OK) == -1 ) continue;

		// obtain information about the file pointed to by filename
		stat(filename,&buf);

		// make sure that the file is a regular file
		if(!S_ISREG(buf.st_mode)) continue;

		// only process files of size filesize
		if( buf.st_size != filesize ) continue;

		// determine the patient ID and name of the file
		// we will look for all files which have the same name and ID
		fp=fopen(filename,"r");
        	if(fp==NULL) file_open_error(filename);

		fseek(fp,exam_hdr_offset+84,0);
		if(fread(patientID2,sizeof(char),12,fp)!= 12) {};
		patientID2[12]='\0';

		if( strcmp(patientID, patientID2) != 0 ) continue;
		        
		fseek(fp,exam_hdr_offset+97,0);
		if( fread(patientname2,sizeof(char),25,fp)!= 25) {};
		patientname2[25]='\0';
		        
		if( strcmp(patientname, patientname2) != 0 ) continue;

		fseek(fp, series_hdr_offset+10, 0);
		if( fread(&seriesnumber2,sizeof(short),1,fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&seriesnumber2, sizeof(short) );

		if( seriesnumber != seriesnumber2 ) continue;

		fseek(fp,im_hdr_offset+202,0);
		if( fread(&echotime2, sizeof(int), 1, fp) != 1) {};
		if( !bigEndian() ) swapByteOrder((char *)&echotime2, sizeof(int) );

		if( echotime != echotime2 )  continue;

		fclose(fp);

		snprintf(filenames[numberoffiles++],512,"%s",filename);
	};

	// close the directory stream and free the structure associated with dp
	closedir(dp);

	nz=numberoffiles;

	if(data!=NULL) free(data);
	data=(short *)calloc(nx*ny*nz,sizeof(short));

	printf("\nReading image data ");
	for(i=0;i<numberoffiles;i++)
	{
		fp=fopen(filenames[i],"r");
                if(fp==NULL) file_open_error(filenames[i]);

		fseek(fp, hdr_size, SEEK_SET);
		if( fread(data+i*np, sizeof(short), np, fp) != (size_t)np) {};
		if( !bigEndian() ) swapN((char *)(data+i*np), np*2);

		fclose(fp);

		printf(". "); fflush(NULL);
	}
	printf("\n\n");

	if(nz==1)
	{
			fp=fopen(filenames[0],"r");
            if(fp==NULL) file_open_error(filenames[0]);

			fseek(fp,im_hdr_offset+130,0);
			if( fread(cntr,sizeof(float),3,fp) != 3) {};
			if(!bigEndian() )
			{
				swapByteOrder( (char *)&cntr[0], sizeof(float));
				swapByteOrder( (char *)&cntr[1], sizeof(float));
				swapByteOrder( (char *)&cntr[2], sizeof(float));
			}
			cntr[0] *= -1.0; cntr[2] *= -1.0; // RAS to FAL system converstion

			centervec[0]=cntr[0];
			centervec[1]=cntr[1];
			centervec[2]=cntr[2];

			fseek(fp,im_hdr_offset+26,0);
			if( fread(&fdum,sizeof(float),1,fp) != 1) {};
			if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
			dz = fdum;
			fclose(fp);

			normalvec[0] = -rowvec[1]*columnvec[2] + rowvec[2]*columnvec[1];
			normalvec[1] = -rowvec[2]*columnvec[0] + rowvec[0]*columnvec[2];
			normalvec[2] = -rowvec[0]*columnvec[1] + rowvec[1]*columnvec[0];

			ddum = sqrt( normalvec[0]*normalvec[0]+normalvec[1]*normalvec[1]+normalvec[2]*normalvec[2]);
			for(int q=0; q<3; q++) normalvec[q]/=ddum;
	} 
	else if( nz>=2)
	{
		float Ca[3];
		float Cb[3];

		for(int i=0;i<numberoffiles;i++)
		{
			if(i==0)
			{
				fp=fopen(filenames[i],"r");
                if(fp==NULL) file_open_error(filenames[i]);

				fseek(fp,im_hdr_offset+130,0);
				if( fread(Ca,sizeof(float),3,fp) != 3) {};
				if(!bigEndian() )
				{
					swapByteOrder( (char *)&Ca[0], sizeof(float));
					swapByteOrder( (char *)&Ca[1], sizeof(float));
					swapByteOrder( (char *)&Ca[2], sizeof(float));
				}
				fclose(fp);

				Ca[0] *= -1.0; Ca[2] *= -1.0; // RAS to FAL system converstion
				break;
			}
		}

		for(int i=0;i<numberoffiles;i++)
		{
			if(i==1)
			{
				fp=fopen(filenames[i],"r");
                if(fp==NULL) file_open_error(filenames[i]);

				fseek(fp,im_hdr_offset+130,0);
				if( fread(Cb,sizeof(float),3,fp) != 3) {};
				if(!bigEndian() )
				{
					swapByteOrder( (char *)&Cb[0], sizeof(float));
					swapByteOrder( (char *)&Cb[1], sizeof(float));
					swapByteOrder( (char *)&Cb[2], sizeof(float));
				}
				fclose(fp);

				Cb[0] *= -1.0; Cb[2] *= -1.0; // RAS to FAL system converstion
				break;
			}
		}

		for(int i=0; i<3; i++) normalvec[i]=Cb[i]-Ca[i];
		
		dz = sqrt( normalvec[0]*normalvec[0] + normalvec[1]*normalvec[1] + normalvec[2]*normalvec[2] );

		for(int i=0; i<3; i++) normalvec[i]/=dz;

		centervec[0]=Ca[0];
		centervec[1]=Ca[1];
		centervec[2]=Ca[2];

		// This is done in the function compute_T() in art.c
		// centervec[0]=Ca[0] + ((nz-1.0)*dz/2.0)*normalvec[0];
		// centervec[1]=Ca[1] + ((nz-1.0)*dz/2.0)*normalvec[1];
		// centervec[2]=Ca[2] + ((nz-1.0)*dz/2.0)*normalvec[2];
	}

	return;
}

void VOLUME::read_gelx(char *pathname)
{
	SUITEDATATYPE   suite_hdr;
	EXAMDATATYPE    exam_hdr;
	SERIESDATATYPE  series_hdr;
	MRIMAGEDATATYPE image_hdr;
	PixHdr          pix_hdr;

	FILE *fp;
	DIR *dp;		  // constant - directory pointer
	struct dirent *dir;	  // variable - directory entry
	struct stat buf;	  // variable - file stats returned by stat()

	char patientID[13];       // constant - patient ID of the input file
	char patientID2[13];      // variable - patient ID's
	char patientname[26];     // constant - patient name in the input file
	char patientname2[26];    // variable - patient names
	char dirname[128];        // directory where files are looked for
	char filename[DEFAULT_STRING_LENGTH];       // variable - full path of files
	char filenames[1024][512]; // names of the files that pass the test

	double ddum;

	int echotime; 	// constant - TE of the input image
	int echotime2;   	// variable - TE's
	int numberoffiles=0;   // number of files that pass the test
	int i;
	int hdr_size;

	float TLHC[3];
	float TRHC[3];
	float BRHC[3];
	float cntr[3];
	float fdum;

	short sdum;
	short seriesnumber;       
	short seriesnumber2;     

	// check to see if the file exist
	if( access(pathname,F_OK) == -1 )
	{
		printf("Error: %s does not exist\n",pathname);
		data=NULL;
		return;
	}

	// check read permission
	if( access(pathname,R_OK) == -1 )
	{
		printf("Error: No read permission for %s\n",pathname);
		data=NULL;
		return;
	}

	// determine the patient ID and name of the input image
	// we will look for all files which have the same name and ID
	fp=fopen(pathname,"r");
    if(fp==NULL) file_open_error(pathname);

	nt=1;
	dt=0;

	if( fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp) != 1) {};
	fseek(fp,116,SEEK_SET);
	if( fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp) != 1) {};
	fseek(fp,116+1040,SEEK_SET);
	if( fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp) != 1) {};
	fseek(fp,116+1040+1028,SEEK_SET);
	if( fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp) != 1) {};

	fdum = image_hdr.pixsize_X;
	if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
	dx = fdum;

	fdum = image_hdr.pixsize_Y;
	if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
	dy = fdum;

	sdum = image_hdr.imatrix_X;
	if(!bigEndian() ) swapByteOrder( (char *)&sdum, sizeof(short));
	nx = sdum;

	sdum = image_hdr.imatrix_Y;
	if(!bigEndian() ) swapByteOrder( (char *)&sdum, sizeof(short));
	ny = sdum;

	np = nx*ny;

	snprintf(patientID,sizeof(patientID),"%s",exam_hdr.patid);
	printf("Patient ID: %s\n",patientID);

	snprintf(patientname,sizeof(patientname),"%s",exam_hdr.patname);
	printf("Patient Name: %s\n",patientname);

	echotime = image_hdr.te;
	if( !bigEndian() ) swapByteOrder((char *)&echotime, sizeof(int) );
	printf("Echo time=%d (microseconds)\n",echotime);

	seriesnumber=series_hdr.se_no;
	if( !bigEndian() ) swapByteOrder((char *)&seriesnumber, sizeof(short) );
	printf("Series number: %d\n",seriesnumber);

	TLHC[0] = image_hdr.tlhc_R;
	TLHC[1] = image_hdr.tlhc_A;
	TLHC[2] = image_hdr.tlhc_S;
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&TLHC[0], sizeof(float));
		swapByteOrder( (char *)&TLHC[1], sizeof(float));
		swapByteOrder( (char *)&TLHC[2], sizeof(float));
	}
	TLHC[0] *= -1.0; TLHC[2] *= -1.0; // RAS to FAL system converstion

	TRHC[0] = image_hdr.trhc_R;
	TRHC[1] = image_hdr.trhc_A;
	TRHC[2] = image_hdr.trhc_S;
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&TRHC[0], sizeof(float));
		swapByteOrder( (char *)&TRHC[1], sizeof(float));
		swapByteOrder( (char *)&TRHC[2], sizeof(float));
	}
	TRHC[0] *= -1.0; TRHC[2] *= -1.0; // RAS to FAL system converstion

	BRHC[0] = image_hdr.brhc_R;
	BRHC[1] = image_hdr.brhc_A;
	BRHC[2] = image_hdr.brhc_S;
	if(!bigEndian() ) 
	{
		swapByteOrder( (char *)&BRHC[0], sizeof(float));
		swapByteOrder( (char *)&BRHC[1], sizeof(float));
		swapByteOrder( (char *)&BRHC[2], sizeof(float));
	}
	BRHC[0] *= -1.0; BRHC[2] *= -1.0; // RAS to FAL system converstion

	for(int i=0; i<3; i++) rowvec[i] = TRHC[i]-TLHC[i];
	for(int i=0; i<3; i++) columnvec[i] = BRHC[i]-TRHC[i];

	ddum = sqrt( columnvec[0]*columnvec[0] + columnvec[1]*columnvec[1] + columnvec[2]*columnvec[2] );
	for(int i=0; i<3; i++) columnvec[i]/=ddum;

	ddum = sqrt( rowvec[0]*rowvec[0] + rowvec[1]*rowvec[1] + rowvec[2]*rowvec[2] );
	for(int i=0; i<3; i++) rowvec[i]/=ddum;

	fclose(fp);

	// determine the image directory
	get_directory(pathname,dirname);

	// open the directory and associate the directory stream dp with it
	dp=opendir(dirname);

	// If opendir() returns a NULL pointer, then dirname cannot be accessed
	// or if it cannot malloc() enough memory to hold the whole thing
	if(dp==NULL)
	{
		printf("Error: Cannot access %s\n", dirname);
		data=NULL;
		return;
	}

	// Get a pointer to the next directory entry using readdir().  
	// If readdir() returns NULL, the end of the directory has been 
	// reached.  
	while( (dir=readdir(dp)) != NULL )
	{
		// create the full path of the file to be examined
		snprintf(filename,sizeof(filename),"%s/%s",dirname,dir->d_name);

		// ensure that the file is readable
		if( access(filename,R_OK) == -1 ) continue;

		// obtain information about the file pointed to by filename
		stat(filename,&buf);

		// make sure that the file is a regular file
		if(!S_ISREG(buf.st_mode)) continue;

		// determine the patient ID and name of the file
		// we will look for all files which have the same name and ID
		fp=fopen(filename,"r");
        if(fp==NULL) file_open_error(filename);

		if( fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp)!=1 ) { fclose(fp); continue; }

		fseek(fp,116,SEEK_SET);
		if( fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp)!=1 )       { fclose(fp); continue; }

		fseek(fp,116+1040,SEEK_SET);
		if( fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp)!=1 )   { fclose(fp); continue; }

		fseek(fp,116+1040+1028,SEEK_SET);
		if( fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp)!= 1)   { fclose(fp); continue; }

		snprintf(patientID2,sizeof(patientID2),"%s",exam_hdr.patid);
		if( strcmp(patientID, patientID2) != 0 ) { fclose(fp); continue; }
		        
		snprintf(patientname2,sizeof(patientname2),"%s",exam_hdr.patname);
		if( strcmp(patientname, patientname2) != 0 ) { fclose(fp); continue; }

		seriesnumber2=series_hdr.se_no;
		if( !bigEndian() ) swapByteOrder((char *)&seriesnumber2, sizeof(short) );
		if( seriesnumber != seriesnumber2 ) { fclose(fp); continue; }

		echotime2 = image_hdr.te;
		if( !bigEndian() ) swapByteOrder((char *)&echotime2, sizeof(int) );
		if( echotime != echotime2 )  { fclose(fp); continue; }

		fclose(fp);

		snprintf(filenames[numberoffiles++],512,"%s",filename);
	};

	// close the directory stream and free the structure associated with dp
	closedir(dp);

	nz=numberoffiles;

	if(data!=NULL) free(data);
	data=(short *)calloc(nx*ny*nz,sizeof(short));

	printf("\nReading image data ");
	for(i=0;i<numberoffiles;i++)
	{
		fp=fopen(filenames[i],"r");
        if(fp==NULL) file_open_error(filenames[i]);

                if( fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp) != 1) {};
		fseek(fp,116,SEEK_SET);
                if( fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp) != 1) {};
		fseek(fp,116+1040,SEEK_SET);
                if( fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp) != 1) {};
		fseek(fp,116+1040+1028,SEEK_SET);
                if( fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp) != 1) {};
		fseek(fp,116+1040+1028+1044,SEEK_SET);
                if( fread(&pix_hdr, sizeof(PixHdr), 1, fp) != 1) {};

                hdr_size=pix_hdr.img_hdr_length;
		if( !bigEndian() ) swapByteOrder( (char *)&hdr_size, sizeof(int));
		hdr_size = hdr_size + sizeof(SUITEDATATYPE) + sizeof(EXAMDATATYPE) + sizeof(SERIESDATATYPE)
                + sizeof(MRIMAGEDATATYPE);

		fseek(fp, hdr_size, SEEK_SET);
		if( fread(data+i*np, sizeof(short), np, fp) != (size_t)np ) {};
		if( !bigEndian() ) swapN((char *)(data+i*np), np*2);

		fclose(fp);

		printf(". "); fflush(NULL);
	}
	printf("\n\n");

	if(nz==1)
	{
		fp=fopen(filenames[0],"r");
        if(fp==NULL) file_open_error(filenames[0]);

                if( fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp) != 1) {};
		fseek(fp,116,SEEK_SET);
                if( fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp) != 1) {};
		fseek(fp,116+1040,SEEK_SET);
                if( fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp) != 1) {};
		fseek(fp,116+1040+1028,SEEK_SET);
                if(fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp) != 1) {};

		cntr[0]=image_hdr.ctr_R;
		cntr[1]=image_hdr.ctr_A;
		cntr[2]=image_hdr.ctr_S;

		if(!bigEndian() )
		{
			swapByteOrder( (char *)&cntr[0], sizeof(float));
			swapByteOrder( (char *)&cntr[1], sizeof(float));
			swapByteOrder( (char *)&cntr[2], sizeof(float));
		}
		cntr[0] *= -1.0; cntr[2] *= -1.0; // RAS to FAL system converstion

		centervec[0]=cntr[0];
		centervec[1]=cntr[1];
		centervec[2]=cntr[2];

		fdum=image_hdr.slthick;
		if(!bigEndian() ) swapByteOrder( (char *)&fdum, sizeof(float));
		dz = fdum;
		fclose(fp);

		normalvec[0] = -rowvec[1]*columnvec[2] + rowvec[2]*columnvec[1];
		normalvec[1] = -rowvec[2]*columnvec[0] + rowvec[0]*columnvec[2];
		normalvec[2] = -rowvec[0]*columnvec[1] + rowvec[1]*columnvec[0];

		ddum = sqrt( normalvec[0]*normalvec[0]+normalvec[1]*normalvec[1]+normalvec[2]*normalvec[2]);
		for(int q=0; q<3; q++) normalvec[q]/=ddum;
	} 

	if( nz>=2)
	{
		float Ca[3];
		float Cb[3];

		for(int i=0;i<numberoffiles;i++)
		{
			if(i==0)
			{
				fp=fopen(filenames[i],"r");
                if(fp==NULL) file_open_error(filenames[i]);

				if(fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp) != 1) {};
				fseek(fp,116,SEEK_SET);
				if(fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp) != 1) {};
				fseek(fp,116+1040,SEEK_SET);
				if(fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp) != 1) {};
				fseek(fp,116+1040+1028,SEEK_SET);
				if(fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp) != 1) {};

				Ca[0]=image_hdr.ctr_R;
				Ca[1]=image_hdr.ctr_A;
				Ca[2]=image_hdr.ctr_S;
				if(!bigEndian() )
				{
					swapByteOrder( (char *)&Ca[0], sizeof(float));
					swapByteOrder( (char *)&Ca[1], sizeof(float));
					swapByteOrder( (char *)&Ca[2], sizeof(float));
				}
				fclose(fp);

				Ca[0] *= -1.0; Ca[2] *= -1.0; // RAS to FAL system converstion
				break;
			}
		}

		for(int i=0;i<numberoffiles;i++)
		{
			if(i==1)
			{
				fp=fopen(filenames[i],"r");
                                if(fp==NULL) file_open_error(filenames[i]);

				if( fread(&suite_hdr, sizeof(SUITEDATATYPE), 1, fp) != 1) {};
				fseek(fp,116,SEEK_SET);
				if( fread(&exam_hdr, sizeof(EXAMDATATYPE), 1, fp) != 1) {};
				fseek(fp,116+1040,SEEK_SET);
				if( fread(&series_hdr, sizeof(SERIESDATATYPE), 1, fp) != 1) {};
				fseek(fp,116+1040+1028,SEEK_SET);
				if( fread(&image_hdr, sizeof(MRIMAGEDATATYPE), 1, fp) != 1) {};

				Cb[0]=image_hdr.ctr_R;
				Cb[1]=image_hdr.ctr_A;
				Cb[2]=image_hdr.ctr_S;
				if(!bigEndian() )
				{
					swapByteOrder( (char *)&Cb[0], sizeof(float));
					swapByteOrder( (char *)&Cb[1], sizeof(float));
					swapByteOrder( (char *)&Cb[2], sizeof(float));
				}
				fclose(fp);

				Cb[0] *= -1.0; Cb[2] *= -1.0; // RAS to FAL system converstion
				break;
			}
		}

		for(int i=0; i<3; i++) normalvec[i]=Cb[i]-Ca[i];
		
		dz = sqrt( normalvec[0]*normalvec[0] + normalvec[1]*normalvec[1] + normalvec[2]*normalvec[2] );

		for(int i=0; i<3; i++) normalvec[i]/=dz;

		centervec[0]=Ca[0];
		centervec[1]=Ca[1];
		centervec[2]=Ca[2];
	}

	return;
}
// Input - full path of an image in DICOM format
// Behaviour - Examines all the files in the same directory as the input
// image.  Selects the readable files that have the same "patient id", 
// "patient name", "study number", TE and size as the input image.
// The images in the files that pass all tests are read into the data member
// of the VOLUME class along with matrix and voxels sizes and other 
// information. The image slices are arranged according to their "image number".
void VOLUME::read_dicom(char *pathname)
{
	// check to see if the file exist
	if( access(pathname,F_OK) == -1 )
	{
		printf("Error: %s does not exist\n",pathname);
		data=NULL;
		return;
	}

	// check read permission
	if( access(pathname,R_OK) == -1 )
	{
		printf("Error: No read permission for %s\n",pathname);
		data=NULL;
		return;
	}
}

int VOLUME::read_analyze(const char *pathname)
{
	char hdrfile[1024];
	char imgfile[1024];
	FILE *fp;

	struct dsr analyzehdr;

	get_analyze_file_names(pathname, hdrfile, imgfile);

	if( !check_F_R_permission(hdrfile) )
	{
		printf("\nError: cannot read %s\n",hdrfile);
		return(0);
	}

	if( !check_F_R_permission(imgfile) )
	{
		printf("\nError: cannot read %s\n",imgfile);
		return(0);
	}

	read_analyze_hdr(&analyzehdr, hdrfile);

    if ( analyzehdr.hk.sizeof_hdr != 348 )
    {
        swapN( (char *) analyzehdr.dime.dim , 16);

        for(int i=0; i<8; i++)
            swapByteOrder( (char *) &analyzehdr.dime.pixdim[i], sizeof(float) );
    }

	dx=analyzehdr.dime.pixdim[1];
	dy=analyzehdr.dime.pixdim[2];
	dz=analyzehdr.dime.pixdim[3];

    ny=analyzehdr.dime.dim[2];
    nx=analyzehdr.dime.dim[1];
    nz=analyzehdr.dime.dim[3];

	np = nx*ny;
	nv = nx*ny*nz;

	data =(short *)calloc(nv, sizeof(short));

	fp = fopen(imgfile,"r");
        if(fp==NULL) file_open_error(imgfile);
	if( fread(data,sizeof(short),nv,fp) != (size_t)nv ) {};
	fclose(fp);

	if ( analyzehdr.hk.sizeof_hdr != 348 )
		swapN( (char *)data, 2*nv);

	return(1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void VOLUME::get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img)
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

void VOLUME::read_analyze_hdr(struct dsr *hdr, char *filename)
{
	FILE *fp;
	int n;

	fp=fopen(filename,"r");
	if(fp==NULL) 
	{
		printf("error: cannot open file %s !\n",filename);
		exit(0);
	}

	n=(int)fread(hdr,sizeof(struct dsr),1,fp);

	if(n!=1) 
	{
		printf("error: reading from file %s !\n",filename);
		fclose(fp);
		exit(0);
	}

	fclose(fp);
}

int VOLUME::check_F_R_permission(char *file)
{

	if( access(file,F_OK) == -1 )
		return(0);

	if( access(file,R_OK) == -1 )
		return(0);

	return(1);
}
