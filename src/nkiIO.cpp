/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  program: nkiIO.c                                          *
*  Copyright 2004 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <nki.h>
#include <babak_lib.h>

// Returns 1 if file can be read and is NKI format, 0 otherwise
int isNKI(char *file)
{
	FILE *fp;
	char id[3];
	int n;

	fp = fopen(file,"r");
	if(fp==NULL) return(0);
	n=fread(id,sizeof(char),3,fp);
	fclose(fp);
	if(n!=3) return(0);

	if(id[0]=='N' && id[1]=='K' && id[2]=='I')
		return(1);
	else
		return(0);
}

// returns 0 on failure, 1 on success
int readNKI(char *filename, nki *image)
{
	int nv;
	FILE *fp;

	// send error message if the file doesn't exist
	if(access(filename,F_OK)==-1)
	{
		printf("\n\nError: %s does not exist!\n\n",filename); 
		return(0);
	}

	// send error message if the file has no read permission   
	if(access(filename,R_OK)==-1)
	{
		printf("\n\nError: read permission to %s denied!\n\n",filename); 
		return(0);
	}  

	if( isNKI(filename)==0 )
	{
		printf("\n\nError: %s does not seem to be in NKI format!\n\n",filename);
		return(0);
	}

	fp=fopen(filename,"r");
	if(fp==NULL)
	{
		printf("\n\nCannot open %s!\n\n",filename);
		return(0);
	}       
                
	(void)fread(image->id,sizeof(char),3,fp);

	(void)fread(&image->flg,sizeof(short),1,fp);
	(void)fread(&image->hs,sizeof(short),1,fp);
	(void)fread(&image->nx,sizeof(int),1,fp);
	(void)fread(&image->ny,sizeof(int),1,fp);
	(void)fread(&image->nz,sizeof(int),1,fp);
	(void)fread(&image->nt,sizeof(int),1,fp);

	if(image->flg!=1)
	{
		swapByteOrder( (char *)&image->hs, sizeof(short));
		swapByteOrder( (char *)&image->nx, sizeof(int));
		swapByteOrder( (char *)&image->ny, sizeof(int));
		swapByteOrder( (char *)&image->nz, sizeof(int));
		swapByteOrder( (char *)&image->nt, sizeof(int));
	}

	(void)fread(&image->dx,sizeof(double),1,fp); 
	(void)fread(&image->dy,sizeof(double),1,fp); 
	(void)fread(&image->dz,sizeof(double),1,fp); 
	(void)fread(&image->dt,sizeof(double),1,fp); 

	if(image->flg!=1) 
	{
		swapByteOrder( (char *)&image->dt, sizeof(double));
		swapByteOrder( (char *)&image->dz, sizeof(double));
		swapByteOrder( (char *)&image->dx, sizeof(double));
		swapByteOrder( (char *)&image->dy, sizeof(double));
	}

	(void)fread(image->cntrv, sizeof(double), 3, fp);
	(void)fread(image->nrmlv, sizeof(double), 3, fp);
	(void)fread(image->rowv, sizeof(double), 3, fp);
	(void)fread(image->colv, sizeof(double), 3, fp);

	if(image->flg!=1)
	{
		swapByteOrder( (char *)&image->cntrv[0], sizeof(double));
		swapByteOrder( (char *)&image->cntrv[1], sizeof(double));
		swapByteOrder( (char *)&image->cntrv[2], sizeof(double));

		swapByteOrder( (char *)&image->nrmlv[0], sizeof(double));
		swapByteOrder( (char *)&image->nrmlv[1], sizeof(double));
		swapByteOrder( (char *)&image->nrmlv[2], sizeof(double));

		swapByteOrder( (char *)&image->rowv[0], sizeof(double));
		swapByteOrder( (char *)&image->rowv[1], sizeof(double));
		swapByteOrder( (char *)&image->rowv[2], sizeof(double));

		swapByteOrder( (char *)&image->colv[0], sizeof(double));
		swapByteOrder( (char *)&image->colv[1], sizeof(double));
		swapByteOrder( (char *)&image->colv[2], sizeof(double));
	}

	(void)fread(&image->datatype, sizeof(char), 1, fp);
	(void)fread(&image->u, sizeof(char), 1, fp);
	(void)fread(&image->c, sizeof(char), 1, fp);

	//printf("Center Vector = (%7.5lf %7.5lf %7.5lf)\n", image->cntrv[0],image->cntrv[1], image->cntrv[2]);
	//printf("Normal Vector = (%7.5lf %7.5lf %7.5lf)\n", image->nrmlv[0],image->nrmlv[1], image->nrmlv[2]);
	//printf("Row Vector = (%7.5lf %7.5lf %7.5lf)\n", image->rowv[0],image->rowv[1], image->rowv[2]);
	//printf("Column Vector = (%7.5lf %7.5lf %7.5lf)\n", image->colv[0],image->colv[1], image->colv[2]);

	if( fseek(fp, image->hs, SEEK_SET) != 0 )
	{
		printf("\n\nError: improper seek operation!\n\n");
		return(0);
	}

	nv = image->nx*image->ny*image->nz*image->nt;

	image->cim = NULL;
	image->sim = NULL;
	image->iim = NULL;
	image->fim = NULL;
	image->dim = NULL;

	switch (image->datatype) {
            case 0:
				image->cim =(char *)calloc(nv,sizeof(char));
				if(image->cim == NULL) { printf("\n\nError: memory allocation.\n\n"); return(0);}
				nv=fread(image->cim, sizeof(char), nv, fp);
                break;
            case 1:
				image->sim =(short *)calloc(nv,sizeof(short));
				if(image->sim == NULL) { printf("\n\nError: memory allocation.\n\n"); return(0);}
				nv=fread(image->sim, sizeof(short), nv, fp);
				if(image->flg!=1) for(int i=0; i<nv; i++) swapByteOrder( (char *)&image->sim[i], sizeof(short));
                break;
            case 2:
				image->iim =(int *)calloc(nv,sizeof(int));
				if(image->iim == NULL) { printf("\n\nError: memory allocation.\n\n"); return(0);}
				nv=fread(image->iim, sizeof(int), nv, fp);
				if(image->flg!=1) for(int i=0; i<nv; i++) swapByteOrder( (char *)&image->iim[i], sizeof(int));
                break;
            case 3:
				image->fim =(float *)calloc(nv,sizeof(float));
				if(image->fim == NULL) { printf("\n\nError: memory allocation.\n\n"); return(0);}
				nv=fread(image->fim, sizeof(float), nv, fp);
				if(image->flg!=1) for(int i=0; i<nv; i++) swapByteOrder( (char *)&image->fim[i], sizeof(float));
                break;
            case 4:
				image->dim =(double *)calloc(nv,sizeof(double));
				if(image->dim == NULL) { printf("\n\nError: memory allocation.\n\n"); return(0);}
				nv=fread(image->dim, sizeof(double), nv, fp);
				if(image->flg!=1) for(int i=0; i<nv; i++) swapByteOrder( (char *)&image->dim[i], sizeof(double));
                break;
	}

	fclose(fp);

	if( nv != image->nx*image->ny*image->nz*image->nt )
	{       
                printf("\n\nProblem reading %s!\n\n",filename);
                return(0);
	}


	//printf("Matrix size = %d x %d x %d (voxels)\n",image->nx, image->ny,image->nz);
	//printf("Voxel size = %8.6f x %8.6f x %8.6f (mm)\n", image->dx,image->dy,image->dz);


	return(1);
}

// returns 0 on failure, 1 on success
int saveNKI(char *filename, nki image)
{
	FILE *fp;
	int nv;

	fp = fopen(filename,"w");

	if(fp==NULL)
	{ 
		printf("\n\nError openning %s.\n\n",filename);
		return(0);
	} 

	fwrite("NKI", sizeof(char), 3, fp);

	fwrite(&image.flg, sizeof(short), 1, fp);
	fwrite(&image.hs, sizeof(short), 1, fp);
	fwrite(&image.nx, sizeof(int), 1, fp); 
	fwrite(&image.ny, sizeof(int), 1, fp);
	fwrite(&image.nz, sizeof(int), 1, fp);
	fwrite(&image.nt, sizeof(int), 1, fp);
	fwrite(&image.dx, sizeof(double), 1, fp);
	fwrite(&image.dy, sizeof(double), 1, fp);
	fwrite(&image.dz, sizeof(double), 1, fp);
	fwrite(&image.dt, sizeof(double), 1, fp);
	fwrite(image.cntrv, sizeof(double), 3, fp);
	fwrite(image.nrmlv, sizeof(double), 3, fp);
	fwrite(image.rowv, sizeof(double), 3, fp);
	fwrite(image.colv, sizeof(double), 3, fp);
	fwrite(&image.datatype, sizeof(char), 1, fp);
	fwrite(&image.u, sizeof(char), 1, fp);
	fwrite(&image.c, sizeof(char), 1, fp);
	fseek(fp,image.hs,SEEK_SET);

	nv = image.nx * image.ny * image.nz * image.nt;

	switch (image.datatype) {
            case 0:
				nv=fwrite(image.cim, sizeof(char), nv, fp);
                break;
            case 1:
				nv=fwrite(image.sim, sizeof(short), nv, fp);
                break;
            case 2:
				nv=fwrite(image.iim, sizeof(int), nv, fp);
                break;
            case 3:
				nv=fwrite(image.fim, sizeof(float), nv, fp);
                break;
            case 4:
				nv=fwrite(image.dim, sizeof(double), nv, fp);
                break;
	}

	fclose(fp);
		
	if( nv != image.nx*image.ny*image.nz*image.nt )
	{
		printf("\n\nError writing to %s.\n\n",filename);
		return(0);
	}

	return(1);
}
