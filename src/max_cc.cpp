#define _max_cc

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <spm_analyze.h>
#include <babak_lib.h>


void Connected_Components( char *im, int nx, int ny, int nz, int *LABEL, int *N, int **Clabel, int **Size);
void adj_top_down(char *im, int *LABEL, int i, int j, int n, int nx, int *L, int **eqtable, int *M);
void adj_bottom_up(int *LABEL, int i, int n, int nx, int *L, int **eqtable, int *M);
void adj_front_back(int *LABEL, int i, int n, int np, int *L, int **eqtable, int *M);
void adj_back_front(int *LABEL, int i, int n, int np, int *L, int **eqtable, int *M);
void eq_table(int *l1, int *l2, int **eqtable, int **M);
void resolve_table(int **eqtable, int M);

void max_Connected_Component(char *im, int nx, int ny, int nz, int *ncc, int *maxsize);
void thr_Connected_Component( char *im, int thr, int nx, int ny, int nz, int *ncc, int *ntcc);
void heightthr_Connected_Component(float *im, float thr, int nx, int ny, int nz, int *ncc, int *ntcc);
void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L);

// Identifies different connected components in the input image and modifies the input image such that,
// all the components are removed from the image except for one component, which has the maximum size
// Also returns ncc (number of connected components) and maxsize (maximum component size) 
void max_Connected_Component( char *im, int nx, int ny, int nz, int *ncc, int *maxsize)
{
	int *labeled_im;
       	int N,L=0,*CClabel,*CCsize;
	int nv;

	nv=nx*ny*nz;	
	labeled_im=(int *)calloc(nv,sizeof(int)); // Allocate memory for the labeled connected component image

	// creates a labeled image 'labeled_im', such that all voxels in the binary image which are zero
	// will remain zero in labeled image and all voxels which are non zero will have a distinct label based on pixel connectivity
	// N: Number of cluster
	// CClable: Lable number of each cluster
	// CCsize: Size of each cluster
	Connected_Components(im,nx,ny,nz,labeled_im,&N,&CClabel,&CCsize);

	*maxsize = 0;
	// Find max cluster size
        for(int i=0; i<N; i++)
        if(CCsize[i] > *maxsize)
	{
		*maxsize = CCsize[i];
		L = CClabel[i]; 
	}

	// Find the connected component with max size
	for(int i=0; i<nv; i++)
	if(im[i])
	{			
		if(labeled_im[i] == L)
		im[i] = 1;
		else im[i] = 0;
	}
	
	*ncc = N;
	free(labeled_im); free(CClabel); free(CCsize); 
}


// Identifies different connected components in the input image and modifies the input image such that,    
// it deletes any components with size lesser than the threshold size
// Also returns ncc (number of connected components) and ntcc (number of connected components with size greater than threshold)
void thr_Connected_Component( char *im, int thr, int nx, int ny, int nz, int *ncc, int *ntcc)                                                                {       
	int *labeled_im;   
	int N,*CClabel,*CCsize;
	int nv;
	nv= nx*ny*nz;   
	labeled_im=(int *)calloc(nv,sizeof(int)); //Allocate memory for the labeled connected component image 
	    
	if(thr <= 0)
	{
		printf("\n Threshold has to be greater than zero\n");
		exit(0);
	}
  
	// creates a labeled image 'labeled_im', such that all voxels in the binary image which are zero
	// will remain zero in labeled image and all voxels which are non zero will have a distinct label based on pixel connectivity
	// N: Number of cluster
	// CClable: Lable number of each cluster
	// CCsize: Size of each cluster
	Connected_Components(im,nx,ny,nz,labeled_im,&N,&CClabel,&CCsize);

        //Find the connected components with size greater than threshold
        for(int i=0; i<nv; i++) im[i]=0;

        *ntcc = 0;
	for(int i=0; i<N; i++)
	if(CCsize[i]>=thr) ++*ntcc;

        int k=0;
        for(int i=0; i<nv; i++)
        {
                if(labeled_im[i])
                {
                        if(labeled_im[i]==CClabel[k] && CCsize[k]>=thr)  im[i]=1;
                        else if(labeled_im[i]<CClabel[N/2])
                        for(int j=0; j<N/2; j++)
                        {
                        	if (labeled_im[i]==CClabel[j] && CCsize[j]>=thr)
                                {
                                	im[i]=1; k=j;
                                        break;
                                }
                        }
                        else
                        for(int j=N/2; j<N; j++)
                        {
                        	if (labeled_im[i]==CClabel[j] && CCsize[j]>=thr)
                                {
                                	im[i]=1; k=j;
                                        break;
                                }
                        }
                }
        }

        *ncc = N;
  	free(labeled_im); free(CClabel); free(CCsize);
}

// Identifies different connected components in the input image and modifies the input image such that,
// it deletes any components with height lesser than the threshold height
// Also returns ncc (number of connected components) and ntcc (number of connected components with height greater than threshold)
void heightthr_Connected_Component(float *im, float thr, int nx, int ny, int nz, int *ncc, int *ntcc)
{
        int *labeled_im;
        char *imtemp;
        int N,*CClabel,*CCsize;
        float *CCheight;
        int nv;

        nv = nx*ny*nz;

        labeled_im=(int *)calloc(nv,sizeof(int)); // Allocate memory for the labeled connected component image
        imtemp=(char *)calloc(nv,sizeof(char));

        for(int i=0; i<nv; i++)
	if(im[i]) imtemp[i]=1;

        // creates a labeled image 'labeled_im', such that all voxels in the binary image which are zero
        // will remain zero in labeled image and all voxels which are non zero will have a distinct label based on pixel connectivity
        // N: Number of cluster
        // CClable: Lable number of each cluster
        // CCsize: Size of each cluster
        Connected_Components(imtemp,nx,ny,nz,labeled_im,&N,&CClabel,&CCsize);

        // Find the height of each connected component
        
	CCheight = (float *)calloc(N,sizeof(float));

        for(int i=0; i<N; i++)
        CCheight[i]=0.0;

        int k=0;
        for(int i=0; i<nv; i++)
        {
                if(labeled_im[i])
                {
                        if(labeled_im[i]==CClabel[k] && im[i]>CCheight[k])
                        CCheight[k]=im[i];
                        else if(labeled_im[i] < CClabel[N/2])
                        for(int j=0; j<N/2; j++)
                        {
                        	if (labeled_im[i]==CClabel[j] && im[i]>CCheight[j])
                                {
                                	CCheight[j]=im[i]; k=j;
                                        break;
                                }
                        }
                        else
                        for(int j=N/2; j<N; j++)
                        {
                        	if (labeled_im[i]==CClabel[j] && im[i]>CCheight[j])
                                {
                                	CCheight[j]=im[i]; k=j;
                                        break;
                                }
                        }
                }
        }

 	// Find the connected components with height greater than threshold
        for(int i=0; i<nv; i++)
        im[i]=0.0;

        *ntcc = 0;
        for(int i=0; i<N; i++)
	if(CCheight[i]>=thr) ++*ntcc;

        k=0;
        for(int i=0; i<nv; i++)
        {
                if(labeled_im[i])
                {
                        if(labeled_im[i]==CClabel[k] && CCheight[k]>=thr)
                        im[i]=1.0;
                        else if(labeled_im[i] < CClabel[N/2])
                        for(int j=0; j<N/2; j++)
                        {
                        	if (labeled_im[i]==CClabel[j] && CCheight[j]>=thr)
                                {
                                	im[i]=1.0; k=j;
                                        break;
                                }
                       	}
                        else
                        for(int j=N/2; j<N; j++)
                        {
                                if (labeled_im[i]==CClabel[j] && CCheight[j]>=thr)
                                {
                                       im[i]=1.0; k=j;
                                       break;
                                }
                        }
                }
        }

        *ncc = N;
        free(labeled_im); free(CClabel); free(CCsize); free(CCheight); free(imtemp);
}

// Identifies connected components in input image and stores the location of voxels in each connected component in array L.
// It also computes array S, where S[i] points to the location in L from where the index of voxels in ith connected component are stored
// Also returns ncc (number of connected components)
void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L)
{
        int *labeled_im;
        int *CClabel,*CCsize;
        int nv;
        int N;
        int nbv;                                                // Number of non-zero voxels in input image
        char *imtemp;

        nv = nx*ny*nz;

        labeled_im=(int *)calloc(nv,sizeof(int));           // Allocate memory for the labeled connected component image

        // creates a labeled image 'labeled_im', such that all voxels in the binary image which are zero
        // will remain zero in labeled image and all voxels which are non zero will have a distinct label based on pixel connectivity
        // N: Number of cluster
        // CClable: Lable number of each cluster
        // CCsize: Size of each cluster
        Connected_Components(im,nx,ny,nz,labeled_im,&N,&CClabel,&CCsize);

        nbv=0;
        for(int i=0; i<N; i++) nbv += CCsize[i];

        (*L)=(int *)calloc(nbv,sizeof(int));                    // Allocate memory to store the voxel index
        (*S)=(int *)calloc(N,sizeof(int));
        imtemp=(char *)calloc(nv,sizeof(char));

        for(int i=0; i<nv; i++) imtemp[i]=im[i];

        (*S)[0]=0;
        for(int i=1; i<N; i++) (*S)[i] = (*S)[i-1] + CCsize[i-1];

        // Store the location of voxels in each connected component
        int ii=0;
        for(int i=0; i<N; i++)
        {
                for(int j=0; j<nv; j++)
                {
                        if(imtemp[j])
                        {
                                if(CClabel[i]==labeled_im[j])
                                {
                                	(*L)[ii++]=j;
                                      	imtemp[j]=0;
                                }
                        }
                }
        }
        *ncc = N;

        free(labeled_im); free(CClabel); free(CCsize); free(imtemp);
}

// Identifies connected components in input image and stores the location of voxels in each connected component in array L.
// It also computes array S, where S[i] points to the location in L from where the index of voxels in ith connected component are stored
// Also returns ncc (number of connected components)
void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L, int **CCsize)
{
        int *labeled_im;
        int *CClabel;
        int nv;
        int N;
        int nbv;                                                // Number of non-zero voxels in input image
        char *imtemp;

        nv = nx*ny*nz;

        labeled_im=(int *)calloc(nv,sizeof(int));           // Allocate memory for the labeled connected component image

        // creates a labeled image 'labeled_im', such that all voxels in the binary image which are zero
        // will remain zero in labeled image and all voxels which are non zero will have a distinct label based on pixel connectivity
        // N: Number of cluster
        // CClable: Lable number of each cluster
        // CCsize: Size of each cluster
        Connected_Components(im,nx,ny,nz,labeled_im,&N,&CClabel,CCsize);

        nbv=0;
        for(int i=0; i<N; i++) nbv += (*CCsize)[i];

        (*L)=(int *)calloc(nbv,sizeof(int));                    // Allocate memory to store the voxel index
        (*S)=(int *)calloc(N,sizeof(int));
        imtemp=(char *)calloc(nv,sizeof(char));

        for(int i=0; i<nv; i++) imtemp[i]=im[i];

        (*S)[0]=0;
        for(int i=1; i<N; i++) (*S)[i] = (*S)[i-1] + (*CCsize)[i-1];

        // Store the location of voxels in each connected component
        int ii=0;
        for(int i=0; i<N; i++)
        {
                for(int j=0; j<nv; j++)
                {
                        if(imtemp[j])
                        {
                                if(CClabel[i]==labeled_im[j])
                                {
                                	(*L)[ii++]=j;
                                      	imtemp[j]=0;
                                }
                        }
                }
        }
        *ncc = N;

        free(labeled_im); free(CClabel); free(imtemp);
}


// Computes image 'LABEL', a labeled image, such that voxels, which are non zero in input image 'im' will be
// assigned a distinct label based on pixel connectivity and voxels which are '0' in input image will remain '0' in labeled image.
// It also returns
// N - number of components
// Clable(Nx1) - Label number of each component
// Size(Nx1) - Size of each component 

void Connected_Components(char *im, int nx, int ny, int nz, int *LABEL,
int *N, int **Clabel, int **Size)
{

   int n,np,nv;
	int NEWLABEL;  // Generates the next chronological label that has not been used
	int L;
	
	int **eqtable; // Label equivalence table
	int M;	       
	//int flag=0;
	
	//float min;
	///////////////////////////////////////////////////////////////////////////////
	// 2-DIMENSIONAL CONNECTED COMPONENTS ON EACH SLICE SEPARATELY 
	
	np=nx*ny;
	nv=np*nz;

   eqtable=(int **)calloc(2,sizeof(int));
   eqtable[0]=(int *)calloc(nx,sizeof(int));
   eqtable[1]=(int *)calloc(nx,sizeof(int));

  	NEWLABEL=0; L=0; 

	for(int  k=0; k<nz ;k++)
	{
		int l=k*np;	
		// PROCESS FROM TOP TO DOWN 
		
	  	for(int j=0 ; j<ny; j++)
	    	{
			int m=j*nx+l;
	      		// Initialize all labels on row j to zero

	      		for(int i=0; i<nx; i++) LABEL[i+m]=0;

			// Pass1 on row j 
	      
			M=0; 
	      		for(int i=0; i<nx; i++)
			{
		  		n=i+m;
		  		if(im[n])
		    		{
					// checks the labels of adjacent voxels (row 'j' and 'j-1') and returns minimum label 'LABEL',
		      			adj_top_down(im,LABEL,i,j,n,nx,&L,eqtable,&M);
   
		      			if(L==0)   LABEL[n]=++NEWLABEL;
		      			else   		  LABEL[n]=L;
		    		}
			}
		
			// Find equivalence classes from 'eqtable' 
			if(M>1)
			resolve_table(eqtable,M);   

			// Pass2 on row j (replace labels with its equivalence class)   
			if(M>=1)
			{		
        			for(int i=0; i<nx; i++)
        			{
                			n=i+m;
                			if(LABEL[n])
                			{
                       				for(int p=1; p<=M; p++)
                       				{
                             				if(LABEL[n]== eqtable[0][p])
                             				{
                                        			LABEL[n]=eqtable[1][p];
                                        			break;
                             				}
                        			}
                			}
         			}		
			}
		}

		// PROCESS FROM BOTTOM TO TOP 
		
		for(int j=ny-2; j>=0; j--)
		{

			int m=j*nx+l;
			//Pass1 on row j

                        M=0;
                        for(int i=0; i<nx; i++)
                        {
                                n=i+m;
                                if(im[n])
                                {
					// checks adjacent pixels (row 'j' and 'j+1') and returns minimum label 'LABEL'
                                        adj_bottom_up(LABEL,i,n,nx,&L,eqtable,&M);
                                        LABEL[n]=L;
                                }
                        }
			
			if(M>1)
			resolve_table(eqtable,M);
			
			// Pass2 on row j (replace labels with its equivalence class)
			if(M>=1)
			{
        			for(int i=0; i<nx; i++)
        			{
                			n=i+ m;
                			if(LABEL[n])
                			{
                       				for(int p=1; p<=M; p++)
                       				{
                             				if(LABEL[n]== eqtable[0][p])
                             				{
                                        			LABEL[n]=eqtable[1][p];
                                        			break;
                             				}
                        			}
                			}
         			}
			}
		}	
	}

	free(eqtable);

	////////////////////////////////////////////////////////////////////////////////////////////
	// 3-DIMENSIONAL CONNECTED COMPONENTS      
	
        eqtable=(int **)calloc(2,sizeof(int));
        eqtable[0]=(int *)calloc(np,sizeof(int));
        eqtable[1]=(int *)calloc(np,sizeof(int));
	
	// PROCESS FROM FIRST TO LAST SLICE  
        
	for(int k=1; k<nz ;k++)
        {
		int l=k*np;
		// Pass1 on slice k
		M=0;
                for(int j=0; j<ny; j++)
                {
			int m=j*nx+l;
                        for(int i=0; i<nx; i++)
                        {
                                n=i+m;
                                if(im[n])
                                {
					// Checks adjacent pixels (slice 'k' and 'k-1') and returns minimum label  
                                        adj_front_back(LABEL,i,n,np,&L,eqtable,&M);
                                        LABEL[n]=L;
                                }

                        }
		}

		// Find equivalence classes from 'eqtable'
		if(M>1)	
		resolve_table(eqtable,M);		
	
		// Pass2 on slice k (replace labels with its equivalence class)
                if(M>=1)
		{
			for(int j=0; j<ny; j++)
        		{
				int m=j*nx+l;
                          	for(int i=0;i<nx;i++)
                		{
                        		n=i+m;
                        		if(LABEL[n])
                        		{
                                		for(int p=1; p<=M; p++)
                                		{
                                       			if(LABEL[n]== eqtable[0][p])
                                       			{
                                               			LABEL[n]=eqtable[1][p];
                                               			break;
                                       			}
                                		}
                        		}
                		}
        		}
		}
	}

	// PROCESS FROM LAST TO FIRST SLICE

	for(int k=nz-2; k>=0; k--)
	{
		int l= k*np;
		// Pass1 on slice k
		M=0;
		for(int j=0; j<ny; j++)
		{
			int m=j*nx+l;
			for(int i=0; i<nx;i++)
			{
				n=i+m;
				if(im[n])
				{
					// Checks adjacent pixels (slice 'k' and 'k+1') and returns minimum label
					adj_back_front(LABEL,i,n,np,&L,eqtable,&M);
					LABEL[n]=L;
				}
			}
		}
	  	
                // Find equivalence classes from 'eqtable'

		if(M>1)
                resolve_table(eqtable,M);

		// Pass2 on slice k (replace labels with its equivalence class)

                if(M>=1)
                {
			for(int j=0; j<ny; j++)
        		{
				int m=j*nx+l;
                		for(int i=0;i<nx;i++)
                		{
                        		n=i+m;
                        		if(LABEL[n])
                        		{
                                		for(int p=1; p<=M; p++)
                                		{
                                        		if(LABEL[n]== eqtable[0][p])
                                        		{
                                                		LABEL[n]=eqtable[1][p];
                                                		break;
                                        		}
                                		}
                        		}
                		}
        		}
		}
	}
	
	free(eqtable);

	// At this point image 'label' has distinct labels for each component
	
	//////////////////////////////////////////////////////////////////////////////////////
	// Calculates total number of distinct component

	(*Clabel) = (int *)calloc(NEWLABEL,sizeof(int));
        for( int i=0; i< NEWLABEL; i++) (*Clabel)[i]=0;

	// Finds the first label number in image 'LABEL' 
	for( int i=0; i< (nx*ny*nz); i++)
	{
		if(im[i])
		{
			(*Clabel)[0]=LABEL[i];
			break;
		}
	}
	
	if( (*Clabel)[0]>0)
	{
		*N=1;
 		for(int i=0; i<nx*ny*nz; i++)
        	{
			if(im[i])
                	{
                      		if(LABEL[i]>(*Clabel)[*N-1])
				(*Clabel)[(*N)++]=LABEL[i];
			}
		}
	}
	else
	{
		*N=0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	// Calculates the size of each component 
	
	(*Size) = (int *)calloc(*N,sizeof(int));

	for( int i=0; i<*N; i++) (*Size)[i]=0;

	int k=0;
        for(int i=0; i<nv; i++)
        {
                if(im[i])
                {
                        if(LABEL[i]==(*Clabel)[k]) ++(*Size)[k];
                        else if(LABEL[i] < (*Clabel)[*N/2])
                        {
                                for(int j=0; j<*N/2; j++)
                                {
                                        if((*Clabel)[j]==LABEL[i])
                                        {
                                                ++(*Size)[j];                  // Size of each component
                                                k=j;
                                                break;
                                        }
                                }
                        }
                        else
                        {
                                for(int j=*N/2; j<*N; j++)
                                {
                                        if((*Clabel)[j]==LABEL[i])
                                        {
                                                ++(*Size)[j];                   // Size of each component
                                                k=j;
                                                break;
                                        }
                             	}
                       	}
                }
        }
	//end//
}

// if adjacent voxels are not connected, then it returns '0'. 
// The equivalent labels are also entered in the equivalence table 'eqtable'and
// 'M' is the number of entries in the eqtable

void adj_top_down(char *im, int *LABEL, int i, int j, int n, int nx, int *L, int **eqtable, int *M)
{
	
	if(j==0 && i==0)   *L=0;

	else if(j==0 && i>0)
	{
      		if(im[n-1]) *L=LABEL[n-1];
		else *L=0;
    	}

	else if(i==0 && j>0)
	{
		if(im[n-nx]) *L=LABEL[n-nx];
		else *L=0;
	}
  
	else if(i>0 && j>0)
	{	
		if(LABEL[n-1]==LABEL[n-nx]) *L=LABEL[n-1];

	  	else if(LABEL[n-1]==0 && LABEL[n-nx]>0)  *L=LABEL[n-nx];

                else if(LABEL[n-nx]==0 && LABEL[n-1]>0)  *L=LABEL[n-1];
 
  		else if(LABEL[n-1]>0 && LABEL[n-nx]>0)
    		{
      			if(LABEL[n-nx]>LABEL[n-1])
			{
				eq_table(&LABEL[n-nx],&LABEL[n-1],eqtable,&M);
			  	*L=LABEL[n-1];
			}
			else 
			{ 
				eq_table(&LABEL[n-1],&LABEL[n-nx],eqtable,&M);
				*L=LABEL[n-nx];
			}
		}
	}
}

void adj_bottom_up(int *LABEL, int i, int n, int nx, int *L, int **eqtable, int *M)
{
        if(!i)
        {
		if(LABEL[n+nx]>0 && LABEL[n+nx]<LABEL[n])
		{
			eq_table(&LABEL[n],&LABEL[n+nx],eqtable,&M);
 			*L=LABEL[n+nx];
		}
                else  *L=LABEL[n];
	}
	else if(LABEL[n-1]==0 && LABEL[n+nx]==0)
	{
		*L=LABEL[n];
	}
	else if(LABEL[n-1]==0 && LABEL[n+nx]>0)
	{
		if(LABEL[n+nx] < LABEL[n])
                {
                	eq_table(&LABEL[n],&LABEL[n+nx],eqtable,&M);
			*L=LABEL[n+nx];
		}
		else *L=LABEL[n];
	}
	else if(LABEL[n-1]>0 && LABEL[n+nx]==0)
	{
		if(LABEL[n-1] < LABEL[n])
		{
			eq_table(&LABEL[n],&LABEL[n-1],eqtable,&M);	
			*L=LABEL[n-1];
		}
		else if(LABEL[n-1] > LABEL[n])
		{
			eq_table(&LABEL[n-1],&LABEL[n],eqtable,&M);
			*L=LABEL[n];
		}
		else  *L=LABEL[n];
	}
	else if(LABEL[n-1]>0 && LABEL[n+nx]>0)
	{
		if(LABEL[n+nx]<LABEL[n-1]) 
		{
			eq_table(&LABEL[n-1],&LABEL[n+nx],eqtable,&M);
			*L=LABEL[n+nx];
		}
		else *L=LABEL[n-1];
	}
}

void adj_front_back(int *LABEL, int i, int n, int np, int *L, int **eqtable, int *M)
{
        if(!i)
	{
		if(LABEL[n-np]>0 && LABEL[n-np]!=LABEL[n])
		{
			if(LABEL[n-np]<LABEL[n])
			{
				eq_table(&LABEL[n],&LABEL[n-np],eqtable,&M);
				*L=LABEL[n-np];
			}
			else if(LABEL[n]<LABEL[n-np])
			{
				eq_table(&LABEL[n-np],&LABEL[n],eqtable,&M);
				*L=LABEL[n];
			}
		}
		else 	*L=LABEL[n];
	}
	else if(LABEL[n-1]==0 && LABEL[n-np]==0) 
        {
		*L=LABEL[n];
	}
	else if(LABEL[n-1]==0 && LABEL[n-np]>0)
	{
		if(LABEL[n-np]<LABEL[n])
                {
                        eq_table(&LABEL[n],&LABEL[n-np],eqtable,&M);
                        *L=LABEL[n-np];
               	}
                else if(LABEL[n]<LABEL[n-np])
                {
                        eq_table(&LABEL[n-np],&LABEL[n],eqtable,&M);
                        *L=LABEL[n];
              	}
                else 	*L=LABEL[n];
	}
	else if(LABEL[n-np]==0 && LABEL[n-1]>0)
	{
		 if(LABEL[n-1]<LABEL[n])
                {
                        eq_table(&LABEL[n],&LABEL[n-1],eqtable,&M);
                        *L=LABEL[n-1];
                }
                else if(LABEL[n]<LABEL[n-1])
                {
                        eq_table(&LABEL[n-1],&LABEL[n],eqtable,&M);
                        *L=LABEL[n];
                }
                else *L=LABEL[n];
	}
	else if(LABEL[n-1]>0 && LABEL[n-np]>0)
        {
                if(LABEL[n-np]<LABEL[n-1])
                {
                        eq_table(&LABEL[n-1],&LABEL[n-np],eqtable,&M);
                        *L=LABEL[n-np];
                }
		else if(LABEL[n-1]<LABEL[n-np])
		{
                        eq_table(&LABEL[n-np],&LABEL[n-1],eqtable,&M);
                        *L=LABEL[n-1];
		}
 	        else *L=LABEL[n-1];
        }
}

void adj_back_front(int *LABEL, int i, int n, int np, int *L, int **eqtable, int *M)
{
	if(!i)
        {
                if(LABEL[n+np]>0 && LABEL[n+np]<LABEL[n])
                {
                	eq_table(&LABEL[n],&LABEL[n+np],eqtable,&M);
                        *L=LABEL[n+np];
                }
                else *L=LABEL[n];
        }
	else if(LABEL[n-1]==0 && LABEL[n+np]==0)
        {
                *L=LABEL[n];
        }
        else if(LABEL[n-1]==0 && LABEL[n+np]>0)
        {
                if(LABEL[n+np]<LABEL[n])
                {
                        eq_table(&LABEL[n],&LABEL[n+np],eqtable,&M);
                        *L=LABEL[n+np];
                }
                else *L=LABEL[n];
        }
	else if(LABEL[n+np]==0 && LABEL[n-1]>0)
        {
                if(LABEL[n-1]<LABEL[n])
                {
                        eq_table(&LABEL[n],&LABEL[n-1],eqtable,&M);
                        *L=LABEL[n-1];
                }
                else if(LABEL[n]<LABEL[n-1])
                {
                        eq_table(&LABEL[n-1],&LABEL[n],eqtable,&M);
                        *L=LABEL[n];
                }
                else *L=LABEL[n];
        }
	else if(LABEL[n-1]>0 && LABEL[n+np]>0)
        {
                if(LABEL[n+np]<LABEL[n-1])
                {
                        eq_table(&LABEL[n-1],&LABEL[n+np],eqtable,&M);
                        *L=LABEL[n+np];
                }
                else *L=LABEL[n-1];
        }
}

void eq_table(int *l1, int *l2, int **eqtable, int **M)
{
	int flag=0;

	for(int i=1; i<=**M; i++)
	{
		if(eqtable[0][i]==*l1 && eqtable[1][i]==*l2)
		{
			flag=1;
			break;
		}
	}
	if (!flag)
	{
		++**M;
		eqtable[0][**M]=*l1; eqtable[1][**M]=*l2;	
	}
}

void resolve_table(int **eqtable, int M)
{
        int flag=0;

        for (int i=1; i<=M; i++)
	{      
        	for(int j=1; j<=M; j++)
                {
                        if(eqtable[0][i]==eqtable[0][j] && eqtable[1][i]!=eqtable[1][j])
                        {
                                if(eqtable[1][j]<eqtable[1][i])
                                {
                                        eqtable[0][j]=eqtable[1][i];
                                        eqtable[1][i]=eqtable[1][j];
                                        flag=1;
					--i;
                                        break;
                                }
                                else
                                {
                                        eqtable[0][j]=eqtable[1][j];
                                        eqtable[1][j]=eqtable[1][i];
					flag=1;
					--i;
                                        break;
                                }
                        }
	         }

		if(flag==1 && i==M)
		{
			i=1;
			flag=0;
		}
        }

	flag=0;
	for (int i=1; i<=M; i++)
	{      
                for(int j=1; j<=M; j++)
                {
			if(eqtable[0][i]==eqtable[1][j])
		    	{
		      		eqtable[1][j]=eqtable[1][i];
		      		flag=1;
		    	}
		}
		if(flag==1 && i==M)
		{
			i=1;
			flag=0;
		}
	}	
}
