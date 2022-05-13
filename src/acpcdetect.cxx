#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <volume.h>
#include <ctype.h>
#include <nifti1_io.h>
#include <niftiimage.h>
#include <babak_lib.h>
#include <minmax.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
   {"-nx",1,'x'},
   {"-ny",1,'y'},
   {"-nz",1,'z'},
   {"-dx",1,'X'},
   {"-dy",1,'Y'},
   {"-dz",1,'Z'},
   {"-i",1,'i'},
   {"-v",0,'v'},
   {"-center-AC", 0, 'M'},
   {"-output-orient",1,'u'},
   {"-oo",1,'u'},
   {"-nn",0,'n'}, // does nearest neighbor interpolation
   {"-standard", 0, 'S'},
   {"-s", 0, 'S'},
   {"-no-tilt-correction", 0, 'R'},
   {"-lm", 1, 'L'},
//   {"-sform",0,'s'},
//   {"-qform",0,'q'},
   {"-noppm",0,'N'},
   {"-nopng",0,'g'},
   {"-notxt",0,'t'},
   {"-version",0,'V'},
//   {"-T2",0,'T'},
//   {"-T2",0,'T'},
   {"-rvsps",1,'0'},
   {"-rac",1,'1'},
   {"-rpc",1,'2'},
//   {"-m",1,'m'},
//   {"-model",1,'m'},
//   {"-model",1,'m'},
   {"-input-orient",1,'O'}, // secret option
   {"-help",0,'h'},
   {0,0,0}
};

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
  printf(
  "\nUsage: acpcdetect [options] -i <input-volume>.nii\n\n"

  "-i <input-volume>.nii\n"
  "\t3D T1W MRI brain volume in NIFTI format of type short or unsigned short\n\n"

  "Options:\n\n"

  "-v\n"
  "\tEnables verbose mode\n\n"

  "-lm <landmarks-file>\n"
  "\tA text file containing manually determined (i, j, k) coordinates\n"
  "\tof the AC, PC and VSPS, respectively. When this file is supplied,\n"
  "\tautomatic detection of these landmarks is suppressed. This is useful\n" 
  "\tin cases when automatic landmark detection fails.\n\n"

  "-no-tilt-correction\n"
  "\tDoes not tilt-correct the output, but the SFORM and QFORM are set\n"
  "\tcorrectly in the output volume header. This is useful for applications\n"
  "\tthat would like to use acpcdetect as a preprocessing tilt-correction\n"
  "\tstep without applying interpolation at this stage.\n\n"

  "-center-AC\n"
  "\tPlaces the output volume's FOV center at AC\n\n"

  "-standard\n"
  "\tTilt-correction is performed without using the Orion landmarks.\n"
  "\tTilt-correction is done using the AC, PC and MSP only. This is the\n"
  "\tmethod used in version 1.0 of acpcdetect.  In the current version, the\n"
  "\t8 orion landmarks are also used to stabilize the standardization of the\n"
  "\torientation. Using this option, therefore, reverts back to the method\n"
  "\tof version 1.0 without using the additional Orion landmarks.\n\n"

  "-output-orient <orientation-code>\n"
  "\tSpecifies the orientation of the output volume (default: RAS).\n"
  "\tIn ART, orientation codes are 3-letter codes consisting of 6 letters:\n"
  "\tA, P, I, S, L, R.  There are 48 possible combinations. For example\n"
  "\tPIL for Posterior-Inferior-Left or RAS for Right-Anterior-Superior.\n\n"

  "-nx <int>\n"
  "\tNumber of voxels in i direction (the fastest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-ny <int>\n"
  "\tNumber of voxels in j direction (the 2nd fastest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-nz <int>\n"
  "\tNumber of voxels in k direction (the slowest varying index) of the\n"
  "\toutput volume. The default value is determined from the input volume.\n\n"

  "-dx <float>\n"
  "\tVoxel dimension of the output volume in i direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-dy <float>\n"
  "\tVoxel dimension of the output volume in j direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-dz <float>\n"
  "\tVoxel dimension of the output volume in k direction. The default\n"
  "\tvalue is determined from the input volume.\n\n"

  "-version\n"
  "\tPrints software version\n\n"

  "-help\n"
  "\tPrints help information\n\n"

  "-noppm\n"
  "\tPrevents outputting *.ppm images\n\n"

  "-nopng\n"
  "\tPrevents outputting *.png images\n\n"

  "-notxt\n"
  "\tPrevents outputting *.txt files\n\n"

  "-rvsps <r>\n"
  "\tSearch radius for VSPS (default = 50 mm)\n\n"

  "-rac <r>\n"
  "\tSearch radius for AC (default = 15 mm)\n\n"

  "-rpc <r>\n"
  "\tSearch radius for PC (default = 15 mm)\n\n"

  "-nn\n"
  "\tUses the nearest neighbor interpolation for tilt-correction.\n\n"

  "Outputs:\n\n"
  "<output-volume>.nii\n"
  "\tWhere the output volume is saved. The default filename is\n"
  "\t<input-volume>_<output-orientation-code> (default\n"
  "\t<output-orientation-code> is RAS). This volume will be the tilt-corrected\n"
  "\tversion of the input volume. However, if the -no-tilt-correction option is\n"
  "\tselected, the output volume will not be resliced (only reoriented). The\n"
  "\ttilt-correction information, however, are still written in the QFORM and\n"
  "\tSFORM entries of the image header as well as in the *.mrx and *.mat files\n"
  "\t(described below).\n\n"

  "<input-volume>.mrx\n"
  "\tTransformation matrix for tilt-correction in ART format\n\n"

  "<input-volume>_FSL.mat\n"
  "\tTransformation matrix for tilt-correction in FSL format\n\n"

  "<input-volume>_ACPC_sagittal.ppm\n"
  "\tSagittal view of the detected AC/PC locations in\n"
  "\tPPM format (output suppressed by -noppm option)\n\n"

  "<input-volume>_ACPC_sagittal.png\n"
  "\tSagittal view of the detected AC/PC locations in\n"
  "\tPNG format (output suppressed by -nopng option)\n\n"

  "<input-volume>_ACPC_axial.ppm\n"
  "\tAxial view of the detected AC/PC locations in PPM\n"
  "\tformat (output suppressed by -noppm option)\n\n"

  "<input-volume>_ACPC_axial.png\n"
  "\tAxial view of the detected AC/PC locations in PNG\n"
  "\tformat (output suppressed by -nopng option)\n\n"

  "<input-volume>_orion.ppm\n"
  "\tMid-sagittal view of the detected Orion landmarks in\n"
  "\tPPM format (output suppressed by -noppm option)\n\n"

  "<input-volume>_orion.png\n"
  "\tMid-sagittal view of the detected Orion landmarks in\n"
  "\tPNG format (output suppressed by -nopng option)\n\n"

  "<input-volume>_ACPC.txt\n"
  "\tStores the detected AC, PC and VSPS (i, j, k) coordinates and the\n"
  "\testimated mid-sagittal plane (output suppressed by -notxt option)\n\n"

  "<input-volume>_orion.txt\n"
  "\tStores (i, j, k) coordinates of the 8 detected Orion\n"
  "\tlandmarks (output suppressed by -notxt option)\n\n"
  );

  exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void computeSiemensVisionOffsets(float *Tmsp, float *AC, float *PC)
{
	double pi;
  double beta;
  double alpha;
	float *invT;
	float ac[4], pc[4];
	float pcac[4];
	float n[4];
	float N[4];
	float TLAI[16];
	float Rx[16], Ry[16], R[16];

	pi=4.0*atan(1.0);

	invT=inv4(Tmsp);

	// assumes original is PIL, N[3]=0.0 is not a mistake, in fact very necessary
	N[0]=0.0; N[1]=0.0; N[2]=1.0; N[3]=0.0; 
	multi(invT,4,4,N,4,1,N); 
	normalizeVector(N,3);

	multi(invT, 4, 4,  AC, 4,  1, ac); // assumes original is PIL
	multi(invT, 4, 4,  PC, 4,  1, pc); // assumes original is PIL

	free(invT);

	pcac[0] = ac[0] - pc[0];
	pcac[1] = ac[1] - pc[1];
	pcac[2] = ac[2] - pc[2];

	normalizeVector(pcac,3);

	crossProduct(pcac,N,n);
	n[3]=1.0;

	inversePILtransform("LAI", TLAI);
	multi(TLAI, 4, 4, n, 4, 1, n);

   	beta = asin( -(float)n[0] );
   	alpha = asin( -n[1]/cos(beta) );

   	beta *= 180/pi;
   	alpha *= 180/pi;

	Rx[0]=1.0; Rx[1]=0.0; Rx[2]=0.0; Rx[3]=0.0;
	Rx[4]=0.0; Rx[5]=(float)cos(alpha*pi/180); Rx[6]=-(float)sin(alpha*pi/180); Rx[7]=0.0;
	Rx[8]=0.0; Rx[9]=(float)sin(alpha*pi/180); Rx[10]=(float)cos(alpha*pi/180); Rx[11]=0.0;
	Rx[12]=0.0; Rx[13]=0.0; Rx[14]=0.0; Rx[15]=1.0;

	Ry[0]=(float)cos(beta*pi/180); Ry[1]=0.0;	Ry[2]=-(float)sin(beta*pi/180); Ry[3]=0.0;
	Ry[4]=0.0; Ry[5]=1.0; Ry[6]=0.0; Ry[7]=0.0;
	Ry[8]=(float)sin(beta*pi/180); Ry[9]=0.0; Ry[10]=(float)cos(beta*pi/180); Ry[11]=0.0;
	Ry[12]=0.0; Ry[13]=0.0; Ry[14]=0.0; Ry[15]=1.0;

	multi(Rx,4,4,Ry,4,4,R);

	multi(TLAI, 4, 4,  ac, 4,  1, ac);
	multi(TLAI, 4, 4,  pc, 4,  1, pc);

	invT=inv4(R);
	multi(invT,4,4,ac,4,1,ac);
	multi(invT,4,4,pc,4,1,pc);
	free(invT);

	printf("\nSiemens Vision FOV offsets (assuming PIL orientation of structural scan):\n");
   	printf("Trans. to coronal = %lf deg.\nto sagittal = %lf deg.\n", alpha, beta);
	printf("Shift = %f mm\n",ac[2]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  int opt_nn=NO;
  float PIL2RAS[16];
  float IPORIENT2PIL[16];
  float PIL2IPORIENT[16];
  float OPORIENT2PIL[16];
  float PIL2OPORIENT[16];
  float Tijk2xyz[16];
  float *invT=NULL;
  short *opimage;
  short *ipimage;
  char opt_tiltcorrect=YES;
  int nPA=0, nLR=0, nIS=0;
  float dPA=0.0, dLR=0.0, dIS=0.0;
  char opt_standard=NO;
  FILE *fp;
  DIM ipdim, opdim;  // input (ip) and output (op) dimension structures (dim)
  opdim.dx = opdim.dy = opdim.dz = 0.0;
  opdim.nx = opdim.ny = opdim.nz = 0;

  nifti_1_header iphdr, ophdr; // 348 bytes
//  int opt_T2=NO;
  char optransformationpath[1024]="";
  char landmarksfilepath[512]="";
//  char modelfile[1024]="";
  char ipimagepath[1024]="";
  char ipimagename[512]="";
  char ipimagedir[512]="";  // important to initialize to ""
  char opimagepath[1024]="";
  char iporient[4]="";
  char oporient[4]="RAS"; // default output orientation is RAS
  float Tout[16]; // transforms ipimage to the specified output orientation
  float TPIL[16]; // transforms ipimage to PIL orientation

  // It is very important to have these initializations.
  int nx=0, ny=0, nz=0;
  float dx=0.0, dy=0.0, dz=0.0;

  // opt_CENTER_AC=NO means that by default the mid-point between AC and PC is set to the FOV
  // center.  If -centerAC is selected, then the AC is made the FOV center.
  // The opt_CENTER_AC variable is defined in PILtransform.cpp and made global in babak_lib.h
  opt_CENTER_AC=NO; 

  while ((opt = getoption(argc, argv, options)) != -1 )
  {
      switch (opt) 
      {
         case 'n':
            opt_nn=YES;
            break;
         case 'M':
            opt_CENTER_AC=YES;
            break;
         case 'S':
            opt_standard=YES;
            break;
         case 't':
            opt_txt = NO;
            break;
         case 'g':
            opt_png = NO;
            break;
         case 'N':
            opt_ppm = NO;
            break;
//         case 's':
//            opt_sform = YES;
//            break;
//         case 'q':
//            opt_qform = YES;
//            break;
         case 'V':
            printf("acpcdetect v2.0, April 19, 2018 release\n");
            exit(0);
         case 'x':
            nx = atoi(optarg);
            break;
         case 'y':
            ny = atoi(optarg);
            break;
         case 'z':
            nz = atoi(optarg);
            break;
         case 'X':
            dx = atof(optarg);
            break;
         case 'Y':
            dy = atof(optarg);
            break;
         case 'Z':
            dz = atof(optarg);
            break;
//         case 'T':
//            opt_T2=YES;
//            break;
         case 'R':
            opt_tiltcorrect=NO;
            break;
         case 'L':
            sprintf(landmarksfilepath,"%s",optarg);
            break;
         case 'i':
            sprintf(ipimagepath,"%s",optarg);
            break;
         case '0':
            searchradius[0] = atof(optarg); // searchradius[0] is for VSPS
            break;
         case '1':
            searchradius[1] = atof(optarg); // searchradius[1] is for AC
            break;
         case '2':
            searchradius[2] = atof(optarg); // searchradius[2] is for PC
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'O':
            sprintf(iporient,"%s",optarg);
            iporient[0]=(char)toupper((int)iporient[0]);
            iporient[1]=(char)toupper((int)iporient[1]);
            iporient[2]=(char)toupper((int)iporient[2]);
            break;
         case 'u':
            sprintf(oporient,"%s",optarg);
            oporient[0]=(char)toupper((int)oporient[0]);
            oporient[1]=(char)toupper((int)oporient[1]);
            oporient[2]=(char)toupper((int)oporient[2]);
            break;
//         case 'm':
//            sprintf(modelfile,"%s",optarg);
//            break;
         case '?':
            print_help_and_exit();
      }
  }

  //////////////////////////////////////////////////////////////////////////////
  // This block of code sets: ipimagepath, ipimagename, ipimagedir
  // iporient, ipimage, iphdr, ipdim, nPA, nIS, nLR, dPA, dIS, dLR 
  // IPORIENT2PIL, PIL2IPORIENT

  if( ipimagepath[0]=='\0' )
  {
    printf("Please specify an input image using: -i <inputimage.nii>\n");
    exit(0);
  }
  if(opt_v) printf("Input image: %s\n",ipimagepath);

  // determine input image filename without the .nii suffix
  if( niftiFilename(ipimagename, ipimagepath)==0 ) { exit(1); }

  // determine input image directory
  getDirectoryName(ipimagepath, ipimagedir);

  // if input orientation is specified using -input-orient, make sure it's valid
  // -input-orient overrides the orientation inferred from the image header
  if(iporient[0]!='\0' && isOrientationCodeValid(iporient)==0)
  {
    printf("%s is not a valid orientation code, aborting ...\n",iporient);
    exit(0);
  }

  // If input orientation is not specified using -input-orient option, 
  // read it from image header. This is almost always going to be the case.
  if(iporient[0]=='\0')
  {
    getNiftiImageOrientation(ipimagepath, iporient);
  }

  if(opt_v) printf("Input image orientation: %s\n",iporient);

  ipimage = (short  *)read_nifti_image(ipimagepath, &iphdr);
  if(ipimage==NULL)
  {
    printf("Error reading %s, aborting ...\n", ipimagepath);
    exit(1);
  }

  if( iphdr.datatype != DT_SIGNED_SHORT && iphdr.datatype != DT_UINT16)
  {
    printf("This program can only handle short and unsigned short data types, aborting ...\n");
    exit(0);
  }

  set_dim(ipdim, iphdr);

  if(opt_v) 
  {
    printf("Input image matrix size: %d x %d x %d\n",ipdim.nx,ipdim.ny,ipdim.nz);
    printf("Input image voxel size: %6.4f x %6.4f x %6.4f\n",ipdim.dx,ipdim.dy,ipdim.dz);
    printf("Input image scl_slope = %f\n",iphdr.scl_slope);
    printf("Input image scl_inter = %f\n",iphdr.scl_inter);
  }

  // set nPA, nIS, nLR, dPA, dIS, dLR depending on iporient
  if(iporient[0]=='P' || iporient[0]=='A') { nPA=ipdim.nx; dPA=ipdim.dx; }
  if(iporient[0]=='I' || iporient[0]=='S') { nIS=ipdim.nx; dIS=ipdim.dx; }
  if(iporient[0]=='L' || iporient[0]=='R') { nLR=ipdim.nx; dLR=ipdim.dx; }
  if(iporient[1]=='P' || iporient[1]=='A') { nPA=ipdim.ny; dPA=ipdim.dy; }
  if(iporient[1]=='I' || iporient[1]=='S') { nIS=ipdim.ny; dIS=ipdim.dy; }
  if(iporient[1]=='L' || iporient[1]=='R') { nLR=ipdim.ny; dLR=ipdim.dy; }
  if(iporient[2]=='P' || iporient[2]=='A') { nPA=ipdim.nz; dPA=ipdim.dz; }
  if(iporient[2]=='I' || iporient[2]=='S') { nIS=ipdim.nz; dIS=ipdim.dz; }
  if(iporient[2]=='L' || iporient[2]=='R') { nLR=ipdim.nz; dLR=ipdim.dz; }

  PILtransform(iporient, IPORIENT2PIL);
  inversePILtransform(iporient, PIL2IPORIENT);
  ///////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  // This block of code sets: oporient, opdim, opimagepath, PIL2OPORIENT
  // OPORIENT2PIL, Tijk2xyz, ophdr

  // if ouput orientation is specified using -oo option, make sure it's valid
  if(isOrientationCodeValid(oporient)==0 )
  {
    printf("%s is not a valid orientation code, aborting ...\n",oporient);
    exit(0);
  }

  // setting opdim
  if(oporient[0]=='P' || oporient[0]=='A') { opdim.nx=nPA; opdim.dx=dPA; }
  if(oporient[0]=='I' || oporient[0]=='S') { opdim.nx=nIS; opdim.dx=dIS; }
  if(oporient[0]=='L' || oporient[0]=='R') { opdim.nx=nLR; opdim.dx=dLR; }
  if(oporient[1]=='P' || oporient[1]=='A') { opdim.ny=nPA; opdim.dy=dPA; }
  if(oporient[1]=='I' || oporient[1]=='S') { opdim.ny=nIS; opdim.dy=dIS; }
  if(oporient[1]=='L' || oporient[1]=='R') { opdim.ny=nLR; opdim.dy=dLR; }
  if(oporient[2]=='P' || oporient[2]=='A') { opdim.nz=nPA; opdim.dz=dPA; }
  if(oporient[2]=='I' || oporient[2]=='S') { opdim.nz=nIS; opdim.dz=dIS; }
  if(oporient[2]=='L' || oporient[2]=='R') { opdim.nz=nLR; opdim.dz=dLR; }
  if(nx > 0) opdim.nx=nx; 
  if(ny > 0) opdim.ny=ny; 
  if(nz > 0) opdim.nz=nz;
  if(dx > 0.0) opdim.dx=dx; 
  if(dy > 0.0) opdim.dy=dy; 
  if(dz > 0.0) opdim.dz=dz;
  opdim.nt=1; 
  opdim.dt=0.0; 
  opdim.np=opdim.nx*opdim.ny; 
  opdim.nv=opdim.np*opdim.nz; 

  sprintf(opimagepath,"%s/%s_%s.nii",ipimagedir,ipimagename,oporient);

  if(opt_v) 
  {
    printf("Output image: %s\n",opimagepath);
    printf("Output image matrix size: %d x %d x %d\n",opdim.nx,opdim.ny,opdim.nz);
    printf("Output image voxel size: %6.4f x %6.4f x %6.4f\n",opdim.dx,opdim.dy,opdim.dz);
    printf("Output image scl_slope = %f\n",1.0);
    printf("Output image scl_inter = %f\n",0.0);
    printf("Output image orientation: %s\n",oporient);
  }

  // PIL2OPORIENT takes points from PIL space to oporient space
  inversePILtransform(oporient, PIL2OPORIENT);

  // OPORIENT2PIL takes points from oporient space to PIL space
  PILtransform(oporient, OPORIENT2PIL);

  // (i,j,k) -> (x,y,z) in oporient
  ijk2xyz(Tijk2xyz, opdim);

  ophdr = iphdr;
  ophdr.pixdim[1]=opdim.dx; 
  ophdr.pixdim[2]=opdim.dy; 
  ophdr.pixdim[3]=opdim.dz;
  ophdr.dim[0] = 4;
  ophdr.dim[1]=opdim.nx; 
  ophdr.dim[2]=opdim.ny; 
  ophdr.dim[3]=opdim.nz;
  ophdr.dim[4] = 1;
  ophdr.scl_slope = 1.0;
  ophdr.scl_inter = 0.0;
  sprintf(ophdr.descrip,"Created by ART acpcdetect");
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // set TPIL and Tout
  
  if(searchradius[0]<=0.0 || searchradius[0]>200.0) searchradius[0]=50.0;
  if(searchradius[1]<=0.0 || searchradius[1]>100.0) searchradius[1]=15.0;
  if(searchradius[2]<=0.0 || searchradius[2]>100.0) searchradius[2]=15.0;

  if(opt_v && landmarksfilepath[0]!='\0') 
  {
    printf("Manually specified landmarks: %s\n",landmarksfilepath);
  }

  // find TPIL which transforms the input image in iporient to tilt-corrected PIL
  if(opt_standard)
  {
    standard_PIL_transformation(ipimagepath, landmarksfilepath, iporient, 0, TPIL);
  }
  else
  {
    new_PIL_transform(ipimagepath, landmarksfilepath, iporient, TPIL, 0);
  }

  //Tout transforms points from the input iporient to tilt-corrected oporient 
  multi(PIL2OPORIENT, 4, 4,  TPIL, 4,  4, Tout);
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //save Tout and Tout_FSL
  {
    float Tout_FSL[16];

    sprintf(optransformationpath,"%s/%s.mrx",ipimagedir,ipimagename);
    if(opt_v) printf("Output transformation matrix: %s\n",optransformationpath);

    fp = fopen(optransformationpath,"w");
    if(fp==NULL) file_open_error(optransformationpath);
    printMatrix(Tout, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
    fclose(fp);

    art_to_fsl(Tout, Tout_FSL, ipdim, opdim);
    sprintf(optransformationpath,"%s/%s_FSL.mat",ipimagedir,ipimagename);
    if(opt_v) printf("Output transformation matrix (FSL format): %s\n",optransformationpath);
    fp=fopen(optransformationpath,"w");
    if(fp==NULL) file_open_error(optransformationpath);
    printMatrix(Tout_FSL,4,4,"",fp);
    fclose(fp);
  }
  //////////////////////////////////////////////////////////////////////////////

  // PIL2RAS takes (x,y,z) points from PIL to RAS space
  inversePILtransform("RAS", PIL2RAS);

  if( opt_tiltcorrect == YES)
  {
    float Ttmp[16];

    ////////////////////////////////////////////////////////////////////
    // setup ophdr sform and qform

    multi(OPORIENT2PIL, 4, 4,  Tijk2xyz, 4,  4, Ttmp);
    multi(PIL2RAS, 4, 4,  Ttmp, 4,  4, Ttmp);

    update_qsform(ophdr, Ttmp);
    ////////////////////////////////////////////////////////////////////

    invT = inv4(Tout);
  }

  if( opt_tiltcorrect == NO)
  {
    float Ttmp[16];
    float IPORIENT2OPORIENT[16];
    float OPORIENT2IPORIENT[16];

    ////////////////////////////////////////////////////////////////////
    // setup ophdr sform and qform

    multi(PIL2IPORIENT, 4, 4,  OPORIENT2PIL, 4,  4, OPORIENT2IPORIENT);
    multi(OPORIENT2IPORIENT, 4, 4,  Tijk2xyz, 4,  4, Ttmp);
    multi(TPIL, 4, 4,  Ttmp, 4,  4, Ttmp);
    multi(PIL2RAS, 4, 4,  Ttmp, 4,  4, Ttmp); 
    update_qsform(ophdr, Ttmp);
    ////////////////////////////////////////////////////////////////////

    multi(PIL2OPORIENT, 4, 4,  IPORIENT2PIL, 4,  4, IPORIENT2OPORIENT);
    invT = inv4(IPORIENT2OPORIENT);
  }

  if(!opt_nn)
  {
    opimage = resliceImage(ipimage,ipdim,opdim,invT,LIN);
  }
  else
  {
    opimage = resliceImage(ipimage,ipdim,opdim,invT,NEARN);
  }
  free(invT);

  save_nifti_image(opimagepath, opimage, &ophdr);

  delete ipimage;
  delete opimage;
}
