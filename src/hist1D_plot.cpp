#include <stdio.h>
#include <stdlib.h>
#include <babak_lib.h>

//*********************************************************************
// "hist1D_plot" plots the histogram of the data1 and data2 arrays based 
// on the bins. The program generates 3 files in the running directory 
// in the following formats:  name.mat 
//                            name.plt 
//                            name.png
//
// Inputs: 
//         name:  name of the output files
//         n:     number of the bins
//         bin:   array of the bin numbers
//         data1: array of the data1 values
//         data2: array of the data2 values
//         max_x: x value of the maximum point
//         thld: threshold value
//


//-----------------------------------------------------------------------------
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2)
{
   FILE *fp;
   char filename[DEFAULT_STRING_LENGTH]="";  // a generic filename for reading/writing stuff
   char command[DEFAULT_STRING_LENGTH]="";   // a generic command for running in the bash

   //////////////////////////////////////////////////////////////////////
   //Writing the data file
   {
      snprintf(filename,sizeof(filename),"%s.dat",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      for(int i=0; i<n; i++)
         fprintf(fp,"%d %f %f\n", bin[i], data1[i], data2[i]);

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////
   //Writing the plot file
   {
      snprintf(filename,sizeof(filename),"%s.plt",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      fprintf(fp,"set terminal png\n");
      fprintf(fp,"set output '%s.png'\n",name);

      fprintf(fp,"set style line 1 linecolor rgb '#0000ff' lt 1 lw 1\n");   //Blue
      fprintf(fp,"set style line 2 linecolor rgb '#FF0000' lt 3 lw 2\n");   //Red

      fprintf(fp,"plot \"%s.dat\" using 1:2 title 'Histogram' with lines ls 1,\\\n",name);
      fprintf(fp,"     \"%s.dat\" using 1:3 title 'EM fit' with lines ls 2\n",name);

      fclose(fp);
   }    

   //Running the gnuplot to generate the .png file.
   snprintf(command,sizeof(command),"gnuplot %s.plt",name);
   //system(command);
   (void)system(command);
}



//----------------------------------------------------------------------------------------
//----------------------------   with maximum line  --------------------------------------
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2, int max_x)
{
   FILE *fp;
   char filename[DEFAULT_STRING_LENGTH]="";  // a generic filename for reading/writing stuff
   char command[DEFAULT_STRING_LENGTH]="";   // a generic command for running in the bash

   //////////////////////////////////////////////////////////////////////
   //Writing the data file
   {
      snprintf(filename,sizeof(filename),"%s.dat",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      for(int i=0; i<n; i++)
         fprintf(fp,"%d %f %f\n", bin[i], data1[i], data2[i]);

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////
   //Writing the plot file
   {
      snprintf(filename,sizeof(filename),"%s.plt",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      fprintf(fp,"set terminal png\n");
      fprintf(fp,"set output '%s.png'\n",name);

      fprintf(fp,"set style line 1 linecolor rgb '#0000ff' lt 1 lw 1\n");   //Blue
      fprintf(fp,"set style line 2 linecolor rgb '#FF0000' lt 3 lw 2\n");   //Red

      fprintf(fp,"set arrow from %d,graph(0,0) to %d,graph(1,1) nohead\n", max_x, max_x);
      fprintf(fp,"set label 3 \"Max\" at %d,graph(1,0.95) right\n", max_x);

      fprintf(fp,"plot \"%s.dat\" using 1:2 title 'Histogram' with lines ls 1,\\\n",name);
      fprintf(fp,"     \"%s.dat\" using 1:3 title 'EM fit' with lines ls 2\n",name);

      fclose(fp);
   }    

   //Running the gnuplot to generate the .png file.
   snprintf(command,sizeof(command),"gnuplot %s.plt",name);
   //system(command);
   (void)system(command);
}

//-------------------------------------------------------------------------------------------------
//-----------------   with maximum and threshold lines  ------------------------------------------
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2, int max_x, int thld)
{
   FILE *fp;
   char filename[DEFAULT_STRING_LENGTH]="";  // a generic filename for reading/writing stuff
   char command[DEFAULT_STRING_LENGTH]="";   // a generic command for running in the bash

   //////////////////////////////////////////////////////////////////////
   //Writing the data file
   {
      snprintf(filename,sizeof(filename),"%s.dat",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      for(int i=0; i<n; i++)
         fprintf(fp,"%d %f %f\n", bin[i], data1[i], data2[i]);

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////
   //Writing the plot file
   {
      snprintf(filename,sizeof(filename),"%s.plt",name);
      fp = fopen(filename,"w");
      if(fp==NULL) file_open_error(filename);

      fprintf(fp,"set terminal png\n");
      fprintf(fp,"set output '%s.png'\n",name);

      fprintf(fp,"set style line 1 linecolor rgb '#0000ff' lt 1 lw 1\n");   //Blue
      fprintf(fp,"set style line 2 linecolor rgb '#FF0000' lt 3 lw 2\n");   //Red

      fprintf(fp,"set arrow 3 from %d,graph(0,0) to %d,graph(1,1) nohead \n", max_x, max_x);
      fprintf(fp,"set label 3 \"Max\" at %d,graph(1,0.95) right\n", max_x);

      fprintf(fp,"set arrow 4 from %d,graph(0,0) to %d,graph(1,1) nohead\n", thld, thld);
      fprintf(fp,"set label 4 \"Thld\" at %d,graph(1,0.95) right \n", thld);

      fprintf(fp,"plot \"%s.dat\" using 1:2 title 'Histogram' with lines ls 1,\\\n",name);
      fprintf(fp,"     \"%s.dat\" using 1:3 title 'EM fit' with lines ls 2\n",name);

      fclose(fp);
   }    

   //Running the gnuplot to generate the .png file.
   snprintf(command,sizeof(command),"gnuplot %s.plt",name);
   //system(command);
   (void)system(command);

   // added back to remove these files afterwards
   snprintf(filename,sizeof(filename),"%s.dat",name);
   //remove(filename);
   snprintf(filename,sizeof(filename),"%s.plt",name);
   //remove(filename);
}
