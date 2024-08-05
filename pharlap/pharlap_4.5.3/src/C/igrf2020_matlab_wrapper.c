/*
 ================================================================
                   igrf2020_matlab_wrapper.c
 ================================================================
  Matlab wrapper for the IGRF subroutines distributed in iri2016
       
  This is a MEX-file for MATLAB.
      
  Change history:
  20/04/23  M. A. Cervera  V1.0  
     Initial version - Based on igrf2016_matlab_wrapper.c (V1.2), 
     updated to use IGRF2016 model from IRI2016

 ================================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int  check_ref_data(char *files_to_check);

void igrf2020_calc_(float *glat, float *glong, float *dec_year, float *height, 
                    float *dipole_moment, float *babs, float *bnorth, 
	   	    float *beast, float *bdown, float *dip, float *dec, 
                    float *dip_lat, float *l_value, int *l_value_code);

int  julday_(int *day, int* month, int *year);
               
/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, ii, jj, idx, idx_mat, arrsize, in_arrsize[6], UT[5];
  int year, month, day, hour, days_in_year, day_of_year, dd1, mm1, dd2, mm2;
  char error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  static int days_in_month[13] = 
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  /* Input parameters to various IGRF routines */
  float glat, glong;              /* Geographic lat/long */
  float dec_year;                 /* decimal year */
  float height;                   /* height range in km */

  /* Outputs from IGRF routines */
  float dipole_moment, bnorth, beast, bdown, babs;
  float dip, dec, dip_lat, l_value;
  int   l_value_code;

  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 5;
  in_arrsize[3] = 1;
 
  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    mexErrMsgTxt("ERROR: 4 inputs required when calling igrf2020.");
    return;
  }

  if (nlhs != 1) {
    mexErrMsgTxt("ERROR: 1 outputs required when calling igrf2020.");
    return;
  }

  /* Check to ensure the inputs are numeric (not strings, NaN or Inf). */
  for (ii = 0; ii < nrhs; ii++) {
    if (mxIsNumeric(prhs[ii]) == 0) {
      snprintf(error_mesg, 200, "ERROR: Input %d must be numeric.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (mxIsNaN(*mxGetPr(prhs[ii])) == 1) {
      snprintf(error_mesg, 200, "ERROR: Input %d must not be NaN.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (mxIsInf(*mxGetPr(prhs[ii])) == 1) {
      snprintf(error_mesg, 200, "ERROR: Input %d must not be Inf.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (mxIsComplex(prhs[ii]) == 1) {
      snprintf(error_mesg, 200, "ERROR: Input %d must not be complex.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Get and check the size of the input arrays. */
  for (ii = 0; ii < nrhs; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Check to see if the reference data files exist. */
  if (check_ref_data("iri2020") == 0) return;

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     ther igrf routines */
  glat = (float) mxGetScalar(prhs[0]);
  glong = (float) mxGetScalar(prhs[1]);
  for (ii = 0; ii < 5; ii++) {
    if (mxIsNaN(*(mxGetPr(prhs[2]) + ii)) == 1) {
      snprintf(error_mesg, 200, "ERROR: Input UT contains a NaN.");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (mxIsInf(*(mxGetPr(prhs[2]) + ii)) == 1) {
      snprintf(error_mesg, 200, "ERROR: Input UT contains an Inf.");
      mexErrMsgTxt(error_mesg);
      return;
    }
    UT[ii] = *(mxGetPr(prhs[2]) + ii);
  }
  height = (float) mxGetScalar(prhs[3]);

  /* calculate decimal year */
  year = UT[0];
  month = UT[1];
  day = UT[2];
  hour = UT[3];      

  dd1 = 0; mm1 = 1; dd2 = 31; mm2 = 12;
  day_of_year  = julday_(&day, &month, &year) - julday_(&dd1, &mm1, &year);
  days_in_year = julday_(&dd2, &mm2, &year) - julday_(&dd1, &mm1, &year);

  dec_year = (float) year + (float) day_of_year / (float) days_in_year + 
             (float) hour / (24 * (float) days_in_year);  

  /* Call the computational subroutines. */
  igrf2020_calc_(&glat, &glong, &dec_year, &height, &dipole_moment, &babs, 
                 &bnorth, &beast, &bdown, &dip, &dec, &dip_lat, &l_value, 
                 &l_value_code);

  /* Create matricies for the return arguments. */
  plhs[0] = mxCreateDoubleMatrix(1, 10, mxREAL);

  /* Load the output data from the igrf routines into the LHS pointers, which 
     are the outputs to MATLAB. Convert Gauss to Tesla. */
  *(mxGetPr(plhs[0]) + 0) = bnorth / 10000.0;
  *(mxGetPr(plhs[0]) + 1) = beast / 10000.0;
  *(mxGetPr(plhs[0]) + 2) = bdown / 10000.0;
  *(mxGetPr(plhs[0]) + 3) = babs / 10000.0;
  *(mxGetPr(plhs[0]) + 4) = dipole_moment;
  *(mxGetPr(plhs[0]) + 5) = l_value;
  *(mxGetPr(plhs[0]) + 6) = l_value_code;
  *(mxGetPr(plhs[0]) + 7) = dip;
  *(mxGetPr(plhs[0]) + 8) = dip_lat;
  *(mxGetPr(plhs[0]) + 9) = dec;

}

