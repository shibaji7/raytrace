/*
  ==============================================================
		irreg_strength_matlab_wrapper.c
  ==============================================================
   Matlab wrapper for irreg_strength_

   This is a MEX-file for MATLAB.

   Change history:
   31/10/05  M.A.Cervera  V1.0  Author
   11/05/07  M.A.Cervera  V2.0  Converted from Fortran
   18/02/10  M.A. Cervera V2.1  Improved error handling

   13/02/22  M. A. Cervera  V2.2 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

  ==============================================================
*/

#include "mex.h"

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

float irreg_strength_(float *glon, float *glat, float *year, int *month,
		      int *day, int *hour, int *minute, float *kp, float *dip,
		      float *dec);
                       

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, arrsize, in_arrsize[4], out_arrsize, ii;
  char error_mesg[200];

  /* Input parameters to irreg_strength_ */
  float glat, glon;                /* Geographic lat/long */
  int   month, day, hour, minute;
  float year, kp;

  /* Output parameters from irreg_strength_ */
  float irregs;
  float dip, dec;

  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 5;
  in_arrsize[3] = 1;
  
  /* Define the expected output array sizes */
  out_arrsize = 1;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    mexErrMsgTxt("ERROR: 4 inputs required when calling irreg_strength.");
    return;
  }
  if (nlhs > 1) {
    mexErrMsgTxt("ERROR: 1 output required when calling irreg_strength.");
    return;
  }

  /* Check to ensure the inputs are numeric (not strings). */
  for (ii = 0; ii < nrhs; ii++) {
    if (mxIsNumeric(prhs[ii]) == 0) {
      snprintf(error_mesg, 200, "ERROR: Input %d must be numeric.", ii+1);
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

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     irreg_strength_ */
  glat = (float) mxGetScalar(prhs[0]);
  glon = (float) mxGetScalar(prhs[1]);
  year = *(mxGetPr(prhs[2]) + 0);
  month = *(mxGetPr(prhs[2]) + 1);
  day = *(mxGetPr(prhs[2]) + 2);
  hour = *(mxGetPr(prhs[2]) + 3);
  minute = *(mxGetPr(prhs[2]) + 4);
  kp = (float) mxGetScalar(prhs[3]);

  /* Create matricies for the return arguments. */
  plhs[0] = mxCreateDoubleMatrix(1, out_arrsize, mxREAL);      
  
  /* Call the computational subroutine. */
  irregs = irreg_strength_(&glon, &glat, &year, &month, &day, &hour, &minute,
		 	   &kp, &dip, &dec);

  /* Load the output data from firic into the LHS pointers, which are the 
  outputs to MATLAB. */
  *mxGetPr(plhs[0]) = irregs;

}
