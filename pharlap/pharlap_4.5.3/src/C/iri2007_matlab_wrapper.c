/*
 ==============================================================
                   iri2007_matlab_wrapper.c
 ==============================================================
  Matlab wrapper for iri_sub_ of iri2007
       
  This is a MEX-file for MATLAB.
      
  Change history:
  28/11/2008  M.A.Cervera  V1.0  Author

  09/12/2009  M.A. Cervera V1.1
      ht_start and ht_step are now optional inputs, if they are not supplied
      then the profiles are not calculated which may speed up the call to IRI

  18/02/2010  M.A. Cervera V1.2  
      Improved error handling

  13/02/2022  M. A. Cervera  V1.3 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

  ==============================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check);

void iri_sub_(int *jf, int *jmag, float *glat, float *glon, int *year,  
              int *mmdd, float *dhour, float *heibeg, float *heiend, 
              float *heistp, float *outf, float *oarr);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, month, ii, jj, idx, idx_mat, arrsize, in_arrsize[6], UT[5];
  float R12, IG12, F107;                               /* solar indices */
  char error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  /* Input parameters to iri_sub_ */
  int   jf[30];                  /* true/false switches for several options */
  int   jmag;                    /* = 0 geographic,  = 1 geomag coordinates */
  float glat, glon;              /* Geographic lat/long */
  int   year, mmdd;              /* year and date */
  float dhour;                   /* decimal hour  */
  float heibeg, heiend, heistp;  /* height range in km, max 100 heights */
  int   num_hts;

  /* Output arrays from iri_sub_ and their sizes*/
  float *outf, *oarr;
  int   outf_M, outf_M_keep, outf_N, oarr_M, oarr_N;

  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 1;
  in_arrsize[3] = 5;
  in_arrsize[4] = 1;
  in_arrsize[5] = 1;
 
  /* Check for proper number of arguments. */
  if (nrhs != 4 && nrhs != 6) {
    mexErrMsgTxt("ERROR: 4 or 6 inputs required when calling iri2007.");
    return;
  }

  if (nlhs != 2) {
    mexErrMsgTxt("ERROR: 2 outputs required when calling iri2007.");
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
 
  /* Check to see if the reference data files exist. */
  if (check_ref_data("iri2007") == 0) return;

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     iri_sub_ . Note, dhour (decimal hour) must have 25 added to its value if 
     iri2007 is to consider it to be UT otherwise iri2007 considers dhour to be
     local time. */
  glat = (float) mxGetScalar(prhs[0]);
  glon = (float) mxGetScalar(prhs[1]);
  R12 = (float) mxGetScalar(prhs[2]);
  for (ii = 0; ii < 5; ii++) {
    UT[ii] = *(mxGetPr(prhs[3]) + ii);
  }
  year = UT[0];
  mmdd = UT[1] * 100 + UT[2];
  dhour = (float) (UT[3] + UT[4] / 60.0 + 25.0);
  if (nrhs == 6) {
    heibeg = (float) mxGetScalar(prhs[4]);
    heistp = (float) mxGetScalar(prhs[5]);   
    num_hts = 100;  /* set num_hts to the maximum allowed by iri_sub */
    heiend = heibeg + (num_hts - 1) * heistp;
  }
  else {
    heibeg = 1;
    heistp = 1;
    heiend = 2;
  }

  /* set up the true/false switches jf array */
  /* 1. IRI standard settings */
  for (ii = 0; ii < 4; ii++) jf[ii] = 1;          
  for (ii = 4; ii < 6; ii++) jf[ii] = 0;
  for (ii = 6; ii < 20; ii++) jf[ii] = 1;
  for (ii = 20; ii < 23; ii++) jf[ii] = 0; 
  for (ii = 23; ii < 27; ii++) jf[ii] = 1;
  for (ii = 27; ii < 30; ii++) jf[ii] = 0;
  /* 2. change some of the standard settings */
  jf[11] = 0;                  /* no messages from IRI */
  jf[20] = 1;                  /* ion drift computed */
  jf[21] = 0;                  /* ion densities in m-3 */
  if (R12 > 0 && R12 <= 200) {
    jf[25] = 0;                  /* storm model off */
    jf[16] = 0;                  /* user input R12 */
    jf[24] = 0;                  /* user input F10.7 */
    jf[26] = 0;                  /* user input IG12 */
  }
  else if (R12 == -1) {
    jf[25] = 1;                  /* storm model on */
    jf[16] = 1;                  /* historical or projected R12 */
    jf[24] = 1;                  /* historical or projected F10.7 */
    jf[26] = 1;                  /* historical of projected IG12 */
  }
  else if (R12 == -2) {
    jf[25] = 0;                  /* storm model off */
    jf[16] = 1;                  /* historical or projected R12 */
    jf[24] = 1;                  /* historical or projected F10.7 */
    jf[26] = 1;                  /* historical of projected IG12 */
  }
  else {
    mexErrMsgTxt("Invalid value for R12");
    return;
  }

  /* malloc some space for the iri_sub output */
  outf_M = 20;
  outf_N = 100;
  oarr_M = 1;
  oarr_N = 50;
  outf = malloc(outf_M * outf_N * 4);   /* 20 by 100 array of floats */
  oarr = malloc(oarr_N * 4);            /* 1 by 50 array of floats */

  /* Calculate IG12 and F10.7 from the user input R12 and put into the oarr 
     array for input into iri_sub_. The F10.7 expression is from Davies, 1990, 
     pp 442. The expession for IG12 was obtained from the irisub.for fortran 
     code (line 721) from IRI2007 and was verified against the data found in 
     ig_rz.dat. */
  if (R12 > 0) {
    F107 = 63.75 + R12 * (0.728 + R12*0.00089);
    IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
    *(oarr + 32) = R12;
    *(oarr + 38) = IG12;
    *(oarr + 40) = F107;
  }

  /* Call the computational subroutine. */
  jmag = 0;                /* geographic coordinates */
  iri_sub_(jf, &jmag, &glat, &glon, &year, &mmdd, &dhour, &heibeg, &heiend,
           &heistp, outf, oarr);

  /* Create matricies for the return arguments. We only want to retain the 
     first 11 outputs of outf */
  outf_M_keep = 11;
  plhs[0] = mxCreateDoubleMatrix(outf_M_keep, outf_N, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(oarr_M, oarr_N, mxREAL);

  /* Load the output data from iri_sub into the LHS pointers, which are the 
  outputs to MATLAB. */
  for (ii = 0; ii < outf_N; ii++) {
    for (jj = 0; jj < outf_M_keep; jj++) {
      idx = ii*outf_M + jj;
      idx_mat = ii*outf_M_keep + jj;
      *(mxGetPr(plhs[0]) + idx_mat) = *(outf + idx);
    }
  }
  
  for (jj = 0; jj < oarr_N; jj++) {
    *(mxGetPr(plhs[1]) + jj) = oarr[jj];
  }

  /* clean up after ourselves - free the malloced space */
  free(outf);
  free(oarr);
}

