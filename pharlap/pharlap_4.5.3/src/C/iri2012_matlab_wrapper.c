/*
 ==============================================================
                   iri2012_matlab_wrapper.c
 ==============================================================
  Matlab wrapper for iri_sub_ of iri2012
       
  This is a MEX-file for MATLAB.
      
  Change history:
  29/08/2012  L.H.Pederick V1.0 
      Initial version - based on iri2007_matlab_wrapper.c (V1.2),
      updated to use IRI2012

  31/10/2013  M.A. Cervera V1.1 
      Fixed bug where some arrays were being addressed out of bounds

  13/12/2013  M.A. Cervera V1.2
      Now accepts user input of foF2, hmF2, foF1, hmF1, foE and hmE
      if required

  25/02/2015  M.A. Cervera V1.3
      Maximum allowable values for hmF2, hmF1 and hmE increased

  14/09/2015  M.A. Cervera V1.4
      Fixed bug where too much space was malloced for the second 
      returned array (was 1000*sizeof(double) now 100*sizeof(double) )

  13/02/2022  M. A. Cervera  V1.5 
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

#define max_num_hts 1000

int check_ref_data(char *files_to_check);

void iri_sub_(int *jf, int *jmag, float *glat, float *glon, int *year,  
              int *mmdd, float *dhour, float *heibeg, float *heiend, 
              float *heistp, float *outf, float *oarr);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, month, ii, jj, idx, idx_mat, arrsize, in_arrsize[10], UT[5];
  float B0B1_model, D_model, storm_flag;   
  float R12, IG12, F107;            /* solar indices */
  float foF2, hmF2, foF1, hmF1, foE, hmE;   /* user specified iono parameters */
  char error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  /* Input parameters to iri_sub_ */
  int   jf[38];                  /* true/false switches for several options */
  int   jmag;                    /* = 0 geographic,  = 1 geomag coordinates */
  float glat, glon;              /* Geographic lat/long */
  int   year, mmdd;              /* year and date */
  float dhour;                   /* decimal hour  */
  float heibeg, heiend, heistp;  /* height range in km, max 100 heights */
  int   num_hts;

  /* Output arrays from iri_sub_ and their sizes*/
  float *outf, *oarr;
  int   outf_M, outf_M_keep, outf_N, outf_N_keep, oarr_M, oarr_N;
  
  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 1;  /* can also be 2 elements */
  in_arrsize[3] = 5;
  in_arrsize[4] = 1;
  in_arrsize[5] = 1;
  in_arrsize[6] = 1;
  in_arrsize[7] = 1;
  in_arrsize[8] = 1;
  in_arrsize[9] = 6;

  /* Check for proper number of arguments. */
  if (nrhs != 4 && nrhs != 7 && nrhs != 8 && nrhs != 9 && nrhs != 10) {
    mexErrMsgTxt("ERROR: 4, 7, 8, 9 or 10 inputs required when calling iri2012.");
    return;
  }

  if (nlhs != 2) {
    mexErrMsgTxt("ERROR: 2 outputs required when calling iri2012.");
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
    if (ii == 2) {
      if (arrsize != 1 && arrsize != 2) {
        snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
        mexErrMsgTxt(error_mesg);
        return;
      }
    }
    else {     
      if (arrsize != in_arrsize[ii]) {
        snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
        mexErrMsgTxt(error_mesg);
        return;
      }
    }
  }
 
  /* Check to see if the reference data files exist. */
  if (check_ref_data("iri2012") == 0) return;

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     iri_sub_ . Note, dhour (decimal hour) must have 25 added to its value if 
     iri2012 is to consider it to be UT otherwise iri2012 considers dhour to be
     local time. */
  glat = (float) mxGetScalar(prhs[0]);
  glon = (float) mxGetScalar(prhs[1]);
  
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  arrsize = mrows*ncols;
  if (arrsize == 1) {
    R12 = (float) mxGetScalar(prhs[2]);
    storm_flag = 0;
  }
  else {
    R12 = (float) *(mxGetPr(prhs[2]) + 0);
    storm_flag = (float) *(mxGetPr(prhs[2]) + 1);
    if (storm_flag != 0 && storm_flag != 1) {
      mexErrMsgTxt("Invalid value for storm flag");
     return;
    }
  }
  
  for (ii = 0; ii < 5; ii++) {
    UT[ii] = *(mxGetPr(prhs[3]) + ii);
  }
  year = UT[0];
  mmdd = UT[1] * 100 + UT[2];
  dhour = (float) (UT[3] + UT[4] / 60.0 + 25.0);
  if (year < 1958 || year > 2018) {
    mexErrMsgTxt("Year must be in the range 1958 - 2018");
    return;
  }
  
  if (nrhs >= 7) {
    heibeg = (float) mxGetScalar(prhs[4]);
    heistp = (float) mxGetScalar(prhs[5]);   
    num_hts = (float) mxGetScalar(prhs[6]);
    /* restrict num_hts to valid range of values */
    if (num_hts < 2) num_hts = 2;
    if (num_hts > max_num_hts) num_hts = max_num_hts;
    heiend = heibeg + (num_hts - 1) * heistp;
  }
  else {
    heibeg = 1;
    heistp = 1;
    heiend = 2;
    num_hts = 0;
  }

  B0B1_model = 2;  /* set the default - override if user input has been supplied */
  if (nrhs >= 8) {
    B0B1_model = (float) mxGetScalar(prhs[7]);
    if (B0B1_model != 1 && B0B1_model != 2 && B0B1_model != 3){
      mexErrMsgTxt("Invalid value for the B0 / B1 model");
      return;
    }
  }

  D_model = 1;  /* set the default - override if user input has been supplied */
  if (nrhs >= 9) {
    D_model = (float) mxGetScalar(prhs[8]);
    if (D_model != 1 && D_model != 2){
      mexErrMsgTxt("Invalid value for the D region model");
      return;
    }
  }

  if (nrhs == 10) {
    /* read in the user defined ionospheric layer parameters */
    foF2 = *(mxGetPr(prhs[9]) + 0);
    hmF2 = *(mxGetPr(prhs[9]) + 1);
    foF1 = *(mxGetPr(prhs[9]) + 2);
    hmF1 = *(mxGetPr(prhs[9]) + 3);
    foE = *(mxGetPr(prhs[9]) + 4);
    hmE = *(mxGetPr(prhs[9]) + 5);

    /* check that these parameters are sensible */
    if ((foF2 < 0.1 || foF2 > 100.0) && (foF2 != -1.0)) {
      snprintf(error_mesg, 200, "ERROR: Invalid value for input foF2");
      mexErrMsgTxt(error_mesg);
      return;
    }

    if ((foF1 < 0.1 || foF1 > 100) && (foF1 != -1)) { 
      snprintf(error_mesg, 200, "ERROR: Invalid value for input foF1");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (foF1 > foF2 && foF2 != -1) {
      snprintf(error_mesg, 200, "ERROR: foF1 larger than foF2");
      mexErrMsgTxt(error_mesg);
      return;
    }

    if ((foE < 0.1 || foE > 100) && (foE != -1)) { 
      snprintf(error_mesg, 200, "ERROR: Invalid value for input foE");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (foE > foF1 && foF1 != -1) {
      snprintf(error_mesg, 200, "ERROR: foE larger than foF1");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (foE > foF2 && foF2 != -1) {
      snprintf(error_mesg, 200, "ERROR: foE larger than foF2");
      mexErrMsgTxt(error_mesg);
      return;
    }

    if ((hmF2 < 50 || hmF2 > 1000) && (hmF2 != -1.0)) { 
      snprintf(error_mesg, 200, "ERROR: Invalid value for input hmF2");
      mexErrMsgTxt(error_mesg);
      return;
    }

    if ((hmF1 < 50 || hmF1 > 1000) && (hmF1 != -1)) { 
      snprintf(error_mesg, 200, "ERROR: Invalid value for input hmF1");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (hmF1 > hmF2 && hmF2 != -1) {
      snprintf(error_mesg, 200, "ERROR: hmF1 larger than hmF2");
      mexErrMsgTxt(error_mesg);
      return;
    }

    if ((hmE < 50 || hmE > 1000) && (hmE != -1)) { 
      snprintf(error_mesg, 200, "ERROR: Invalid value for input hmE");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (hmE > hmF1 && hmF1 != -1) {
      snprintf(error_mesg, 200, "ERROR: hmE larger than hmF1");
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (hmE > hmF2 && hmF2 != -1) {
      snprintf(error_mesg, 200, "ERROR: hmE larger than hmF2");
      mexErrMsgTxt(error_mesg);
      return;
    }
  }
  else {
    /* set the ionospheric layer parameters to dummy value so the appropriate */
    /* true/false switches in the jf array are set appropriately */
    foF2 = -1;   
    hmF2 = -1;
    foF1 = -1;
    hmF1 = -1;
    foE = -1;
    hmE = -1;
  }

  /* malloc some space for the iri_sub input/output */
  outf_M = 20;
  outf_N = max_num_hts;
  oarr_M = 1;
  oarr_N = 100;
  outf = mxMalloc(outf_M * outf_N * sizeof(float));
  oarr = mxMalloc(oarr_N * sizeof(float));

  /* set up the true/false switches jf array and input values of oarr array */
  /* 1. IRI standard settings */
  for (ii = 0; ii < 3; ii++) jf[ii] = 1;  
  for (ii = 3; ii < 6; ii++) jf[ii] = 0;
  for (ii = 6; ii < 20; ii++) jf[ii] = 1;
  for (ii = 20; ii < 23; ii++) jf[ii] = 0; 
  for (ii = 23; ii < 27; ii++) jf[ii] = 1;
  for (ii = 27; ii < 30; ii++) jf[ii] = 0;
  for (ii = 30; ii < 32; ii++) jf[ii] = 1;
  jf[32] = 0;
  jf[33] = 1;
  jf[34] = 0;
  jf[35] = 1;
  jf[36] = 1;
  jf[37] = 1;

  /* 2. change some of the standard settings */
  if (B0B1_model == 1) {
    jf[3]  = 0;                  /* B0/B1 - APT-2009 option */
    jf[30] = 1;
  }
  else if (B0B1_model == 2) {
    jf[3]  = 1;                  
    jf[30] = 1;                  /* B0/B1 - Bil-2000 option */
  }
  else {            /* OK B0B1_model must = 3 */
    jf[3]  = 0;                  
    jf[30] = 0;                  /* B0/B1 - Gul-1987 option */
  }
  
  if (D_model == 1) {
    jf[23] = 1;        /*  D-region: IRI-1990 */
  }
  else {    /* OK B0B1_model must = 2 */
    jf[23] = 0;        /*  D-region: FT-2001 and DRS-1995 */
  }
  
  jf[33] = 0;                  /* no messages from IRI */
  jf[20] = 1;                  /* ion drift computed */
  jf[21] = 0;                  /* ion densities in m-3 */
  if (R12 > 0 && R12 <= 200) {
    if (storm_flag == 0) {
      jf[25] = 0;                  /* foF2 storm model off */
    }
    else {
      jf[25] = 1;                  /* foF2 storm model on */
    }      
    jf[16] = 0;                  /* user input R12 */
    jf[24] = 0;                  /* user input F10.7 */
    jf[26] = 0;                  /* user input IG12 */
    jf[31] = 0;                  /* user input F10.7_81 */
  }
  else if (R12 == -1) {
    if (storm_flag == 0) {
      jf[25] = 0;                  /* foF2 storm model off */
    }
    else {
      jf[25] = 1;                  /* foF2 storm model on */
    }      
    jf[16] = 1;                  /* historical or projected R12 */
    jf[24] = 1;                  /* historical or projected F10.7 */
    jf[26] = 1;                  /* historical or projected IG12 */
    jf[31] = 1;                  /* historical or projected F10.7_81 */
  }
  else {
    mexErrMsgTxt("Invalid value for R12");
    return;
  }
  if (foF2 != -1) {
    jf[7] = 0;
    *(oarr + 0) = foF2;
  }
  if (hmF2 != -1) {
    jf[8] = 0;
    *(oarr + 1) = hmF2;
  }
  if (foF1 != -1) {
    jf[12] = 0;
    *(oarr + 2) = foF1;
  }
  if (hmF1 != -1) {
    jf[13] = 0;
    *(oarr + 3) = hmF1;
  }
  if (foE != -1) {
    jf[14] = 0;
    *(oarr + 4) = foE;
  }
  if (hmE != -1) {
    jf[15] = 0;
    *(oarr + 5) = hmE;
  }

  /* Calculate IG12, F10.7 and F10.7_81 from the user input R12 and put into the  
     oarr array for input into iri_sub_. The F10.7 expression is from Davies, 1990, 
     pp 442. The expession for IG12 was obtained from the irisub.for fortran 
     code (line 851) from IRI2012 and was verified against the data found in 
     ig_rz.dat. */
  if (R12 > 0) {
    F107 = 63.75 + R12 * (0.728 + R12*0.00089);
    IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
    *(oarr + 32) = R12;
    *(oarr + 38) = IG12;
    *(oarr + 40) = F107;
    *(oarr + 45) = F107;     /* F10.7_81 is set to F10.7 */
  }

  /* Call the computational subroutine. */
  jmag = 0;                /* geographic coordinates */
  iri_sub_(jf, &jmag, &glat, &glon, &year, &mmdd, &dhour, &heibeg, &heiend,
           &heistp, outf, oarr);

  /* Create matricies for the return arguments. We only want to retain the 
     first 14 x num_hts outputs of outf */
  outf_M_keep = 14;
  outf_N_keep = num_hts;
  plhs[0] = mxCreateDoubleMatrix(outf_M_keep, outf_N_keep, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(oarr_M, oarr_N, mxREAL);

  /* Load the output data from iri_sub into the LHS pointers, which are the 
  outputs to MATLAB. */
  for (ii = 0; ii < outf_N_keep; ii++) {
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
  mxFree(outf);
  mxFree(oarr);
}

