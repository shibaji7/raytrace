/*
 ======================================================================
                           iri2016_matlab_wrapper.c
 ======================================================================
  Matlab wrapper for iri2016_calc_ of iri2016
       
  This is a MEX-file for MATLAB.
      
  Change history:
  26/02/2016  M. A. Cervera V1.0 
      Initial version - based on iri2012_matlab_wrapper.c (V1.4),
      updated to use IRI2016
  
  13/07/2017  M. A. Cervera V1.1
      Updated to call the latest version of IRI2016 released 23/02/2017

  13/02/2022  M. A. Cervera  V1.2 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

  ======================================================================
*/

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define max_num_hts 1000   /* limitation of IRI2016 */

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
  #define fscanf fscanf_s       /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check);

void iri2016_calc_(int *jf, int *jmag, float *glat, float *glon, int *year,  
                   int *mmdd, float *dhour, float *heibeg, float *heiend, 
                   float *heistp, float *outf, float *oarr);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, month, ii, jj, idx, idx_mat, arrsize, in_arrsize[11], UT[5];
  float B0B1_model, D_model, hmf2_model, storm_flag;   
  float R12, IG12, F107;            /* solar indices */
  float foF2, hmF2, foF1, hmF1, foE, hmE;   /* user specified iono parameters */
  float B0, B1;               /* user specified iono profile shape parameters */
  float HNEA, HNEE;           /* user specified lower/upper Ne bounderies */
  char error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  /* Input parameters to iri2016_ */
  int   jf[50];                  /* true/false switches for several options */
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
  in_arrsize[2] = 1;  
  in_arrsize[3] = 5;
  in_arrsize[4] = 1;
  in_arrsize[5] = 1;
  in_arrsize[6] = 1;
  in_arrsize[7] = 1;
  
  /* Check for proper number of arguments. */
  if (nrhs != 4 && nrhs != 7 && nrhs != 8 ) {
    mexErrMsgTxt("ERROR: 4, 7, or 8 inputs required when calling iri2016.");
    return;
  }

  if (nlhs != 2) {
    mexErrMsgTxt("ERROR: 2 outputs required when calling iri2016.");
    return;
  }

  /* Check to ensure the inputs are numeric (not strings) or a structure if it
     is the 8th input */
  for (ii = 0; ii < nrhs; ii++) {
    if (ii < 7) {
      if (mxIsNumeric(prhs[ii]) == 0) {
        snprintf(error_mesg, 200, "ERROR: Input %d must be numeric.", ii+1);
        mexErrMsgTxt(error_mesg);
        return;
      }
    }
    else {
      if(!mxIsStruct(prhs[ii])) {
          mexErrMsgIdAndTxt( "MATLAB:iri_options:inputNotStruct",
	  	  	     "ERROR: Input %d must be a structure.", ii+1);
        return;
      }
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
  if (check_ref_data("iri2016") == 0) return;

  /* read the date extents of the id_rz.dat file */
#if !defined(WIN64)
  int dd, mm, yyyy, igrz_mnth_start, igrz_mnth_end, igrz_year_start, 
      igrz_year_end;
  char dummy[10];
  refdata_dir = getenv("DIR_MODELS_REF_DAT");

  snprintf(filename, 200, "%s/iri2016/ig_rz.dat", refdata_dir); 
  FILE *ig_rz_file = fopen(filename, "r");

  fscanf(ig_rz_file, "%d,%d,%d", &dd, &mm, &yyyy);
  fscanf(ig_rz_file, "%s", dummy);
  fscanf(ig_rz_file, "%d,%d,%d,%d\n", &igrz_mnth_start, &igrz_year_start, 
         &igrz_mnth_end, &igrz_year_end);
  fclose(ig_rz_file);
#endif
  
  /* malloc some space for the iri_sub input/output */
  outf_M = 20;
  outf_N = max_num_hts;
  oarr_M = 1;
  oarr_N = 100;
  outf = mxMalloc(outf_M * outf_N * sizeof(float));
  oarr = mxMalloc(oarr_N * sizeof(float));

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     iri_sub_ . Note, dhour (decimal hour) must have 25 added to its value if 
     iri2016 is to consider it to be UT otherwise iri2016 considers dhour to be
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
#if !defined(WIN64)  
  if (year < igrz_year_start || year > igrz_year_end) {
    snprintf(error_mesg, 200, "Year must be in the range %d - %d\nUpdate ig_rz.dat and apf107.dat with the latest avaliable from irimodel.org", igrz_year_start, igrz_year_end);
    mexErrMsgTxt(error_mesg);
    return;
  }
#endif  
  if (nrhs >= 5) {
    heibeg = (float) mxGetScalar(prhs[4]);
    heistp = (float) mxGetScalar(prhs[5]);   
    num_hts = (float) mxGetScalar(prhs[6]);
    /* restrict num_hts to valid range of values for IRI2016*/
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
  jf[38] = 0;
  jf[39] = 1;
  jf[40] = 1;
  jf[41] = 1;
  jf[42] = 1;
  jf[43] = 1;
  jf[44] = 1;
  jf[45] = 1;
  for (ii = 46; ii <= 49; ii++) jf[ii] = 0;   /* these are currently unused */

  /* 2. change some of the standard settings */ 
  jf[33] = 0;                  /* no messages from IRI */
  jf[20] = 1;                  /* ion drift computed */
  jf[21] = 0;                  /* ion densities in m-3 */
  jf[27] = 1;                  /* spread-F probability computed */


  /* 3. override settings based on user input */
  if (R12 > 0 && R12 <= 200) {
    jf[16] = 0;                /* user input R12 */
    jf[24] = 0;                /* user input F10.7 */
    jf[26] = 0;                /* user input IG12 */
    jf[31] = 0;                /* user input F10.7_81 */
    jf[25] = 0;                /* default setting for foF2 storm model is off */
    
    /* Calculate IG12, F10.7 and F10.7_81 from the user input R12 and put into
     the oarr array for input into iri_sub_. The F10.7 expression is from 
     Davies, 1990, pp 442. The expession for IG12 was obtained from the 
     irisub.for fortran code from IRI2016 and was verified against 
     the data found in ig_rz.dat. */
    F107 = 63.75 + R12 * (0.728 + R12*0.00089);
    IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
    *(oarr + 32) = R12;
    *(oarr + 38) = IG12;
    *(oarr + 40) = F107;
    *(oarr + 45) = F107;     /* F10.7_81 is set to F10.7 */
  
  }
  else if (R12 == -1) {
    jf[16] = 1;               /* historical or projected R12 */
    jf[24] = 1;               /* historical or projected F10.7 */
    jf[26] = 1;               /* historical or projected IG12 */
    jf[31] = 1;               /* historical or projected F10.7_81 */
    jf[25] = 1;               /* default setting for foF2 storm model is on */
  }
  else {
    mexErrMsgTxt("Invalid value for R12");
    return;
  }

  /* If it has been input, read in the iri_options structure and change the IRI 
     standard settings if required */
  mxArray *tmp_ptr;
  int invalid_field_value_flag, num_fields, num_valid_fields;
  int foF2_valid, foF1_valid, foE_valid, hmF2_valid, hmF1_valid, hmE_valid;
  int   B0_valid, B1_valid, HNEA_valid, HNEE_valid, result, sizeM, sizeN;
  char tmp_str[50];
  
  invalid_field_value_flag = 0;
  num_valid_fields = 0;
  num_fields = 0;
  foF2_valid = 0; foF2 = 0; hmF2_valid = 0; hmF2 = 0;
  foF1_valid = 0; foF1 = 0; hmF1_valid = 0; hmF1 = 0;
  foE_valid = 0;  foE = 0;  hmE_valid = 0;  hmE = 0;
  B0_valid = 0; B0 = 0;
  B1_valid = 0; B1 = 0;
  HNEA_valid = 0; HNEE_valid  = 0; 
  
  if (nrhs == 8) {

    num_fields = mxGetNumberOfFields(prhs[nrhs-1]);

    /* read in iri_messages field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "iri_messages");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[33] = 0; 
	else if (strcmp(tmp_str, "on") == 0) jf[33] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in foF2 field */   
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foF2");    
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	foF2 = (float) mxGetScalar(tmp_ptr);
	foF2_valid = 1; 
	jf[7] = 0;
	*(oarr + 0) = foF2;
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in hmF2 field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "hmF2");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	hmF2 = (float) mxGetScalar(tmp_ptr);
	hmF2_valid = 1;
	jf[8] = 0;
        *(oarr + 1) = hmF2;
      }
      else invalid_field_value_flag = 1;
    }
 
    /* read in foF1 field */   
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foF1");    
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	foF1 = (float) mxGetScalar(tmp_ptr);
	foF1_valid = 1; 
        jf[12] = 0;
        *(oarr + 2) = foF1;
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in hmF1 field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "hmF1");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	hmF1 = (float) mxGetScalar(tmp_ptr);
	hmF1_valid = 1;
        jf[13] = 0;
        *(oarr + 3) = hmF1;
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in foE field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foE");    
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	foE = (float) mxGetScalar(tmp_ptr);
	foE_valid = 1; 
        jf[14] = 0;
        *(oarr + 4) = foE;
      }
      else invalid_field_value_flag = 1;
    }

    /* read in hmE field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "hmE");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	hmE = (float) mxGetScalar(tmp_ptr);
	hmE_valid = 1;  
        jf[15] = 0;
        *(oarr + 5) = hmE;
      }
      else invalid_field_value_flag = 1;
    }

    /* read in B0 field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "B0");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	B0 = (float) mxGetScalar(tmp_ptr);
	B0_valid = 1;  
        jf[42] = 0;
        *(oarr + 9) = B0;
      }
      else invalid_field_value_flag = 1;
    }

    /* read in B1 field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "B1");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	B1 = (float) mxGetScalar(tmp_ptr);
	B1_valid = 1;  
        jf[43] = 0;
        *(oarr + 34) = B1;
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in HNEA field (Ne lower boundary) */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "HNEA");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	HNEA = (float) mxGetScalar(tmp_ptr);
	HNEA_valid = 1;  
        jf[44] = 0;
        *(oarr + 88) = HNEA;
      }
      else invalid_field_value_flag = 1;
    }

    /* read in HNEE field (Ne upper boundary) */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "HNEE");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsNumeric(tmp_ptr) && !mxIsComplex(tmp_ptr)) {
	HNEE = (float) mxGetScalar(tmp_ptr);
	HNEE_valid = 1;  
        jf[45] = 0;
        *(oarr + 89) = HNEE;
      }
      else invalid_field_value_flag = 1;
    }
    /* read in foF2_coeffs field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foF2_coeffs");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "URSI") == 0) jf[4] = 0;
	else if (strcmp(tmp_str, "CCIR") == 0) jf[4] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in Ni_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Ni_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "RBV-2010 & TTS-2005") == 0) jf[5] = 0; 
	else if (strcmp(tmp_str, "DS-1995 & DY-1985") == 0) jf[5] = 1; 
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in Te_profile field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Te_profile");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "standard") == 0) jf[9] = 1;
	else if (strcmp(tmp_str, "Te/Ne correlation") == 0) jf[9] = 0;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in Te_topside_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Te_topside_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "TBT-2012") == 0) jf[22] = 0;
	else if (strcmp(tmp_str, "Bil-1985") == 0) jf[22] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in Te_PF107_dependance field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Te_PF107_dependance");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[41] = 0;
	else if (strcmp(tmp_str, "on") == 0) jf[41] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }    
    
    /* read in Ne_tops_limited field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Ne_tops_limited");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "f10.7 unlimited") == 0) jf[6] = 0;
	else if (strcmp(tmp_str, "f10.7 limited") == 0) jf[6] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in  Ne_profile_calc field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Ne_profile_calc");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "Lay-function") == 0) jf[10] = 0;
	else if (strcmp(tmp_str, "standard") == 0) jf[10] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in Ne_B0B1_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Ne_B0B1_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "ABT-2009") == 0) {jf[3] = 0; jf[30] = 1;}
	else if (strcmp(tmp_str, "Bil-2000") == 0) {jf[3] = 1; jf[30] = 1;}
	else if (strcmp(tmp_str, "Gul-1987") == 0) {jf[3] = 0; jf[30] = 0;}
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in Ne_topside_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "Ne_topside_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "IRI-2001") == 0) {jf[28] = 1; jf[29] = 1;}
	else if (strcmp(tmp_str, "IRI-2001 corrected") == 0) {
	  jf[28] = 0; jf[29] = 1;}
	else if (strcmp(tmp_str, "NeQuick") == 0) {jf[28] = 0; jf[29] = 0;}
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in F1_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "F1_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "Scotto-1997 no L") == 0) {
	  jf[18] = 1;
	  jf[19] = 1;
	}
	else if (strcmp(tmp_str, "Scotto-1997 with L") == 0) {
	  jf[18] = 1;
	  jf[19] = 0;
	}
	else if (strcmp(tmp_str, "solar zenith") == 0) {
          jf[18] = 0;
	  jf[19] = 1;
	}
	else if (strcmp(tmp_str, "none") == 0) {
          jf[18] = 0;
	  jf[19] = 0;
	}        
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in D_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "D_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "IRI-1990") == 0) jf[23] = 1; 
	else if (strcmp(tmp_str, "FT-2001") == 0) jf[23] = 0;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in hmF2_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "hmF2_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "AMTB") == 0) {jf[38] = 0; jf[39] = 1;}
	else if (strcmp(tmp_str, "Shubin-COSMIC") == 0) {jf[38] = 0; jf[39] =0;}
	else if (strcmp(tmp_str, "M3000F2") == 0) jf[38] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in foF2_storm field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foF2_storm");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[25] = 0; 
	else if (strcmp(tmp_str, "on") == 0) jf[25] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    if (R12 > 0 && R12 <= 200) jf[25] = 0;   /* R12 supplied so turn it off */
    
    /* read in hmF2_storm field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "hmF2_storm");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[35] = 1; 
	else if (strcmp(tmp_str, "on") == 0) jf[35] = 0;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in foE_storm field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "foE_storm");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[34] = 0; 
	else if (strcmp(tmp_str, "on") == 0) jf[34] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in topside_storm field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "topside_storm");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[36] = 1; 
	else if (strcmp(tmp_str, "on") == 0) jf[36] = 0;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }

    /* read in auroral_boundary_model field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "auroral_boundary_model");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "off") == 0) jf[32] = 0; 
	else if (strcmp(tmp_str, "on") == 0) jf[32] = 1;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    /* read in covington field */
    tmp_ptr = mxGetField(prhs[nrhs-1], (mwIndex) 0, "covington");
    if (tmp_ptr != NULL) {
      num_valid_fields = num_valid_fields + 1;
      if (mxIsChar(tmp_ptr)) {
	mxGetString(tmp_ptr, tmp_str, (mwSize) 50);
	if (strcmp(tmp_str, "F10.7_365") == 0) jf[40] = 1; 
	else if (strcmp(tmp_str, "IG12") == 0) jf[40] = 0;
	else invalid_field_value_flag = 1; 
      }
      else invalid_field_value_flag = 1;
    }
    
    
  }

  /* check for invalid fields in the iri_options structure */
  if (num_valid_fields != num_fields) {
    mexPrintf("Warning: IRI2016 - Some of the fields of the supplied iri_options structure\n");
    mexPrintf("         are not valid fields. These fields have been ignored.\n\n");
  }

  /* check for invalid field values in the iri_options structure */
  if (invalid_field_value_flag) {
    mexPrintf("Warning: IRI2016 - Some of the fields of the supplied iri_options structure\n");
    mexPrintf("         have invalid values. These fields have been ignored.\n\n");
  }

  /* check that the user defined ionospheric layer parameters are sensible */
  if (foF2_valid && (foF2 < 0.1 || foF2 > 100.0)) {
    snprintf(error_mesg, 200, "ERROR: Invalid value for input foF2");
    mexErrMsgTxt(error_mesg);
    return;
  }

  if (foF1_valid && (foF1 < 0.1 || foF1 > 100)) { 
    snprintf(error_mesg, 200, "ERROR: Invalid value for input foF1");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (foF2_valid && foF1_valid && foF1 > foF2) {
    snprintf(error_mesg, 200, "ERROR: foF1 larger than foF2");
    mexErrMsgTxt(error_mesg);
    return;
  }

  if (foE_valid && (foE < 0.1 || foE > 100)) { 
    snprintf(error_mesg, 200, "ERROR: Invalid value for input foE");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (foF1_valid && foE_valid && foE > foF1) {
    snprintf(error_mesg, 200, "ERROR: foE larger than foF1");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (foF2_valid && foE_valid && foE > foF2) {
    snprintf(error_mesg, 200, "ERROR: foE larger than foF2");
    mexErrMsgTxt(error_mesg);
    return;
  }

  if (hmF2_valid && (hmF2 < 50 || hmF2 > 1000)) { 
    snprintf(error_mesg, 200, "ERROR: Invalid value for input hmF2");
    mexErrMsgTxt(error_mesg);
    return;
  }

  if (hmF1_valid && (hmF1 < 50 || hmF1 > 1000)) { 
    snprintf(error_mesg, 200, "ERROR: Invalid value for input hmF1");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (hmF2_valid && hmF1_valid && hmF1 > hmF2 ) {
    snprintf(error_mesg, 200, "ERROR: hmF1 larger than hmF2");
    mexErrMsgTxt(error_mesg);
    return;
  }

  if (hmE_valid && (hmE < 50 || hmE > 1000)) { 
    snprintf(error_mesg, 200, "ERROR: Invalid value for input hmE");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (hmF1_valid && hmE_valid && hmE > hmF1) {
    snprintf(error_mesg, 200, "ERROR: hmE larger than hmF1");
    mexErrMsgTxt(error_mesg);
    return;
  }
  if (hmF2_valid && hmE_valid && hmE > hmF2) {
    snprintf(error_mesg, 200, "ERROR: hmE larger than hmF2");
    mexErrMsgTxt(error_mesg);
    return;
  }
   

  /* Call the computational subroutine. */
  jmag = 0;                /* geographic coordinates */
  iri2016_calc_(jf, &jmag, &glat, &glon, &year, &mmdd, &dhour, &heibeg, &heiend,
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

