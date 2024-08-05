/*
 ==============================================================
                        check_ref_data.c
 ==============================================================
  This function checks the existance of the reference data files 
  required by the various propagation, ionopsheric and geomagnetic
  models.

  Inputs: 
    char *files_to_check :  string indicating which files to check the existence
                            "iri2007"  -  IRI2007 reference data files
                            "iri2012"  -  IRI2012 reference data files
                            "iri2016"  -  IRI2016 reference data files
                            "land_sea" -  Global Land/Sea mask data file

  Return Value:
    1 if files exist, otherwise 0

  Change history:
  18/02/2010  M. A. Cervera  V1.0  Author
  
  16/05/2014  M. A. Cervera  V1.1
    Added input string to control which files will be checked

  25/02/2016  M. A. Cervera  V1.2
    Updated to also check for IRI2016 reference data files

  13/02/2022  M. A. Cervera  V1.3 
     Now using sprintf_s (under Windows) instead of _snprintf 
     for generating strings

 ==============================================================
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mex.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check) 
{

  int ii, month, return_val;
  char error_mesg[250], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  return_val = 1;
  
  /* Obtain the reference data directory fom the relevant environment 
     variable.  Check to see if the reference data files exist. */
  refdata_dir = getenv("DIR_MODELS_REF_DAT");

  if (refdata_dir == NULL) {
    snprintf(error_mesg, 250, "Environment variable DIR_MODELS_REF_DAT is not set");
    mexErrMsgTxt(error_mesg);
    return_val = 0;
  }      


  /* check for existence of IRI2007 data files */
  if (strncmp(files_to_check, "iri2007", 7) == 0) {

    snprintf(filename, 200, "%s/iri2007/dgrf00.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (ii = 45; ii <= 95; ii = ii+5) {
      snprintf(filename, 200, "%s/iri2007/dgrf%d.dat", refdata_dir, ii); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);
    }  

    snprintf(filename, 200, "%s/iri2007/igrf05.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2007/igrf05s.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2007/ap.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2007/ig_rz.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (month = 1; month <= 12; month++) {
      snprintf(filename, 200, "%s/iri2007/ccir%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2007/ursi%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;

      }
      fclose(file_ptr);
    }
  }    /* end of check for existence of iri2007 reference data files */

  
  /* check for existence of IRI2012 data files */
  if (strncmp(files_to_check, "iri2012", 7) == 0) {

    snprintf(filename, 200, "%s/iri2012/igrf2015.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2012/igrf2015s.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (ii = 1945; ii <= 2010; ii = ii+5) {
      snprintf(filename, 200, "%s/iri2012/dgrf%d.dat", refdata_dir, ii); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);
    }  

    snprintf(filename, 200, "%s/iri2012/apf107.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2012/ig_rz.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (month = 1; month <= 12; month++) {
      snprintf(filename, 200, "%s/iri2012/ccir%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2012/ursi%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;

      }
      fclose(file_ptr);
    }  

  }    /* end of check for existence of iri2012 reference data files */

  
  /* check for existence of IRI2016 data files */
  if (strncmp(files_to_check, "iri2016", 7) == 0) {

    snprintf(filename, 200, "%s/iri2016/igrf2015.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2016/igrf2015s.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (ii = 1945; ii <= 2010; ii = ii+5) {
      snprintf(filename, 200, "%s/iri2016/dgrf%d.dat", refdata_dir, ii); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);
    }  

    snprintf(filename, 200, "%s/iri2016/apf107.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2016/ig_rz.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (month = 1; month <= 12; month++) {
      snprintf(filename, 200, "%s/iri2016/mcsat%d.dat", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2016/ccir%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2016/ursi%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;

      }
      fclose(file_ptr);
    }  

  }    /* end of check for existence of iri2016 reference data files */


  /* check for existence of IRI2020 data files */
  if (strncmp(files_to_check, "iri2020", 7) == 0) {

    snprintf(filename, 200, "%s/iri2020/igrf2020.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2020/igrf2020s.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (ii = 1945; ii <= 2015; ii = ii+5) {
      snprintf(filename, 200, "%s/iri2020/dgrf%d.dat", refdata_dir, ii); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);
    }  

    snprintf(filename, 200, "%s/iri2020/apf107.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    snprintf(filename, 200, "%s/iri2020/ig_rz.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0;
    }  
    fclose(file_ptr);

    for (month = 1; month <= 12; month++) {
      snprintf(filename, 200, "%s/iri2020/mcsat%d.dat", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2020/ccir%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;
      }
      fclose(file_ptr);

      snprintf(filename, 200, "%s/iri2020/ursi%d.asc", refdata_dir, month+10); 
      file_ptr = fopen(filename, "r");
      if (file_ptr == NULL) {
	snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
	mexErrMsgTxt(error_mesg);
	return_val = 0;

      }
      fclose(file_ptr);
    }  

  }    /* end of check for existence of iri2020 reference data files */


  /* check for existence of global land-sea map data file */
  if (strncmp(files_to_check, "land_sea", 7) == 0) {

    snprintf(filename, 200, "%s/global_land_Mask_3600_by_1800.dat", refdata_dir); 
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL) {
      snprintf(error_mesg, 250, "Reference data file, %s, does not exist or is unreadable. Check that the environment variable DIR_MODELS_REF_DAT is set to the directory containing the reference data.", filename);
      mexErrMsgTxt(error_mesg);
      return_val = 0; 
    }
    fclose(file_ptr);

  }   /* end of check for existence of global land-sea map data file */

  return return_val;
}
