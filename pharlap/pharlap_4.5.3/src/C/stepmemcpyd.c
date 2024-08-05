/*
 ==============================================================
                        stepmemcpyd.c
 ==============================================================
  This function performs a stepped memory copy. Compare with
  the standard C memcpy. Instead of copying a contiguous chunk
  of memory, a skip is performed. The size of the skip is 
  determined by the input parameter step.    

  Inputs: 
    double *dest     :  pointer to memory destination
    double *src      :  pointer to source memory
    int     step     :  step size to use
    int     num_vals :  number of values (doubles) to copy

  Change history:
  03/12/2015  M.A.Cervera  V1.0  Author
  
 ==============================================================
*/

void stepmemcpyd(double *dest, double *src, int step, int num_vals)
{
  int idx, srcidx;
  
  srcidx = 0;
  for (idx = 0; idx < num_vals; idx++) {
    dest[idx] = src[srcidx];
    srcidx += step;
  }

}
