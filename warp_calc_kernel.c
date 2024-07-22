/******************************************************************************
FUNCTION:	WARP_CALC_KERNEL

PURPOSE: Warp a byte image with table based 4x4 resampling.  
   	Use warp_calc_kernel() to calculate the 4x4 mapping kernel to 1/32 pixel
   	increments.  Only cubic convolution is currently supported.
*******************************************************************************/
#include "worgen.h"

FUNCTION warp_calc_kernel(resamp_tbl,alpha)

float **resamp_tbl;	/* 4 x 4 resampling weight table 		*/
float alpha;		/* parametric cubic convolution alpha parameter */
{
float *wtab;		/* ptr to resampling weight table 		*/
float delta, tmp;	/* temp vars for table calculation 		*/
float wght[4][33];	/* one dimensional resampling weight table 	*/
long i,j,k,l;		/* loop counters 				*/
long nonfatal=1;	/* nonfatal error flag				*/

/* Allocate table memory
  ---------------------*/
*resamp_tbl = (float *)malloc(17424*sizeof(float));
if (*resamp_tbl == 0) 
    {
    c_errmsg("Error allocating dynamic memory","avhrrchips-alloc",&nonfatal);
    return(E_FAIL);
    }
wtab = *resamp_tbl;

/* Calculate resampling table
  --------------------------*/
for (i = 0; i < 33; i++) 
   {
   /* delta = i * 0.03125; */
   delta = i / 32.0; 
   tmp = 1.0 - delta;
   wght[0][i] = alpha * delta * tmp * tmp;
   wght[1][i] = 1.0-((alpha+3.0)*delta*delta)+((alpha+2.0)*delta*delta*delta);
   wght[2][i] = 1.0-((alpha+3.0)*tmp*tmp) + ((alpha+2.0)*tmp*tmp*tmp);
   wght[3][i] = alpha * tmp * delta * delta;
   }

for (i=0; i<33; i++)
   for (j=0; j<33; j++)
      for (k=0; k<4; k++)
         for (l=0; l < 4; l++)
            *wtab++ = wght[k][i] * wght[l][j];

return(E_SUCC);
}
