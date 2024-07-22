/******************************************************************************
FUNCTION:	WARP_TB_BYTE

PURPOSE: Warp a byte image with table based 4x4 resampling.  This inverse 
	 mapping function maps in the form:

     insamp = a[0] + (outsamp*a[1]) + (outline*a[2]) + (outline*outsamp*a[3]);
     inline = b[0] + (outsamp*b[1]) + (outline*b[2]) + (outline*outsamp*b[3]);

   Any point (insamp,inline) which lies outside the input image buffer is set
   to the background fill value given in the "fill" variable.
*******************************************************************************/
#include "las.h"

FUNCTION warp_tb_byte(inbuf,outbuf,nl,ns,inwind,a,b,resamp_tbl,fill)

unsigned char *inbuf;		/* input image buffer 			*/
unsigned char *outbuf;		/* output image buffer 			*/
long nl,ns;			/* number of output buffer lines/samples*/
long inwind[];			/* number of input buffer lines/samples */
double *a,*b;			/* bilinear mapping coefficients 	*/
float *resamp_tbl;		/* resampling weight table 		*/
unsigned char fill;		/* fill value 				*/
{
long i,j,k,l;			/* loop counters 			*/
long x,y;			/* sample & line values			*/
long idp,idl;
float temp;			/* temporary value			*/
float *tblptr;			/* value in table			*/
double a2,b2,a3,b3;		/* accumulators for coefficients 	*/
double dline,dsamp;		/* increment in line/samp coefficients 	*/
double line,samp;		/* input line/sample location 		*/
unsigned char *ptr;		/* value in inbuf			*/

for (a2=a3=b2=b3=0.0, i=nl; i--; a2+=a[2],a3+=a[3],b2+=b[2],b3+=b[3]) 
   {
   samp  = a[0] + a2; 
   line  = b[0] + b2;
   dsamp = a[1] + a3; 
   dline = b[1] + b3;
   for (j=ns; j--; samp+=dsamp, line+=dline) 
      {
      y = (long)line; 
      x = (long)samp;
      idp = (long)(32 * (0.015625 + samp - x)); /* 0.015625 = 1.0/64.0 */
      idl = (long)(32 * (0.015625 + line - y));
      x--; 
      y--;
      if ((x < inwind[SS]) || (x > inwind[SS]+inwind[NS]-4) || 
	  (y < inwind[SL]) || (y > inwind[SL]+inwind[NL]-4)) 
	 *outbuf++=fill;
      else 
	 {
         tblptr = resamp_tbl + ((idl * 33 + idp) * 16);
         for (temp = 0.0, k = 4; k--; y++)
            for (ptr=inbuf+(inwind[NS]*(y-inwind[SL]))+(x-inwind[SS]),l=4;l--;)
               temp += *ptr++ * *tblptr++;
         /* *outbuf++ = (temp > 254.5) ? 255 : (long)(temp + 0.5); */
	 if (temp < 0.0) temp = 0.0;
	 if (temp > 255.0) temp = 255.0;
         *outbuf++ = (long)(temp + 0.5);
         }
      }
   }
}
