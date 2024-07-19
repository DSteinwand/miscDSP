/*****************************************************************************
NAME:           kernel

PURPOSE:        
The kernel module contains the routines needed to create and manage the 
resampling kernel.  The kernel is the weight table used for spatially 
resampling pixels.

ROUTINES:
    Kernel_Setup
    Cleanup_Kernel
    Create_Resampling_Kernel
    Get_Resample_Weight_Table_Ptr
    get_lines_in_kernel
    get_samples_in_kernel
    num_left_kernel_samples
    num_right_kernel_samples
    num_top_kernel_lines
    num_bottom_kernel_lines
    get_kernel_step_size

DEVELOPMENT HISTORY:
  Baselined from revision 1.14 of kernel.c

NOTES:
- The kernel module will build the resampling weight table for either 
  cubic convolution or mtf resampling.  The kernel is built for the current 
  band being resampled.
- The Kernel_Setup routine also handles the bilinear resampling option by
  setting the kernel size correctly for bilinear resampling.  However, no
  resampling weight table is built for bilinear resampling.


ALGORITHM REFERENCES:
    The algorithm for constructing the kernel from the five
    coefficients is documented in:
       "Cubic Convolution for One-Pass Restoration and Resampling"
       by Stephen E. Reichenbach and Kevin Haake
       Presented at IGARSS '96
    Note: this reference is actually for MTF resampling which is currently
    disabled (although MTF is a generic form of cubic convolution).


*****************************************************************************/

#include <stdio.h>         /* i/o routines                    */
#include <stdlib.h>        /* memory management routines      */
#include <math.h>          /* fabs routine                    */
#include "resample.h"
#include "debug.h"         /* ASSERT macro                    */

/*
**
** MACROS
**
*/

#define SUBPIXEL_STEPS 32 /* defines the number of subpixel steps to
                             include in the kernel weights.                 */
#define SUBPIXEL_TABLE_ENTRIES (SUBPIXEL_STEPS + 1)
                          /* defines the number of entries in the kernel
                             weight table for the subpixel steps.  This is
                             one more than the number of SUBPIXEL steps since
                             the fractional line number is rounded to the
                             nearest 1/32 subpixel in the resampler.  To
                             avoid having to check for rounding up to the
                             next whole line number, the extra subpixel step
                             is included in the table.                      */
#define LINES_IN_KERNEL 4 /* defines the number of lines in the resampling
                             kernel. This is the normal cubic convolution
                             kernel size.                                   */
#define SAMPLES_IN_KERNEL 4
                          /* defines the number of samples in the resampling
                             kernel.  This is two more than the normal
                             convolution kernel size to account for the
                             detector delays that may be present in the sample
                             direction.                                     */ 
#define BI_LINES_IN_KERNEL 2
                          /* defines the number of lines in the resampling
                             kernel for bilinear resampling.                */
#define BI_SAMPLES_IN_KERNEL 2
                          /* defines the number of samples in the resampling
                             kernel for bilinear resampling.                */
#define DISABLE_MTF       /* disable the MTF compiling for this version     */

/*
**
** MODULE VARIABLES
**
*/

static double *kernel_weights_ptr = NULL;
                          /* kernel_weights_ptr is a pointer to the table of
                             resampling weights for cubic convolution and mtf
                             (general form of cubic convolution) resampling.

                             The table is laid out as a 6-dimensional array
                             with the following indices: table[line subpixel
                             offset (size SUBPIXEL_TABLE_ENTRIES)][sample
                             subpixel offset (size SUBPIXEL_TABLE_ENTRIES)]
                             [kernel line (size LINES_IN_KERNEL)][kernel
                             sample (size SAMPLES_IN_KERNEL)]               */
static int lines_in_kernel = 1;
                          /* number of lines in the kernel.  Its default
                             initialization is one line (assuming the nearest
                             neighbor kernel).  This is set to LINES_IN_KERNEL
                             when Kernel_Setup is called.                   */
static int samples_in_kernel = 1;
                          /* number of samples in the kernel.  Its default
                             initialization is one sample (assuming the
                             nearest neighbor kernel).  This is set to
                             SAMPLES_IN_KERNEL when Kernel_Setup is called. */
static int kernel_size_accessed;
                          /* the kernel_size_accessed is set to a non-zero
                             value when either of the functions to return a
                             kernel size parameter is called.  This allows a
                             debugging check to be put into the Kernel_Setup
                             routine to make sure the kernel size isn't
                             accessed before Kernel_Setup is called and
                             changes the kernel size.  kernel_size_access bit
                             0 means the kernel size has been accessed, and
                             bit 1 means the kernel has been setup.         */


/*****************************************************************************
NAME:           Create_Resampling_Kernel

PURPOSE:
Create_Resampling_Kernel builds the weight table for the resampling kernel.
The weights are setup properly for mtf or cubic convolution resampling, with
or without detector delays taken into account.

RETURN VALUE:
Either E_SUCC or E_FAIL are returned, depending on whether the kernel 
creation is successful or not.

NOTES:
The MTF parameters allow for having separate forward and reverse parameters.
This is not currently used, but is included for future growth.

ALGORITHM REFERENCES:
none

*****************************************************************************/

static int Create_Resampling_Kernel
(
    double *weight_ptr,            /* I/O: pointer to weights table         */
    RESAMPLE_TYPE resampling_type, /* I: resampling type                    */
    double alpha                   /* I: cc alpha parameter                 */
)
{
    int    kernel_line;            /* Counters                              */
    int    kernel_sample;
    int    subpixel_line;
    int    subpixel_sample;
    int    left_samp;              /* size of the kernel, relative to the   */
    int    right_samp;             /* center of the kernel                  */
    int    top_line;
    int    bottom_line;
    double dx, dy;                 /* Sub-pixel spacing for cubic           */
    double x, y;                   /* X and Y location of value to calc     */
    double *w_ptr;                 /* pointer to current weight table       */
    double ccwy;                   /* weight for current y (line) component */


#ifndef DISABLE_MTF
    const MTF_PARAMETERS_TYPE *mtf_along_param_ptr 
        = get_forward_MTF_weights_along(band);
    const MTF_PARAMETERS_TYPE *mtf_across_param_ptr
        = get_forward_MTF_weights_across(band);
#endif


    /* Verify the resampling kernel is only set up for CC since MTF has 
       been disabled for this version of the resampler.
       ----------------------------------------------------------------*/
    ASSERT (resampling_type == CC_RESAMPLING);

    /* Get the kernel dimensions.
       --------------------------*/
    left_samp = -num_left_kernel_samples();
    right_samp = num_right_kernel_samples();
    top_line = -num_top_kernel_lines();
    bottom_line = num_bottom_kernel_lines();

    /* Set temporary pointer to the weight table.
       ------------------------------------------*/
    w_ptr = weight_ptr;

    /* Build the kernel for either cubic convolution or mtf.
       -----------------------------------------------------*/
    if (resampling_type == CC_RESAMPLING)
    {
        for (subpixel_line = 0; subpixel_line <= SUBPIXEL_STEPS; 
                subpixel_line++)
        {
            /* Calculate subpixel line step.
               -----------------------------*/
            dy = subpixel_line / (float)(SUBPIXEL_STEPS);

            for (subpixel_sample = 0; subpixel_sample <= SUBPIXEL_STEPS; 
                    subpixel_sample++)
            {
                /* Calculate subpixel sample step.
                   -------------------------------*/ 
                dx = subpixel_sample / (float)(SUBPIXEL_STEPS);

                for (kernel_line = top_line; kernel_line <= bottom_line;
                        kernel_line++)
                {
                    /* Adjust the kernel line for the current subpixel step.  
                       The subpixel step is subtracted since since a
                       subpixel step down makes the weight at the current
                       line position act like it moved up.
                       -----------------------------------------------------*/
                    y = (double)kernel_line - dy;

                    /* Get the weight for the line.
                       ----------------------------*/
                    ccwy = cubic_convolution(alpha, y);

                    for (kernel_sample=left_samp; kernel_sample<=right_samp;
                            kernel_sample++)
                    {
                        /* Adjust the kernel sample location for the current 
                           subpixel step.
                           -------------------------------------------------*/
                        x = (double)kernel_sample - dx;

                        /* Get the weight for the current pixel.
                           -------------------------------------*/
                        *w_ptr++ = ccwy * cubic_convolution(alpha, x);

                    } /* kernel samples */
                } /* kernel lines */
            } /* subpixel samples */
        } /* subpixel lines */
    } /* CC weights */

#ifndef DISABLE_MTF
    else
    {
        /* MTF weights */
        for (subpixel_line=0; subpixel_line<=SUBPIXEL_STEPS; subpixel_line++)
        {
            /* Calculate subpixel line step.
               -----------------------------*/
            dy = subpixel_line / (float)(SUBPIXEL_STEPS);

            for (subpixel_sample=0; subpixel_sample<=SUBPIXEL_STEPS; 
                    subpixel_sample++)
            {
                /* Calculate subpixel sample step.
                   -------------------------------*/ 
                dx = subpixel_sample / (float)(SUBPIXEL_STEPS);

                for (kernel_line=top_line; kernel_line<=bottom_line;
                        kernel_line++)
                {
                    /* Adjust the kernel line for the current subpixel step.
                       The subpixel step is subtracted since since a
                       subpixel step down makes the weight at the current
                       line position act like it moved up.
                       -----------------------------------------------------*/
                    y = (double)kernel_line - dy;

                    /* Get the weight for the line.
                       ----------------------------*/
                    ccwy = mtf_convolution(mtf_across_param_ptr, y);

                    for (kernel_sample=left_samp; kernel_sample<=right_samp; 
                            kernel_sample++)
                    {
                        /* Adjust the kernel sample location for the current 
                           subpixel step.
                           -------------------------------------------------*/
                        x = (double)kernel_sample - dx;

                        /* Get the weight for the current pixel.
                           -------------------------------------*/
                        *w_ptr++ = ccwy*mtf_convolution(mtf_along_param_ptr,
                                x);

                    } /* kernel samples */
                } /* kernel lines */
            } /* subpixel samples */
        } /* subpixel lines */
    } /* MTF weights */
#endif /* disable MTF */

    else                  /* not CC_RESAMPLING */
        return (E_FAIL);

    return (E_SUCC);
}


/*****************************************************************************
NAME:           Kernel_Setup

PURPOSE:                
Kernel_Setup initializes the row weights for creating the extended image and
creates the resampling weight table for the current resampling method.

RETURN VALUE:
Either E_SUCC or E_FAIL are returned, depending on whether the kernel 
setup is successful or not.

NOTES:
Kernel_Setup sets up the resampling weight table for the current band being
processed.


ALGORITHM REFERENCES:
none

*****************************************************************************/

int Kernel_Setup
(
    RESAMPLE_TYPE resampling_type, /* I: resampling type                    */
    double cc_alpha                /* I: cubic convolution alpha parameter  */
)
{
    struct RWTAB table;     /* User entered resampling weight table         */
    int kernel_line, kernel_sample;    /* counters                          */
    int subpixel_line, subpixel_sample;
    int kernel_table_size;  /* number of weights in the kernel weight table */
    double *w_ptr;          /* temporary pointer to kernel weights          */


    /* Verify the kernel size hasn't been accessed before the Kernel_Setup
       is done the first time.
       -------------------------------------------------------------------*/
    ASSERT(kernel_size_accessed != 1);
    kernel_size_accessed |= 2;            /* set bit 1 */

    if (resampling_type == BI_RESAMPLING)
    {
        /* Setup the kernel size for bilinear resampling.
           ----------------------------------------------*/
        lines_in_kernel = BI_LINES_IN_KERNEL;
        samples_in_kernel = BI_SAMPLES_IN_KERNEL;
    }
    else if ((resampling_type == CC_RESAMPLING) || 
            /* presently disabled (resampling_type == MTF_RESAMPLING) || */
            (resampling_type == TABLE_RESAMPLING))
    {
        if (resampling_type == TABLE_RESAMPLING)
        {
            if (c_getrwt((char *)get_table_filename(), &table) != E_SUCC)
                return E_FAIL;

            /* Verify the resampling table is defined in two dimensions.
               ---------------------------------------------------------*/
            if ((table.dim[0] <= 0) || (table.dim[1] <= 0) ||
                    (table.dim[2] > 0))
                return (E_FAIL);                

            /* Initialize the lines and samples in the kernel.
               -----------------------------------------------*/
            lines_in_kernel = table.dim[0];
            samples_in_kernel = table.dim[1];

            /* If the kernel_weight_ptr has already been allocated from a 
               previous band, free it since the size may change for the 
               current band.
               ----------------------------------------------------------*/
            if (kernel_weights_ptr != NULL)
            {
                free(kernel_weights_ptr);
                kernel_weights_ptr = NULL;
            }

            /* Calculate the size of the kernel weight table.
               ----------------------------------------------*/
            kernel_table_size = SUBPIXEL_TABLE_ENTRIES*SUBPIXEL_TABLE_ENTRIES
                *lines_in_kernel*samples_in_kernel; 

            kernel_weights_ptr = (double *)malloc(sizeof(double) * 
                    kernel_table_size);
            if (kernel_weights_ptr == NULL)
                return (E_FAIL);

            /* Create the resampling table by multiplying the two 1-d
               arrays.
               ------------------------------------------------------*/
            w_ptr = kernel_weights_ptr;
            for (subpixel_line = 0; subpixel_line < SUBPIXEL_TABLE_ENTRIES;
                    subpixel_line++)
            {
                for (subpixel_sample = 0; 
                        subpixel_sample < SUBPIXEL_TABLE_ENTRIES; 
                        subpixel_sample++)
                {
                    for (kernel_line = 0; kernel_line < lines_in_kernel;
                            kernel_line++)
                    {
                        for (kernel_sample = 0;
                                kernel_sample < samples_in_kernel;
                                kernel_sample++)
                        {
                            *w_ptr++ = table.table1[subpixel_line][kernel_line]
                                *table.table2[subpixel_sample][kernel_sample];
                        }
                    }
                }
            }
        }
        else
        {
            /* Initialize the lines and samples in the kernel.
               -----------------------------------------------*/
            lines_in_kernel = LINES_IN_KERNEL;
            samples_in_kernel = SAMPLES_IN_KERNEL;

            /* If the kernel_weight_ptr has already been allocated from a 
               previous band, free it since the size may change for the 
               current band.
               ----------------------------------------------------------*/
            if (kernel_weights_ptr != NULL)
            {
                free(kernel_weights_ptr);
                kernel_weights_ptr = NULL;
            }

            /* Calculate the size of the kernel weight table.
               ----------------------------------------------*/
            kernel_table_size = SUBPIXEL_TABLE_ENTRIES*SUBPIXEL_TABLE_ENTRIES
                *LINES_IN_KERNEL*SAMPLES_IN_KERNEL; 

            kernel_weights_ptr = (double *)malloc(sizeof(double) * 
                    kernel_table_size);
            if (kernel_weights_ptr == NULL)
                return (E_FAIL);

            if (Create_Resampling_Kernel(kernel_weights_ptr,
                        resampling_type, cc_alpha) != E_SUCC)
                return (E_FAIL);
        }
    }

    return (E_SUCC);
}



/*****************************************************************************
NAME:           Cleanup_Kernel

PURPOSE:                
Cleanup_Kernel cleans up the effects of the Kernel_Setup routine.

RETURN VALUE:
nothing is returned

NOTES:
Kernel_Setup sets up the resampling weight table for the current band being
processed.


ALGORITHM REFERENCES:
none

*****************************************************************************/

void Cleanup_Kernel(void)
{
    /* Free the space used for the resampling weights.
       -----------------------------------------------*/
    if (kernel_weights_ptr != NULL)
    {
        free (kernel_weights_ptr);
        kernel_weights_ptr = NULL;
    }
}



/*****************************************************************************
NAME:           Get_Resample_Weight_Table_Ptr

PURPOSE:
Get_Resample_Weight_Table_Ptr provides access to the resampling kernel weight
table.

RETURN VALUE:
The pointer to the resampling weight table is returned.

NOTES:
The pointer to the weight table could be passed back to the calling routine
when it is created, but this is a little more efficient.

ALGORITHM REFERENCES:
none

*****************************************************************************/

const double *Get_Resample_Weight_Table_Ptr(void)
{
    ASSERT (kernel_weights_ptr != NULL);
    return (kernel_weights_ptr);
}


/*****************************************************************************
NAME:           get_lines_in_kernel, get_samples_in_kernel

PURPOSE:
The get_lines_in_kernel and get_samples_in_kernel routines provide access 
to the number of lines and samples in the resampling kernel.

RETURN VALUE:
Type = int
The number of lines/samples in the resampling kernel is returned

NOTES:
none

ALGORITHM REFERENCES:
none

*****************************************************************************/

int get_lines_in_kernel(void)
{
    /* Set the accessed flag to allow debugging check.
       -----------------------------------------------*/
    kernel_size_accessed |= 1;

    return (lines_in_kernel);
}


int get_samples_in_kernel(void)
{
    /* Set the accessed flag to allow debugging check.
       -----------------------------------------------*/
    kernel_size_accessed |= 1;

    return (samples_in_kernel);
}

/* These routines return the number of pixels to the left, right, top, and
   bottom of the kernel center.
--------------------------------------------------------------------------*/

int num_left_kernel_samples(void)
{
    return ((samples_in_kernel) - ((samples_in_kernel)/2) - 1);
}

int num_right_kernel_samples(void)
{
    return (samples_in_kernel/2);
}

int num_top_kernel_lines(void)
{
    return ((lines_in_kernel) - ((lines_in_kernel)/2) - 1);
}

int num_bottom_kernel_lines(void)
{
    return (lines_in_kernel/2);
}

double get_kernel_step_size(void)
{
    return (1.0/(float)SUBPIXEL_STEPS);
}

void get_kernel_table_info
(
    int *subpixel_steps, /* O: number of subpixel steps in table            */
    int *table_entries,  /* O: number of table entries for a whole pixel
                               step (one more than subpixel_steps           */
    int *kernel_width,   /* O: width of kernel (samples across kernel)      */
    int *kernel_height   /* O: height of kernel (lines across kernel)       */
)
{
    *subpixel_steps = SUBPIXEL_STEPS;
    *table_entries = SUBPIXEL_TABLE_ENTRIES;
    *kernel_width = samples_in_kernel;
    *kernel_height = lines_in_kernel;
}

