/*****************************************************************************
NAME: cc_resample

PURPOSE:        
The cc_resample module implements the routines needed to perform cubic 
convolution resampling.


ROUTINES:
    
    initialize_cc_resample
    cc_resample
    cc_resample_nogrid
    

DEVELOPMENT HISTORY:


NOTES:
- The user must first initialize the module by calling initialize_cc_resample
for each new band resampled.
- The static variables located in this module are cached values that should
not be written to by the actual resampling routine.  If they are written to by
the resampling routine, it will cause unpredictable results when
multithreading.  If this caching was not used, there would be a performance
penalty.

ALGORITHM REFERENCES:
None

*****************************************************************************/

/*
**
** INCLUDE FILES 
**
*/

#include "resample.h"
#include "debug.h"       /* ASSERT definition                              */

/*
**
** LOCAL STATIC VARIABLES
**
*/

/* NOTE:  These variables are initialized by the initialize_cc_resample
          routine and should be variables that are simply READ (and not
          written) by each resampling thread.  Any variables that are written
          by a resampling thread should have a separate copy for each thread.
*/

static const gxx_GEOM_GRID_BAND_TYPE *grid_band_ptr;
                                /* pointer to the current band of the grid  */
static struct PROJ *inproj, *outproj; /* input and output projection params */
static int *sample_reached_ptr; /* pointer to the array that tracks what the
                                   last sample number reached was for each
                                   line                                     */
static int input_image_samples; /* number of samples across the input image */
static int do_terrain;          /* flag that indicates terrain correction
                                   should be performed                      */
static const SCAN_MAN_TYPE *scan_man_ptr;
                                /* pointer to te scan manager structure for
                                   the input image                          */
static double fill;             /* pixel fill value                         */
static int output_window_last_sample;
                                /* last sample number contained in the output
                                   window                                   */
static double round;            /* rounding error in resampling kernel      */
static double min_range,        /* local versions of the minimum and maximum*/
              max_range;        /* output range, adjusted for the currently
                                   selected output type.                    */
static const double *kernel_weight_table_ptr;
                                /* pointer to kernel weight table           */
static int kernel_samples;      /* resampling kernel information            */
static int kernel_lines;
static int kernel_subpixel_steps;
static double doub_kernel_subpixel_steps;
static int kernel_table_steps;
static int kernel_size;         /* kernel samples * lines                   */
static int left_kernel_samps;   /* samples to the left of the kernel center */
static int top_kernel_lines;    /* lines to the top of the kernel center    */
static double min_in_line,      /* minimum and maximum valid line/sample    */
              max_in_line,      /* numbers                                  */
              min_in_sample,
              max_in_sample;
static const double *ipixsize,  /* pixel size array ptr(line, samp)(in, out)*/
                    *opixsize;                      
static const double *iulcoor,   /* UL corner coor ptr(lat, long) (in, out)  */
                    *oulcoor;

/*****************************************************************************
FUNCTION NAME:  initialize_cc_resample

PURPOSE:        
initialize_cc_resample initializes the cached static variable in the module.

RETURNS:
Nothing is returned by this routine.

NOTES:
See notes above.

ALGORITHM REFERENCES:
none

*****************************************************************************/

void initialize_cc_resample
(
    const gxx_GEOM_GRID_BAND_TYPE *i_grid_band_ptr,
                               /* I: current band of grid                   */
    const SCAN_MAN_TYPE *i_scan_man_ptr, /* I: input image scan manager     */
    int *i_sample_reached_ptr, /* I: ptr to the last sample reached array   */
    int i_input_image_samples, /* I: samples across the input image         */
    int i_output_window_last_sample, /* I: last sample (not included) in
                                           output window                    */
    double i_min_range,        /* I: minimum output range                   */
    double i_max_range         /* I: maximum output range                   */
)
{
    const WINDOW_TYPE *input_window_ptr; /* input window                    */
    long prtprm[2]={FALSE, FALSE};     /* print dest. flags for c_transinit */


    /* Copy the parameters passed in to their local versions.
       ------------------------------------------------------*/
    if (i_grid_band_ptr)
        grid_band_ptr = i_grid_band_ptr;
    else
    {
        /* No grid, so initialize the projection transformation.
           -----------------------------------------------------*/
        opixsize = get_opixsize();
        ipixsize = get_ipixsize();
        oulcoor = get_oulcoor();
        iulcoor = get_iulcoor();
        inproj = (struct PROJ *)get_inproj();
        outproj = (struct PROJ *)get_outproj();
        c_transinit(&outproj->proj_code, &outproj->units, &outproj->zone_code,
                &outproj->datum_code, outproj->proj_coef,
                &inproj->proj_code, &inproj->units, &inproj->zone_code,
                &inproj->datum_code, inproj->proj_coef, prtprm, NULL);
    }
    sample_reached_ptr = i_sample_reached_ptr;
    input_image_samples = i_input_image_samples;
    scan_man_ptr = i_scan_man_ptr;
    min_range = i_min_range;
    max_range = i_max_range;
    output_window_last_sample = i_output_window_last_sample;

    /* Cache the other local values.
       -----------------------------*/
    fill = get_fill_pixel();
    do_terrain = do_terrain_correction();
    round = 1.0/64.0;

    /* Setup the kernel weight table.
       ------------------------------*/
    kernel_weight_table_ptr = Get_Resample_Weight_Table_Ptr();
    get_kernel_table_info(&kernel_subpixel_steps, &kernel_table_steps,
            &kernel_samples, &kernel_lines);
    left_kernel_samps = num_left_kernel_samples();
    top_kernel_lines = num_top_kernel_lines();
    doub_kernel_subpixel_steps = (double)kernel_subpixel_steps;
    kernel_size = kernel_samples * kernel_lines;

    /* Set the input image limits.
       ---------------------------*/
    input_window_ptr = get_input_window();
    min_in_line = input_window_ptr->start_line + top_kernel_lines;
    max_in_line = input_window_ptr->end_line - num_bottom_kernel_lines();
    min_in_sample = input_window_ptr->start_sample + left_kernel_samps;
    max_in_sample = input_window_ptr->end_sample - num_right_kernel_samples();

    /* Verify the kernel size is correct for the loop that 
       actually does the resampling, which is hard coded to use
       a 4 line by 4 sample kernel for speed.
       --------------------------------------------------------*/
    ASSERT((kernel_lines == 4) && (kernel_samples == 4));

}


/*****************************************************************************
FUNCTION NAME:  cc_resample

PURPOSE:        
cc_resample is the actual cubic convolution resampling routine.  The routine 
attempts to resample an entire output line for each call.  The routine is 
re-entrant so multiple threads can call it.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
0               no input image data used
non-zero        input image data used
The resampled data and the last sample created are also returned.

NOTES:
See notes above.

ALGORITHM REFERENCES:
none

*****************************************************************************/

int cc_resample
(
    int output_line,         /* I: output space line number to create       */
    double *buffer_ptr,      /* I/O: buffer for resampled data              */
    int *ending_sample_ptr,  /* O: last sample reached in the line          */
    double *minpix,          /* I/O: minimum pixel value for band           */
    double *maxpix           /* I/O: maximum pixel value for band           */
)
{
    COEFFICIENTS *coeff_ptr;/* pointer to reverse coefficients in grid      */

    int grid_row_index;     /* index to start of row in grid coefficients   */
    int grid_cell_index;    /* index to cell in row for grid coefficients   */
    int input_line;         /* input space line number                      */
    int input_sample;       /* input space sample number                    */
    int kernel_line, kernel_sample; /* index into the kernel table for the 
                                       fractional input line/sample         */
    int line_counter;       /* resampling loop counter                      */ 
    int output_sample;      /* output space sample number                   */

    double d_input_line;    /* double input space line number               */
    double d_input_sample;  /* double input space sample number             */
    double d_output_line;   /* double output space line number              */
    double d_output_sample; /* double output space sample number            */
    double ds, dl;          /* fractional input space line and sample       */
    double minz = *minpix;  /* local variables for min/max pixel values     */
    double maxz = *maxpix;
    double total;           /* resampling total                             */
    const double *w_ptr;    /* pointer to current resampling weights        */

    const float *data_ptr;  /* Input data pointer                           */


    /* Convert the output line to a double.
       ------------------------------------*/
    d_output_line = output_line;

    /* Convert the output line into a grid row.
       ----------------------------------------*/
    grid_row_index = output_line/grid_band_ptr->cell_lines * 
        grid_band_ptr->grid_cols;

    /* Set the line counter to zero so it can be used to detect when new 
       data is resampled.  Normally a separate flag would have to be
       created, but since the line counter is only modified when new data is
       created, it can be used instead.
       ---------------------------------------------------------------------*/
    line_counter = 0;

    /* Loop through the samples in the output window.
       ----------------------------------------------*/
    for (output_sample = sample_reached_ptr[output_line]; 
            output_sample < output_window_last_sample; 
            output_sample++)
    {
        d_output_sample = output_sample;

        /* Convert the output sample into a grid column.
           ---------------------------------------------*/
        grid_cell_index = grid_row_index + 
            output_sample/grid_band_ptr->cell_samps;

        /* Get the pointer to the grid coefficients.
           -----------------------------------------*/
        coeff_ptr = &grid_band_ptr->reverse_coeffs[grid_cell_index];

        /* Using the grid's reverse coefficients, calculate the input
           space line/sample.
           ----------------------------------------------------------*/
        d_input_sample = coeff_ptr->a[0] +
            coeff_ptr->a[1]*d_output_sample + 
            coeff_ptr->a[2]*d_output_line + 
            coeff_ptr->a[3]*d_output_sample*d_output_line;
        d_input_line = coeff_ptr->b[0] +
            coeff_ptr->b[1]*d_output_sample + 
            coeff_ptr->b[2]*d_output_line + 
            coeff_ptr->b[3]*d_output_sample*d_output_line;

        /* Add terrain effects if needed.
           ------------------------------*/
        if (do_terrain)
            d_input_sample += get_terrain_adjustment(d_input_line, 
                    d_input_sample, output_line, output_sample);

        /* If the input line/sample are within the input image bounds...
           -------------------------------------------------------------*/
        if ((d_input_line >= min_in_line) && (d_input_line < max_in_line) &&
                (d_input_sample >= min_in_sample) &&
                (d_input_sample < max_in_sample))
        {
            /* Convert the input line/sample into whole and fractional 
               line/sample.
               -------------------------------------------------------*/
            input_line = (int)d_input_line;
            input_sample = (int)d_input_sample;
            dl = d_input_line - input_line;
            ds = d_input_sample - input_sample;

            /* Get the resampling kernel index for the fractional
               line/sample.
               --------------------------------------------------*/
            kernel_line = (int)(doub_kernel_subpixel_steps * (round + dl));
            kernel_sample = (int)(doub_kernel_subpixel_steps * (round + ds));

            /* Get pointer to weight table location to use.
               --------------------------------------------*/ 
            w_ptr = kernel_weight_table_ptr + (kernel_line*kernel_table_steps
                    + kernel_sample)*kernel_size;

            /* Get the address of the data at input line/sample.
               -------------------------------------------------*/
            data_ptr = Get_Sample_Data(scan_man_ptr, 
                    input_line - top_kernel_lines,
                    input_sample - left_kernel_samps);

            /* If the data pointer is NULL (i.e. data not in memory),
               break out of the loop.
               ------------------------------------------------------*/
            if (data_ptr == NULL)
                break;

            /*-------------------
              -------------------
              Perform resampling.
              -------------------
              ---------------------*/

            total = 0.0;     /* initialize total to zero */

            /* Loop through the kernel lines for resampling.
               ---------------------------------------------*/
            for (line_counter = -1; line_counter < 3; line_counter++)
            {
                /* Unroll the inner loop of kernel(from samples -2 to < 4).
                   Assumes a 4 x 4 kernel.
                   --------------------------------------------------------*/
                total += w_ptr[0] * data_ptr[0] + w_ptr[1] * data_ptr[1] +
                    w_ptr[2] * data_ptr[2] + w_ptr[3] * data_ptr[3]; 

                /* Advance weight pointer and input sample pointer.
                   ------------------------------------------------*/
                w_ptr = &w_ptr[4];
                data_ptr += input_image_samples;
            }

            /* Put new resampled value in intermediate buffer, limiting it 
               to the output range.
               -----------------------------------------------------------*/
            if (total < min_range)
                total = min_range;
            else if (total > max_range)
                total = max_range;

            /* Update the min/max pixel values.  (Global min/max pixel
               values are initialized to the fill value in main.)
               -------------------------------------------------------*/
            if (total < minz)
                minz = total;
            else if (total > maxz)
                maxz = total;

        } /* in image */
        else
        {
            /* Set the output pixel to the fill value.
               ---------------------------------------*/
            total = fill;
        }

        /* Update the buffer pointer and increment to fill next sample.
           ------------------------------------------------------------*/
        *buffer_ptr = total;
        buffer_ptr++;

    } /* end of line loop */

    /* Return the ending sample and min/max pixel values.
       --------------------------------------------------*/
    *ending_sample_ptr = output_sample;
    *minpix = minz;
    *maxpix = maxz;

    /* Return the data used indication.  (Line counter not zero means
       data used.)
       --------------------------------------------------------------*/
    if (line_counter)
        return (1);
    else
        return (0);
}


/*****************************************************************************
FUNCTION NAME:  cc_resample_nogrid

PURPOSE:        
cc_resample_nogrid is the actual cubic convolution resampling routine.  The
routine attempts to resample an entire output line for each call.  The routine
is re-entrant so multiple threads can call it.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
0               no input image data used
non-zero        input image data used
The resampled data and the last sample created are also returned.

NOTES:
See notes above.

ALGORITHM REFERENCES:
none

*****************************************************************************/

int cc_resample_nogrid
(
    int output_line,         /* I: output space line number to create       */
    double *buffer_ptr,      /* I/O: buffer for resampled data              */
    int *ending_sample_ptr,  /* O: last sample reached in the line          */
    double *minpix,          /* I/O: minimum pixel value for band           */
    double *maxpix           /* I/O: maximum pixel value for band           */
)
{
    int input_line;         /* input space line number                      */
    int input_sample;       /* input space sample number                    */
    int kernel_line, kernel_sample; /* index into the kernel table for the 
                                       fractional input line/sample         */
    int line_counter;       /* resampling loop counter                      */ 
    int output_sample;      /* output space sample number                   */
    int status;             /* function return status                       */

    double d_input_line;    /* double input space line number               */
    double d_input_sample;  /* double input space sample number             */
    double input_x;         /* input space projection x coordinate          */
    double input_y;         /* input space projection y coordinate          */
    double d_output_line;   /* double output space line number              */
    double d_output_sample; /* double output space sample number            */
    double output_x;        /* output space projection x coordinate         */
    double output_y;        /* output space projection y coordinate         */
    double ds, dl;          /* fractional input space line and sample       */
    double total;           /* resampling total                             */
    double minz = *minpix;  /* local variables for min/max pixel values     */
    double maxz = *maxpix;
    const double *w_ptr;    /* pointer to current resampling weights        */

    const float *data_ptr;  /* Input data pointer                           */


    /* Convert the output line to a double.
       ------------------------------------*/
    d_output_line = output_line;

    output_y = oulcoor[0] - d_output_line*opixsize[0];

    /* Set the line counter to zero so it can be used to detect when new 
       data is resampled.  Normally a separate flag would have to be
       created, but since the line counter is only modified when new data is
       created, it can be used instead.
       ---------------------------------------------------------------------*/
    line_counter = 0;

    /* Loop through the samples in the output window.
       ----------------------------------------------*/
    for (output_sample = sample_reached_ptr[output_line]; 
            output_sample < output_window_last_sample; 
            output_sample++)
    {
        d_output_sample = output_sample;

        output_x = oulcoor[1] + d_output_sample*opixsize[1];

        status = c_trans(&outproj->proj_code, &outproj->units,
                &inproj->proj_code, &inproj->units,
                &output_x, &output_y, &input_x, &input_y);

        /* If transformation not possible (possibly in break region),
           use fill value for pixel.
           ----------------------------------------------------------*/
        if (status != E_SUCC)
        {
            *buffer_ptr = fill;

            /* Increment the buffer pointer to fill next sample.
               -------------------------------------------------*/
            buffer_ptr++;

            continue;
        }

        d_input_line = (iulcoor[0] - input_y)/ipixsize[0];
        d_input_sample = (input_x - iulcoor[1])/ipixsize[1];

        /* Add terrain effects if needed.
           ------------------------------*/
        if (do_terrain)
            d_input_sample += get_terrain_adjustment(d_input_line, 
                    d_input_sample, output_line, output_sample);

        /* If the input line/sample are within the input image bounds...
           -------------------------------------------------------------*/
        if ((d_input_line >= min_in_line) && (d_input_line < max_in_line) &&
                (d_input_sample >= min_in_sample) &&
                (d_input_sample < max_in_sample))
        {
            /* Convert the input line/sample into whole and fractional 
               line/sample.
               -------------------------------------------------------*/
            input_line = (int)d_input_line;
            input_sample = (int)d_input_sample;
            dl = d_input_line - input_line;
            ds = d_input_sample - input_sample;

            /* Get the resampling kernel index for the fractional
               line/sample.
               --------------------------------------------------*/
            kernel_line = (int)(doub_kernel_subpixel_steps * (round + dl));
            kernel_sample = (int)(doub_kernel_subpixel_steps * (round + ds));

            /* Get pointer to weight table location to use.
               --------------------------------------------*/ 
            w_ptr = kernel_weight_table_ptr + (kernel_line*kernel_table_steps
                    + kernel_sample)*kernel_size;

            /* Get the address of the data at input line/sample.
               -------------------------------------------------*/
            data_ptr = Get_Sample_Data(scan_man_ptr, 
                    input_line - top_kernel_lines,
                    input_sample - left_kernel_samps);

            /* If the data pointer is NULL (i.e. data not in memory),
               break out of the loop.
               ------------------------------------------------------*/
            if (data_ptr == NULL)
                break;

            /*-------------------
              -------------------
              Perform resampling.
              -------------------
              ---------------------*/

            total = 0.0;     /* initialize total to zero */

            /* Loop through the kernel lines for resampling.
               ---------------------------------------------*/
            for (line_counter = -1; line_counter < 3; line_counter++)
            {
                /* Unroll the inner loop of kernel(from samples -2 to < 4).
                   --------------------------------------------------------*/
                total += w_ptr[0]*data_ptr[0] + w_ptr[1]*data_ptr[1] +
                    w_ptr[2]*data_ptr[2] + w_ptr[3]*data_ptr[3]; 

                /* Advance weight pointer and input sample pointer.
                   ------------------------------------------------*/
                w_ptr = &w_ptr[4];
                data_ptr += input_image_samples;
            }

            /* Put new resampled value in intermediate buffer, limiting it 
               to the output range.
               -----------------------------------------------------------*/
            if (total < min_range)
                total = min_range;
            else if (total > max_range)
                total = max_range;

            /* Update the min/max pixel values.  (Global min/max pixel
               values are initialized to the fill value in main.)
               -------------------------------------------------------*/
            if (total < minz)
                minz = total;
            else if (total > maxz)
                maxz = total;

        } /* in image */
        else
        {
            /* Set the output pixel to the fill value.
               ---------------------------------------*/
            total = fill;
        }

        /* Update the buffer pointer and increment to fill next sample.
           ------------------------------------------------------------*/
        *buffer_ptr = total;
        buffer_ptr++;

    } /* end of line loop */

    /* Return the ending sample and min/max pixel values.
       --------------------------------------------------*/
    *ending_sample_ptr = output_sample;
    *minpix = minz;
    *maxpix = maxz;

    /* Return the data used indication.  (Line counter not zero means
       data used.)
       --------------------------------------------------------------*/
    if (line_counter)
        return (1);
    else
        return (0);
}
