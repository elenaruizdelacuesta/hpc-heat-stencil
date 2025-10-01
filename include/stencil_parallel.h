/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <omp.h>
#include <mpi.h>


#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

// tags for MPI messages
#define TAG_NORTH 0
#define TAG_SOUTH 1
#define TAG_EAST  2
#define TAG_WEST  3

typedef unsigned int uint;

typedef uint    vec2_t[2];
typedef double *restrict buffers_t[4];

typedef struct {
    double   * restrict data;
    vec2_t     size;
} plane_t;



extern int inject_energy ( const int      ,
                           const int      ,
			   const vec2_t  *,
			   const double   ,
                                 plane_t *,
                           const vec2_t   );


extern int update_plane ( const int      ,
                          const vec2_t   ,
                          const plane_t *,
                                plane_t * );


extern int get_total_energy( plane_t *,
                             double  * );

int initialize ( MPI_Comm *,
                 int       ,
		 int       ,
		 int       ,
		 char    **,
                 vec2_t   *,
                 vec2_t   *,                 
		 int      *,
                 int      *,
		 int      *,
		 int      *,
		 int      *,
		 int      *,
                 vec2_t  **,
                 double   *,
                 plane_t  *,
                 buffers_t * );


int memory_release (plane_t *, buffers_t * );


int output_energy_stat ( int      ,
                         plane_t *,
                         double   ,
                         int      ,
                         MPI_Comm *);


/* -----------------------------------------------------------*/
/* inline functions (parallel versions)                       */
/* -----------------------------------------------------------*/

inline int inject_energy ( const int      periodic,
                           const int      Nsources,
			   const vec2_t  *Sources,
			   const double   energy,
                                 plane_t *plane,
                           const vec2_t   N
                           )
{
    const uint register sizex = plane->size[_x_]+2;
    double * restrict data = plane->data;
    
   #define IDX( i, j ) ( (j)*sizex + (i) )
    for (int s = 0; s < Nsources; s++)
        {
            int x = Sources[s][_x_];
            int y = Sources[s][_y_];
            
            data[ IDX(x,y) ] += energy;
            
            if ( periodic )
                {

                    if ( (N[_x_] == 1)  )
                        {
                            // propagate the boundaries if needed
                            // check the serial version
                            if (x == 1)
                                data[ IDX( plane->size[_x_]+1, y) ] += energy;
                            if (x == plane->size[_x_])
                                data[ IDX( 0, y) ] += energy;
                        }
                    
                    if ( (N[_y_] == 1) )
                        {
                            // propagate the boundaries if needed
                            // check the serial version
                            if (y == 1)
                                data[ IDX( x, plane->size[_y_]+1) ] += energy;
                            if (y == plane->size[_y_])
                                data[ IDX( x, 0) ] += energy;
                        }
                }                
        }
 #undef IDX
    
  return 0;
}





inline int update_plane ( const int      periodic, 
                          const vec2_t   N,         // the grid of MPI tasks
                          const plane_t *oldplane,
                                plane_t *newplane
                          )
    
{
    
    uint register fxsize = oldplane->size[_x_]+2;
    
    // internal sizes (excluding ghost cells)
    uint register xsize = oldplane->size[_x_];
    uint register ysize = oldplane->size[_y_];
    
   #define IDX( i, j ) ( (j)*fxsize + (i) ) // index in the flattened array
    
    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization

    double * restrict old = oldplane->data; // actual data
    double * restrict new = newplane->data; // new data

    const double alpha = 0.6;
    const double constant = (1.0 - alpha) / 4.0;

    uint i, j;
    
    #pragma omp parallel for schedule(static)
    for (j = 1; j <= ysize; j++) {
        for (i = 1; i <= xsize; i++)
            {
                
                // NOTE: (i-1,j), (i+1,j), (i,j-1) and (i,j+1) always exist even
                //       if this patch is at some border without periodic conditions;
                //       in that case it is assumed that the +-1 points are outside the
                //       plate and always have a value of 0, i.e. they are an
                //       "infinite sink" of heat
                
                // five-points stencil formula
                //
                // HINT : check the serial version for some optimization
                //
                const double center = old[ IDX(i,j)];
                const double neighbors = old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                                         old[IDX(i, j-1)] + old[IDX(i, j+1)];
                new[ IDX(i,j) ] = center * alpha + neighbors * constant; 
            }
        }

    if ( periodic )
        {
            // if there is only a column of tasks
            if ( N[_x_] == 1 )
                { 
                    // propagate the boundaries as needed
                    // check the serial version
                // copy the values of the first column to the right ghost column (xsize+1)
                // and the values of the last column to the left ghost column (0)
                for (j = 1; j <= ysize; j++ )
                    {
                        new[ IDX(0, j) ] = new[ IDX(xsize, j) ];
                        new[ IDX(xsize+1, j) ] = new[ IDX(1, j) ];
                    }
                }

            // if there is only a row of tasks
            if ( N[_y_] == 1 ) 
                {
                    // propagate the boundaries as needed
                    // check the serial version
                for (i = 1; i <= xsize; i++ )
                    {
                        new[ IDX(i, 0) ] = new[ IDX(i, ysize) ];
                        new[ IDX(i, ysize+1) ] = new[ IDX(i, 1) ];
                    }
                }
        }

    
 #undef IDX
  return 0;
}

inline int update_internal(const plane_t *oldplane, plane_t *newplane) {

    const uint xsize  = oldplane->size[_x_];
    const uint ysize  = oldplane->size[_y_];

    const uint fxsize = xsize + 2;

    #define IDX(i,j) ((j)*fxsize + (i))

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    const double alpha = 0.6;     
    const double constant =  (1.0 - alpha) / 4.0; 

    uint i, j;

    #pragma omp parallel for schedule(static)
    for (j = 2; j <= ysize-1; j++) {
        for (i = 2; i <= xsize-1; i++) {
            const double center    = old[IDX(i,j)];
            const double neighbors = old[IDX(i-1,j)] + old[IDX(i+1,j)] +
                                     old[IDX(i,j-1)] + old[IDX(i,j+1)];
            new[IDX(i,j)] = center * alpha + neighbors * constant;
        }
    }

    #undef IDX
    return 0;
}

inline int update_border(const int periodic, const vec2_t N,
                         const plane_t *oldplane, plane_t *newplane) {
    const uint xsize  = oldplane->size[_x_];
    const uint ysize  = oldplane->size[_y_];

    const uint fxsize = xsize + 2;

    #define IDX(i,j) ((j)*fxsize + (i))

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    const double alpha = 0.6;
    const double constant = (1.0 - alpha) / 4.0;

    uint i, j;

    // Top y bottom rows
    #pragma omp parallel for schedule(static)
    for (i = 2; i <= xsize-1; i++) { // excluding corners
        // top row (j=1)
        double center    = old[IDX(i,1)];
        double neighbors = old[IDX(i-1,1)] + old[IDX(i+1,1)] +
                           old[IDX(i,0)]   + old[IDX(i,2)];
        new[IDX(i,1)] = center * alpha + neighbors * constant;

        // bottom row (j=ysize)
        center    = old[IDX(i,ysize)];
        neighbors = old[IDX(i-1,ysize)] + old[IDX(i+1,ysize)] +
                    old[IDX(i,ysize-1)] + old[IDX(i,ysize+1)];
        new[IDX(i,ysize)] = center * alpha + neighbors * constant;
    }

    // Left y right columns
    #pragma omp parallel for schedule(static)
    for (j = 1; j <= ysize; j++) {
        // left column (i=1)
        double center    = old[IDX(1,j)];
        double neighbors = old[IDX(0,j)] + old[IDX(2,j)] +
                           old[IDX(1,j-1)] + old[IDX(1,j+1)];
        new[IDX(1,j)] = center * alpha + neighbors * constant;

        // right column (i=xsize)
        center    = old[IDX(xsize,j)];
        neighbors = old[IDX(xsize-1,j)] + old[IDX(xsize+1,j)] +
                    old[IDX(xsize,j-1)] + old[IDX(xsize,j+1)];
        new[IDX(xsize,j)] = center * alpha + neighbors * constant;
    }

    // Local periodicity 
    if (periodic) {
        if (N[_x_] == 1) {
            for (j = 1; j <= ysize; j++) {
                new[IDX(0,j)]       = new[IDX(xsize,j)];
                new[IDX(xsize+1,j)] = new[IDX(1,j)];
            }
        }
        if (N[_y_] == 1) {
            for (i = 1; i <= xsize; i++) {
                new[IDX(i,0)]       = new[IDX(i,ysize)];
                new[IDX(i,ysize+1)] = new[IDX(i,1)];
            }
        }
    }

    #undef IDX
    return 0;
}



inline int get_total_energy( plane_t *plane,
                             double  *energy )
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    const int register xsize = plane->size[_x_];
    const int register ysize = plane->size[_y_];
    const int register fsize = xsize+2;

    double * restrict data = plane->data;
    
   #define IDX( i, j ) ( (j)*fsize + (i) )

   #if defined(LONG_ACCURACY)    
    long double totenergy = 0;
   #else
    double totenergy = 0;    
   #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    #pragma omp parallel for reduction(+:totenergy) schedule(static)
    for ( int j = 1; j <= ysize; j++ )
        for ( int i = 1; i <= xsize; i++ )
            totenergy += data[ IDX(i, j) ];

    
   #undef IDX

    *energy = (double)totenergy;
    return 0;
}
