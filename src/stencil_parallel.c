/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_parallel.h"

int test = 0;
int seed = 0;

// Function declaration
double* merged_data(int iter, plane_t *plane, int Rank, int Ntasks, uint *N, const uint [2], MPI_Comm *myCOMM_WORLD);

int dump (double *, const uint [2], const char *, double *, double *);

// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  
  // comm_time accumulates the time spent in communication phases.
  // This includes:
  //   (1) posting non-blocking communications (MPI_Isend / MPI_Irecv),
  //   (2) waiting for their completion (MPI_Waitall).

  // comp_time accumulates the time spent in computation phases.
  // This includes:
  //   (1) updating the internal grid points (update_internal),
  //   (2) updating the border grid points after halo exchange (update_border).
  double start_time_comm, start_time_comp, init_time, total_time;
  double comm_time = 0.0, comp_time = 0.0;

  MPI_Comm myCOMM_WORLD; //
  int  Rank, Ntasks; // my rank and number of tasks
  uint neighbours[4]; // my neighbours in the 4 directions

  int  Niterations;
  int  periodic;
  vec2_t S, N; // size of the plate, and grid of MPI tasks
  
  int      Nsources;
  int      Nsources_local; // number of local sources
  vec2_t  *Sources_local; // coordinates of local sources
  double   energy_per_source;

  plane_t   planes[2];  // two planes, old and new
  buffers_t buffers[2]; // two buffers, send and receive
  
  int output_energy_stat_perstep;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    //
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) {
      printf("MPI_thread level obtained is %d instead of %d\n",
	     level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  }
  
  
  /* argument checking and setting */
  init_time = MPI_Wtime();
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0] );

  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      
      MPI_Finalize();
      return 0;
    }
  
  
  int current = OLD;

  uint ysize = planes[current].size[_y_];
  uint xsize = planes[current].size[_x_];
  uint fxsize = xsize + 2;

  init_time = MPI_Wtime() - init_time;
  total_time = MPI_Wtime();   /* take wall-clock time */
  
  for (int iter = 0; iter < Niterations; ++iter)
    
    {
      double * restrict old = planes[current].data;
      double * restrict new = planes[!current].data;

      MPI_Request reqs[8];
      int nreqs = 0;
      uint j;
      
      /* new energy from sources */
      inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );


      /* -------------------------------------- */



      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
      // pack columns (E/W) and set row pointers (N/S)

      // for EAST and WEST
      if (neighbours[WEST] != MPI_PROC_NULL && buffers[SEND][WEST] != NULL) {
          for (j = 0; j < ysize; j++) {
              // left column (i=1)
              buffers[SEND][WEST][j] = old[(j+1) * fxsize + 1];
          }
      }
      if (neighbours[EAST] != MPI_PROC_NULL && buffers[SEND][EAST] != NULL) {
          for (j = 0; j < ysize; j++) {
              // right column (i=xsize)
              buffers[SEND][EAST][j] = old[(j+1) * fxsize + xsize];
          }
      }

      // for NORTH and SOUTH, we will use direct pointers to contiguous data
      if (neighbours[NORTH] != MPI_PROC_NULL) {
        buffers[RECV][NORTH] = &old[0 * fxsize + 1]; // first ghost row
        buffers[SEND][NORTH] = &old[1 * fxsize + 1]; // first real row
      }

      if (neighbours[SOUTH] != MPI_PROC_NULL) {
        buffers[RECV][SOUTH] = &old[(ysize + 1) * fxsize + 1]; // last ghost row
        buffers[SEND][SOUTH] = &old[ysize * fxsize + 1];       // last real row
      }

      start_time_comm = MPI_Wtime(); 
      // [B] perform the halo communications
      //     (1) use Send / Recv
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      
      // for EAST and WEST, we use the pre-filled buffers
      if (neighbours[EAST] != MPI_PROC_NULL) {
        if (neighbours[EAST] == Rank) {
          //neighbor is myself --> just copy the data
          for (j = 0; j < ysize; j++) {
            buffers[RECV][EAST][j] = buffers[SEND][EAST][j];
          }
        } else {
          MPI_Isend(buffers[SEND][EAST], ysize, MPI_DOUBLE, neighbours[EAST], TAG_EAST, myCOMM_WORLD, &reqs[nreqs++]);
          MPI_Irecv(buffers[RECV][EAST], ysize, MPI_DOUBLE, neighbours[EAST], TAG_WEST, myCOMM_WORLD, &reqs[nreqs++]);
        }
      }

      if (neighbours[WEST] != MPI_PROC_NULL) {
        if (neighbours[WEST] == Rank) {
          for (j = 0; j < ysize; j++) {
            buffers[RECV][WEST][j] = buffers[SEND][WEST][j];
          }
        } else {
          MPI_Isend(buffers[SEND][WEST], ysize, MPI_DOUBLE, neighbours[WEST], TAG_WEST, myCOMM_WORLD, &reqs[nreqs++]);
          MPI_Irecv(buffers[RECV][WEST], ysize, MPI_DOUBLE, neighbours[WEST], TAG_EAST, myCOMM_WORLD, &reqs[nreqs++]);
        }
      }


      if (neighbours[NORTH] != MPI_PROC_NULL) {
        if (neighbours[NORTH] == Rank) {
          for (j = 0; j < xsize; j++) {
            buffers[RECV][NORTH][j] = buffers[SEND][NORTH][j];
          }
        } else {
          MPI_Isend(buffers[SEND][NORTH], xsize, MPI_DOUBLE, neighbours[NORTH], TAG_NORTH, myCOMM_WORLD, &reqs[nreqs++]);
          MPI_Irecv(buffers[RECV][NORTH], xsize, MPI_DOUBLE, neighbours[NORTH], TAG_SOUTH, myCOMM_WORLD, &reqs[nreqs++]);
        }
      }

      if (neighbours[SOUTH] != MPI_PROC_NULL) {
        if (neighbours[SOUTH] == Rank) {
          for (j = 0; j < xsize; j++) {
            buffers[RECV][SOUTH][j] = buffers[SEND][SOUTH][j];
          }
        } else {
          MPI_Isend(buffers[SEND][SOUTH], xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_SOUTH, myCOMM_WORLD, &reqs[nreqs++]);
          MPI_Irecv(buffers[RECV][SOUTH], xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_NORTH, myCOMM_WORLD, &reqs[nreqs++]);
        }
      }
      comm_time += MPI_Wtime() - start_time_comm;
      
      start_time_comp = MPI_Wtime(); 
      update_internal(&planes[current], &planes[!current]); // compute internal grid points (no dependency on halo data yet)
      comp_time += MPI_Wtime() - start_time_comp;

      // wait for all the non-blocking operations to complete
      start_time_comm = MPI_Wtime();
      if (nreqs > 0) {
          MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
      }
      comm_time += MPI_Wtime() - start_time_comm; 



      // [C] copy the haloes data
      //  EAST and WEST
      if (neighbours[WEST] != MPI_PROC_NULL && buffers[RECV][WEST] != NULL) {
          for (j = 0; j < ysize; j++) {
              // left ghost column (i=0)
              old[(j+1) * fxsize + 0] = buffers[RECV][WEST][j];
          }
      }
      if (neighbours[EAST] != MPI_PROC_NULL && buffers[RECV][EAST] != NULL) {
          for (j = 0; j < ysize; j++) {
              // right ghost column (i=xsize+1)
              old[(j+1) * fxsize + (xsize + 1)] = buffers[RECV][EAST][j];
          }
      }


      /* --------------------------------------  */
      
      /* update grid points */
      start_time_comp = MPI_Wtime();
      //update_plane( periodic, N, &planes[current], &planes[!current] );
      update_border( periodic, N, &planes[current], &planes[!current] ); // compute border grid points (requires received halo data)
      comp_time += MPI_Wtime() - start_time_comp;

      /* output if needed */
      #if defined(VERBOSE_LEVEL) && VERBOSE_LEVEL >= 1
              if ( output_energy_stat_perstep )
                  output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
      #endif

      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  total_time = MPI_Wtime() - total_time;

#if defined(VERBOSE_LEVEL) && VERBOSE_LEVEL >= 1
    output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
#endif

  /* release allocated memory */ //MIRAR
  memory_release(planes, buffers );

/* gather maxima across ranks (we report the slowest process values) */
double max_comp_time = 0.0, max_comm_time = 0.0, max_total_time = 0.0;
MPI_Reduce(&comp_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD); 

/* free duplicated communicator (tidy) */
MPI_Comm_free(&myCOMM_WORLD);

/* print only on rank 0 using the reduced (max) values */
if (Rank == 0) {
    /* safety: avoid division by zero */
    double safe_total = (max_total_time > 0.0 ? max_total_time : 1e-30);

    double overhead = max_total_time - (max_comp_time + max_comm_time);
    if (overhead < 0.0) overhead = 0.0; // numerical safety

    printf("=== PERFORMANCE SUMMARY ===\n");
    printf("Total time (max over ranks):          %.6f s\n", max_total_time);
    printf("Computation time (max over ranks):    %.6f s\n", max_comp_time);
    printf("Communication time (max over ranks):  %.6f s\n", max_comm_time);
    printf("Overhead (total - (comp + comm)):     %.6f s\n", overhead);
    printf("Comp/Total ratio:   %.2f%%\n", (max_comp_time / safe_total) * 100.0);
    printf("Comm/Total ratio:   %.2f%%\n", (max_comm_time / safe_total) * 100.0);
    printf("=============================\n");
}

/* optional CSV output for automated tests (only rank 0) */
if (Rank == 0 && test) {
    const char* test_type = getenv("TEST_TYPE"); // must be set by caller
    if (test_type == NULL) {
        fprintf(stderr, "TEST_TYPE not set; skipping CSV output\n");
    } else {
        char filename[256];
        snprintf(filename, sizeof(filename), "data/%s_results.csv", test_type);
        #ifdef _WIN32
          system("mkdir data 2>nul");
        #else
          system("mkdir -p data 2>/dev/null");
        #endif

        FILE *results_file = fopen(filename, "a");
        if (results_file != NULL) {
            /* write header if file empty */
            fseek(results_file, 0, SEEK_END);
            long size = ftell(results_file);
            if (size == 0) {
                fprintf(results_file, "TestType,Nodes,TotalTasks,TasksPerNode,ThreadsPerTask,XDim,YDim,Iterations,TotalTime,ComputationTime,CommunicationTime\n");
            }

            const char* nodes_str = getenv("SLURM_NNODES");
            const char* total_tasks_str = getenv("SLURM_NTASKS");
            const char* tasks_per_node_str = getenv("SLURM_NTASKS_PER_NODE");
            const char* threads_per_task_str = getenv("OMP_NUM_THREADS");

            int nodes = nodes_str ? atoi(nodes_str) : 0;
            int total_tasks = total_tasks_str ? atoi(total_tasks_str) : 0;
            int tasks_per_node = tasks_per_node_str ? atoi(tasks_per_node_str) : 0;
            int threads_per_task = threads_per_task_str ? atoi(threads_per_task_str) : 0;

            fprintf(results_file, "%s,%d,%d,%d,%d,%u,%u,%d,%.6f,%.6f,%.6f\n",
                test_type,
                nodes,
                total_tasks,
                tasks_per_node,
                threads_per_task,
                S[_x_],
                S[_y_],
                Niterations,
                max_total_time,
                max_comp_time,
                max_comm_time );

            fclose(results_file);
        } else {
            fprintf(stderr, "ERROR: cannot open results file %s for writing\n", filename);
        }
    }
}

  MPI_Finalize();
  return 0;
}


/* ==========================================================================
   =                                                                        =
   =   routines called within the integration loop                          =
   ========================================================================== */





/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */


int simple_factorization( uint, int *, uint ** );

int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
			int       ,
			int      *,
			vec2_t  ** );


int memory_allocate ( const int       *,
		      const vec2_t     ,
		            buffers_t *,
		            plane_t   * );
		      

int initialize ( MPI_Comm *Comm,
		 int      Me,                  // the rank of the calling process
		 int      Ntasks,              // the total number of MPI ranks
		 int      argc,                // the argc from command line
		 char   **argv,                // the argv from command line
		 vec2_t  *S,                   // the size of the plane
		 vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
		 int     *periodic,            // periodic-boundary tag
		 int     *output_energy_stat,
		 int     *neighbours,          // four-int array that gives back the neighbours of the calling task
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int     *Nsources_local,
		 vec2_t **Sources_local,
		 double  *energy_per_source,   // how much heat per source
		 plane_t *planes,
		 buffers_t *buffers
		 )
{
  int halt = 0;
  int ret;
  int verbose = 0;
  
  // ··································································
  // set deffault values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 4;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 1000;
  *energy_per_source = 1.0;

  if ( planes == NULL ) {
    // manage the situation
  }

  // initialize planes and buffers structures
  planes[OLD].size[0] = planes[OLD].size[1] = 0;
  planes[NEW].size[0] = planes[NEW].size[1] = 0;

  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL;
  
  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'v': verbose = atoi(optarg);
	    break;

	  case 'h': {
	    if ( Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [1000]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n\n"
		      );
	    halt = 1; }
	    break;
	    
	    
	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;
	    
	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;
  
  
  // ··································································
  /*
   * here we should check for all the parms being meaningful
   *
   */

  // ...

  if (Ntasks <= 0) {
    printf("Error: the number of MPI tasks must be positive\n");
    return 1;
  } else if (Ntasks > (*S)[_x_] * (*S)[_y_]) {
    printf("Error: the number of MPI tasks must be less than or equal to the number of grid points\n");
    return 1;
  }

  if (*Nsources <= 0) {
    printf("Error: the number of heat sources must be positive\n");
    return 1;
  } else if (*Nsources > (*S)[_x_] * (*S)[_y_]) {
    printf("Error: the number of heat sources must be less than or equal to the number of grid points\n");
    return 1;
  }

  if (*Niterations <= 0) {
    printf("Error: the number of iterations must be positive\n");
    return 1;
  } else if (*Niterations > 1000000) {
    printf("Warning: the number of iterations must be less than or equal to 1000000\n");
  }

  if (*energy_per_source <= 0.0) {
    printf("Error: the energy per source must be positive\n");
    return 1;
  }

  if (*periodic != 0 && *periodic != 1) {
    printf("Error: periodic must be either 0 or 1\n");
    return 1;
  }

  if (*output_energy_stat != 0 && *output_energy_stat != 1) {
    printf("Error: output_energy_stat must be either 0 or 1\n");
    return 1;
  }

  // ··································································
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) ); // 1 or 2

  
  if ( dimensions == 1 ) { // means that we can split the plane in a single row or column of tasks
      if ( (*S)[_x_] >= (*S)[_y_] )
	        Grid[_x_] = Ntasks, Grid[_y_] = 1;
      else
	        Grid[_x_] = 1, Grid[_y_] = Ntasks;
  } else { // dimensions == 2
      int   Nf;         // number of factors
      uint *factors;    // the factors
      uint  first = 1;
      ret = simple_factorization( Ntasks, &Nf, &factors );

      if ( ret != 0 ) {
        printf("Error: factorization failed\n");
        return 1;
      }

      for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
          first *= factors[i];

      uint f1 = first;
      uint f2 = Ntasks / first;

      if ( (*S)[_x_] > (*S)[_y_] ) {
          if (f1 >= f2) 
              Grid[_x_] = f1, Grid[_y_] = f2;
          else
              Grid[_x_] = f2, Grid[_y_] = f1;
      } else {
          if (f1 >= f2)
              Grid[_y_] = f1, Grid[_x_] = f2;
          else 
              Grid[_y_] = f2, Grid[_x_] = f1;
      }
      
      if ( factors != NULL ) {
          free(factors);
      }
  }


  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  

  // ··································································
  // my cooridnates in the grid of processors
  // Me = rank MPI (0, .., Ntasks-1)
  // Grid[_x_] = number of tasks in x direction
  // Grid[_y_] = number of tasks in y direction
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  // ··································································
  // find my neighbours
  //

  if ( Grid[_x_] > 1 ) // if there is more than one task in the x direction
    {  
      if ( *periodic ) {    
        neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_];
        neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1);
      } else {
        neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
        neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL );
      }
    }

  if ( Grid[_y_] > 1 ) // if there is more than one task in the y direction
    {
      if ( *periodic ) {      
        neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
        neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks;
      } else {  
        neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
        neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL );
      }
    }

  // ··································································
  // the size of my patch
  //

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_]; 
  mysize[_x_] = s + (X < r); // distribute the remainder
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r); // distribute the remainder

  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  

  if ( verbose > 0 )
    {
      if ( Me == 0 ) {
          printf("Tasks are decomposed in a grid %d x %d\n\n",
            Grid[_x_], Grid[_y_] );
          fflush(stdout);
      }

      MPI_Barrier(*Comm);
      
      for ( int t = 0; t < Ntasks; t++ ) {
          if ( t == Me ) {
              printf("Task %4d :: "
              "\tgrid coordinates : %3d, %3d\n"
              "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
              Me, X, Y,
              neighbours[NORTH], neighbours[EAST],
              neighbours[SOUTH], neighbours[WEST] );
              fflush(stdout);
            }
            MPI_Barrier(*Comm);
      }
      
    }

  
  // ··································································
  // allocate the needed memory
  //
  ret = memory_allocate( neighbours, mysize, buffers, planes );
  if (ret != 0) {
    // an error has occurred during memory allocation
      fprintf(stderr, "Task %d: memory allocation failed\n", Me );
      return 1;
    }

  // ··································································
  // allocate the heat sources
  //
  ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
  if ( ret != 0 ) {
    printf("Error: initialize_sources failed\n");
    return 1;
  }
  
  return 0;  
}


int simple_factorization( uint Ntasks, int *Nfactors, uint **factors )
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 #
 */
{
  int N = 0;
  uint f = 2;
  uint _A_ = Ntasks;

  // first pass: count the number of factors
  uint temp_A = Ntasks;
  while (f*f <= temp_A) {
      while (temp_A % f == 0) {
          N++;
          temp_A /= f;
      }
      f++;
  }
  // last factor
  if (temp_A > 1) {
      N++;
  }
  
  *Nfactors = N;
  if (N == 0) {
      *factors = NULL;
      return 0;
  }
  // allocate memory for the factors
  uint *_factors_ = (uint*)malloc( (size_t)N * sizeof(uint) );
  if (_factors_ == NULL) {
      printf("Error: memory allocation failed\n");
      return 1;
  }
  // second pass: store the factors
  N   = 0;
  f   = 2;
  _A_ = Ntasks;

  while (f*f <= _A_) {
      while (_A_ % f == 0) {
          _factors_[N++] = f;
          _A_ /= f;
      }
      f++;
  }
  if (_A_ > 1) {
      _factors_[N++] = _A_;
  }

  *factors = _factors_;
  return 0;
}

int initialize_sources( int       Me,
			int       Ntasks,
			MPI_Comm *Comm,
			vec2_t    mysize,
			int       Nsources,
			int      *Nsources_local,
			vec2_t  **Sources )

{

  srand48(time(NULL) ^ Me);
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
	tasks_with_sources[i] = (int)lrand48() % Ntasks;
    }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 )
    {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
	{
	  helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	  helper[s][_y_] = 1 + lrand48() % mysize[_y_];
	}

      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}



int memory_allocate ( const int       *neighbours  ,
		      const vec2_t     N           ,
		            buffers_t *buffers_ptr ,
		            plane_t   *planes_ptr
		      )

{
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

      The memory layout that I propose to you is as follows:

      (i) --- calculations
      you need 2 memory regions: the "OLD" one that contains the
      results for the step (i-1)th, and the "NEW" one that will contain
      the updated results from the step ith.

      Then, the "NEW" will be treated as "OLD" and viceversa.

      These two memory regions are indexed by *plate_ptr:

      planew_ptr[0] ==> the "OLD" region
      plamew_ptr[1] ==> the "NEW" region


      (ii) --- communications

      you may need two buffers (one for sending and one for receiving)
      for each one of your neighnours, that are at most 4:
      north, south, east amd west.      

      To them you need to communicate at most mysizex or mysizey
      daouble data.

      These buffers are indexed by the buffer_ptr pointer so
      that

      (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
      (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
      
      --->> Of course you can change this layout as you prefer
      
     */

  if (planes_ptr == NULL ) {
    // an invalid pointer has been passed
    printf("Error: an invalid pointer has been passed to memory_allocate\n");
    return 1;
  }


  if (buffers_ptr == NULL ) {
    // an invalid pointer has been passed
    printf("Error: an invalid pointer has been passed to memory_allocate\n");
    return 1;
  }

  // ··················································
  // allocate memory for data
  // we allocate the space needed for the plane plus a contour frame
  // that will contains data form neighbouring MPI tasks


  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[OLD].data == NULL ) {
    // manage the malloc fail
    printf("Error: memory allocation failed for OLD plane\n");
    return 2;
  }
  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );

  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[NEW].data == NULL ) {
    printf("Error: memory allocation failed for NEW plane\n");
    return 2;
  }
  
  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );


  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you preer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers
  //

  // for north and south I can just use pointers to the correct positions
  (buffers_ptr)[SEND][NORTH] = NULL;
  (buffers_ptr)[SEND][SOUTH] = NULL;

  // for east and west I need to allocate memory
  if ( neighbours[EAST] != MPI_PROC_NULL )
    (buffers_ptr)[SEND][EAST] = (double*)calloc( planes_ptr[OLD].size[_y_], sizeof(double) );
  else
    (buffers_ptr)[SEND][EAST] = NULL;

  if ( neighbours[WEST] != MPI_PROC_NULL )
    (buffers_ptr)[SEND][WEST] = (double*)calloc( planes_ptr[OLD].size[_y_], sizeof(double) );
  else
    (buffers_ptr)[SEND][WEST] = NULL;
  
  (buffers_ptr)[RECV][NORTH] = NULL;
  (buffers_ptr)[RECV][SOUTH] = NULL;

  if ( neighbours[EAST] != MPI_PROC_NULL )
    (buffers_ptr)[RECV][EAST] = (double*)calloc( planes_ptr[OLD].size[_y_], sizeof(double) );
  else
    (buffers_ptr)[RECV][EAST] = NULL;

  if ( neighbours[WEST] != MPI_PROC_NULL )
    (buffers_ptr)[RECV][WEST] = (double*)calloc( planes_ptr[OLD].size[_y_], sizeof(double) );
  else
    (buffers_ptr)[RECV][WEST] = NULL;
  

  // ··················································

  
  return 0;
}



 int memory_release ( plane_t   *planes, buffers_t  *buffers)
  
{

  if ( planes != NULL ) {
      if ( planes[OLD].data != NULL )
          free (planes[OLD].data);
      
      if ( planes[NEW].data != NULL )
          free (planes[NEW].data);
  }

  // free the buffers
  if ( buffers != NULL ) {
    if ( buffers[SEND][EAST] != NULL )
      free ( buffers[SEND][EAST] );

    if ( buffers[SEND][WEST] != NULL )
      free ( buffers[SEND][WEST] );

    if ( buffers[RECV][EAST] != NULL )
      free ( buffers[RECV][EAST] );

    if ( buffers[RECV][WEST] != NULL )
      free ( buffers[RECV][WEST] );
  }
  

  return 0;
}



int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
    {
      if ( step >= 0 ) {
	        printf(" [ step %4d ] ", step ); fflush(stdout);
      }

      printf( "total injected energy is %g, "
	      "system energy is %g "
	      "( in avg %g per grid point)\n",
	      budget,
	      tot_system_energy,
	      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }
  
  return 0;
}

// Merged data
double* merged_data(int iter, plane_t *plane, int Rank, int Ntasks, uint *N, const uint S[2], MPI_Comm *myCOMM_WORLD) {
  // Gather all data to rank 0 for output
  int xsize = plane->size[_x_];
  int ysize = plane->size[_y_];
  int fxsize = xsize + 2;
  size_t local_size = xsize * ysize;
  double *global_grid = NULL;

  // Get the global grid decomposition
  int Grid_x = N[_x_];
  int Grid_y = N[_y_];
  // Get the global grid dimensions
  int global_xsize = S[_x_];
  int global_ysize = S[_y_];

  // Get local data without halos
  double *localbuf = (double*)malloc(local_size * sizeof(double));
  for (int j = 0; j < ysize; j++) {
      for (int i = 0; i < xsize; i++) {
          localbuf[j * xsize + i] = plane->data[(j + 1) * fxsize + (i + 1)];
      }
  }

  // If Ntasks == 1 just return local data as global data
  if (Ntasks == 1) {
      return localbuf;
  }

  double *gathered = NULL;
  int *recvcounts = NULL;
  int *displs = NULL;
  recvcounts = (int*)calloc((size_t)Ntasks, sizeof(int));
  displs = (int*)calloc((size_t)Ntasks, sizeof(int));

  // Rank 0 allocates arrays to gather sizes and data
  if (Rank == 0) {
    gathered = (double*)malloc((size_t)(global_xsize * global_ysize) * sizeof(double));
    // Compute the sizes of each rank's patch
    int offset = 0;
    for (int r = 0; r < Ntasks; r++) {
        int rx = r % Grid_x;
        int ry = r / Grid_x;
        int sx = global_xsize / Grid_x + (rx < (global_xsize % Grid_x));
        int sy = global_ysize / Grid_y + (ry < (global_ysize % Grid_y));
        // Compute the displacements and counts for Gatherv
        displs[r] = offset;
        recvcounts[r] = sx * sy;
        offset += sx * sy;
    }
  }
  MPI_Gatherv(localbuf, (int)local_size, MPI_DOUBLE,
              gathered, recvcounts, displs, MPI_DOUBLE,
              0, *myCOMM_WORLD);
  free(localbuf);

  // Rank 0 rearranges the gathered data into the correct global grid
  if (Rank == 0) {
      global_grid = (double*)malloc((size_t)(global_xsize * global_ysize) * sizeof(double));
      
      // Rearrange each rank's block into its correct position
      for (int r = 0; r < Ntasks; r++) {
          int rx = r % Grid_x;
          int ry = r / Grid_x;
          int sx = global_xsize / Grid_x + (rx < (global_xsize % Grid_x));
          int sy = global_ysize / Grid_y + (ry < (global_ysize % Grid_y));
          
          // Compute the starting indices in the global grid
          int startx = 0;
          for (int i = 0; i < rx; i++) {
              startx += global_xsize / Grid_x + (i < (global_xsize % Grid_x));
          }
          int starty = 0;
          for (int j = 0; j < ry; j++) {
              starty += global_ysize / Grid_y + (j < (global_ysize % Grid_y));
          }

          // Copy the patch into the correct position in the global grid
          for (int j = 0; j < sy; j++) {
              for (int i = 0; i < sx; i++) {
                  global_grid[(starty + j) * global_xsize + (startx + i)] = gathered[displs[r] + j * sx + i];
              }
          }
      }
      free(gathered);
      free(recvcounts);
      free(displs);
  }
  MPI_Barrier(*myCOMM_WORLD);
  return global_grid;
}

int dump (double *data, const uint size[2], const char *filename, double *min, double *max){
    if ( (filename != NULL) && (filename[0] != '\0') ) {
        FILE *f = fopen(filename, "wb");
        if (f == NULL) {
          return 2;
        }
        fwrite(size, sizeof(uint), 2, f);

        float *array = (float*)malloc( size[0] * sizeof(float));

        double _min_ = DBL_MAX;
        double _max_ = -DBL_MAX;

        for ( int j = 0; j < size[1]; j++ ) {
          const double * restrict line = data + j*size[0];
          for ( int i = 0; i < size[0]; i++ ) {
            array[i] = (float)line[i];
            _min_ = (line[i] < _min_ ? line[i] : _min_);
            _max_ = (line[i] > _max_ ? line[i] : _max_);
          }
          fwrite(array, sizeof(float), size[0], f);
        }
        free(array);
        if (min != NULL) {
          *min = _min_;
        }
        if (max != NULL) {
          *max = _max_;
        }
        fclose(f);
        return 0;
      }
      else return 1;
    }