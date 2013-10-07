//#include "functions.h"
#include "Defaults.h" // For Life's constants

#include <time.h>     // For seeding random
#include <stdlib.h>   // For malloc et al.
#include <stdbool.h>  // For true/false
#include <getopt.h>   // For argument processing
#include <stdio.h>    // For file i/o


int         initialize (struct data * vars, int * c, char *** v);
void        eval_grid (struct data * vars);
void        send_receive (struct data * vars);
void        copy_grid (struct data * vars);
void    	allocate (struct data * vars);
void        init_grids (struct data * vars);
void        write_alive_cells (struct data * vars);
void        free_grids (struct data * vars);
double      rand_double ();
void    	rand_grid (struct data * vars, double prob);
void        seed_random (int my_rank);
void        cleanup (struct data * vars);
void        parse_args (struct data * vars, int argc, char ** argv);
void        help ();
void get_surrounding_ranks(struct data * vars);
int terminate(struct data *vars);

int main(int argc, char ** argv) {

	int generations;
	int foo = 0;
//	int threshold = 0;
	struct data vars;
//	int i;
	vars.print = 0;
	if( initialize(&vars, &argc, &argv) == -1 )
		exit(EXIT_FAILURE);
	MPI_Barrier(vars.my_comm);
	double time =  MPI_Wtime();
	for (generations = 0; generations < vars.generations; generations++) {

		send_receive(&vars);

		eval_grid(&vars);
		copy_grid(&vars);
	  foo = terminate(&vars);
	  if(foo){
			printf("Dead grid %d exiting\n",vars.my_rank);
			break;
		}
		MPI_Barrier(vars.my_comm);
	}
	MPI_Barrier(vars.my_comm);
	time = MPI_Wtime() - time;
	printf("Elapsed time = %10f\n",time);

	//vars.changes = generations;
	cleanup(&vars);
	exit(EXIT_SUCCESS);
}

int terminate(struct data *vars){
	int i,error;
	int terminate = 0;
//	int msg_exists = 0;
	int ter_count = 0;
	MPI_Request d_req;

	if(vars->changes != 0 && vars->my_rank == MY_MPI_ROOT_PROC){//if there were changes and you are root
		//then just discard all inbox
			for( i =1; i < vars->mpi_size; i++){
					error =  MPI_Irecv(&error, 1, MPI_INT, i, FROM_MPI,
										vars->my_comm, &d_req);
					error =  MPI_Isend(&terminate, 1, MPI_INT, i, FROM_MPI,vars->my_comm,&d_req);
			}
	}else if(vars->my_rank != MY_MPI_ROOT_PROC){
		if(vars->changes != 0)
			terminate = 0;
		else
			terminate = 1;
		error = MPI_Isend(&terminate, 1, MPI_INT, MY_MPI_ROOT_PROC, FROM_MPI, vars->my_comm,&d_req);

		error =  MPI_Irecv(&terminate, 1, MPI_INT, MY_MPI_ROOT_PROC, FROM_MPI, vars->my_comm, &d_req);

	}else if(vars->changes == 0 && vars->my_rank == MY_MPI_ROOT_PROC){ //if no changes
		ter_count = 1;
		for( i =1; i < vars->mpi_size; i++){
			terminate = 0;
			error =  MPI_Irecv(&terminate, 1, MPI_INT, i, FROM_MPI, vars->my_comm, &d_req);
			if(terminate)
				ter_count++;
		}
		if( ter_count == vars->mpi_size){// if everyone sent no changes then
			terminate = 1;
		}else{
			terminate = 0;
		}
		for( i =1; i < vars->mpi_size; i++){
			error =  MPI_Isend(&terminate, 1, MPI_INT, i, FROM_MPI,vars->my_comm,&d_req);
		}
	}
	vars->changes = 0;
	return terminate;
}
/*
	initialize()
		Initialize runtime environment and initializes MPI.
*/
int initialize (struct data * vars, int * c, char *** v) {
	int *argc          = c;
	char *** argv      = v;
	int dim_size[2];
	int periods[2];
	int error = 0;
	MPI_Comm tmp_comm;

	vars->colno       = DEFAULT_SIZE;
	vars->rowno       = DEFAULT_SIZE;
	vars->generations = DEFAULT_GENS;
	vars->my_rank        = -1;
	vars->mpi_size       = -1;
	vars->infile      = NULL;
	vars->outfile     = NULL;

	error = MPI_Init(argc, argv);
	if(error != 0){
		fprintf(stderr, "\nError in MPI_Init\n");
		return -1;
	}
	error = MPI_Comm_rank(MPI_COMM_WORLD, &vars->my_rank);
	if(error != 0){
		fprintf(stderr, "Error in MPI_Comm_rank\n");
		return -1;
	}
	error = MPI_Comm_size(MPI_COMM_WORLD, &vars->mpi_size);
	if(error != 0){
		fprintf(stderr, "Error in MPI_Comm_size\n");
		return -1;
	}
	int size = vars->mpi_grid_size = sqrt(vars->mpi_size);
	int ndims = MY_MPI_DIMS;				/* 2D matrix */
	dim_size[0] = size;					/* N rows */
	dim_size[1] = size;					/* N columns */
	periods[0] = 1;						/* row periodic */
	periods[1] = 1;						/* column periodic */
	int reorder = 1;					/* allows processes reordered for efficiency */

	error = MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &tmp_comm);
	if(error != 0)
		fprintf(stderr, "Error in MPI_Cart_create");
	vars->my_comm = tmp_comm;

	get_surrounding_ranks(vars);

	seed_random(vars->my_rank);
	parse_args(vars, *argc, *argv);
	init_grids(vars);
	return 0;
}

void get_surrounding_ranks(struct data * vars){
	int error;
	int coords[2];
	int left;
	int right;
	int top;
	int bottom;
	int t_left;
	int t_right;
	int b_left;
	int b_right;

	MPI_Comm my_comm = vars->my_comm;

	/* top - bottom */
	error =  MPI_Cart_shift(my_comm, 0, 1, &top, &bottom);
	if( error != 0)
		fprintf(stderr,"Error in MPI_Cart_shift\n");
	/* left - right */
	error =	 MPI_Cart_shift(my_comm, 1, 1, &left, &right);
	if( error != 0)
		fprintf(stderr,"Error in MPI_Cart_shift\n");
	/*top - right*/
	error = MPI_Cart_coords(my_comm, top , MY_MPI_DIMS, coords);
	if(coords[1] == vars->mpi_grid_size - 1)
		coords[1] = 0;
	else
		coords[1]++;
	error = MPI_Cart_rank( my_comm, coords, &t_right );

	/* top - left*/
	error = MPI_Cart_coords(my_comm, top , MY_MPI_DIMS, coords);
	if(coords[1] == 0)
		coords[1] = vars->mpi_grid_size - 1;
	else
		coords[1]--;
	error = MPI_Cart_rank( my_comm, coords, &t_left );

	/*bottom -left*/
	error = MPI_Cart_coords(my_comm, bottom , MY_MPI_DIMS, coords);
	if(coords[1] == 0)
		coords[1] = vars->mpi_grid_size - 1;
	else
		coords[1]--;
	error = MPI_Cart_rank( my_comm, coords, &b_left );

	/* bottom - right*/
	error = MPI_Cart_coords(my_comm, bottom , MY_MPI_DIMS, coords);
	if(coords[1] == vars->mpi_grid_size - 1)
		coords[1] = 0;
	else
		coords[1]++;
	error = MPI_Cart_rank( my_comm, coords, &b_right );

	vars->left = left;
	vars->right = right;
	vars->top = top;
	vars->bot = bottom;
	vars->t_left = t_left;
	vars->t_right = t_right;
	vars->b_left = b_left;
	vars->b_right = b_right;

	/*if(vars->my_rank == 17)
		printf(" top: %d bottom: %d left: %d right: %d t_l: %d t_r: %d b_l: %d b_r: %d\n",
		top,bottom,left,right,t_left,t_right,b_left,b_right);

	if(vars->my_rank == 30){
		int my_coords[2];
		error = MPI_Cart_coords(my_comm, vars->my_rank , MY_MPI_DIMS, my_coords);
		printf("err = %d  mpi_size = %d X: %d, Y: %d , rank = %d my bottom = %d, i'm %d my coords are X = %d Y = %d\n",
		error, vars->mpi_size, coords[0],coords[1], b_right, bottom ,vars->my_rank, my_coords[0], my_coords[1]);
	//MPI_Cart_rank( my_comm, coords, ranks[1] );
	}*/
}

/*
	Evaluate the rules of life for each cell count
	neighbors and update current state accordingly.
*/
void eval_grid (struct data * vars) {
	int i,j,k,l,neighbors;

	int colno = vars->colno;
	//int rowno = vars->rowno;

	int ** grid      = vars->c_grid;
	int ** f_grid = vars->f_grid;
	vars->changes = 0;

	for (i = 1; i <= colno; i++) {
		for (j = 1; j <= colno; j++) {
			neighbors = 0;

			// count neighbors
			for (k = i-1; k <= i+1; k++) {
				for (l = j-1; l <= j+1; l++) {
					if (!(k == i && l == j) && grid[k][l] != DEAD)
						neighbors++;
				}
			}

			// update state
			if (neighbors < LOWER_THRESH || neighbors > UPPER_THRESH){
				f_grid[i][j] = DEAD;
				if(grid[i][j] != DEAD)
					vars->changes++;
			}
			else if (grid[i][j] != DEAD || neighbors == SPAWN_THRESH){
				f_grid[i][j] = ALIVE;
				if(grid[i][j] != ALIVE)
					vars->changes++;
			}
		}
	}
}

void print_grid(struct data * vars, char *msg){
int i,j;
	fflush(stdout);
	printf("%s grid size %d x %d\n", msg, vars->rowno, vars->colno);
	for(i = 0; i < vars->rowno; i++) {
		for(j = 0; j < vars->colno; j++){
			printf("%2d",vars->c_grid[i][j]);
		}
		printf("\n");
	}
	printf("\n");fflush(stdout);
}

/*
	send_receive()
		Copies sides, top, and bottom to their respective locations.
		All boundaries are considered periodic.

		In the MPI model, processes are aligned side-by-side.
		Left and right sides are sent to neighboring processes.
		Top and bottom are copied from the process's own grid.
*/
void send_receive (struct data * vars) {
	int i;//,j;

	int error;
	//int my_rank  = vars->my_rank;
	int mpi_size  = vars->mpi_size;
	int colno = vars->colno;
//	int rowno = vars->rowno;
	MPI_Comm my_comm = vars->my_comm;
	int ** grid = vars->c_grid;
	int *l_send_tmp = (int*)vars->l_send_tmp;
	int *l_rcv_tmp = (int*)vars->l_rcv_tmp;
	int *r_send_tmp = (int*)vars->r_send_tmp;
	int *r_rcv_tmp = (int*)vars->r_rcv_tmp;
//	MPI_Status status;
	MPI_Request requests[9];
	int reqs = 0;
	/*
		int xpos = my_rank /vars->colno;
		int ypos = my_rank % vars-> colno;
		printf("My rank = %d. my xpos = %d  my ypos = %d and mpisize = %d\n", vars->my_rank, xpos, ypos, mpi_size);
	*/

	enum TAGS {
		TOLEFT = 0,
		TORIGHT = 1,
		FROMLEFT = 2,
		FROMRIGHT = 3,
		TOTOP = 4,
		TOBOTTOM = 5,
		FROMTOP = 6,
		FROMBOTTOM = 7,
		TOTOPLEFT = 8,
		FROMTOPLEFT = 9,
		TOTOPRIGHT = 10,
		FROMTOPRIGHT = 11,
		TOBOTTOMLEFT =12,
		FROMBOTTOMLEFT = 13,
		TOBOTTOMRIGHT = 14,
		FROMBOTTOMRIGHT= 15
	};

	//	Some MPIs deadlock if a single process tries
	//to communicate with itself
	if (mpi_size != 1) {

		/* Copy your top to the top proccess*/
//		error = MPI_Send( &(grid[1][1]), vars->rowno, MPI_INT, vars->top, TOTOP,my_comm);
		error = MPI_Isend(&grid[1][1], vars->rowno, MPI_INT, vars->top, TOTOP, my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");
        /* Copy your bottom to the bottom proccess*/
		error = MPI_Isend(&(grid[vars->rowno][1]), vars->rowno, MPI_INT, vars->bot, TOBOTTOM, my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");
      /* Copy your top left to the top left proccess*/
		error = MPI_Isend( &(grid[1][1]), 1, MPI_INT, vars->t_left, TOTOPLEFT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");

      /* Copy your top rigth to the top right proccess*/
		error = MPI_Isend( &(grid[1][vars->colno]), 1, MPI_INT, vars->t_right, TOTOPRIGHT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");

      /* Copy your bottom left to the bottom left proccess*/
		error = MPI_Isend( &(grid[vars->rowno][1]), 1, MPI_INT, vars->b_left, TOBOTTOMLEFT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");

      /* Copy your bottom right to the bottom right proccess*/
		error = MPI_Isend( &(grid[vars->rowno][vars->colno]), 1, MPI_INT, vars->b_right, TOBOTTOMRIGHT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");

      /* Fill buffers for sending */
      for( i = 0; i< vars-> rowno; i++){
		l_send_tmp[i] = grid[i][1];
		r_send_tmp[i] = grid[i][vars->colno];
		//printf("%d == %d\n",l_send_tmp[i], grid[i][1]);
	  }
	  /* Copy your left to the left proccess*/
		error = MPI_Isend( l_send_tmp, colno, MPI_INT, vars->left, TOLEFT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");

      /* Copy your right to the rigth proccess*/
		error = MPI_Isend( r_send_tmp, colno, MPI_INT, vars->right, TORIGHT,
             my_comm, &(requests[0]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Isend\n");


        /* Receive top process' bottom as your top shadow line */
        error =  MPI_Irecv(&(grid[0][0]), vars->rowno, MPI_INT, vars->top, TOBOTTOM ,
							my_comm, &(requests[++reqs]));
        if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");

        /* Receive bottom process' top as your bottom shadow line */
        //error =  MPI_Recv(&(grid[vars->rowno + 1][0]), vars->rowno, MPI_INT, vars->bot, TOTOP ,
	//						my_comm, MPI_STATUS_IGNORE);
	error =  MPI_Irecv(&(grid[vars->rowno + 1][0]), vars->rowno, MPI_INT, vars->bot, TOTOP ,
                                                      my_comm, &(requests[++reqs]));

		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive top - left process' bottom - right as your top - left shadow cell */
        error =  MPI_Irecv(&(grid[0][0]), 1, MPI_INT, vars->t_left, TOBOTTOMRIGHT ,
							my_comm, &(requests[++reqs]));
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive top - right process' bottom - left as your top - right shadow cell */
        error =  MPI_Irecv(&(grid[0][vars->colno + 1]), 1, MPI_INT, vars->t_right, TOBOTTOMLEFT,
							my_comm, &(requests[++reqs]));
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive bottom - left process' top - right as your bottom - left shadow cell */
        error =  MPI_Irecv(&(grid[vars->rowno + 1][0]), 1, MPI_INT, vars->b_left, TOTOPRIGHT,
							my_comm, &(requests[++reqs]));
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive bottom - right process' top - left as your bottom - right shadow cell */
        error =  MPI_Irecv(&(grid[vars->rowno + 1][vars->colno + 1]), 1, MPI_INT, vars->b_right, TOTOPLEFT,
							my_comm, &(requests[++reqs]));
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive left process' right as your left shadow column */
        error =  MPI_Irecv(l_rcv_tmp, colno, MPI_INT, vars->left, TORIGHT ,
							my_comm, &(requests[++reqs]));
		if(error != 0)
			fprintf(stderr,"Error in MPI_Receive\n");
		/* Receive right process' left as your right shadow column */
        error =  MPI_Irecv(r_rcv_tmp, colno, MPI_INT, vars->right, TOLEFT ,
							my_comm, &(requests[++reqs]));
		for( i = 1; i< 9; i++){

			error = MPI_Wait(&(requests[i]), MPI_STATUS_IGNORE);
			if(error != 0)
				fprintf(stderr,"Error in MPI_Receive\n");
		}
		 /* Update table columns */
      for( i = 0; i< vars-> rowno; i++){
		 grid[i][1] = l_rcv_tmp[i];
		 grid[i][vars->colno] = r_rcv_tmp[i];
	  }

	}
		return;
}

/*
	copy_grid()
		Copies temporary values from f_grid into c_grid.
*/
void copy_grid (struct data * vars) {
	int i,j;
	int colno = vars->colno;
	int rowno = vars->rowno;
	int ** grid      = vars->c_grid;
	int ** f_grid = vars->f_grid;

	for (i = 0; i < rowno+2; i++)
		for (j = 0; j < colno+2; j++)
			grid[i][j] = f_grid[i][j];
}// END copy_grid()

/*
	allocate()
		Allocates memory for a 2D array of integers.
*/
void allocate (struct data * vars) {
	int i;//,j;
	int colno = vars->colno;
	int rowno = vars->rowno;

	vars->c_grid      = (int **) malloc(sizeof(int *) * (rowno+2));
	vars->f_grid = (int **) malloc(sizeof(int *) * (rowno+2));
	vars->l_send_tmp = (int **) malloc(sizeof(int) * rowno);
	vars->r_send_tmp = (int **) malloc(sizeof(int) * rowno);
	vars->l_rcv_tmp =  (int **) malloc(sizeof(int) * rowno);
	vars->r_rcv_tmp =  (int **) malloc(sizeof(int) * rowno);

	for (i = 0; i < colno+2; i++) {
		vars->c_grid[i]      = (int *) malloc(sizeof(int) * (colno+2));
		vars->f_grid[i] = (int *) malloc(sizeof(int) * (colno+2));
	}
}

/*
		Initialize cells based on input file, otherwise all cells
		are DEAD.
*/
void init_grids (struct data * vars) {
	FILE * fd = NULL;
	int i,j;

	if (vars->infile != NULL) {
		if ((fd = fopen(vars->infile, "r")) == NULL) {
			perror("Failed to open file for input");
			exit(EXIT_FAILURE);
		}

		if (fscanf(fd, "%d %d\n", &vars->colno, &vars->rowno) == EOF) {
			printf("File must at least define grid dimensions!\nExiting.\n");
			exit(EXIT_FAILURE);
		}
	}

	allocate(vars);

	for (i = 0; i < vars->colno+2; i++) {
		for (j = 0; j < vars->rowno+2; j++) {
			vars->c_grid[i][j]      = DEAD;
			vars->f_grid[i][j] = DEAD;
		}
	}

	if (vars->infile != NULL) {
		while (fscanf(fd, "%d %d\n", &i, &j) != EOF) {
			vars->c_grid[i][j]      = ALIVE;
			vars->f_grid[i][j] = ALIVE;
		}

		fclose(fd);
	} else {
		rand_grid(vars, INIT_PROB);
	}
}

/*
		Dumps the current state of life.grid to life.outfile.
		Only outputs the coordinates of !DEAD cells.
*/
void write_alive_cells(struct data * vars) {
	FILE * fd;
	int i,j;
	int colno   = vars->colno;
	int rowno   = vars->rowno;
	int ** grid = vars->c_grid;

	if (vars->outfile != NULL) {
		if ((fd = fopen(vars->outfile, "w")) == NULL) {
			perror("Failed to open file for output");
			exit(EXIT_FAILURE);
		}

		fprintf(fd, "%d %d\n", colno, rowno);

		for (i = 1; i <= rowno; i++) {
			for (j = 1; j <= colno; j++) {
				if (grid[i][j] != DEAD)
					fprintf(fd, "Process %d-> %d %d\n",vars->my_rank, i, j);
			}
		}

		fclose(fd);
	}
	if(vars->print == 1){
		char str[5];
		sprintf(str, "rank: %d", vars->my_rank);
		print_grid(vars,str);
	}
}

void free_grids (struct data * vars) {
	int i;
	int colno = vars->colno;

	for (i = 0; i < colno+2; i++) {
		free(vars->c_grid[i]);
		free(vars->f_grid[i]);
	}
	free(vars->l_send_tmp);
	free(vars->r_send_tmp);
	free(vars->l_rcv_tmp);
	free(vars->r_rcv_tmp);
	//printf("received:\n");
	free(vars->c_grid);
	free(vars->f_grid);

}

double rand_double() {
	return (double)random()/(double)RAND_MAX;
}

void rand_grid (struct data * vars, double prob) {
	int i,j;
	int colno = vars->colno;
	int rowno = vars->rowno;

	for (i = 1; i <= rowno; i++) {
		for (j = 1; j <= colno; j++) {
			if (rand_double() < prob)
				vars->c_grid[i][j] = ALIVE;
		}
	}
}

void seed_random (int my_rank) {
	srandom(time(NULL) + 100*my_rank);
}
void cleanup (struct data * vars) {
	write_alive_cells(vars);
	free_grids(vars);

	MPI_Finalize();
}

void help () {
	printf("\nUsage: Life [options]\n");
	printf("  -c   Number of columns in grid. Default: %d\n", DEFAULT_SIZE);
	printf("  -r   Number of rows in grid. Default: %d\n", DEFAULT_SIZE);
	printf("  -g   Number of generations to run. Default: %d\n", DEFAULT_GENS);
	printf("  -i   Input file. See README for format. Default: none.\n");
	printf("  -o   Output file. Default: none.\n");
	printf("  -p   Print result of each process to stdout. Default: none.\n");
	printf("  -h   This help page.\n");
	printf("\nSee README for more information.\n\n");

	exit(EXIT_FAILURE);
}

void parse_args (struct data * vars, int argc, char ** argv) {
	int opt       = 0;
	int opt_index = 0;
//	int i;

	for (;;) {
		opt = getopt_long(argc, argv, opts, long_opts, &opt_index);

		if (opt == -1) break;

		switch (opt) {
			case 'c':
				vars->colno = strtol(optarg, (char**) NULL, 10);
				break;
			case 'r':
				vars->rowno = strtol(optarg, (char**) NULL, 10);
				break;
			case 'g':
				vars->generations = strtol(optarg, (char**) NULL, 10);
				break;
			case 'i':
				vars->infile = optarg;
				break;
			case 'o':
				vars->outfile = optarg;
				break;
			case 'h':
			case '?':
				help();
				break;
			case 'p':
				vars->print = 1;
				break;
			default:
				break;
		}
	}

	// Backwards compatible argument parsing
	if (optind == 1) {
		if (argc > 1)
			vars->rowno       = strtol(argv[1], (char**) NULL, 10);
		if (argc > 2)
			vars->colno       = strtol(argv[2], (char**) NULL, 10);
		if (argc > 3)
			vars->generations = strtol(argv[3], (char**) NULL, 10);
	}
}
