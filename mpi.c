#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

void play_game(){
}
void get_surrounding_ranks(int *left, int *right, int *upper, int *lower,
						   int *top_left, int *top_right, int *bottom_left,
						   int *bottom_right,
						   MPI_Comm my_comm, int ndims, int *coords){
	
	/* Get the rank of the left and right elements of the process */
	MPI_Cart_shift(my_comm, 0, 1, left, right);

	/* Get the rank of the upper and lower elements of the process */
	MPI_Cart_shift(my_comm, 1, 1, upper, lower);

	/* Get the rank of the top-left, top-right, lower-left and lower-right elements of the process */
	MPI_Cart_coords(my_comm, *upper, ndims, coords/*must be 1-dimentional int array of 2 elements*/);
	coords[1] -= 1;
	MPI_Cart_rank( my_comm, coords, top_left );

	MPI_Cart_coords(my_comm, *upper, ndims, coords);
	coords[1] += 1;
	MPI_Cart_rank( my_comm, coords, top_right );

	MPI_Cart_coords(my_comm, *lower, ndims, coords);
	coords[1] -= 1;
	MPI_Cart_rank( my_comm, coords, bottom_left );

	MPI_Cart_coords(my_comm, *lower, ndims, coords);
	coords[1] += 1;
	MPI_Cart_rank( my_comm, coords, bottom_right );
}

int main(int argc, char *argv[])
{
	MPI_Comm my_comm;
	MPI_Request upper_request_send, upper_request_rcv,
				lower_request_send, lower_request_rcv,
				left_request_send, left_request_rcv,
				right_request_send, right_request_rcv,
				upper_left_request_send, upper_left_request_rcv,
				upper_right_request_send, upper_right_request_rcv,
				bottom_right_request_send, bottom_right_request_rcv,
				bottom_left_request_send, bottom_left_request_rcv;
	int com_size;
	int my_rank;
	int num_loops, grid_size;						/* Dimension size of the whole array */
	int ndims, reorder, periods[2], dim_size[2];
	int pr_grid_size;
	int *my_ext_grid;
	int *tmp_buff_s_l, *tmp_buff_r_l, *tmp_buff_s_r, *tmp_buff_r_r;
	int my_grid_size;
	int my_ext_grid_size;
	int i,j,generations;
	int error = 0;
	int left,right,upper,lower,top_left,top_right,bottom_left,bottom_right;
	int coords[2];

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &com_size);

	/* Check the validity of the given arguments */
	if (argc == 3)
	{
		if (atoi(argv[1]) <= 0)
		{
			if (my_rank == 0)
			{
				printf("1st argument (number of loops) must be greater than zero\n");
			}
			error = 1;
		}
		else
		{
			num_loops = atoi(argv[1]);
		}

		if (atoi(argv[2]) <= 0)
		{
			if (my_rank == 0)
			{
				printf("2nd argument (dimension of array) must be greater than zero\n");
			}
			error = 1;
		}
		else
		{
			grid_size = atoi(argv[2]);
		}
	}
	else
	{
		if (my_rank == 0)
		{
			printf("Must give exactly 2 arguments\n");
		}
		error = 1;
	}

	if (error)
	{
		MPI_Finalize();
		exit(1);
	}


	pr_grid_size = sqrt(com_size);				/* Dimension size of the Cartesian Topology */
	my_grid_size = grid_size / pr_grid_size;	/* Dimension size of the subArray of every process */

	ndims = 2;							/* 2D matrix */
	dim_size[0] = pr_grid_size;			/* N rows */
	dim_size[1] = pr_grid_size;			/* N columns */
	periods[0] = 1;						/* row periodic */
	periods[1] = 1;						/* column periodic */
	reorder = 1;						/* allows processes reordered for efficiency */

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &my_comm);

	/* Dimension size + the ghost elements of the subArray of every process */
	my_ext_grid_size = my_grid_size + 2;

	/* Use of 1-dimensional array to facilitate message sending */
	my_ext_grid = (int *)malloc(my_ext_grid_size * my_ext_grid_size * sizeof(int));

	tmp_buff_s_l = malloc(my_grid_size * sizeof(int));
	tmp_buff_r_l = malloc(my_grid_size * sizeof(int));
	tmp_buff_s_r = malloc(my_grid_size * sizeof(int));
	tmp_buff_r_r = malloc(my_grid_size * sizeof(int));
	
	/* initialize random seed: */
	srand (time(NULL));

	/* Randomize the position of organizations (non-ghost elements) */
	for(j=1;j<my_ext_grid_size-1;j++){
		for (i=j*my_ext_grid_size+1; i<j*(my_ext_grid_size-1); i++){
			my_ext_grid[i] = rand()%2;
		}
	}

get_surrounding_ranks(&left, &right, &upper, &lower, &top_left,
					  &top_right, &bottom_left, &bottom_right,
					  my_comm, ndims, coords);

	/* Begin receiving other processes elements (for fillin ghost-elements positions */
		 /*Protocol:
		  * ka8e process stelnei tags:
		  * aristera: 1
		  * panw-aristera: 2
		  * panw: 3
		  * panw-de3ia: 4
		  * de3ia: 5
		  * katw-de3ia:6
		  * katw: 7
		  *katw-aristera: 8
		  */
		for(generations = 0; generations < num_loops; generations++){

			/*Send _your_ upper line*/
			MPI_Isend( (void *)&my_ext_grid[my_ext_grid_size+1],my_grid_size, MPI_INT, upper, 3, my_comm, &upper_request_send);
			/* Receive top process' bottom line */
			MPI_Irecv( (void *)&my_ext_grid[1], my_grid_size, MPI_INT, upper, 7, my_comm, &upper_request_rcv);

			/*Send _your_ bottom line*/
			MPI_Isend( (void *)&my_ext_grid[(my_ext_grid_size-2)*my_ext_grid_size+1],my_grid_size, MPI_INT, lower, 7, my_comm, &lower_request_send);
			/* Receive bottom process' top line */
			MPI_Irecv( (void *)&my_ext_grid[(my_ext_grid_size-1)*my_ext_grid_size+1], my_grid_size, MPI_INT, upper, 3, my_comm, &lower_request_rcv);

			/* Send your first element to Top_left*/
			MPI_Isend( (void *)&my_ext_grid[my_ext_grid_size+1], 1, MPI_INT, top_left, 2, my_comm, &upper_left_request_send);
			/*Receive top_left process' bottom_right*/
			MPI_Irecv( (void *)&my_ext_grid[0], 1, MPI_INT, top_left, 6, my_comm, &upper_left_request_rcv);

			/* Send your top_right element to Top_right*/
			MPI_Isend( (void *)&my_ext_grid[2*my_ext_grid_size-2], 1, MPI_INT, top_right, 4, my_comm, &upper_right_request_send);
			/*Receive top_right process' bottom_left*/
			MPI_Irecv( (void *)&my_ext_grid[my_ext_grid_size-1], 1, MPI_INT, top_right, 8, my_comm, &upper_right_request_rcv);

			/* Send your bottom_right element to bottom_right*/
			MPI_Isend( (void *)&my_ext_grid[(my_ext_grid_size-1)*my_ext_grid_size-1], 1, MPI_INT, bottom_right, 6, my_comm, &bottom_right_request_send);
			/*Receive bottom_right process' top_left*/
			MPI_Irecv( (void *)&my_ext_grid[(my_ext_grid_size)*my_ext_grid_size-1], 1, MPI_INT, bottom_right, 2, my_comm, &bottom_right_request_rcv);

			/* Send your bottom_left element to bottom_left*/
			MPI_Isend( (void *)&my_ext_grid[(my_ext_grid_size-1)*my_ext_grid_size+1], 1, MPI_INT, bottom_left, 8, my_comm, &bottom_left_request_send);
			/*Receive bottom_left process' top_right*/
			MPI_Irecv( (void *)&my_ext_grid[(my_ext_grid_size-1)*my_ext_grid_size], 1, MPI_INT, bottom_left, 4, my_comm, &bottom_left_request_rcv);
			
			for(i=1;i<my_ext_grid_size;i++){
				tmp_buff_s_l[i] = my_ext_grid[i*my_ext_grid_size];
				tmp_buff_s_r[i] = my_ext_grid[(i+1)*my_ext_grid_size-1];
			}
			/*Send _your_ left column*/
			MPI_Isend( (void *)&tmp_buff_s_l, my_grid_size, MPI_INT, left, 1, my_comm, &left_request_send);
			/* Receive left process' right column */
			MPI_Irecv( (void *)&tmp_buff_r_l, my_grid_size, MPI_INT, left, 5, my_comm, &left_request_rcv);

			/*Send _your_ right column*/
			MPI_Isend( (void *)&tmp_buff_s_r, my_grid_size, MPI_INT, right, 5, my_comm, &right_request_send);
			/* Receive left process' right column */
			MPI_Irecv( (void *)&tmp_buff_r_r, my_grid_size, MPI_INT, right, 1, my_comm, &right_request_rcv);

			/*Twra pou exoume to grid*/
			play_game(&my_ext_grid, "inner");


			/* ena wait gia ka8e rcv kai meta*/
			for(i=1;i<my_ext_grid_size;i++){
				my_ext_grid[i*my_ext_grid_size] = tmp_buff_r_l[i];
				my_ext_grid[(i+1)*my_ext_grid_size-1] = tmp_buff_r_r[i];
			}
			play_game(&my_ext_grid, "outer");
		}
	free(my_ext_grid);
	MPI_Finalize();
return 0;
}

/* Na dw an borw se define na exw akyres metavlites*/
/*#define access_two2(width,xpos,ypos)
#define twodim[x][y] = onedim[x*size + y]


/* telika exoume 2dim array kai paizoume me deiktes gia na steiloume */
