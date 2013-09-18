#ifndef BCCD_LIFE_DEFAULTS_H
#define BCCD_LIFE_DEFAULTS_H

#include <stddef.h>
#include <stdbool.h>
#include <getopt.h>
#include <mpi.h>
#include <math.h>

static const char * opts = "c:r:g:i:o:t::xh?";
static const struct option long_opts[] = {
	{ "columns", required_argument, NULL, 'c' },
	{ "rows", required_argument, NULL, 'r' },
	{ "gens", required_argument, NULL, 'g' },
	{ "output", required_argument, NULL, 'o' },
	{ "input", required_argument, NULL, 'i' },
	{ "throttle", optional_argument, NULL, 't' },
	{ "help", no_argument, NULL, 'h' },
	{ NULL, no_argument, NULL, 0 }
};

// Default parameters for the simulation
const int DEFAULT_THROTTLE = 60;
const int     DEFAULT_SIZE = 4;//5;
const int     DEFAULT_GENS = 10;
const double     INIT_PROB = 0.25;

// All the data needed by an instance of Life
struct data {
	int  my_rank;
	int  mpi_size;
	int  mpi_grid_size;
	int  colno;
	int  rowno;
	int  ** c_grid;
	int  ** f_grid;
	int  ** l_send_tmp;
	int  ** l_rcv_tmp;
	int  ** r_send_tmp;
	int  ** r_rcv_tmp;
	int  generations;
	char * infile;
	char * outfile;
	MPI_Comm my_comm;
/* Ranks of neighboring processes*/
	int left;
	int right;
	int top;
	int bot;
	int t_left;
	int t_right;
	int b_left;
	int b_right;

	int changes;
};

enum CELL_STATES {
	DEAD = 0,
	ALIVE = 1
};

// Cells become DEAD with more than UPPER_THRESH 
// or fewer than LOWER_THRESH neighbors
const int UPPER_THRESH = 3;
const int LOWER_THRESH = 2;

// Cells with exactly SPAWN_THRESH neighbors become ALIVE
const int SPAWN_THRESH = 3;

const int MPI_DIMS = 2;
#endif
