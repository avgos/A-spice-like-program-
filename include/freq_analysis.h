#define _GNU_SOURCE 
#include <complex.h>
#include <math.h>
#include "element_structs.h"
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include <string.h>
#include <gsl/gsl_complex.h>
#include "dc_analysis.h"
#include <gsl/gsl_linalg.h> //this requires gsl
#include "cs.h"
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <errno.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>


typedef struct solver_info_{
	unsigned int sparse_size;
	cs_ci * alpha_sparse;
	cs_ci * alpha_sparse_cbar;
	//--------------
	int analysis_type;
	double complex * alpha_matrix;
	unsigned int vector_size;
	double complex * b_vector;
	//--------------
	gsl_permutation *p;
	gsl_vector *x;
	double complex * x_vector;
	gsl_matrix_view m;
} solver_info;


double complex * g_bar_matrix;
//int run_freq_analysis(void);
//int run_freq_analysis(int command_id);
int run_freq_analysis(int command_id, command_type * cmnd);
int g_bar_construction_im(solver_info * sol_info,double freq);
int e_vector_construction_im(double complex * e_vector, double freq);
double complex * complex_bicg_solver(solver_info * sol_info);
double complex * complex_cg_solver(solver_info * sol_info);
int init_freq_output(char *** filenames_, FILE *** fp_, int command_id, plot_print * plist);
void write_freq_output(double complex * x_vector, int is_log, char ** filenames, double freq, plot_print * plist);