#ifndef DC_ANALYSIS_H
#define DC_ANALYSIS_H
#include <gsl/gsl_linalg.h> //this requires gsl
//#include "csparse.h"
#include "cs.h"

#include "analysis_info.h"
cs_di * alpha_matrix_sparse;
cs_di * alpha_matrix_compressed;

cs_dis * csS;
cs_din * csN; 

double * alpha_matrix;
double * alpha_matrix_initial;
double * x_vector;
double * x2_vector; //this is for transient analysis vectorization dual core i7 socket north bridge
double * b_vector;
unsigned int vector_size;
double * z_vector;
double * r_vector;
double * p_vector;
double * q_vector;
double itol;
unsigned int k_index;

int alpha_matrix_size;
int init_analysis(unsigned int n, unsigned int m2);
void printf_matrix();
int mna_construction();
//void mna_solver();
gsl_permutation *  mna_solver();
//void dc_sweep_analysis(gsl_permutation * p);
//void dc_sweep_analysis(gsl_permutation * p, command_type * sol_info);
void dc_sweep_analysis(gsl_permutation * p, command_type * sol_info, int command_id);
void cg_iter_solver(void);
void bicg_iter_solver(void);

#endif
