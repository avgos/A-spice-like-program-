/*
 * Transient analysis declaration, info and options:
 *	.TRAN <time_step> <fin_time>
 *			<time_step>  --> The sampling step
 *			<fin_time>	 --> The finishing time of transient analysis(starting from time=0)
 *
 *	.OPTIONS METHOD=TR --> Approximate of derivative with trapezoidal method
 *	.OPTIONS METHOD=BE --> Approximate of derivative with Backward Euler(BE) method
 *
 *
 *
 *
 *
 *
 * Transient analysis types in Voltage/Current sources:
 *	EXP (i1 i2 td1 tc1 td2 tc2)
 *	SIN (i1 ia fr td df ph)
 *	PULSE (i1 i2 td tr tf pw per)
 *	PWL (t1 i1) (t2 i2) ... (tn in)
 * 	PWL (t1 v1) (t2 v2) ... (tn vn)
 *
 *
 *
 *
 *
 *
 * Trasient analysis system:
 *	
 *  _          _   dx(t)
 *  G * x(t) + C * ----  = e
 *                  dt
 *
 */



#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <stdlib.h>
#include "dc_analysis.h"
#include "element_structs.h"
#include "node_to_int.h"
#include <string.h>
#include "util.h"
#include <stdio.h>
#include <math.h>
#include "transient_analysis.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>



#define BICG_FPS 1e-14
#define M_PI 3.14159265358979323846

#define ALLOC_CHECK(a) if(a == NULL) { \
		fprintf(stderr, "[ERR] MALLOC FAILED FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__ );\
       	exit(-1);\
    }








typedef struct solver_info_{
	int analysis_type;
	cs * alpha_matrix_sparse;
	double *alpha_matrix;
	unsigned int alpha_matrix_size;
	double * b_vector;
	//--------------
	gsl_permutation *p;
	gsl_vector *x;
	double *x_vector;
	gsl_matrix_view m;
} solver_info;




double t_spec_calc(transient_spec * tmp, int step,command_type * cmnd);
//double * _cg_iter_solver(double * alpha_matrix, double *b_vector, unsigned int vector_size);
double * _cg_iter_solver(solver_info * sol_info);
//double * _bicg_iter_solver(double * alpha_matrix, double *b_vector, unsigned int vector_size) ;
double * _bicg_iter_solver(solver_info * sol_info);

/*
 *	
 *  _          _   dx(t)
 *  G * x(t) + C * ----  = e
 *                  dt     |
 *						   |
 *							--> e_vector is the RHS vector from the above system for the current step(t).
 */

int e_vector_construction(double * e_vector, int step,command_type * cmnd){
	current_source * temp2 = (current_source *)crnt_list;
	while(temp2 != NULL){
		//satiny check 2
		if(temp2->positive_node == temp2->negative_node){
			temp2 = temp2->next;
			continue;
		}else if(temp2->positive_node == 0){
			if(temp2->transient_info != NULL)
				e_vector[temp2->negative_node-1] += t_spec_calc(temp2->transient_info, step,cmnd);
			else
				e_vector[temp2->negative_node-1] += temp2->value;
		}else if(temp2->negative_node == 0){
			if(temp2->transient_info != NULL)
				e_vector[temp2->positive_node-1] -= t_spec_calc(temp2->transient_info, step,cmnd);
			else
				e_vector[temp2->positive_node-1] -= temp2->value;
		}else{
			if(temp2->transient_info != NULL){	
				e_vector[temp2->negative_node-1] += t_spec_calc(temp2->transient_info, step,cmnd);
				e_vector[temp2->positive_node-1] -= t_spec_calc(temp2->transient_info, step,cmnd);
			}else{
				e_vector[temp2->negative_node-1] += temp2->value;
				e_vector[temp2->positive_node-1] -= temp2->value;
			}
		}
		temp2 = temp2->next;
	}
	unsigned int e_vector_index = n; // n+1 - hope so
	voltage_source * temp3 = (voltage_source *)vltg_list;
	while(temp3 != NULL){
		//sanity check
		if(temp3->positive_node == temp3->negative_node){
			temp3 = temp3->next;
			continue;
		}else if(temp3->negative_node == 0){
			if(temp3->transient_info != NULL)
				e_vector[e_vector_index] += t_spec_calc(temp3->transient_info, step,cmnd);
			else
				e_vector[e_vector_index] += temp3->value;
		}else if(temp3->positive_node == 0){
			if(temp3->transient_info != NULL)
				e_vector[e_vector_index] += t_spec_calc(temp3->transient_info, step,cmnd);
			else
				e_vector[e_vector_index] += temp3->value;
		}else{
			if(temp3->transient_info != NULL)
				e_vector[e_vector_index] += t_spec_calc(temp3->transient_info, step,cmnd);
			else
				e_vector[e_vector_index] += temp3->value;
		}
		e_vector_index++;
		temp3 = temp3->next;
	}
	return 1;
}




double t_spec_calc(transient_spec * tmp, int step,command_type * cmnd){
	//tmp->spec.sin_info.i1 = atof(t1);s
	if(tmp == NULL) return 0.0;
	double tmp_calc = 0.0;
	double time_point = cmnd->d.t_info.tran_time_step*(double)step;
	if(tmp->type == 'e'){
		//exp
		if(time_point <= tmp->spec.exp_info.td1)
			tmp_calc = tmp->spec.exp_info.i1;
		else if(tmp->spec.exp_info.td1 < time_point && time_point <= tmp->spec.exp_info.td2){
			tmp_calc = tmp->spec.exp_info.i1 + 
			(tmp->spec.exp_info.i2 - tmp->spec.exp_info.i1)*(1 - exp((tmp->spec.exp_info.td1-time_point)/(tmp->spec.exp_info.tc1)) );
		}else{
			tmp_calc = tmp->spec.exp_info.i1 + 
			(tmp->spec.exp_info.i2 - tmp->spec.exp_info.i1)*(
				exp((tmp->spec.exp_info.td2-time_point)/(tmp->spec.exp_info.tc2)) -
				 exp((tmp->spec.exp_info.td1-time_point)/(tmp->spec.exp_info.tc1)) );
		}		
	}else if(tmp->type == 's'){
		//sin
		//tran_spec_info->spec.sin_info.i1 = atof(t1);
		if(time_point <= tmp->spec.sin_info.td){
			tmp_calc = tmp->spec.sin_info.i1 + tmp->spec.sin_info.ia*sin(2*M_PI*tmp->spec.sin_info.ph/360);
		}else{
			tmp_calc = tmp->spec.sin_info.i1 + tmp->spec.sin_info.ia*sin(
				2*M_PI*tmp->spec.sin_info.fr*(time_point - tmp->spec.sin_info.td)
			  + 2*M_PI*tmp->spec.sin_info.ph/360)*exp((tmp->spec.sin_info.td - time_point)*tmp->spec.sin_info.df);
		}
	}else if(tmp->type == 'p'){
		//pulse

		double where = fmod(cmnd->d.t_info.tran_time_step*step, tmp->spec.pulse_info.per);
		//tran_spec_info->spec.pulse_info.i1 = atof(t1);
		if(where <= tmp->spec.pulse_info.td){
			tmp_calc = tmp->spec.pulse_info.i1;
		}else if(where >= tmp->spec.pulse_info.td && where <= tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr){
			//lin
			tmp_calc = ((tmp->spec.pulse_info.i2 - tmp->spec.pulse_info.i1)/(tmp->spec.pulse_info.td + 
				tmp->spec.pulse_info.tr - tmp->spec.pulse_info.td))*(where - tmp->spec.pulse_info.td) + 
			tmp->spec.pulse_info.i1;

		}else if(where >= tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr && where <= tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr + tmp->spec.pulse_info.pw){
			tmp_calc = tmp->spec.pulse_info.i2;
		}else if(where >= tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr + tmp->spec.pulse_info.pw && where <= 
			tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr + tmp->spec.pulse_info.pw + tmp->spec.pulse_info.tf){
			//lin
			tmp_calc = ((tmp->spec.pulse_info.i1 - tmp->spec.pulse_info.i2)/
				((tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr + tmp->spec.pulse_info.pw + tmp->spec.pulse_info.tf)
					-(tmp->spec.pulse_info.td + tmp->spec.pulse_info.tr + tmp->spec.pulse_info.pw)
					))*(where - tmp->spec.pulse_info.td - tmp->spec.pulse_info.tr - tmp->spec.pulse_info.pw) 
			+ tmp->spec.pulse_info.i2;
		}else{
			tmp_calc = tmp->spec.pulse_info.i1;
		}


	}else if(tmp->type == 'w'){
		//pwl
		double y1,y2;
		double t1,t2;
		int step_idx  = 0;
		while(tmp->spec.pwl_info.t1[step_idx+1] <= time_point ) step_idx++;
		y1 = tmp->spec.pwl_info.i1[step_idx];
		y2 = tmp->spec.pwl_info.i1[step_idx+1];
		t1 = tmp->spec.pwl_info.t1[step_idx];
		t2 = tmp->spec.pwl_info.t1[step_idx+1];
		double m  = (y2 - y1) / (t2 - t1);
		tmp_calc = m*(time_point - t1) + y1;
	}

	return tmp_calc; 
}




/*
 * get_analysis_type() == 0 --> LU solver
 * get_analysis_type() == 1 --> Cholesky solver
 * get_analysis_type() == 2 --> Iterative solver using Conjugate Gradient(CG)
 * get_analysis_type() == 3 --> Iterative solver using Bi-Conjugate Gradient(Bi-CG)
 * 
 * get_sparse() == 0 --> Small Scale System Solve (no use of sparse matrices)
 * get_sparse() == 1 --> Large Scale System Solve (use of sparse matrices)
 *
 */

void init_solver(solver_info * sol_info ){
	gsl_matrix_view m;
	sol_info->x = gsl_vector_alloc (sol_info->alpha_matrix_size);
	sol_info->p = gsl_permutation_alloc (sol_info->alpha_matrix_size);
	if(get_sparse() == 0){
		if(get_analysis_type() == 0){
			int s;
			m = gsl_matrix_view_array(sol_info->alpha_matrix, sol_info->alpha_matrix_size, sol_info->alpha_matrix_size);
			gsl_linalg_LU_decomp (&m.matrix, sol_info->p, &s);
		}else if(get_analysis_type() == 1){
			m = gsl_matrix_view_array(sol_info->alpha_matrix, sol_info->alpha_matrix_size, sol_info->alpha_matrix_size);
			gsl_linalg_cholesky_decomp(&m.matrix);
		}
	}else{
		if(get_analysis_type() == 0){
			//LU
			/*csS = cs_sqr(2,alpha_matrix_compressed,0);
			csN = cs_lu(alpha_matrix_compressed, csS,1);			
			b_initial = (double *)malloc(sizeof(double)*vector_size);
			memcpy(b_initial,b_vector, sizeof(double)*vector_size);
			*/

			cbar_csS = cs_sqr(2,sol_info->alpha_matrix_sparse,0);
			sol_info->alpha_matrix_sparse->nzmax = sol_info->alpha_matrix_sparse->n;
			cbar_csN = cs_lu(sol_info->alpha_matrix_sparse, cbar_csS,1);
			//cs_spfree(alpha_matrix_compressed);
			//end of alloc
			//x_size = vector_size;
		}else if(get_analysis_type() == 1){
			//cholesky
			cbar_csS = cs_schol(1, sol_info->alpha_matrix_sparse);
			cbar_csN = cs_chol(sol_info->alpha_matrix_sparse, cbar_csS);
			cs_spfree(sol_info->alpha_matrix_sparse);
			//b_initial = (double *)malloc(sizeof(double)*vector_size);
			//memcpy(b_initial,b_vector, sizeof(double)*vector_size);
		}
	}
}


double * solver_step( solver_info * sol_info){
	//printf("get_analysis_type: %d get_sparse: %d\n", get_analysis_type(), get_sparse());
	//exit(-1);

	gsl_matrix_view m;
	if(get_analysis_type() == 0){
		if(get_sparse() == 0){
		m = gsl_matrix_view_array(sol_info->alpha_matrix, sol_info->alpha_matrix_size, sol_info->alpha_matrix_size);
		gsl_vector_view b = gsl_vector_view_array(sol_info->b_vector, sol_info->alpha_matrix_size);
		gsl_linalg_LU_solve (&(m).matrix, sol_info->p, &b.vector, sol_info->x);
		sol_info->x_vector = sol_info->x->data;
	}else{
		cs_ipvec(cbar_csN->pinv,sol_info->b_vector,sol_info->x_vector,sol_info->alpha_matrix_size);
		cs_lsolve(cbar_csN->L,sol_info->x_vector);
		cs_usolve(cbar_csN->U,sol_info->x_vector);
		cs_ipvec(cbar_csS->q,sol_info->x_vector,sol_info->b_vector,sol_info->alpha_matrix_size);
		//sol_info->x_vector = sol_info->b_vector;
		memcpy(sol_info->x_vector, sol_info->b_vector, sizeof(double)*vector_size);
/*

		memcpy(b_vector,b_initial, sizeof(double)*vector_size);
		cs_ipvec(csN->pinv,b_vector,x_vector,vector_size);
		cs_lsolve(csN->L,x_vector);
		cs_usolve(csN->U,x_vector);
		cs_ipvec(csS->q,x_vector,b_vector,vector_size);
		x_data = b_vector;*/

	}
	//printf("Hello\n");
	}else if(get_analysis_type() == 1){
		if(get_sparse() == 0){
			m = gsl_matrix_view_array(sol_info->alpha_matrix, sol_info->alpha_matrix_size, sol_info->alpha_matrix_size);
			gsl_vector_view b = gsl_vector_view_array(sol_info->b_vector, sol_info->alpha_matrix_size);
			gsl_linalg_cholesky_solve(&(m).matrix, &b.vector, sol_info->x);
			sol_info->x_vector = sol_info->x->data;
		}else{
			// sparsing..
//				memcpy(b_vector,b_initial, sizeof(double)*vector_size);
			cs_ipvec(cbar_csS->pinv,sol_info->b_vector,sol_info->x_vector,sol_info->alpha_matrix_size);
			cs_lsolve(cbar_csN->L,sol_info->x_vector);
			cs_ltsolve(cbar_csN->L,sol_info->x_vector);
			cs_pvec(cbar_csS->pinv,sol_info->x_vector,sol_info->b_vector,sol_info->alpha_matrix_size);
			//x_data = b_vector;
			memcpy(sol_info->x_vector, sol_info->b_vector, sizeof(double)*(sol_info->alpha_matrix_size));
		}
	}else if(get_analysis_type() == 2){	
		//if(get_sparse() != 0) exit(-1);
		//Conjugate Gradient - SPD Matrices
		//cg_iter_solver();
		sol_info->x_vector = _cg_iter_solver(sol_info);
	}else if(get_analysis_type() == 3){
		//if(get_sparse() != 0) exit(-1);
	//Bi-Conjugate Gradient - Genaral Matrices
		//bicg_iter_solver();
		//sol_info->x_vector = _bicg_iter_solver(sol_info->alpha_matrix, sol_info->b_vector, sol_info->alpha_matrix_size);
		sol_info->x_vector = _bicg_iter_solver(sol_info);
		
	}
	

	return NULL;
}



/* print a sparse matrix */
int mycs_print (const cs *A, int brief)
{
	if(!get_sparse()) return -1;
    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (nz < 0)
    {
	printf ("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax,
		Ap [n], cs_norm (A)) ;
	for (j = 0 ; j < n ; j++)
	{
	    printf ("    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		printf ("      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
		if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
	    }
	}
    }
    else
    {
	printf ("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
	for (p = 0 ; p < nz ; p++)
	{
	    printf ("    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
	    if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
	}
    }
    return (1) ;
}




/*
 *                                        dx(t)
 * get_approx_type() == 1 --> Approximate ----   with trapezoidal(TR) method
 *										   dt
 *
 *                                        dx(t)
 * get_approx_type() == 2 --> Approximate ----   with Backward Euler(BE) method
 *										   dt
 *
 */

void trnsient_analysis(int command_id, command_type *cmnd){
/*
	double * matrix = malloc(sizeof(double)*16);
	double * vector = malloc(sizeof(double)*4);
	double * sol;
	matrix[0] = 9;
	matrix[1] = -3;
	matrix[2] = 3;
	matrix[3] = 9;
	//---------------
	matrix[4] = -3;
	matrix[5] = 17;
	matrix[6] = -1;
	matrix[7] = -7;
	//-------------
	matrix[8] = 3;
	matrix[9] = -1;
	matrix[10] = 17;
	matrix[11] = 15;
	//---------------
	matrix[12] = 9;
	matrix[13] = -7;
	matrix[14] = 15;
	matrix[15] = 44;


	//--------------
	vector[0] = 1;
	vector[1] = 2;
	vector[2] = 3;
	vector[3] = 4;
 	sol = _cg_iter_solver(matrix, vector, 4);
	for(int i=0;i<4;i++){
		printf("%20.16f\n",sol[i]);
	}
	exit(-1);*/


	solver_info * sol_info = NULL;
	sol_info = malloc(sizeof(solver_info));
	unsigned int step_num;
	step_num = (int)ceil(cmnd->d.t_info.tran_fin_time/cmnd->d.t_info.tran_time_step);


	plot_print * pptemp;
	pptemp = cmnd->plot_list;
	int plot_counter = 0;
	while(pptemp != NULL){
		if(pptemp->index != -1)
			plot_counter++;
		pptemp = pptemp->next;
	}
	if(plot_counter == 0){
    	fprintf(stderr, "ERR: Nothing to print, printing all the results in stdout\n" );
    	exit(-1);
    }

	char ** filenames;
	FILE **fp = NULL;
	fp = malloc(sizeof(FILE *)*(plot_counter+1));
	filenames = malloc(sizeof(char *)*(plot_counter+1));
	int i = 0;


	//time_t timestamp = time(NULL);
	pptemp = cmnd->plot_list;//plot_print_list;
	while(pptemp != NULL){
		if(pptemp->index == -1){
			pptemp = pptemp->next;
			continue;
		}
		if(!fopen("sol","r")){
			mkdir("sol", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
		if(get_sparse()){
			asprintf(&filenames[i],"%s/%d_TRANSIENT_sparse_node_%s",cs_info->sol_dir,command_id, &pptemp->name[2]);
		}else{
			asprintf(&filenames[i],"%s/%d_TRANSIENT_node_%s",cs_info->sol_dir, command_id, &pptemp->name[2]);
		}

		fp[i] = fopen(filenames[i], "w+");
		if(!fp[i]){
			fprintf(stderr, "failed to open: %s\n", filenames[i]);
		}else{
			fprintf(fp[i], "#*********************************************************************\n" );
			fprintf(fp[i], "#Time step:%lf to time: %lf steps: %d node name: %s\n", cmnd->d.t_info.tran_time_step, 
			cmnd->d.t_info.tran_fin_time, step_num, pptemp->name);
			fprintf(fp[i], "#*********************************************************************\n" );	
			fclose(fp[i]);
		}
		pptemp = pptemp->next;
		i++;
	}











	if(get_approx_type() == 1){
		//trapezoidal
		printf("trapezoidal\n");
		printf("-------------------------------------------------------------\n");
		int step;
		int i,j;
		double * alpha_temp;
		alpha_temp = malloc(sizeof(double)*alpha_matrix_size*alpha_matrix_size);
		double * alpha_temp_b;
		double * x_vector_prev = NULL;
		double * e_vector_prev, *e_vector;
		x_vector_prev = malloc(sizeof(double)*vector_size);
		e_vector_prev = malloc(sizeof(double)*vector_size);
		e_vector = malloc(sizeof(double)*vector_size);
		
		for(i=0;i<vector_size;i++){
			e_vector[i] = e_vector_prev[i] = 0.0;
			x_vector_prev[i] = x2_vector[i];
		}

		e_vector_construction(e_vector_prev, 0,cmnd);
		if(plot_counter == 0){
			printf("E_VECTOR START-------------------------------------------\n");
			//for(int i=0;i<vector_size;i++) printf("%8.5f\t%8.5f\t%8.5f\n", 0.0,0.0,e_vector_prev[i]);
			for(int i=0;i<vector_size;i++) printf("%8.5f\t%8.5f\t%8.5f\n", 0.0,0.0,e_vector_prev[i]);
			printf("E_VECTOR END-------------------------------------------\n");
		}
			
		alpha_temp_b = malloc(sizeof(double)*alpha_matrix_size*alpha_matrix_size); 
		// A + (2*c_bar)/h
		cs * alpha_temp_sparse_;
		cs * alpha_temp_sparse;
		if(get_sparse() == 0){
			for(i=0;i<alpha_matrix_size*alpha_matrix_size; i++){
				alpha_temp[i]   = alpha_matrix_initial[i] + (2*c_bar_matrix[i]) / cmnd->d.t_info.tran_time_step;
				alpha_temp_b[i] = alpha_matrix_initial[i] - (2*c_bar_matrix[i]) / cmnd->d.t_info.tran_time_step;
			}
			/*if(plot_counter == 0){
				printf("\ncbar\n");
				for(i=0;i<alpha_matrix_size*alpha_matrix_size;i++)
					printf("%8.5f", c_bar_matrix[i]);
				printf("\n");
			}*/
		}else{
			alpha_temp_sparse = cs_add(alpha_matrix_compressed, cbar_compressed, 1.0, 2.0/cmnd->d.t_info.tran_time_step);
			alpha_temp_sparse_ = cs_add(alpha_matrix_compressed, cbar_compressed, 1.0, -2.0/cmnd->d.t_info.tran_time_step);
			sol_info->alpha_matrix_sparse = alpha_temp_sparse;
		}




		sol_info->alpha_matrix_size = alpha_matrix_size;
		sol_info->alpha_matrix = alpha_temp;
		init_solver(sol_info);
		sol_info->x_vector = malloc(sizeof(double)*vector_size);
		for(step=1;step<step_num;step++){

			
			e_vector_construction(e_vector,step,cmnd);
			double temp_mult;
			double * y_temp = NULL;
			if(get_sparse() == 0){
				for(i=0;i<alpha_matrix_size;i++){
					temp_mult = 0.0;
					for(j=0;j<alpha_matrix_size;j++){
						temp_mult += alpha_temp_b[i*alpha_matrix_size+j]* x_vector_prev[j];
					}
					double t = e_vector_prev[i];
					e_vector_prev[i] = e_vector[i] + t - temp_mult;
				}
			}else{
				y_temp = malloc(sizeof(double)*vector_size);
				for(int i=0;i<vector_size;i++) y_temp[i] = 0.0;

				cs_gaxpy(alpha_temp_sparse_, x_vector_prev, y_temp);
				for(i=0;i<alpha_matrix_size;i++){
					double t = e_vector_prev[i];
					e_vector_prev[i] = e_vector[i] + t - y_temp[i];
				}
				free(y_temp);
			}
			sol_info->b_vector = e_vector_prev; //movable
			//mycs_print(alpha_temp_sparse_,0);
			if(plot_counter == 0){				
				printf("E_VECTOR START-------------------------------------------\n");
				printf("b_vector\te_vector\te_vector_prev\tyt\tx_vector_prev\n");
				for(int i=0;i<vector_size;i++){


					double yt = 0.0;
					if(get_sparse() == 1){
						yt = y_temp[i];
					}
				printf("%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n", 
					sol_info->b_vector[i],
					e_vector[i],
					e_vector_prev[i],
					yt,
					x_vector_prev[i]);
				}	
				printf("E_VECTOR END-------------------------------------------\n");
			}
			solver_step(sol_info);
			//exit(-1);
			//if(get_sparse() == 1)
			//	;//mycs_print(sol_info->alpha_matrix_sparse,0);
			/*if(get_sparse() != 1){
				for(int i=0;i<vector_size;i++){
					for(j=0;j<vector_size;j++){
						printf("%8.5f ", alpha_matrix_initial[i*vector_size+j]);
					}
					printf("\n");
				}
			}*/

			i = 0;
			pptemp = cmnd->plot_list;//plot_print_list;
			while(pptemp != NULL){
				if(pptemp->index == -1){
					pptemp = pptemp->next;
					continue;
				}
				fp[i] = fopen(filenames[i], "a+");
				if(fp[i]){
				fprintf(fp[i], "%.20e %.20e\n",((double)step)*cmnd->d.t_info.tran_time_step,sol_info->x_vector[pptemp->index-1] );
				fflush(fp[i]);
				fclose(fp[i]);
			}else{
				printf("failed to write to file %s\n",filenames[i] );
			}
				pptemp = pptemp->next;
				i++;
			}
			//afou exw lusei giati exw lusei
			if(plot_counter == 0){
				for(i=0;i<alpha_matrix_size;i++){
					//printf("x_data: %f\n",skata[i]); //sol_info->x_vector[i]);
					printf("x_data: %f\n",sol_info->x_vector[i]);
					//x_vector_prev[i] = sol_info->x_vector[i];// x_vector[i];
					//e_vector_prev[i] = e_vector[i];
					//e_vector[i] = 0.0;
				}
				printf("\n");
			}
			for(i=0;i<alpha_matrix_size;i++){
					//printf("x_data: %f\n",skata[i]); //sol_info->x_vector[i]);
					//printf("x_data: %f\n",sol_info->x_vector[i]);
					x_vector_prev[i] = sol_info->x_vector[i];// x_vector[i];
					e_vector_prev[i] = e_vector[i];
					e_vector[i] = 0.0;
			  
			}
			

		}
		free(x_vector_prev);
		printf("--------------------------------------------------------------\n");
		printf("trapezoidal\n");
		printf("-------------------------------------------------------------\n");
	}else if(get_approx_type() == 2){
		//backward euler
		int step;
		int i,j;
		double * alpha_temp;
		alpha_temp = malloc(sizeof(double)*alpha_matrix_size*alpha_matrix_size);
		double * x_vector_prev = NULL;
		double * e_vector;
		if(vector_size == 0){
			fprintf(stderr, "[ERROR] vector_size too small\n");
			exit(-1);
		}
		x_vector_prev = malloc(sizeof(double)*vector_size);
		for(i=0;i<vector_size;i++) x_vector_prev[i] = x2_vector[i];
		e_vector = malloc(sizeof(double)*vector_size);
		for(i=0;i<vector_size;i++) e_vector[i] = 0.0;


		//cs * alpha_temp_sparse_;
		cs * alpha_temp_sparse;
		if(get_sparse() == 0){
			// A + (1*c_bar)/h
			for(i=0;i<alpha_matrix_size*alpha_matrix_size;i++){
				alpha_temp[i] = alpha_matrix_initial[i] + ((1*c_bar_matrix[i]) / cmnd->d.t_info.tran_time_step);
			}
		}else{

			//printf("alpha_matrix_compressed ------------------------\n");
			//mycs_print(alpha_matrix_compressed,0);
			//printf("alpha_temp_sparse ------------------------\n");
			alpha_temp_sparse = cs_add(alpha_matrix_compressed, cbar_compressed, 1.0, 1.0/cmnd->d.t_info.tran_time_step);
			//alpha_temp_sparse_ = cs_compress(alpha_temp_sparseb);
			//cs_dupl(alpha_temp_sparse_);
						//cs_dupl()
			sol_info->alpha_matrix_sparse = alpha_temp_sparse;

			//cs * alpha_temp_sparseb = cs_add(alpha_matrix_compressed, cbar_compressed, 1.0, -2.0/tran_time_step);
			//sol_info->alpha_matrix_sparse = alpha_temp_sparse;
			//sol_info->x_vector
		}





		sol_info->alpha_matrix_size = alpha_matrix_size;
		sol_info->alpha_matrix = alpha_temp;
		sol_info->b_vector = e_vector;

		init_solver(sol_info);
		sol_info->x_vector = malloc(sizeof(double)*vector_size);

		for(step=1;step<step_num;step++){
			/*e_vector_construction(e_vector, step);

			double temp_mult;
			for(i=0;i<alpha_matrix_size;i++){
				temp_mult = 0.0;
				for(j=0;j<alpha_matrix_size;j++){
					temp_mult += c_bar_matrix[i*alpha_matrix_size+j]* x_vector_prev[j];
				}
				e_vector[i] += temp_mult/tran_time_step;
			}*/
			for(i=0;i<vector_size;i++) e_vector[i] = 0.0;
			e_vector_construction(e_vector,step,cmnd);
			double temp_mult;
			double * y_temp;
			if(get_sparse() == 0){
				for(i=0;i<alpha_matrix_size;i++){
					temp_mult = 0.0;
					for(j=0;j<alpha_matrix_size;j++){
						temp_mult += c_bar_matrix[i*alpha_matrix_size+j]* x_vector_prev[j];
					}
					e_vector[i] += temp_mult/cmnd->d.t_info.tran_time_step;
				}
			}else{
				y_temp = malloc(sizeof(double)*vector_size);
				for(int i=0;i<vector_size;i++) y_temp[i] = 0.0;

				cs_gaxpy(cbar_compressed, x_vector_prev, y_temp);
				for(i=0;i<alpha_matrix_size;i++){
					e_vector[i] += y_temp[i]/cmnd->d.t_info.tran_time_step;
				}
				free(y_temp);
			}
			sol_info->b_vector = e_vector; //movable
			if(plot_counter == 0){				
				printf("E_VECTOR START-------------------------------------------\n");
				printf("b_vector\n");
				for(int i=0;i<vector_size;i++){
					printf("%8.5f\n",e_vector[i] );
				}
				printf("E_VECTOR END-------------------------------------------\n");
			}


			solver_step(sol_info);

			i = 0;
			pptemp = cmnd->plot_list;//plot_print_list;
			while(pptemp != NULL){
				if(pptemp->index == -1){
					pptemp = pptemp->next;
					continue;
				}
				fp[i] = fopen(filenames[i], "a+");
				if(fp[i]){

				fprintf(fp[i], "%.20e %.20e\n",step*cmnd->d.t_info.tran_time_step, sol_info->x_vector[pptemp->index-1] );
				fflush(fp[i]);
				fclose(fp[i]);
			}else{
				printf("failed to write to file %s\n",filenames[i] );
			
			}
				pptemp = pptemp->next;
				i++;
			}

			if(plot_counter == 0){			
				//afou exw lusei giati exw lusei
				for(i=0;i<alpha_matrix_size;i++){
					printf("x_data: %f\n",sol_info->x_vector[i]);
					//x_vector_prev_bicg[i] = x_res[i];
					//x_vector_prev[i] = sol_info->x_vector[i];// x_vector[i];
					//e_vector[i] = 0.0;
				}
				printf("\n");
			}
			for(i=0;i<alpha_matrix_size;i++){
					x_vector_prev[i] = sol_info->x_vector[i];// x_vector[i];
					e_vector[i] = 0.0;
			  
			}
		}
		free(x_vector_prev);
	}
	if(plot_counter > 0){
		free(fp);
		//for(int i=0;i<plot_counter;i++){
		//	fclose(fp[i]);
		//}
		FILE * gp_pipe = NULL;

		char * command = "gnuplot --persist";	
		char * mode = "w";
	    int fd[2];
	    int pid;
	    pipe(fd);
	    pid = vfork(); //DON'T COPY THE WHOLE FUCKING PAGE TABLE
	    if(pid == 0){
	    	//in the child
	        dup2(fd[0], STDIN_FILENO);
    		close(fd[0]);   /* close other gp_pipe fds */
    		close(fd[1]);
    		execl("/bin/sh", "sh", "-c", command, (char *) NULL);
    		_exit(1);
	    }else if(pid != -1){
            close(fd[0]);
	        if (!(gp_pipe = fdopen(fd[1], mode))) {
    	        close(fd[1]);
        	}
    	}
	    
		if(gp_pipe == NULL){
			printf( "The error generated was %d\n", errno );
	 	   	printf( "That means: %s\n", strerror( errno ) );
			fprintf(stderr, "Install gnuplot for graphs..\n");
			//exit(-1);
		}
		printf("------------------------------------------\n");
		printf("RESULT FILES -----------------------------\n");
		printf("------------------------------------------\n");
		for(int i=0;i<plot_counter;i++){
			printf("%s\n", filenames[i]);
		}

		printf("gnuplot command start: --------------------------------------\n");
		fprintf(stdout, "plot ");

		if(gp_pipe){
			printf("|%s|\n",cmnd->raw_cmnd );
			fprintf(gp_pipe, "set title \"%s\" font \",14\" noenhanced\n",cmnd->raw_cmnd );
			fprintf(gp_pipe, "plot ");
		}
		for(int i=0;i<plot_counter;i++){
			char * node_name;
			node_name = strstr(filenames[i],"node");
			if(node_name == NULL) node_name = filenames[i];
			else node_name += 5;
			
			if(gp_pipe){
				fprintf(gp_pipe, "\"%s\" using 1:2 title \"node: %s\" noenhanced smooth unique %s", filenames[i], node_name,
				(i != plot_counter-1)?", ": " \n");
			}
			fprintf(stdout, "\"%s\" using 1:2 title \"node: %s\" smooth unique %s", filenames[i], node_name,
				(i != plot_counter-1)?", ": " \n");

		}
		
		printf("gnuplot command end: --------------------------------------\n");
	    if(gp_pipe){	
		    fflush(gp_pipe); 
			pclose(gp_pipe);
	    }
		
		printf("------------------------------------------\n");
	}
	free(sol_info);
}





//double * _cg_iter_solver(double * alpha_matrix, double *b_vector, unsigned int vector_size){

/*
 * CG iterative solver 
 *
 */
double * _cg_iter_solver(solver_info * sol_info){
	
	double * alpha_matrix = sol_info->alpha_matrix;
	double *b_vector = sol_info->b_vector;
	unsigned int vector_size = sol_info->alpha_matrix_size;	

	int i,j;
	printf("------------------------------------------------------\n");
	printf("cg_iter_solver\n");
	printf("------------------------------------------------------\n");
	printf("alpha_matrix\n");
	printf("------------------------------------------------------\n");
	for(i=0;i<vector_size;i++){
		for(j=0;j<vector_size;j++){
			printf("%8.5f ", alpha_matrix[i*vector_size+j]);
		}
		printf("\n");
	}

	printf("------------------------------------------------------\n");
	printf("\n");
	printf("b_vector\n");
	for(i=0;i<vector_size;i++){
		printf("%7.5f\n", b_vector[i]);
	}
	printf("------------------------------------------------------\n");
	
	if(vector_size == 0){
		fprintf(stderr, "[ERROR] vector_size too small\n");
		exit(-1);
	}


	double * r_vector;
	double * z_vector;
	double * q_vector;
	double * p_vector;
	double * x_vector;

	r_vector = malloc(sizeof(double)*vector_size);
	z_vector = malloc(sizeof(double)*vector_size);
	q_vector = malloc(sizeof(double)*vector_size);
	p_vector = malloc(sizeof(double)*vector_size);
	x_vector = malloc(sizeof(double)*vector_size);
	
	for(i=0;i<vector_size;i++){
		r_vector[i] = b_vector[i];
		z_vector[i] = 0.0;
		q_vector[i] = 0.0;
		p_vector[i] = 0.0;
	}
	unsigned int iter = 0;
	double * M;
	//double * M2;
	double norm, norm_b = 0.0, norm_r;
	

	M = malloc(sizeof(double)*(vector_size));
	memset(M,0, sizeof(double)*vector_size);
	//M2 = malloc(sizeof(double)*(vector_size));
	if(get_sparse() == 0){
		for(i=0;i<(vector_size);i++){
			x_vector[i] = 0.0;
			M[i] = (alpha_matrix[(i*(vector_size))+i] == 0.0) ? 1.0 : alpha_matrix[(i*(vector_size))+i];
			norm_b += pow(b_vector[i],2);
		}
	}else{
    	//sparse
    	norm_b = 0.0;
 		int column = 0;
 		M[column] = 1.0;
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = sol_info->alpha_matrix_sparse->p[column+1] - sol_info->alpha_matrix_sparse->p[column];
 		for(i=0;i<non_zero;i++){
 			if(sol_info->alpha_matrix_sparse->i[i] == column){
 				M[column] = (sol_info->alpha_matrix_sparse->x[i] == 0) ? 1.0 : sol_info->alpha_matrix_sparse->x[i];
 				norm_b += pow(b_vector[column],2);
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				M[column] = 1.0;
 				column_count = sol_info->alpha_matrix_sparse->p[column + 1] - sol_info->alpha_matrix_sparse->p[column];
 			}
 		}
	
	}
	
	// Calculate ||b||
	norm_b = sqrt(norm_b);
	if (!norm_b)
		norm_b = 1;

	memset(x_vector,0,sizeof(double)*vector_size);

	double rho,rho1 = 1.0,beta,alpha;
	while(iter < vector_size){
		iter += 1;
		
		for(i=0;i<vector_size;i++)
			z_vector[i] = r_vector[i] / M[i];
		
		rho = 0.0;
		for(i=0;i<vector_size;i++)
			rho += r_vector[i]*z_vector[i];

		if(iter == 1){
			for(i=0;i<vector_size;i++)
				p_vector[i] = z_vector[i];
		}else{
			beta = rho / rho1;
			for(i=0;i<vector_size;i++)
				p_vector[i] = z_vector[i] + beta*p_vector[i];
		}
		rho1 = rho;
		if(get_sparse() == 0){
			for(i=0;i<vector_size;i++){
				q_vector[i] = 0;
				for(j=0;j<vector_size;j++){
					q_vector[i] += alpha_matrix[i*(vector_size)+j]*p_vector[j];
				}
			}
		}else{
			for(i=0;i<vector_size;i++) q_vector[i] = 0.0;
			cs_gaxpy (sol_info->alpha_matrix_sparse, p_vector, q_vector);
		}

		double temp = 0.0;
		for(i=0;i<vector_size;i++)
			temp += p_vector[i]*q_vector[i];
		alpha = rho / temp;
		norm_r = 0.0;
		for(i=0;i<vector_size;i++){
			x_vector[i] = x_vector[i] + alpha*p_vector[i]; 			
			r_vector[i] = r_vector[i] - alpha*q_vector[i];
			norm_r += pow(r_vector[i],2.0);
		}
		norm_r = sqrt(norm_r);
		norm = norm_r / norm_b;
		if(norm < itol){
			break;
		}

	}
	free(M); M = NULL;
	free(p_vector); p_vector = NULL;
	free(q_vector);	q_vector = NULL;
	free(r_vector);	r_vector = NULL;
	free(z_vector); z_vector = NULL;
	/*printf("\n");
	printf("x_vector\n");
	for(i=0;i<vector_size;i++){
		printf("%7.5f\n", x_vector[i]);
	}
	printf("------------------------------------------------------\n");*/
	return x_vector;
}












//double * _bicg_iter_solver(double * alpha_matrix, double *b_vector, unsigned int vector_size) {

/*
 * Bi-CG iterative solver 
 *
 */
double * _bicg_iter_solver(solver_info * sol_info){

	double * alpha_matrix = sol_info->alpha_matrix;
	double *b_vector = sol_info->b_vector;
	unsigned int vector_size = sol_info->alpha_matrix_size;
	//printf("bicg_iter_solver\n");

	if(vector_size == 0){
		fprintf(stderr, "[ERROR] vector_size too small\n");
		exit(-1);
	}
	double * r_vector_tilda = (double*)calloc(vector_size,sizeof(double));
	double * z_vector_tilda = (double*)malloc(vector_size*sizeof(double));
	double * p_vector_tilda = (double*)malloc(vector_size*sizeof(double));
	double * q_vector_tilda = (double*)malloc(vector_size*sizeof(double));
	double * M = (double*)malloc(vector_size*sizeof(double));
	double * x_vector = (double *) malloc(sizeof(double)*vector_size);
	//initial values:
	//                 _x[] = 0,...,0
	//                 _r = b

	double tol = 0.0001;
	int i = 0;
	for(int i =0;i<vector_size;i++) M[i] = 1.0;
	if(get_sparse() == 0){
		for(i=0;i<vector_size;i++){
			r_vector_tilda[i] = r_vector[i] = b_vector[i];
			M[i] = (alpha_matrix[i*vector_size+i] != 0) ? alpha_matrix[i*vector_size+i] : 1.0;
		}
	}else{
		int column = 0;
 		M[column] = 1.0;
		r_vector_tilda[column] = r_vector[column] = b_vector[column];
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = sol_info->alpha_matrix_sparse->p[column+1] - sol_info->alpha_matrix_sparse->p[column];
 		/*
 		XXX: nzmax ? non_zero
 		*/
 		for(i=0;i<sol_info->alpha_matrix_sparse->nzmax;i++){
 			if(sol_info->alpha_matrix_sparse->i[i] == column){
 				M[column] = (sol_info->alpha_matrix_sparse->x[i] == 0) ? 1.0 : sol_info->alpha_matrix_sparse->x[i];
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				//printf("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE %d %d %d %d\n", column, vector_size, i, non_zero);
 				r_vector_tilda[column] = r_vector[column] = b_vector[column];
 				M[column] = 1.0;
 				column_count = sol_info->alpha_matrix_sparse->p[column + 1] - sol_info->alpha_matrix_sparse->p[column];
 			}
 		}

	} 

	double temp = 0.0;
	for(i=0;i<vector_size;i++) temp += r_vector[i]*r_vector[i];
	double rho1 = temp; 
	
	// Calculate ||b||
	temp = 0.0;
	for(i=0;i<vector_size;i++) temp += b_vector[i]*b_vector[i];
	double norm_b = sqrt(temp);
	if (norm_b == 0) norm_b = 1;

	memset(x_vector,0,sizeof(double)*vector_size);
	const int max_iter = 2*vector_size; // (vector_size > 20) ? vector_size : 20; 
	int iter = 0;
	while(iter<max_iter) {

		for(i=0;i<vector_size;i++){
			z_vector[i] = r_vector[i]/M[i];
			z_vector_tilda[i] = r_vector_tilda[i]/M[i];
		}
		double rho = 0; 
		for(i=0;i<vector_size;i++) rho += r_vector_tilda[i]*z_vector[i];
		double abs_rho = fabs(rho);
		if(abs_rho < BICG_FPS){
			fprintf(stderr, "%s failed (abs(rho) < BICG_FPS) - bye bye now\n",__FUNCTION__);
			exit(-1);
		}
		if(iter == 0){
			for(i=0;i<vector_size;i++){
				p_vector[i] = z_vector[i];
				p_vector_tilda[i] = z_vector_tilda[i];
			}
		}else{
			double beta = rho/rho1;
			for(i=0;i<vector_size;i++){            	
				p_vector[i] =  z_vector[i] + beta*p_vector[i];
				p_vector_tilda[i] = z_vector_tilda[i] + beta*p_vector_tilda[i]; 
			}
		}
		rho1 = rho;
		int j;
		if(get_sparse() == 0){	
			for(i=0;i<vector_size;i++){
				q_vector[i] = 0;
				for(j=0;j<vector_size;j++){
					q_vector[i] += alpha_matrix[i*(vector_size)+j]*p_vector[j];
				}
			}

			for(i=0; i<vector_size; ++i) {
				double *col = &alpha_matrix[i];
				q_vector_tilda[i] = 0.0; 
				for(j=0;j<vector_size;j++){
					q_vector_tilda[i] += col[j*vector_size]*p_vector_tilda[j];
				}
			}
		}else{
			for(i=0;i<vector_size;i++) q_vector[i] = 0.0;
			cs_gaxpy (sol_info->alpha_matrix_sparse, p_vector, q_vector);
			int ap;
			for(i=0;i<vector_size;i++){
				q_vector_tilda[i] = 0.0;
				for(ap = sol_info->alpha_matrix_sparse->p[i]; ap < sol_info->alpha_matrix_sparse->p[i+1]; ap++)
					q_vector_tilda[i] = q_vector_tilda[i] + sol_info->alpha_matrix_sparse->x[ap] * p_vector_tilda[sol_info->alpha_matrix_sparse->i[ap]];
			}
		}
		

		double omega = 0.0;
		for(i=0;i<vector_size;i++) omega += p_vector_tilda[i] * q_vector[i];
		
		double abs_omega = fabs(rho);
		if(abs_omega < BICG_FPS){
			fprintf(stderr,"%s: failed (abs(omega) < BICG_FPS) - bye bye\n",__FUNCTION__);
			exit(-1);
		}
		double alpha = rho/omega;
		for(i=0;i<vector_size;i++){        	
			x_vector[i]  = x_vector[i]  + alpha*p_vector[i]; 
			r_vector[i]  = r_vector[i]  - alpha*q_vector[i]; 
			r_vector_tilda[i] = r_vector_tilda[i] - alpha*q_vector_tilda[i]; 
		}

		// Calculate ||r||
		temp = 0.0;
		for(i=0;i<vector_size;i++) temp += r_vector[i] * r_vector[i];

		double cond = sqrt(temp)/norm_b;
		if (cond < tol) break;

		iter++;
	}
	
	//for(i=0;i<vector_size;i++){
	//	printf(":%d :%lf\n",i,x_vector[i] );
	//}

	free(r_vector_tilda);
	free(z_vector_tilda);
	free(p_vector_tilda);
	free(q_vector_tilda);
	free(M); //I6
	return x_vector;
}
