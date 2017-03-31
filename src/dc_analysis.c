/*
 * DC analysis declaration, info and options:
 *	.DC <input_variable> <start_value> <end_value> <increment>
 *			<input_variable>  --> Current/Voltage source
 *			<start_value>	  --> dc starting value
 *			<end_value>       --> dc final value
 *			<increment>       --> step for DC sweep analysis
 *
 *
 *
 *	.OPTIONS                 --> LU
 *	.OPTIONS SPD             --> Cholesky
 *	.OPTIONS ITER            --> Bi-DG
 *	.OPTIONS ITER SPD        --> CG
 *	.OPTIONS SPARSE          --> Sparse LU
 *	.OPTIONS SPARSE SPD      --> Sparse Cholesky
 *	.OPTIONS SPARSE ITER     --> Sparse Bi-CG
 *	.OPTIONS SPARSE ITER SPD --> Sparse CG
 *
 *
 *
 *
 *
 * DC analysis system:
 * 
 *  G * x(t) = e
 *  
 *
 */


#define _GNU_SOURCE
#include <stdlib.h>
#include "dc_analysis.h"
//#include "element_structs.h"
#include "node_to_int.h"
#include "transient_analysis.h"
#include <string.h>
#include "util.h"
#include "freq_analysis.h"
#include <stdio.h>
#include <math.h>


#define BICG_FPS 1e-14
#define M_PI 3.14159265358979323846


void run_gnuplot(char **filenames, int plot_counter, command_type *cmnd){

	if(plot_counter > 0){

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
			fprintf(gp_pipe, "set title \"%s\" font \",14\"\n",cmnd->raw_cmnd );
			fprintf(gp_pipe, "plot ");
		}
		for(int i=0;i<plot_counter;i++){
			char *node_name;
			node_name = strstr(filenames[i],"node");
			if(node_name == NULL) node_name = filenames[i];
			else node_name += 5;
			//printf("%s\n", filenames[i]);
				if(gp_pipe){
					fprintf(gp_pipe, "\"%s\" using 1:2 title \"node: %s\" noenhanced smooth unique %s", filenames[i], 
						node_name,
						//&filenames[i][35],
					(i != plot_counter-1)?", ": " \n");
			//		 fflush(gp_pipe); 
				}
				fprintf(stdout, "\"%s\" using 1:2 title \"node: %s\" smooth unique %s", filenames[i], 
					node_name,
					//&filenames[i][35],
				(i != plot_counter-1)?", ": " \n");
			}
		
		printf("gnuplot command end: --------------------------------------\n");
	    if(gp_pipe)
		    fflush(gp_pipe);   // flush the gp_pipe to update the plot
		
		if(gp_pipe)
			pclose(gp_pipe);
		
		printf("------------------------------------------\n");
	}else{
		printf("plot_counter 0\n");
	}
}


int init_analysis(unsigned int _n, unsigned int m2){
	int i;
	if(_n == 0 && m2 == 0)
		return -1;
	vector_size = _n + m2;
	//n Num of nodes - GND
	alpha_matrix_size = _n+m2;





	if(get_sparse() == 0){
		alpha_matrix = malloc(sizeof(double)*(_n+m2)*(_n+m2));
		alpha_matrix_initial = malloc(sizeof(double)*(_n+m2)*(_n+m2));
		if(get_tran()){
			printf("ALLOCATED c_bar_matrix\n");
			c_bar_matrix = malloc(sizeof(double)*(_n+m2)*(_n+m2));
		}
		if(get_freq_type()){
			g_bar_matrix = malloc(sizeof(double complex)*(_n+m2)*(_n+m2));
		}

		for(i=0;i<(_n+m2)*(_n+m2);i++){
			alpha_matrix[i] = 0.0;
			if(get_tran()){
				alpha_matrix_initial[i] = 0.0;
				c_bar_matrix[i] = 0.0;
			}
			if(get_freq_type()){
				g_bar_matrix[i] = 0.0;
			}
		}
	}else{
		alpha_matrix_sparse = cs_spalloc(alpha_matrix_size, alpha_matrix_size, non_zero, 1, 1);
		alpha_matrix_sparse->nzmax = non_zero;
		alpha_matrix_sparse->nz = non_zero;
		printf("non_zero: %d\n", non_zero);
		printf("c_bar_non_zero: %d\n", c_bar_non_zero);
		//alpha_matrix_sparse->nz = _n + m2;
		if(get_tran()){
			cbar_sparse = cs_spalloc(alpha_matrix_size, alpha_matrix_size, c_bar_non_zero, 1, 1);
			cbar_sparse->nzmax = c_bar_non_zero;
			cbar_sparse->nz = c_bar_non_zero;
		}
	}
	
	if(get_tran())
		x2_vector = malloc(sizeof(double)*(_n+m2));

	x_vector = malloc(sizeof(double)*(_n+m2));
	for(i=0;i<(vector_size);i++){
		x_vector[i] = 0.0;
		if(get_tran())
			x2_vector[i] = 0.0;
	}
	b_vector = malloc(sizeof(double)*(_n+m2));
	z_vector = malloc(sizeof(double)*(_n+m2));
	r_vector = malloc(sizeof(double)*(_n+m2));
	if(!r_vector) exit(-1);
	printf("NM2: %d\n", _n+m2);
	p_vector = malloc(sizeof(double)*(_n+m2));
	q_vector = malloc(sizeof(double)*(_n+m2));
	n = _n;
	for(i=0;i<(_n+m2);i++){ x_vector[i] = b_vector[i] = z_vector[i] = p_vector[i] = r_vector[i] = q_vector[i] = 0.0;}
	return 1;
}




/* print a sparse matrix */
int mycs_print_DZ (const cs *A, int brief)
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


int mna_construction(){
	//printf("c_bar_non_zero: %d\n", c_bar_non_zero);
	//exit(-1);
	printf("get_tran(): %d\n", get_tran());

	resistor * temp;
	k_index = 0;
	unsigned int cbar_index = 0;
	temp = (resistor *)resistor_list;
	while(temp != NULL){
		//sanity check
		if(temp->positive_node == temp->negative_node){
			temp = temp->next;
			continue;
		}else if(temp->positive_node == 0 && temp->negative_node > 0){
			if(get_sparse() == 0)
				alpha_matrix[(temp->negative_node-1)*(n+m2)+temp->negative_node-1] += 1/(temp->value);
			else{	
				alpha_matrix_sparse->i[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->p[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->x[k_index] = 1/(temp->value);
				k_index++;
			}
		}else if(temp->positive_node  > 0 && temp->negative_node == 0){
			if(get_sparse() == 0)
				alpha_matrix[(temp->positive_node-1)*(n+m2)+temp->positive_node-1] += 1/(temp->value);
			else{
				alpha_matrix_sparse->i[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->p[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1/(temp->value);
				k_index++;
			}
		}else{
			if(get_sparse() == 0){
				alpha_matrix[(temp->negative_node-1)*(n+m2)+temp->negative_node-1] += 1/(temp->value);
				alpha_matrix[(temp->positive_node-1)*(n+m2)+temp->positive_node-1] += 1/(temp->value);
				alpha_matrix[(temp->positive_node-1)*(n+m2)+temp->negative_node-1] -= 1/(temp->value); 
				alpha_matrix[(temp->negative_node-1)*(n+m2)+temp->positive_node-1] -= 1/(temp->value); 
			}else{
				alpha_matrix_sparse->i[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->p[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->x[k_index] = 1/(temp->value);
				k_index++;
				alpha_matrix_sparse->i[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->p[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1/(temp->value);
				k_index++;
				alpha_matrix_sparse->i[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->p[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->x[k_index] = (-1)/(temp->value);
				k_index++;
				alpha_matrix_sparse->i[k_index] = temp->negative_node-1;
				alpha_matrix_sparse->p[k_index] = temp->positive_node-1;
				alpha_matrix_sparse->x[k_index] = (-1)/(temp->value);
				k_index++;
			}
		}
		temp = temp->next;
	}

	current_source * temp2 = (current_source *)crnt_list;
	while(temp2 != NULL){
		//satiny check 2
		if(temp2->positive_node == temp2->negative_node){
			temp2 = temp2->next;
			continue;
		}else if(temp2->positive_node == 0){
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp2->is_ac))
				b_vector[temp2->negative_node-1] += temp2->value;
		}else if(temp2->negative_node == 0){
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp2->is_ac))
				b_vector[temp2->positive_node-1] -= temp2->value;
		}else{
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp2->is_ac)){
				b_vector[temp2->negative_node-1] += temp2->value;
				b_vector[temp2->positive_node-1] -= temp2->value;
			}
		}
		temp2 = temp2->next;
	}
	unsigned int b_vector_index = n; // n+1 - hope so
	voltage_source * temp3 = (voltage_source *)vltg_list;
	while(temp3 != NULL){
		//printf("NAME: %s %d\n",temp3->name , temp3->negative_node-1 );
		//printf("vlt:%lf nn:%d pn:%d\n", temp3->value, temp3->negative_node, temp3->positive_node);
		//sanity check
		if(temp3->positive_node == temp3->negative_node){
			temp3 = temp3->next;
			continue;
		}else if(temp3->negative_node == 0){
			//printf("b_vector_index:%d \n", b_vector_index);
			//printf("positive_node-1 : %d\n",temp3->positive_node-1 );
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp3->is_ac))
				b_vector[b_vector_index] += temp3->value;

			if(get_sparse() == 0){
				alpha_matrix[(temp3->positive_node-1)*(n+m2)+b_vector_index] += 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp3->positive_node-1] += 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp3->positive_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp3->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
			}
		}else if(temp3->positive_node == 0){
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp3->is_ac))
				b_vector[b_vector_index] += temp3->value;

			if(get_sparse() == 0){
				alpha_matrix[(temp3->negative_node-1)*(n+m2)+b_vector_index] -= 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp3->negative_node-1] -= 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp3->negative_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp3->negative_node-1;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++; 
			}
		}else{
			if(get_analysis_type() < 4 || (get_analysis_type() > 3 && temp3->is_ac))
				b_vector[b_vector_index] += temp3->value;
		
			if(get_sparse() == 0){
				alpha_matrix[(temp3->positive_node-1)*(n+m2)+b_vector_index] += 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp3->positive_node-1] += 1.0;
				alpha_matrix[(temp3->negative_node-1)*(n+m2)+b_vector_index] -= 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp3->negative_node-1] -= 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp3->positive_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp3->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = temp3->negative_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index; 
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp3->negative_node-1;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
			}
		}
		temp3->b_vector_pos = b_vector_index;
		b_vector_index++;
		temp3 = temp3->next;
	}


	coil * temp4 = (coil *)coil_list;
	while(temp4 != NULL){
		if(get_tran() &&  temp4->positive_node != temp4->negative_node){
			if(get_sparse() == 0){
				c_bar_matrix[b_vector_index*(n+m2)+b_vector_index] -= temp4->value;
			}else{
				//SParse..
				cbar_sparse->x[cbar_index] = -temp4->value;
				cbar_sparse->p[cbar_index] = b_vector_index;
				cbar_sparse->i[cbar_index] = b_vector_index;
				cbar_index++;
			}
		}
		if(temp4->positive_node == temp4->negative_node){
			temp4 = temp4->next;
			continue;
		}else if(temp4->negative_node == 0){
			if(get_sparse() == 0){
				alpha_matrix[(temp4->positive_node-1)*(n+m2)+b_vector_index] += 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp4->positive_node-1] += 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp4->positive_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp4->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
			}
		}else if(temp4->positive_node == 0){
			if(get_sparse() == 0){
				alpha_matrix[(temp4->negative_node-1)*(n+m2)+b_vector_index] -= 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp4->negative_node-1] -= 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp4->negative_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp4->negative_node-1;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
 			}
		}else {
			if(get_sparse() == 0){	
				alpha_matrix[(temp4->positive_node-1)*(n+m2)+b_vector_index] += 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp4->positive_node-1] += 1.0;
				alpha_matrix[(temp4->negative_node-1)*(n+m2)+b_vector_index] -= 1.0;
				alpha_matrix[(b_vector_index)*(n+m2)+temp4->negative_node-1] -= 1.0;
			}else{
				alpha_matrix_sparse->i[k_index] = temp4->positive_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp4->positive_node-1;
				alpha_matrix_sparse->x[k_index] = 1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = temp4->negative_node-1;
				alpha_matrix_sparse->p[k_index] = b_vector_index;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
				alpha_matrix_sparse->i[k_index] = b_vector_index;
				alpha_matrix_sparse->p[k_index] = temp4->negative_node-1;
				alpha_matrix_sparse->x[k_index] = -1.0;
				k_index++;
			}
		}
		b_vector_index++;
		temp4 = temp4->next;
	}

	capacitor * temp5 = (capacitor *)capacitor_list;
	while(temp5 != NULL && get_tran()){
		//sanity check
		if(temp5->positive_node == temp5->negative_node){
			temp5 = temp5->next;
			continue;
		}else if(temp5->positive_node == 0 && temp5->negative_node > 0){
			if(get_sparse() == 0)
				c_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->negative_node-1] += (temp5->value);
			else{
				cbar_sparse->x[cbar_index] = temp5->value;
				cbar_sparse->i[cbar_index] = temp5->negative_node-1;
				cbar_sparse->p[cbar_index] = temp5->negative_node-1;
				cbar_index++;
			}
		}else if(temp5->positive_node  > 0 && temp5->negative_node == 0){
			if(get_sparse() == 0)
				c_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->positive_node-1] += (temp5->value);
			else{
				cbar_sparse->i[cbar_index] = temp5->positive_node-1;
				cbar_sparse->p[cbar_index] = temp5->positive_node-1;
				cbar_sparse->x[cbar_index] = temp5->value;
				cbar_index++;
			}
		}else{
			if(get_sparse() == 0){
				c_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->negative_node-1] += (temp5->value);
				c_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->positive_node-1] += (temp5->value);
				c_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->negative_node-1] -= (temp5->value); 
				c_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->positive_node-1] -= (temp5->value); 
			}else{
				cbar_sparse->i[cbar_index] = temp5->negative_node-1;
				cbar_sparse->p[cbar_index] = temp5->negative_node-1;
				cbar_sparse->x[cbar_index] = temp5->value;
				cbar_index++;
				cbar_sparse->i[cbar_index] = temp5->positive_node-1;
				cbar_sparse->p[cbar_index] = temp5->positive_node-1;
				cbar_sparse->x[cbar_index] = temp5->value;
				cbar_index++;
				cbar_sparse->i[cbar_index] = temp5->positive_node-1;
				cbar_sparse->p[cbar_index] = temp5->negative_node-1;
				cbar_sparse->x[cbar_index] = -temp5->value;
				cbar_index++;
				cbar_sparse->i[cbar_index] = temp5->negative_node-1;
				cbar_sparse->p[cbar_index] = temp5->positive_node-1;
				cbar_sparse->x[cbar_index] = -temp5->value;
				cbar_index++;
			}
		}
		temp5 = temp5->next;
	}


	if(get_tran() && get_sparse() == 0){
		for(int i=0;i<alpha_matrix_size*alpha_matrix_size;i++)
			alpha_matrix_initial[i] = alpha_matrix[i];	
	}

	if(get_sparse() == 1){
		//printf("QQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n\n\n");
		alpha_matrix_compressed = cs_compress(alpha_matrix_sparse);
		//cs_spfree(alpha_matrix_sparse);
		cs_dupl(alpha_matrix_compressed);

		cbar_compressed = cs_compress(cbar_sparse);
		cs_dupl(cbar_compressed);
	}

	return 1;
}



void printf_matrix(){
	int i,j;
	if(get_sparse() == 0){

		n = table_counter;	
		printf("m2: %d n: %d\n",m2,n);
		for(i=0;i<(n+m2);i++){
			for(j=0;j<(n+m2);j++){ 
				printf("%10lf ",alpha_matrix[i*(n+m2)+j]);		
			}
			printf(" | ");
			printf("%9lf  =", x_vector[i]);
			printf(" %9lf\n", b_vector[i]);
		}
	}else{
		//for(i=0;i<non_zero;i++){
		//printf("%d %lf\n",i,alpha_matrix_sparse->x[i]);
		//}
		//printf("NZ: %d\n", alpha_matrix_sparse->nz );
		//cs_print(alpha_matrix_compressed, "alpha_matrix_compressed",0);
	//	cs_print(alpha_matrix_sparse, "alpha_matrix_sparse",0);

	}
}


gsl_permutation *  mna_solver(){
	char * filenames[] = {"LU_sol","cholesky_sol", "cg_sol", "bicg_sol"};
	char * filenames_sp[] = {"LU_sol_sparse","cholesky_sol_sparse", "cg_sol_sparse", "bicg_sol_sparse"};
	gsl_vector *x = NULL;
	gsl_permutation * p = NULL;
	double *x_data = NULL;
	unsigned int x_size = 0;
	// LU - Genaral Matrices
	if(get_analysis_type() == 0){
		if(get_sparse() == 0){
			int s;
			printf_matrix();
			x = gsl_vector_alloc (alpha_matrix_size);
			gsl_matrix_view m = gsl_matrix_view_array(alpha_matrix, alpha_matrix_size, alpha_matrix_size);
			gsl_vector_view b = gsl_vector_view_array(b_vector, alpha_matrix_size);
			p = gsl_permutation_alloc (alpha_matrix_size);
			gsl_linalg_LU_decomp (&m.matrix, p, &s);
			gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
			x_data = x->data;
			x_size = x->size;
		}else{
			//printf_matrix();
			//printf("MATRIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			//mycs_print_DZ(alpha_matrix_compressed,0);
			csS = cs_di_sqr(2,alpha_matrix_compressed,0);
			alpha_matrix_compressed->nzmax = alpha_matrix_compressed->n;
			csN = cs_di_lu(alpha_matrix_compressed, csS,1);
			//cs_spfree(alpha_matrix_compressed);
			//end of alloc
			cs_di_ipvec(csN->pinv,b_vector,x_vector,vector_size);
			double * b_vector_init = malloc(sizeof(double)*vector_size);
			memcpy(b_vector_init,b_vector,sizeof(double)*vector_size);
			//for(int i=0;i<vector_size;i++){
			//	printf("%d %f\n",i, b_vector[i] );
			//}
			cs_di_lsolve(csN->L,x_vector);
			cs_di_usolve(csN->U,x_vector);
			cs_di_ipvec(csS->q,x_vector,b_vector,vector_size);
			x_data = b_vector;
			b_vector = b_vector_init;
			x_size = vector_size;
		}
	}
	//Cholesky - SPD Matrices
	else if(get_analysis_type() == 1){
		if(get_sparse() == 0){
			gsl_matrix_view m2 = gsl_matrix_view_array(alpha_matrix, alpha_matrix_size, alpha_matrix_size);
			gsl_vector_view b2 = gsl_vector_view_array(b_vector, alpha_matrix_size);
			gsl_linalg_cholesky_decomp(&m2.matrix);
			x = gsl_vector_alloc (alpha_matrix_size);
			gsl_linalg_cholesky_solve(&m2.matrix, &b2.vector, x);
			x_data = x->data;
			x_size = x->size;
		}else{
			csS = cs_schol(1, alpha_matrix_compressed);
			csN = cs_chol(alpha_matrix_compressed, csS);
			//cs_spfree(alpha_matrix_compressed);
			//end of alloc
			cs_ipvec(csS->pinv,b_vector,x_vector,vector_size);
			double * b_vector_init = malloc(sizeof(double)*vector_size);
			memcpy(b_vector_init,b_vector,sizeof(double)*vector_size);
			cs_lsolve(csN->L,x_vector);
			cs_ltsolve(csN->L,x_vector);
			cs_pvec(csS->pinv,x_vector,b_vector,vector_size);
			x_data = b_vector;
			b_vector = b_vector_init;
			x_size = vector_size;

		}
	}
	//Conjugate Gradient - SPD Matrices
	else if(get_analysis_type() == 2){	
		cg_iter_solver();
		x_data = x_vector;
		x_size = vector_size;
	}
	//Bi-Conjugate Gradient - Genaral Matrices
	else if(get_analysis_type() == 3){
		bicg_iter_solver();
		x_data = x_vector;
		x_size = vector_size;
	}
	//copy 
	int i;
	if(get_tran())
		for(i=0;i<x_size;i++)
			x2_vector[i] = x_data[i];



	// write output to file
	FILE *fp;
	if(get_sparse() == 0)
		fp = fopen(filenames[get_analysis_type()], "w+");
	else
		fp = fopen(filenames_sp[get_analysis_type()], "w+");
	if(fp && x_data != NULL){
		int i;	
		for(i=0;i<x_size;i++)
			if(get_nodename(i) == NULL){
				/*
				if(temp2 != NULL && temp2->name != NULL){
					fprintf(fp, "[I%s]: %g\n",temp2->name, x_data[i] );
					temp2 = temp2->next;
				}else{
					fprintf(fp, "[I]: %g\n",x_data[i] );
				}*/
			}else
				//fprintf(fp, "[V%s]: %g\n",get_nodename(i), x_data[i] );
				fprintf(fp, "%s  %.5e\n",get_nodename(i), x_data[i] );
		fclose(fp);
	}
	if(x != NULL)
		gsl_vector_free(x);
	if(get_analysis_type() == 0)
		return p;
	else 
		return NULL;
}

int init_dc_sweep_output(char *** filenames_, int command_id, command_type * sol_info){
	printf("init_dc_sweep_output\n");
	plot_print * pplist;
	dc_options * sweep_node;
	sweep_node = sol_info->d.d_info;
	plot_print * pptemp;
	pplist = sol_info->plot_list;
	pptemp = pplist;
	int plot_counter = 0;
	while(pptemp != NULL){
		if(pptemp->index != -1){
			printf("node_name: %s\n",pptemp->name );
			plot_counter++;
		}
		pptemp = pptemp->next;
	}


	if(plot_counter == 0) return 0;
	char ** filenames;
	FILE * fp;
	filenames = malloc(sizeof(char *)*(plot_counter+1));
	if(!filenames){
		fprintf(stderr, "ERROR: failed to allocate memory FILE:%s FUNCTION:%s LINE:%d\n", __FILE__, __FUNCTION__, __LINE__);
	} 
	int i = 0;
	pptemp = pplist;
	//time_t timestamp = time(NULL);
	if(!fopen("sol","r")){
		mkdir("sol", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	char * cmnd_names[] = {"LU","CHOLESKY","CG","BiCG"};
	while(pptemp != NULL){
		if(pptemp->index == -1){
			pptemp = pptemp->next;
			continue;
		}

		char *node_name_ptr;
		node_name_ptr = pptemp->name;
		if(node_name_ptr[0] == 'V' || node_name_ptr[0] == 'v' ||
			node_name_ptr[0] == 'I' || node_name_ptr[0] == 'i')
			node_name_ptr++;
		if(node_name_ptr[0] == '(') node_name_ptr++;
		char *node_name_ptr2 = node_name_ptr;
		while(*node_name_ptr2 != ')' && *node_name_ptr2 != '\0') node_name_ptr2++;
		if(*node_name_ptr2 == ')') *node_name_ptr2 = '\0';

		if(get_sparse()){

			asprintf(&filenames[i],"%s/%d_DC_SWEEP_%s_SPARSE_%s_node_%s",cs_info->sol_dir, command_id,sweep_node->node_name,
				cmnd_names[get_analysis_type()],node_name_ptr); // &pptemp->name[2]);
		}else{
			asprintf(&filenames[i],"%s/%d_DC_SWEEP_%s___%s_node_%s",cs_info->sol_dir, command_id,sweep_node->node_name,
				cmnd_names[get_analysis_type()], node_name_ptr); //&pptemp->name[2]);
		}

		printf("filename: %s\n", filenames[i]);
		fp = fopen(filenames[i], "w+");
		if(fp){
			fprintf(fp, "#*********************************************************************\n" );
			fprintf(fp, "#* %s from_voltage: %2.5f to_voltage:%2.5f step:%2.5f\n", sweep_node->node_name,sweep_node->from_voltage, 
				sweep_node->to_voltage, sweep_node->step);
			fprintf(fp, "#*********************************************************************\n" );
			fclose(fp);
		}else{
			fprintf(stderr, "Failed to open: %s\n",filenames[i] );
		}

		pptemp = pptemp->next;
		i++;
	}
	*filenames_ = filenames;
	return plot_counter;
}



void write_dc_sweep_output(double * x_vector, char ** filenames, double step, command_type * cs_info){
	if(x_vector == NULL) return;
	int i = 0;
	plot_print * pptemp;
	pptemp = cs_info->plot_list;
	FILE * tmp_fp;
	while(pptemp != NULL){
		if(pptemp->index == -1){
			pptemp = pptemp->next;
			continue;
		}

		tmp_fp = fopen(filenames[i], "a+");
		if(!tmp_fp){
			fprintf(stderr, "failed to open: %s\n", filenames[2*i]);
			pptemp = pptemp->next;
			i++;
			continue;
		}
		if(pptemp->index > 0)
			fprintf(tmp_fp, "%.20e %.20e\n",step,(x_vector[pptemp->index-1]));
    	else
			fprintf(tmp_fp, "%.20e %.20e\n",step,(x_vector[0]));

    	fclose(tmp_fp);
		pptemp = pptemp->next;
		i++;
	}
}



void dc_sweep_analysis(gsl_permutation * p, command_type * sol_info, int command_id){

	plot_print * pp_tmp;

		pp_tmp = sol_info->plot_list;
		while(pp_tmp != NULL){
			printf("\t---PLOT: node_name: %s l:%d index:%d\n",pp_tmp->name, pp_tmp->line, pp_tmp->index );
			pp_tmp = pp_tmp->next;
		}	




	//char * filenames[] = {"LU_sol_nodes","cholesky_sol_nodes","CG_sol_nodes","BiCG_sol_nodes"};
	//char * filenames_sp[] = {"LU_sol_nodes_sparse","cholesky_sol_nodes_sparse","CG_sol_nodes_sparse","BiCG_sol_nodes_sparse"};
	gsl_vector *x = NULL;
	gsl_matrix_view m2,m;
	//gsl_permutation * p;
	int i;	
	double *b_initial = NULL;
	for(i=0;i<(vector_size);i++){
		x_vector[i] = 0.0;
	}

	
	//---------------------------------------------------------------------------------------------
	// SOLVER INITIALIZATION BEGINS
	//---------------------------------------------------------------------------------------------
	// Cholesky - SPD Matrices
	if(get_analysis_type() == 1){
		if(get_sparse() == 0){
			m2 = gsl_matrix_view_array(alpha_matrix, alpha_matrix_size, alpha_matrix_size);
			x = gsl_vector_alloc (alpha_matrix_size);
		}else{
			csS = cs_schol(1, alpha_matrix_compressed);
			csN = cs_chol(alpha_matrix_compressed, csS);
			cs_spfree(alpha_matrix_compressed);

				b_initial = (double *)malloc(sizeof(double)*vector_size);
				memcpy(b_initial,b_vector, sizeof(double)*vector_size);
		}
	}
	// LU - Genaral Matrices
	else if (get_analysis_type() == 0){
		if(get_sparse() == 0){
			x = gsl_vector_alloc (alpha_matrix_size);
			m = gsl_matrix_view_array(alpha_matrix, alpha_matrix_size, alpha_matrix_size);
		}else{
			csS = cs_sqr(2,alpha_matrix_compressed,0);
			csN = cs_lu(alpha_matrix_compressed, csS,1);

				b_initial = (double *)malloc(sizeof(double)*vector_size);
				memcpy(b_initial,b_vector, sizeof(double)*vector_size);
		}
		printf("b_vector-----------------------------\n");
		for(int i=0;i<vector_size;i++){
			printf("%d %f\n",i,b_vector[i] );
		}
		printf("---------------------------------------\n");
	}
	//---------------------------------------------------------------------------------------------
	// SOLVER INITIALIZATION ENDS
	//---------------------------------------------------------------------------------------------
	
	
	voltage_source * temp2;// = vltg_list;
	dc_options * temp;
	//temp = dc_option_list;
	temp = sol_info->d.d_info;
	printf("dc_option_list 0x%x \n",temp );
	int pos = -1;
	double *x_data = NULL;
	double voltage_step; 

	int plot_counter;
	char **filenames;

	
	// FOR VOLTAGE SOURCE SWEEP
	if(temp->node_name[0] == 'V' || temp->node_name[0] == 'v'){
		pos = -1;
		temp2 = vltg_list;
		while(temp2 != NULL){
			if(strcasecmp(temp2->name,temp->node_name) == 0){
				pos = temp2->b_vector_pos; break;
			}
			temp2 = temp2->next;
		}
	}

	//while(temp != NULL){
	if(pos != -1){	
		printf("IN DC_SWEEP: %s\n",temp->node_name );
		/*printf("SWEEP V SRC: %s\n",temp->node_name );
		temp = temp->next;
		continue;
		temp2 = vltg_list;
		pos = -1;
		while(temp2 != NULL){
			if(strcmp(temp2->name,temp->node_name) == 0){
				pos = temp2->b_vector_pos; break;
			}
			temp2 = temp2->next;
		}
		if(pos == -1){
			temp = temp->next;
			continue;
		}*/
		double b_vector_pos_initial_value = b_vector[pos];
		b_vector[pos] -= temp2->value;
		double b_vector_temp = b_vector[pos];


/*		char * filename;

		FILE *fp;
		plot_print * pptemp;
		if(get_sparse() == 0){
//			asprintf(&filename, "%s/%d_DC_SWEEP_%s_SPARSE_");
			fp = fopen(filenames[get_analysis_type()], "w+");
		}else{
			fp = fopen(filenames_sp[get_analysis_type()], "w+");
		}
		fprintf(fp, "*********************************************************************\n" );
		fprintf(fp, "* %s from_voltage: %2.5f to_voltage:%2.5f step:%2.5f\n", temp->node_name,temp->from_voltage, 
			temp->to_voltage, temp->step);
		fprintf(fp, "* Watching Nodes: \n" );
		pptemp = plot_print_list;
		while(pptemp != NULL){

			if(pptemp->index != -1)
				fprintf(fp, "* %s\n", &pptemp->name[2] );
			pptemp = pptemp->next;
		}
		fprintf(fp, "*********************************************************************\n\n" );
		*/

		printf("INITE init_dc_sweep_output\n");
		plot_counter = init_dc_sweep_output(&filenames, command_id, sol_info);





		for(voltage_step = temp->from_voltage; voltage_step< temp->to_voltage; voltage_step += temp->step){
			//memset
			for(i=0;i<(vector_size);i++){
				x_vector[i] = 0.0;
			}

			b_vector[pos] = b_vector_temp + voltage_step;

			//---------------------------------------------------------------------------------------------
			// SOLVER BEGINS
			//---------------------------------------------------------------------------------------------
			//// SPD matrices /
			// Cholesky
			if(get_analysis_type() == 1){
				if(get_sparse() == 0){
					gsl_vector_view b2 = gsl_vector_view_array(b_vector, alpha_matrix_size);
					gsl_linalg_cholesky_solve(&m2.matrix, &b2.vector, x);
					x_data = x->data;
				}else{
					// sparsing..
					memcpy(b_vector,b_initial, sizeof(double)*vector_size);
					b_vector[pos] = b_vector_temp + voltage_step;
					cs_ipvec(csS->pinv,b_vector,x_vector,vector_size);
					cs_lsolve(csN->L,x_vector);
					cs_ltsolve(csN->L,x_vector);
					cs_pvec(csS->pinv,x_vector,b_vector,vector_size);
					x_data = b_vector;
				}
			}
			// Conjugate Gradient
			else if (get_analysis_type() == 2){
				cg_iter_solver();
				x_data = x_vector;
			}
			/// General matrices //
			// LU
			else if (get_analysis_type() == 0){
				//printf("DC_SWEEP WITH LU\n");
				if(get_sparse() == 0){
					gsl_vector_view b = gsl_vector_view_array(b_vector, alpha_matrix_size);
					gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
					x_data = x->data;
				}else{
					// sparsing ..
					memcpy(b_vector,b_initial, sizeof(double)*vector_size);
					b_vector[pos] = b_vector_temp + voltage_step;

					cs_ipvec(csN->pinv,b_vector,x_vector,vector_size);
					cs_lsolve(csN->L,x_vector);
					cs_usolve(csN->U,x_vector);
					cs_ipvec(csS->q,x_vector,b_vector,vector_size);
					x_data = b_vector;
				}

			}
			// Bi-Conjugate Gradient
			else if (get_analysis_type() == 3){					
				bicg_iter_solver();
				x_data = x_vector;
			}
			//---------------------------------------------------------------------------------------------
			// SOLVER END
			//---------------------------------------------------------------------------------------------

			write_dc_sweep_output(x_data, filenames, voltage_step, sol_info);

	
/*			if(fp){
				pptemp = plot_print_list;
				fprintf(fp, "cur_step: %2.5f\n", voltage_step );
				while(pptemp != NULL){
					if(pptemp->index != -1)
						fprintf(fp, "[%s]: %g - %d\n",&pptemp->name[2], x_data[pptemp->index-1], pptemp->index );
					pptemp = pptemp->next;
				}
				fprintf(fp, "\n" );
			}*/

		}
		
		run_gnuplot(filenames, plot_counter,sol_info);
		if(plot_counter > 0){

		for(int i=0;i<plot_counter;i++){
			free(filenames[i]);
		}
			free(filenames);
		}

		//fclose(fp);
		b_vector[pos] = b_vector_pos_initial_value;


	}

















	// FOR CURRENT SOURCE SWEEP
	current_source * temp3;// = current_list;
	current_source * curr_src = NULL; 
	if(temp->node_name[0] == 'I' || temp->node_name[0] == 'i'){
		//pos = -1;
		temp3 = crnt_list;
		while(temp3 != NULL){
			if(strcasecmp(temp3->name,temp->node_name) == 0){
				curr_src = temp3;
				//pos = temp3->b_vector_pos; break;
			}
			temp3 = temp3->next;
		}
	}

	//while(temp != NULL){
	if(curr_src != NULL){	
	

	//
	// FOR CURRENT SOURCE SWEEP
	/*temp = dc_option_list;

	while(temp != NULL){
		temp3 = crnt_list;
		while(temp3 != NULL){
			if(strcmp(temp3->name,temp->node_name) == 0){
				curr_src = temp3;
				break;
			}
			temp3 = temp3->next;
		}
		if(curr_src == NULL){
			temp = temp->next;
			continue;
		}*/


		double b_vector_pos_initial_value = b_vector[curr_src->positive_node];
		double b_vector_neg_initial_value = b_vector[curr_src->negative_node];

		printf("pos %g, neg %g\n", curr_src->value, curr_src->value);

		b_vector[curr_src->positive_node] += curr_src->value;
		b_vector[curr_src->negative_node] -= curr_src->value;
		double b_vector_temp_pos = b_vector[curr_src->positive_node];
		double b_vector_temp_neg = b_vector[curr_src->negative_node];
		/*
		FILE *fp;
		plot_print * pptemp;
		if(get_sparse() == 0)
			fp = fopen(filenames[get_analysis_type()], "w+");
		else
			fp = fopen(filenames_sp[get_analysis_type()], "w+");
		fprintf(fp, "*********************************************************************\n" );
		fprintf(fp, "* %s from_current: %2.5f to_current:%2.5f step:%2.5f\n", temp->node_name,temp->from_voltage, 
			temp->to_voltage, temp->step);
		fprintf(fp, "* Watching Nodes: \n" );
		pptemp = plot_print_list;
		while(pptemp != NULL){

			if(pptemp->index != -1)
				fprintf(fp, "* %s\n", &pptemp->name[2] );
			pptemp = pptemp->next;
		}
		fprintf(fp, "*********************************************************************\n\n" );
		*/


		plot_counter = init_dc_sweep_output(&filenames, command_id, sol_info);



		for(voltage_step = temp->from_voltage; voltage_step< temp->to_voltage; voltage_step += temp->step){
			for(i=0;i<(vector_size);i++){
				x_vector[i] = 0.0;
			}
			printf("current_step: %f\n", voltage_step);
			
			b_vector[curr_src->negative_node] = b_vector_temp_neg + voltage_step;
			b_vector[curr_src->positive_node] = b_vector_temp_pos - voltage_step;



			//----------------------------------------------------------------------------------------
			// SOLVER BEGINS
			//----------------------------------------------------------------------------------------
			/// SPD matrices /
			// Cholesky
			if(get_analysis_type() == 1){
				if(get_sparse() == 0){
					gsl_vector_view b2 = gsl_vector_view_array(b_vector, alpha_matrix_size);
					gsl_linalg_cholesky_solve(&m2.matrix, &b2.vector, x);
					if(x == NULL){
						fprintf(stderr, "%s %s %d\n",__FILE__,__FUNCTION__,__LINE__ );
						exit(-1);
					}
						x_data =  x->data;		
				}else{ // sparse..
					memcpy(b_vector,b_initial, sizeof(double)*vector_size);
					cs_ipvec(csS->pinv,b_vector,x_vector,vector_size);
					cs_lsolve(csN->L,x_vector);
					cs_ltsolve(csN->L,x_vector);
					cs_pvec(csS->pinv,x_vector,b_vector,vector_size);
					x_data = b_vector;
				}
			}
			// Conjugate Gradient
			else if (get_analysis_type() == 2){
				cg_iter_solver();
				x_data =  x_vector;
			}
			/// General matrices //
			// LU
			else if (get_analysis_type() == 0){
				if(get_sparse() == 0){
					gsl_vector_view b = gsl_vector_view_array(b_vector, alpha_matrix_size);
					gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
					if(x == NULL){
						fprintf(stderr, "%s %s %d\n",__FILE__,__FUNCTION__,__LINE__ );
						exit(-1);
					}
					x_data =  x->data;
				}else{ // sparse ..
					memcpy(b_vector,b_initial, sizeof(double)*vector_size);
					cs_ipvec(csN->pinv,b_vector,x_vector,vector_size);
					cs_lsolve(csN->L,x_vector);
					cs_usolve(csN->U,x_vector);
					cs_ipvec(csS->q,x_vector,b_vector,vector_size);
					x_data = b_vector;
				}

			}
			// Bi-Conjugate Gradient
			else if (get_analysis_type() == 3){					
				bicg_iter_solver();
				x_data =  x_vector;
			}
			//----------------------------------------------------------------------------------------
			// SOLVER ENDS
			//----------------------------------------------------------------------------------------
			write_dc_sweep_output(x_data, filenames, voltage_step, sol_info);

	
			/*if(fp && x_data != NULL){
				pptemp = plot_print_list;
				fprintf(fp, "cur_step: %2.5f\n", voltage_step );
				while(pptemp != NULL){
					if(pptemp->index != -1)
						fprintf(fp, "[%s]: %g - %d\n",&pptemp->name[2], x_data[pptemp->index-1], pptemp->index );
					pptemp = pptemp->next;
				}
				fprintf(fp, "\n" );
			}*/

		}
		printf("running run_gnuplot\n");
		run_gnuplot(filenames, plot_counter,sol_info);

		//fclose(fp);
		b_vector[curr_src->positive_node] = b_vector_pos_initial_value;
		b_vector[curr_src->negative_node] = b_vector_neg_initial_value;

		// temp = temp->next;	
	}

	//

	//if(get_analysis_type() == 0 ) 	    
		//gsl_permutation_free (p);
}


/*
 * CG iterative solver
 */
void cg_iter_solver(void){
	printf("vector_size: %d\n", vector_size);
	int i,j;
	for(i=0;i<vector_size;i++)
		r_vector[i] = b_vector[i];
	unsigned int iter = 0;
	double * M;
	//double * M2;
	double norm, norm_b = 1.0, norm_r;
	if(vector_size == 0){
		fprintf(stderr, "[ERROR] vector_size too small\n");
		exit(-1);
	}
	M = malloc(sizeof(double)*(vector_size));
	memset(M,0, sizeof(double)*vector_size);
	//M2 = malloc(sizeof(double)*(vector_size));
	if(get_sparse() == 0){
		for(i=0;i<(vector_size);i++){
			M[i] = (alpha_matrix[(i*(vector_size))+i] == 0.0) ? 1.0 : alpha_matrix[(i*(vector_size))+i];
			norm_b += pow(b_vector[i],2);
		}
	}else{
    	norm_b = 0.0;
 		int column = 0;
 		M[column] = 1.0;
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = alpha_matrix_compressed->p[column+1] - alpha_matrix_compressed->p[column];
 		for(i=0;i<non_zero;i++){
 			if(alpha_matrix_compressed->i[i] == column){
 				M[column] = (alpha_matrix_compressed->x[i] == 0) ? 1.0 : alpha_matrix_compressed->x[i];
 				norm_b += pow(b_vector[column],2);
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				M[column] = 1.0;
 				column_count = alpha_matrix_compressed->p[column + 1] - alpha_matrix_compressed->p[column];
 			}
 		}

	}
	memset(x_vector,0,sizeof(double)*vector_size);
	
	// Calculate ||b||
	norm_b = sqrt(norm_b);
	if (!norm_b)
		norm_b = 1;


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
			cs_gaxpy (alpha_matrix_compressed, p_vector, q_vector);
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
		//printf("x_vector cg iter: %d norm: %f\n",iter, norm);
		//for(i=0;i<vector_size;i++){
		//	printf("%lf\n",x_vector[i]);
		//}
		//printf("-----------------------------------------------------------\n");
	}
	free(M);
}


/*
 * Bi-CG iterative solver
 */
void bicg_iter_solver(void) {
	printf("bicg_iter_solver vector_size: %d\n", vector_size);
	if(vector_size == 0){
		fprintf(stderr, "[ERROR] vector_size too small\n");
		exit(-1);
	}
	double * r_vector_tilda = 	(double*)calloc(vector_size,sizeof(double));
	double * z_vector_tilda = 	(double*)malloc(vector_size*sizeof(double));
	double * p_vector_tilda = 	(double*)malloc(vector_size*sizeof(double));
	double * q_vector_tilda = 	(double*)malloc(vector_size*sizeof(double));
	double * M = 				(double*)malloc(vector_size*sizeof(double));

	//initial values:
	//                 _x[] = 0,...,0
	//                 _r = b

	double tol = 0.0001;
	int i = 0;
	for(int i =0;i<vector_size;i++) M[i] = 1.0;
	if(get_sparse() == 0){
		for(i=0;i<vector_size;i++){
			r_vector_tilda[i] = r_vector[i] = b_vector[i];
			z_vector[i] = 0.0;
			M[i] = (alpha_matrix[i*vector_size+i] != 0) ? alpha_matrix[i*vector_size+i] : 1.0;
		}
	}else{
		int column = 0;
 		M[column] = 1.0;
		r_vector_tilda[column] = r_vector[column] = b_vector[column];
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = alpha_matrix_compressed->p[column+1] - alpha_matrix_compressed->p[column];
 		/*
 		XXX: nzmax ? non_zero
 		*/
 		for(i=0;i<alpha_matrix_compressed->nzmax;i++){
 			if(alpha_matrix_compressed->i[i] == column){
 				M[column] = (alpha_matrix_compressed->x[i] == 0) ? 1.0 : alpha_matrix_compressed->x[i];
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				//printf("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE %d %d %d %d\n", column, vector_size, i, non_zero);
 				r_vector_tilda[column] = r_vector[column] = b_vector[column];
 				M[column] = 1.0;
 				column_count = alpha_matrix_compressed->p[column + 1] - alpha_matrix_compressed->p[column];
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
	const int max_iter = vector_size; // (vector_size > 20) ? vector_size : 20; 
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
			cs_gaxpy (alpha_matrix_compressed, p_vector, q_vector);
			int ap;
			for(i=0;i<vector_size;i++){
				q_vector_tilda[i] = 0.0;
				for(ap = alpha_matrix_compressed->p[i]; ap < alpha_matrix_compressed->p[i+1]; ap++)
					q_vector_tilda[i] = q_vector_tilda[i] + alpha_matrix_compressed->x[ap] * p_vector_tilda[alpha_matrix_compressed->i[ap]];
			}
		}
		

		double omega = 0.0;
		for(i=0;i<vector_size;i++) omega += p_vector_tilda[i] * q_vector[i];
		
		double abs_omega = fabs(omega);
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
}
