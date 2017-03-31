/*
 * AC analysis declaration and info:
 *	.AC <sweep> <points> <start_freq> <end_freq>
 *			<sweep>      --> can be LIN or LOG, depending on how the user wants to sweep the space [<start_freq>,<end_freq>]
 *			<points>	 --> number of points in the space [<start_freq>,<end_freq>]		
 *			<start_freq> --> starting frequency
 *			<end_freq>   --> ending frequency
 *
 *
 *
 *
 *
 *
 * AC analysis system:
 *	_
 *  G * x = e
 *
 */

#include "freq_analysis.h"



#define BICG_FPS 1e-14
#define M_PI 3.14159265358979323846




/* x=A\b where A is unsymmetric; b overwritten with solution */
int cs_ci_lusol_ (int order, const cs_ci *A, cs_complex_t *b, cs_complex_t * x, double tol)
{
    cs_cis *S ;
    cs_cin *N ;
    int n, ok ;
   if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;
    S = cs_ci_sqr (order, A, 0) ;              /* ordering and symbolic analysis */
    N = cs_ci_lu (A, S, tol) ;                 /* numeric LU factorization */
   // x = cs_malloc (n, sizeof (cs_complex_t)) ;    /* get workspace */
    ok = (S && N) ;
    if (ok)
    {
        cs_ci_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
        cs_ci_lsolve (N->L, x) ;               /* x = L\x */
        cs_ci_usolve (N->U, x) ;               /* x = U\x */
        cs_ci_ipvec (S->q, x, b, n) ;          /* b(q) = x */
    }
    cs_ci_sfree (S) ;
    cs_ci_nfree (N) ;
    return (ok) ;
}


/* x=A\b where A is symmetric positive definite; b overwritten with solution */
int cs_cholsol_ (int order, const cs_ci *A, cs_complex_t *b, cs_complex_t *x)
{
    cs_cis *S ;
    cs_cin *N ;
    int n, ok ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;
    S = cs_ci_schol (order, A) ;               /* ordering and symbolic analysis */
    N = cs_ci_chol (A, S) ;                    /* numeric Cholesky factorization */
    ok = (S && N) ;
    if (ok)
    {
        cs_ci_ipvec (S->pinv, b, x, n) ;   /* x = P*b */
        cs_ci_lsolve (N->L, x) ;           /* x = L\x */
        cs_ci_ltsolve (N->L, x) ;          /* x = L'\x */
        cs_ci_pvec (S->pinv, x, b, n) ;    /* b = P'*x */
    }
    cs_ci_sfree (S) ;
    cs_ci_nfree (N) ;
    return (ok) ;
}




double * complex_solver_step( solver_info * sol_info){
	//printf("get_analysis_type: %d get_sparse: %d\n", get_analysis_type(), get_sparse());
	
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

	 //LU solver
	if(get_analysis_type() == 0){
		// No use of sparse
		if(get_sparse() == 0){

			gsl_complex b_temp;
			gsl_matrix_complex * a;
  			gsl_vector_complex * b;
  			gsl_vector_complex * x;
			gsl_permutation * p;
 
  			a = gsl_matrix_complex_alloc(sol_info->vector_size,sol_info->vector_size);
  			b = gsl_vector_complex_alloc(sol_info->vector_size);
  			x = gsl_vector_complex_alloc(sol_info->vector_size);
		  	p = gsl_permutation_alloc(sol_info->vector_size);

	 		for(int i=0;i<sol_info->vector_size;i++){
  				b_temp.dat[0] = creal(sol_info->b_vector[i]);
  				b_temp.dat[1] = cimag(sol_info->b_vector[i]);
  				gsl_vector_complex_set(b,i,b_temp);
			  	for(int j=0;j<sol_info->vector_size;j++){
			  		b_temp.dat[0] = creal(g_bar_matrix[i*sol_info->vector_size+j]);
			  		b_temp.dat[1] = cimag(g_bar_matrix[i*sol_info->vector_size+j]);
			  		gsl_matrix_complex_set(a, i, j, b_temp);
			  	}
  			}
  			//solve
		  	int temp;
		  	gsl_linalg_complex_LU_decomp(a, p, &temp);
		  	gsl_linalg_complex_LU_solve (a, p, b, x);
		  	//return results x_vector must be allocated
		  	for(int i=0;i<sol_info->vector_size;i++){
	  			gsl_complex x_val = gsl_vector_complex_get (x, i);
	  			sol_info->x_vector[i] = CMPLX(GSL_REAL(x_val), GSL_IMAG(x_val));
		  	}
		  	gsl_matrix_complex_free(a);
		  	gsl_vector_complex_free(b);
		  	gsl_vector_complex_free(x);
		  	gsl_permutation_free(p);

		// Use of sparse
		}else{
			int res = cs_ci_lusol_ (0,  sol_info->alpha_sparse, (cs_complex_t *) sol_info->b_vector, (cs_complex_t *) sol_info->x_vector, 0.1);
			if(!res){
				fprintf(stderr, "DUDEEEE sparse LU failed goodbye!\n");
				exit(-1);
			}
		}
	//Cholesky solver
	}else if(get_analysis_type() == 1){
		// No use of sparse
		if(get_sparse() == 0){

			gsl_complex b_temp;
			gsl_matrix_complex * a;
  			gsl_vector_complex * b;
  			gsl_vector_complex * x;
 
  			a = gsl_matrix_complex_alloc(sol_info->vector_size,sol_info->vector_size);
  			b = gsl_vector_complex_alloc(sol_info->vector_size);
  			x = gsl_vector_complex_alloc(sol_info->vector_size);

	 		for(int i=0;i<sol_info->vector_size;i++){
  				b_temp.dat[0] = creal(sol_info->b_vector[i]);
  				b_temp.dat[1] = cimag(sol_info->b_vector[i]);
  				gsl_vector_complex_set(b,i,b_temp);
			  	for(int j=0;j<sol_info->vector_size;j++){
			  		b_temp.dat[0] = creal(g_bar_matrix[i*sol_info->vector_size+j]);
			  		b_temp.dat[1] = cimag(g_bar_matrix[i*sol_info->vector_size+j]);
			  		gsl_matrix_complex_set(a, i, j, b_temp);
			  	}
  			}
  			//solve
			gsl_linalg_complex_cholesky_decomp(a);
			gsl_linalg_complex_cholesky_solve(a, b, x);
		  	//return results x_vector must be allocated
		  	for(int i=0;i<sol_info->vector_size;i++){
	  			gsl_complex x_val = gsl_vector_complex_get (x, i);
	  			sol_info->x_vector[i] = CMPLX(GSL_REAL(x_val), GSL_IMAG(x_val));
		  	}
		  	gsl_matrix_complex_free(a);
		  	gsl_vector_complex_free(b);
		  	gsl_vector_complex_free(x);

		// Use of sparse
		}else{
			int res = cs_cholsol_ (0,  sol_info->alpha_sparse, (cs_complex_t *) sol_info->b_vector, (cs_complex_t *) sol_info->x_vector);
			if(!res){
				fprintf(stderr, "DUDEEEE sparse Cholesky failed goodbye!\n");
				exit(-1);
			}
		}
	// Conjugate Gradient(CG) solver
	}else if(get_analysis_type() == 2){	
		sol_info->x_vector = complex_cg_solver(sol_info);
	// Bi-Conjugate Gradient(Bi-CG) solver
	}else if(get_analysis_type() == 3){
		sol_info->x_vector = complex_bicg_solver(sol_info);			
	}
	return NULL;
}







int run_freq_analysis(int command_id, command_type * cmnd){
	//ac_end_freq = ac_end_freq*2.0*M_PI;
	//ac_start_freq = ac_start_freq*2.0*M_PI;
	printf("ac_end_freq: %lf\n",ac_end_freq );
	printf("ac_start_freq: %lf\n",ac_start_freq );
	//exit(-1);

    cs_ci * temp;
	solver_info sol_info;
	int vector_size = (n+m2);

    if(get_sparse()){
    	sol_info.alpha_sparse = cs_ci_spalloc(vector_size,vector_size,c_bar_non_zero,1,1);
    	if(sol_info.alpha_sparse == NULL){
    		fprintf(stderr, "[ERROR] %s %s %d\n",__FILE__,__FUNCTION__,__LINE__ );
    		exit(-1);
    	}
    	temp = cs_ci_spalloc(vector_size,vector_size,non_zero,1,1);
    	for(int i=0;i<k_index;i++){
    		int x = alpha_matrix_sparse->i[i];
    		int y = alpha_matrix_sparse->p[i];
    		double value =  alpha_matrix_sparse->x[i];
    		cs_ci_entry(temp, (int) x, (int) y, value + 0.0*I);
    	}
    	cs_spfree((cs_di *)alpha_matrix_sparse);


 		cs_ci * alpha_sparse_cbar_t = cs_ci_compress(temp);
 		if(alpha_sparse_cbar_t){
 			cs_ci_free(temp);
 			temp = alpha_sparse_cbar_t;
 		}
 		if(!cs_ci_dupl(temp)){
    		cs_ci_print(temp,0);
    		printf("cs_ci_dupl failed sol_info->alpha_sparse_cbar\n");
    		exit(-1);
    	}


    	sol_info.alpha_sparse_cbar = temp;
    	sol_info.sparse_size = c_bar_non_zero;
    }

    char ** filenames;
    FILE ** fp;
    int plot_counter = 0;
    plot_counter = init_freq_output(&filenames, &fp, command_id, cmnd->plot_list);
    printf("%d\n",plot_counter );
    if(plot_counter == 0){
    	fprintf(stderr, "ERR: Nothing to print, printing all the results in stdout\n" );
    	//exit(-1);
    }

    for(int i=0;i<plot_counter;i++){
    	printf("a: %s b:%s\n",filenames[2*i], filenames[2*i+1] );
    }
    //exit(-1);

	double complex * b_vector;
	double complex * x_vector;
	double complex * b_vector_initial = NULL;
	b_vector = malloc(sizeof(double complex)*vector_size);
	x_vector = malloc(sizeof(double complex)*vector_size);
	memset(b_vector,0,sizeof(double complex)*vector_size);
	memset(x_vector,0,sizeof(double complex)*vector_size);
	e_vector_construction_im(b_vector, 0.0);
	if(get_sparse()){
		b_vector_initial = malloc(sizeof(double complex)*vector_size);
		memcpy(b_vector_initial,b_vector, sizeof(double complex)*vector_size);
	}

	sol_info.alpha_matrix = g_bar_matrix;
	sol_info.b_vector = b_vector;
	sol_info.x_vector = x_vector;
	sol_info.vector_size = vector_size;

	double step_delta = 0.0;
	if(cmnd->d.a_info.ac_sweep == 0)
  		step_delta = (cmnd->d.a_info.ac_end_freq - cmnd->d.a_info.ac_start_freq) / cmnd->d.a_info.ac_steps;
    else
    	step_delta = exp(log(10.0)/cmnd->d.a_info.ac_steps);
    	
    //for(double freq_i = ac_start_freq; freq_i <= ac_end_freq-step; freq_i += step){
    printf("step_delta:%f %d\n",step_delta,ac_steps );

    //exit(-1);
 //   step_delta = step_delta*ac_start_freq;

   	double freq_i = ac_start_freq;
   	unsigned int step = 0;
   	//for(int step=0;step<8*ac_steps;step++){
   	while(freq_i < cmnd->d.a_info.ac_end_freq){
   		//*exp(step*log(step_delta));
   		//,freq_i *= step_delta
   		if(cmnd->d.a_info.ac_sweep == 0)
   			freq_i = cmnd->d.a_info.ac_start_freq + step*step_delta;
   		else
	   		freq_i = cmnd->d.a_info.ac_start_freq*exp(step*log(step_delta));


    	if(get_sparse()){
    		memcpy(b_vector,b_vector_initial, sizeof(double complex)*vector_size);
		}
		g_bar_construction_im(&sol_info, 2*M_PI*freq_i);
		complex_solver_step(&sol_info);
		write_freq_output(sol_info.x_vector, ac_sweep, filenames, freq_i,cmnd->plot_list);
    
    	//hopefully here's the end!
	    if(plot_counter == 0){

    		//printf("freq_[%d] : %.10e %.10e\n",step,freq_i,step*step_delta );
		    printf("FREQ: %f --------------------------------------------\n", freq_i);
		  	for(int i =0;i<vector_size;i++){
		  		for(int j=0;j<vector_size;j++){
		  			printf("%06.4f %06.4fi\t",creal(g_bar_matrix[i*vector_size+j]), cimag(g_bar_matrix[i*vector_size+j]));
		  		}
		  		printf(" || %6.e %6.ei \t\t=\t %e %ei",
		  		//GSL_REAL(b_i),GSL_IMAG(b_i),
		  		creal(b_vector[i]), cimag(b_vector[i]),
		  		//GSL_REAL(x_i),GSL_IMAG(x_i)
		  		creal(sol_info.x_vector[i]),cimag(sol_info.x_vector[i])
		  		);
		  		printf("\n");
		  	}
		  	printf("------------------------------------------------------------\n");
	    }
	    step++;
  	}

/*  	if(plot_counter > 0){
		printf("------------------------------------------\n");
		printf("RESULT FILES -----------------------------\n");
		printf("------------------------------------------\n");
		for(int i=0;i<2*plot_counter;i++){
			printf("%s\n", filenames[i]);
		}
		printf("------------------------------------------\n");
  	}*/
 
	if(plot_counter > 0){

		FILE * gp_pipe;

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

		printf("gnuplot command start: ------------------------------------%d--\n",plot_counter);
		fprintf(stdout, "plot ");


		plot_print * pptemp;
		pptemp = cmnd->plot_list;
	
		for(int i=0;i<2*plot_counter;i++){
			//printf("%s\n", filenames[i]);
				if(gp_pipe){
						char * node_name;
						node_name = strstr(filenames[i],"node");
						if(node_name == NULL) node_name = filenames[i];
						else node_name += 5;

					if(i%2 == 0){
						if(i>0)
							fprintf(gp_pipe, "set terminal qt %d\n",i );
						if(cmnd->d.a_info.ac_sweep == 1)
							fprintf(gp_pipe, "set logscale x\n");
						if(pptemp->is_db)
							fprintf(gp_pipe, "set multiplot layout 1,2 title \"%s - (db)\" font \",14\"\n",cmnd->raw_cmnd);
						else
							fprintf(gp_pipe, "set multiplot layout 1,2 title \"%s\" font \",14\"\n",cmnd->raw_cmnd);
						fprintf(gp_pipe, "set xtics rotate\n");
						fprintf(gp_pipe, "set bmargin 2\n");

						fprintf(gp_pipe, "set title \"Magnitude of node: %s\" noenhanced\n",node_name);
						fprintf(gp_pipe, "unset key\n");
						fprintf(gp_pipe, "plot \"%s\" using 1:2 smooth unique\n",filenames[i] );
						pptemp = pptemp->next;
					}else{
						fprintf(gp_pipe, "set title \"Argument of node: %s\" noenhanced\n",node_name);
						fprintf(gp_pipe, "unset key\n");
						fprintf(gp_pipe, "plot \"%s\" using 1:2 smooth unique\n",filenames[i] );	
						fprintf(gp_pipe, "unset multiplot\n");					
					}
				}
			}
		
		printf("gnuplot command end: --------------------------------------\n");
	    if(gp_pipe)
		    fflush(gp_pipe);   // flush the gp_pipe to update the plot
		
		if(gp_pipe)
			pclose(gp_pipe);
		
		printf("------------------------------------------\n");
	}

  	return 1;
}



void write_freq_output(double complex * x_vector, int is_log, char ** filenames, double freq, plot_print * plist){
	int i = 0;
	plot_print * pptemp;
	pptemp = plist;//plot_print_list;
	FILE * tmp_fp;
	while(pptemp != NULL){
		if(pptemp->index == -1){
			pptemp = pptemp->next;
			continue;
		}

		tmp_fp = fopen(filenames[2*i], "a+");
		if(!tmp_fp){
			fprintf(stderr, "failed to open: %s\n", filenames[2*i]);
			pptemp = pptemp->next;
			i++;
			continue;
		}
		if(pptemp->is_db){ //log	
			fprintf(tmp_fp, "%.20e %.20e\n",
				freq,20 * log10(cabs(x_vector[pptemp->index-1])));
    	}else{ //lin
    		//printf("LIN LIN LIN %.20e %.20e %.20e\n", cabs(x_vector[pptemp->index-1]), creal(x_vector[pptemp->index-1]), cimag(x_vector[pptemp->index-1]));
    		//printf("LIN LIN LIN %.20e %.20e %.20e\n", cabs(x_vector[pptemp->index]), creal(x_vector[pptemp->index]), cimag(x_vector[pptemp->index]));

			fprintf(tmp_fp, "%.20e %.20e\n", freq,creal(x_vector[pptemp->index-1]));
    		//exit(-1);
    	}
    	fclose(tmp_fp);


		tmp_fp = fopen(filenames[2*i+1], "a+");
		if(!tmp_fp){
			fprintf(stderr, "failed to open: %s\n", filenames[2*i+1]);
			pptemp = pptemp->next;
			i++;
			continue;
		}
		if(is_log){ //log	
			fprintf(tmp_fp, "%.20e %.20e\n",
				freq,180*carg(x_vector[pptemp->index-1])/M_PI);
    	}else{ //lin
			fprintf(tmp_fp, "%.20e %.20e\n",
				freq,180*carg(x_vector[pptemp->index-1])/M_PI);
    	}
    	fclose(tmp_fp);




		pptemp = pptemp->next;
		i++;
	}
}

int init_freq_output(char *** filenames_, FILE *** fp_, int command_id, plot_print * plist){
	plot_print * pptemp;
	pptemp = plist;
	int plot_counter = 0;
	while(pptemp != NULL){
		if(pptemp->index != -1){
			printf("node name: %s\n",pptemp->name );
			plot_counter++;
		}
		pptemp = pptemp->next;
	}
	printf("plot_counter %d\n",plot_counter );
	//exit(-1);
	char ** filenames;
	FILE * fp;
//	fp = malloc(sizeof(FILE *)*(plot_counter+1));
//	if(!fp) return 0;
	filenames = malloc(sizeof(char *)*(2*plot_counter+1));
	if(!filenames) return 0;
	int i = 0;
	pptemp = plist;
	//time_t timestamp = time(NULL);
	if(!fopen("sol","r")){
		mkdir("sol", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	int index = 0;
	while(pptemp != NULL){
		if(pptemp->index == -1){
			pptemp = pptemp->next;
			continue;
		}

		int pos = 1;
		if(pptemp->name[1] == '(')
			pos = 2;


		if(get_sparse()){
			asprintf(&filenames[2*i  ],"%s/%d_AC_%d_sparse_MAG_node_%s",cs_info->sol_dir, command_id,index,&pptemp->name[pos]);
			asprintf(&filenames[2*i+1],"%s/%d_AC_%d_sparse_ARG_node_%s",cs_info->sol_dir, command_id,index,&pptemp->name[pos]);
		}else{
			asprintf(&filenames[2*i  ],"%s/%d_AC_%d_MAG_node_%s",cs_info->sol_dir, command_id,index,&pptemp->name[pos]);
			asprintf(&filenames[2*i+1],"%s/%d_AC_%d_ARG_node_%s",cs_info->sol_dir, command_id,index,&pptemp->name[pos]);
		}
		index++;

		for(int j=0;j<2;j++){
			fp = fopen(filenames[2*i+j], "w+");
			if(fp){
				fprintf(fp, "#*********************************************************************\n" );
				fprintf(fp, "#Starting Frequency:%lf Stop Frequency: %lf steps: %d node name: %s\n", ac_start_freq, 
					ac_end_freq, ac_steps, pptemp->name);
				fprintf(fp, "#*********************************************************************\n" );	
				fclose(fp);
			}else{
				fprintf(stderr, "Failed to open: %s\n",filenames[i] );
			}
		}

		pptemp = pptemp->next;
		i++;
	}
	*filenames_ = filenames;
	*fp_ = NULL;
	return plot_counter;
}




 /*
 *	_
 *  G * x = e
 *  |
 *  --> g_bar_matrix is the LHS matrix from the above system.
 *
 */

int g_bar_construction_im(solver_info * sol_info,double freq){
	//memset(g_bar_matrix,0,sizeof(double complex)*(n+m2)*(n+m2));
	if(!get_sparse())
		for(int i=0;i<(n+m2)*(n+m2);i++) g_bar_matrix[i] = 0.0 + _Complex_I*0.0;
	else{
		// 18+ XXX 
		sol_info->alpha_sparse->nz = 0;
	}
	capacitor * temp5 = (capacitor *)capacitor_list;
	while(temp5 != NULL && get_freq_type()){
		//sanity check
		if(temp5->positive_node == temp5->negative_node){
			temp5 = temp5->next;
			continue;
		}else if(temp5->positive_node == 0 && temp5->negative_node > 0){
			if(get_sparse() == 0)
				g_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->negative_node-1] += CMPLX(0.0,freq*(temp5->value));
			else{				
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->negative_node-1), 
					(int) (temp5->negative_node-1), 0.0 + (freq*temp5->value)*I);
			}
		}else if(temp5->positive_node  > 0 && temp5->negative_node == 0){
			if(get_sparse() == 0){
				g_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->positive_node-1] += CMPLX(0.0,freq*(temp5->value));
			}
			else{
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->positive_node-1), 
					(int) (temp5->positive_node-1), 0.0 + (freq*temp5->value)*I);
			}
		}else{
			if(get_sparse() == 0){
				g_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->negative_node-1] += CMPLX(0.0,freq*(temp5->value));
				g_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->positive_node-1] += CMPLX(0.0,freq*(temp5->value));
				g_bar_matrix[(temp5->positive_node-1)*(n+m2)+temp5->negative_node-1] -= CMPLX(0.0,freq*(temp5->value)); 
				g_bar_matrix[(temp5->negative_node-1)*(n+m2)+temp5->positive_node-1] -= CMPLX(0.0,freq*(temp5->value)); 
			}else{
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->negative_node-1), 
					(int) (temp5->negative_node-1), 0.0 + (freq*temp5->value)*I);
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->positive_node-1), 
					(int) (temp5->positive_node-1), 0.0 + (freq*temp5->value)*I);
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->positive_node-1), 
					(int) (temp5->negative_node-1), 0.0 - (freq*temp5->value)*I);	
				cs_ci_entry(sol_info->alpha_sparse, (int) (temp5->negative_node-1), 
					(int) (temp5->positive_node-1), 0.0 - (freq*temp5->value)*I);													
			}
		}
		temp5 = temp5->next;
	}

	int g_bar_index = n+1;
	coil * temp3 = (coil *)coil_list;
	while(temp3 != NULL){
		//sanity check
		if(temp3->positive_node == temp3->negative_node){
			temp3 = temp3->next;
			continue;
		}

		if(get_sparse()){
			cs_ci * matrix = sol_info->alpha_sparse; 
			int x = g_bar_index;
			int y = g_bar_index;
			double complex val = 0.0 - (freq*(temp3->value)*I);
			cs_ci_entry(matrix, x, y, val);
		}else{
			g_bar_matrix[g_bar_index*(n+m2)+g_bar_index] -= _Complex_I*freq*(temp3->value);
		}
		g_bar_index++;
		temp3 = temp3->next;
	}

	if(get_sparse()){
 		cs_complex_t one;
 		one = 1.0 + 0.0*I;
 		cs_ci * res;
 		cs_ci * alpha_sparse_t = cs_ci_compress(sol_info->alpha_sparse);
 		if(alpha_sparse_t){
 			cs_ci_free(sol_info->alpha_sparse);
 			sol_info->alpha_sparse = alpha_sparse_t;
 		}
    	if(!cs_ci_dupl(sol_info->alpha_sparse)){
    		cs_ci_print(sol_info->alpha_sparse,0);
    		printf("cs_ci_dupl failed sol_info->alpha_sparse\n");
    		exit(-1);
    	}
 		res = cs_ci_add(sol_info->alpha_sparse, sol_info->alpha_sparse_cbar, one, one);
 		cs_ci_free(sol_info->alpha_sparse);
 		sol_info->alpha_sparse = res;
 	}else{
		for(int i=0;i<(n+m2)*(n+m2);i++){
			g_bar_matrix[i] = g_bar_matrix[i] + alpha_matrix[i];
 		}
 	}

	return 0;
}



 /*
 *	_
 *  G * x = e
 *          |
 *          --> e_vector is the RHS vector from the above system.
 *
 */

int e_vector_construction_im(double complex * e_vector, double freq){
	current_source * temp2 = (current_source *)crnt_list;
	while(temp2 != NULL){
		//satiny check 2
		if(temp2->positive_node == temp2->negative_node){
			temp2 = temp2->next;
			continue;
		}else if(temp2->positive_node == 0){
			e_vector[temp2->negative_node-1] += temp2->ac_mag*CMPLX(cos(temp2->ac_phase*M_PI / 180.0), sin(temp2->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp2->ac_phase));
			//eulering
			//e_vector[temp2->negative_node-1] += temp2->ac_mag*CMPLX(cos(2*M_PI*freq + temp2->ac_phase), sin(2*M_PI*freq + temp2->ac_phase));
		}else if(temp2->negative_node == 0){
			e_vector[temp2->positive_node-1] -= temp2->ac_mag*CMPLX(cos(temp2->ac_phase*M_PI / 180.0), sin(temp2->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp2->ac_phase));
			//e_vector[temp2->positive_node-1] -= temp2->ac_mag*CMPLX(cos(2*M_PI*freq + temp2->ac_phase), sin(2*M_PI*freq + temp2->ac_phase));
		}else{
			e_vector[temp2->negative_node-1] += temp2->ac_mag*CMPLX(cos(temp2->ac_phase*M_PI / 180.0), sin(temp2->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp2->ac_phase));
			e_vector[temp2->positive_node-1] -= temp2->ac_mag*CMPLX(cos(temp2->ac_phase*M_PI / 180.0), sin(temp2->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp2->ac_phase));
			//e_vector[temp2->negative_node-1] += temp2->ac_mag*CMPLX(cos(2*M_PI*freq + temp2->ac_phase), sin(2*M_PI*freq + temp2->ac_phase));
			//e_vector[temp2->positive_node-1] -= temp2->ac_mag*CMPLX(cos(2*M_PI*freq + temp2->ac_phase), sin(2*M_PI*freq + temp2->ac_phase));
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
			//exit(1);
			e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(temp3->ac_phase*M_PI / 180.0), sin(temp3->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp3->ac_phase));
			//e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(2*M_PI*freq + temp3->ac_phase), sin(2*M_PI*freq + temp3->ac_phase));
		}else if(temp3->positive_node == 0){
			e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(temp3->ac_phase*M_PI / 180.0), sin(temp3->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp3->ac_phase));
			//e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(2*M_PI*freq + temp3->ac_phase), sin(2*M_PI*freq + temp3->ac_phase));
		}else{
			e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(temp3->ac_phase*M_PI / 180.0), sin(temp3->ac_phase*M_PI/180.0));
			//exp(_Complex_I*(temp3->ac_phase));
			//e_vector[e_vector_index] += temp3->ac_mag*CMPLX(cos(2*M_PI*freq + temp3->ac_phase), sin(2*M_PI*freq + temp3->ac_phase));
		}
		e_vector_index++;
		temp3 = temp3->next;
	}
	return 1;
}









 /*
 * Bi-CG iterative solver 
 * for complex numbers
 *
 */
double complex * complex_bicg_solver(solver_info * sol_info){

	double complex * alpha_matrix = sol_info->alpha_matrix;
	double complex * b_vector = sol_info->b_vector;
	double complex * x_vector = sol_info->x_vector;
	unsigned int vector_size = sol_info->vector_size;


	double complex * r_vector = 	(double complex *)calloc(vector_size,sizeof(double complex));
	double complex * z_vector = 	(double complex *)malloc(vector_size*sizeof(double complex));
	double complex * p_vector = 	(double complex *)malloc(vector_size*sizeof(double complex));
	double complex * q_vector = 	(double complex *)malloc(vector_size*sizeof(double complex));


	double complex * r_vector_tilda = 	(double complex *)calloc(vector_size,sizeof(double complex));
	double complex * z_vector_tilda = 	(double complex *)malloc(vector_size*sizeof(double complex));
	double complex * p_vector_tilda = 	(double complex *)malloc(vector_size*sizeof(double complex));
	double complex * q_vector_tilda = 	(double complex *)malloc(vector_size*sizeof(double complex));
	double complex * M = 				(double complex *)malloc(vector_size*sizeof(double complex));
	//double complex * x_vector = 		(double complex *)malloc(sizeof(double complex)*vector_size);
	//initial values:
	//                 _x[] = 0,...,0
	//                 _r = b

	double tol = 0.0001;
	int i = 0;
	if(get_sparse() == 0){
		for(i=0;i<vector_size;i++){
			r_vector_tilda[i] = r_vector[i] = b_vector[i];
			z_vector[i] = CMPLX(0.0,0.0);
			M[i] = (creal(alpha_matrix[i*vector_size+i]) != 0 && cimag(alpha_matrix[i*vector_size+i]) != 0) ? 
			alpha_matrix[i*vector_size+i] : CMPLX(1.0,0.0);
		}
	}else{
		int column = 0;
 		M[column] = 1.0;
		r_vector_tilda[column] = r_vector[column] = b_vector[column];
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = sol_info->alpha_sparse->p[column+1] - sol_info->alpha_sparse->p[column];
 		/*
 		XXX: nzmax ? non_zero
 		*/
 		for(i=0;i<sol_info->alpha_sparse->nzmax;i++){
 			if(sol_info->alpha_sparse->i[i] == column){
 				M[column] = (sol_info->alpha_sparse->x[i] == 0) ? 1.0 : sol_info->alpha_sparse->x[i];
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				r_vector_tilda[column] = r_vector[column] = b_vector[column];
 				M[column] = 1.0;
 				column_count = sol_info->alpha_sparse->p[column + 1] - sol_info->alpha_sparse->p[column];
 			}
 		}
	} 

	double temp = 0.0;
	for(i=0;i<vector_size;i++) temp += r_vector[i]*r_vector[i];
	
	double complex rho1 = CMPLX(temp,0.0); 
	
	// Calculate ||b||
	temp = 0.0;
	for(i=0;i<vector_size;i++) temp += cabs(b_vector[i])*cabs(b_vector[i]);
	double norm_b = sqrt(temp);
	if (norm_b == 0) norm_b = 1;

	memset(x_vector,0,sizeof(double complex)*vector_size);
	const int max_iter = 2 *vector_size; // (vector_size > 20) ? vector_size : 20; 
	int iter = 0;
	while(iter<max_iter) {

		for(i=0;i<vector_size;i++){
			z_vector[i] = r_vector[i]/ M[i];
			z_vector_tilda[i] = r_vector_tilda[i]/conj(M[i]);
		}

		double complex rho = CMPLX(0.0,0.0); 
		for(i=0;i<vector_size;i++) rho += conj(r_vector_tilda[i])*z_vector[i];
		
		double abs_rho = cabs(rho);
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
			double complex beta = rho/rho1;
			for(i=0;i<vector_size;i++){            	
				p_vector[i] =  z_vector[i] + beta*p_vector[i];
				p_vector_tilda[i] = z_vector_tilda[i] + conj(beta)*p_vector_tilda[i]; 
			}
		}
		rho1 = rho;
		int j;
		if(get_sparse() == 0){	
			for(i=0;i<vector_size;i++){
				q_vector[i] = CMPLX(0.0,0.0);
				for(j=0;j<vector_size;j++){
					q_vector[i] += alpha_matrix[i*(vector_size)+j]*p_vector[j];
				}
			}	

			for(j=0;j<vector_size;j++){
				q_vector_tilda[j] = 0;
				for(i=0;i<vector_size;i++){
					q_vector_tilda[j] += conj(alpha_matrix[i*(vector_size)+j])*p_vector_tilda[i];
				}
			}

			/*for(i=0; i<vector_size; ++i) {
				double *col = &alpha_matrix[i];
				q_vector_tilda[i] = 0.0; 
				for(j=0;j<vector_size;j++){
					q_vector_tilda[i] += col[j*vector_size]*p_vector_tilda[j];
				}
			}*/
		}else{

			for(i=0;i<vector_size;i++) q_vector[i] = 0.0;
			cs_ci_gaxpy (sol_info->alpha_sparse, p_vector, q_vector);
			int ap;
			for(i=0;i<vector_size;i++){
				q_vector_tilda[i] = 0.0;
				for(ap = sol_info->alpha_sparse->p[i]; ap < sol_info->alpha_sparse->p[i+1]; ap++)
					q_vector_tilda[i] = q_vector_tilda[i] + conj(sol_info->alpha_sparse->x[ap]) * p_vector_tilda[sol_info->alpha_sparse->i[ap]];
			}
		}
		

		double complex omega = CMPLX(0.0,0.0);
		for(i=0;i<vector_size;i++){
			omega += conj(p_vector_tilda[i]) * q_vector[i];
		}
		
		double abs_omega = cabs(rho); //sqrt(rho*conj(rho)); //fabs(rho);
		if(abs_omega < BICG_FPS){
			fprintf(stderr,"%s: failed (abs(omega) < BICG_FPS) - bye bye\n",__FUNCTION__);
			exit(-1);
		}

		//if(creal(omega) == 0.0 && cimag(omega) == 0.0) omega = CMPLX(1.0,1.0);
		double complex alpha = rho/omega;
		for(i=0;i<vector_size;i++){        	
			x_vector[i]  = x_vector[i]  + alpha*p_vector[i]; 
			r_vector[i]  = r_vector[i]  - alpha*q_vector[i]; 
			r_vector_tilda[i] = r_vector_tilda[i] - conj(alpha)*q_vector_tilda[i]; 
		}

		// Calculate ||r||
		temp = 0.0;
		for(i=0;i<vector_size;i++) temp += cabs(r_vector[i]) * cabs(r_vector[i]);

		double cond = sqrt(temp)/norm_b;
		if (cond < tol) break;

		iter++;
	}

	free(r_vector_tilda);
	free(z_vector_tilda);
	free(p_vector_tilda);
	free(q_vector_tilda);
	free(r_vector);
	free(z_vector);
	free(p_vector);
	free(q_vector); 
	free(M); //I6
	return x_vector;
}









//double * _cg_iter_solver(double * alpha_matrix, double *b_vector, unsigned int vector_size){

/*
 * CG iterative solver 
 * for complex numbers.
 *
 */
double complex * complex_cg_solver(solver_info * sol_info){
	
	double complex * alpha_matrix = sol_info->alpha_matrix;
	double complex * b_vector = sol_info->b_vector;
	double complex * x_vector = sol_info->x_vector;
	unsigned int vector_size = sol_info->vector_size;	

	int i,j;

	if(vector_size == 0){
		fprintf(stderr, "[ERROR] vector_size too small\n");
		exit(-1);
	}
	double complex * r_vector;
	double complex * z_vector;
	double complex * q_vector;
	double complex * p_vector;
	//double complex * x_vector;

	r_vector = malloc(sizeof(double complex)*vector_size);
	z_vector = malloc(sizeof(double complex)*vector_size);
	q_vector = malloc(sizeof(double complex)*vector_size);
	p_vector = malloc(sizeof(double complex)*vector_size);
	//x_vector = malloc(sizeof(double complex)*vector_size);
	
	for(i=0;i<vector_size;i++){
		r_vector[i] = b_vector[i];
		z_vector[i] = CMPLX(0.0,0.0);
		q_vector[i] = CMPLX(0.0,0.0);
		p_vector[i] = CMPLX(0.0,0.0);
	}
	unsigned int iter = 0;
	double complex * M;
	//double * M2;
	double norm, norm_b = 0.0, norm_r;
	

	M = malloc(sizeof(double complex)*(vector_size));
	memset(M,0, sizeof(double complex)*vector_size);
	//M2 = malloc(sizeof(double)*(vector_size));
	if(get_sparse() == 0){
		for(i=0;i<(vector_size);i++){
			x_vector[i] = CMPLX(0.0,0.0);

			M[i] = (creal(alpha_matrix[i*vector_size+i]) != 0 && cimag(alpha_matrix[i*vector_size+i]) != 0) ? 
			alpha_matrix[i*vector_size+i] : CMPLX(1.0,0.0);
			norm_b += cabs(b_vector[i])*cabs(b_vector[i]);
		}
	}else{
    	//sparse
    	norm_b = 0.0;
 		int column = 0;
 		M[column] = 1.0;
 		int column_count = 0;
 		if(vector_size > 0) 
 			column_count = sol_info->alpha_sparse->p[column+1] - sol_info->alpha_sparse->p[column];
 		for(i=0;i<non_zero;i++){
 			if(sol_info->alpha_sparse->i[i] == column){
 				M[column] = (sol_info->alpha_sparse->x[i] == 0) ? 1.0 : sol_info->alpha_sparse->x[i];
 				norm_b += cabs(b_vector[i])*cabs(b_vector[i]); //pow(b_vector[column],2);
 			}
 			column_count--;
 			if(column_count == 0){
 				column++;
 				if(column == vector_size) break;
 				M[column] = 1.0;
 				column_count = sol_info->alpha_sparse->p[column + 1] - sol_info->alpha_sparse->p[column];
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
			rho += r_vector[i]*conj(z_vector[i]);

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
			cs_ci_gaxpy (sol_info->alpha_sparse, p_vector, q_vector);
		}

		double temp = 0.0;
		for(i=0;i<vector_size;i++)
			temp += p_vector[i]*conj(q_vector[i]);
		alpha = rho / temp;
		norm_r = 0.0;
		for(i=0;i<vector_size;i++){
			x_vector[i] = x_vector[i] + alpha*p_vector[i]; 			
			r_vector[i] = r_vector[i] - alpha*q_vector[i];
			norm_r += cabs(r_vector[i])*cabs(r_vector[i]);
		}
		norm_r = sqrt(norm_r);
		norm = norm_r / norm_b;
		if(norm < itol){
			break;
		}

	}
	free(M); M = NULL;
	free(p_vector); p_vector = NULL;
	free(q_vector); q_vector = NULL;
	free(r_vector); r_vector = NULL;
	free(z_vector); z_vector = NULL;
	/*printf("\n");
	printf("x_vector\n");
	for(i=0;i<vector_size;i++){
		printf("%7.5f\n", x_vector[i]);
	}
	printf("------------------------------------------------------\n");*/
	return x_vector;
}


