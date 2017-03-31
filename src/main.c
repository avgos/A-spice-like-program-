#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include "node_to_int.h"
#include "util.h"
#include "dc_analysis.h"
#include "element_structs.h"
#include <time.h>
#include "transient_analysis.h"
#include "freq_analysis.h"
#include "analysis_info.h"
#include <gsl/gsl_linalg.h> //this requires gsl

//remove gcc warning
extern int yylex(void);

char * command_names[] = {"DC", "DC_SWEEP", "TRANSIENT", "AC"};
char * solver_names[]  = {"LU", "CHOLESKY", "CG", "BICG"};


void mksoldir(crct_sim_spec * cs_info){
	int res;
	if(!fopen("sol","r")){
		res = mkdir("sol", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if(res == -1){
			fprintf(stderr, "ERROR: Failed to create sol dir. Check application permissions\n");
			exit(-1);
		}
	}
	char * dir = NULL;
	int folder_sub_id = 1;
	while(1){
		if(dir == NULL)
			res = asprintf(&dir, "sol/%lld", (long long)cs_info->timestamp);
		if(res == -1){
			fprintf(stderr, "ERROR: asprintf failed FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__);
			exit(-1);
		}
		if(fopen(dir,"r")){ //shit folder already exists (corner case)
			free(dir);
			res = asprintf(&dir, "sol/%lld_%d", (long long)cs_info->timestamp, folder_sub_id);
			if(res == -1){
				fprintf(stderr, "ERROR: asprintf failed FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__);
				exit(-1);
			}
			if(folder_sub_id > 100){
				fprintf(stderr, "ERROR: tried %d to make a dir and failed try again later\n",folder_sub_id);
				exit(-1);
			}
			folder_sub_id++;
		}else{
			res = mkdir(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if(res == -1){
				fprintf(stderr, "ERROR: Failed to create dir: (%s). Check application permissions\n",dir);
				exit(-1);
			}
			cs_info->sol_dir = dir;
			break;
		}
	}
}


crct_sim_spec * init_simulator(){
	crct_sim_spec * cs_info;
	cs_info = malloc(sizeof(crct_sim_spec));
	if(cs_info == NULL){
		fprintf(stderr, "Simulator initialization failed %s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
		exit(-1);
	}

	cs_info->commands_count = 0 		;
	cs_info->timestamp 		= time(NULL);
	cs_info->sol_dir		= NULL 		;
	cs_info->cmnd 			= NULL		;
	cs_info->end 			= NULL  	;
	cs_info->info_count 	= 0			;
	cs_info->info_size  	= 0			;
	mksoldir(cs_info);
	init_lists();
	init_hash();
	return cs_info;
}



void attach_plots_to_command(plot_print * print_list_,command_type * command_list){
	plot_print * plot_iter;
	unsigned int line_1;

	if(command_list == NULL){
		return;
	}else if(command_list->next == NULL){
		//one node
		command_type * cmnd_node = command_list;
		//command_list->plot_list = print_list;
		for(plot_iter = print_list_; plot_iter != NULL;){
			if(plot_iter->name == NULL){
				plot_iter = plot_iter->next;
				continue;
			}
			int plot_name_len = strlen(plot_iter->name);
			
			if(plot_name_len > 1 && plot_iter->name[plot_name_len-1] == ')')
				plot_iter->name[plot_name_len-1] = '\0';
			
			
			int name_pos = 1;
			if(plot_name_len > 1 && plot_iter->name[1] == '(')
				name_pos = 2;

			plot_iter->index = get_key(&(plot_iter->name[name_pos]));



			//line_1 = cmnd_node->line;
			if(cmnd_node->next == NULL){
				plot_print * tmp;
				tmp = plot_iter->next;
				plot_iter->next = cmnd_node->plot_list;
				cmnd_node->plot_list = plot_iter;
				plot_iter = tmp;
			}
		}		


	}else{
		//at least one node
		command_type * cmnd_node = command_list;
		for(plot_iter = print_list_; plot_iter != NULL;){
			if(plot_iter->name == NULL){
				plot_iter = plot_iter->next;
				continue;
			}
			int plot_name_len = strlen(plot_iter->name);
			if(plot_name_len > 1 && plot_iter->name[plot_name_len-1] == ')')
				plot_iter->name[plot_name_len-1] = '\0';
			
			if(plot_name_len > 1 && plot_iter->name[1] == '(')
				plot_iter->index = get_key(&(plot_iter->name[2]));
			else
				plot_iter->index = get_key(&(plot_iter->name[1]));


			line_1 = cmnd_node->line;
			if(cmnd_node->next == NULL){
				plot_print * tmp;
				tmp = plot_iter->next;
				plot_iter->next = cmnd_node->plot_list;
				cmnd_node->plot_list = plot_iter;
				plot_iter = tmp;
				
				//cmnd_node->plot_list = plot_iter;
				//break;
			}else{
				//printf("plot_line: %d line_1:%d line_2:%d\n",plot_iter->line,line_1,
				//cmnd_node->next->line );
				if( (line_1 < plot_iter->line && cmnd_node->next->line > plot_iter->line)
					|| line_1 == plot_iter->line){
					plot_print * tmp;
					tmp = plot_iter->next;
					plot_iter->next = cmnd_node->plot_list;
					cmnd_node->plot_list = plot_iter;
					plot_iter = tmp;
				}else{
					cmnd_node = cmnd_node->next;
				}
			} 
		}
	}
}


int main(int argc, char* argv[]){

	struct timespec  tv1, tv2;
	clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);

	//init_lists();
	//init_hash();
	cs_info = init_simulator();


	extern FILE *yyin;
	//fclose(yyin);
	 if ( argc > 1 )
		yyin = fopen( argv[1], "r" );
	 else
		yyin = stdin;

	yylex();
	if(argc > 1) fclose(yyin);
	table_counter--;

	clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);

	printf ("parsing time = %10g seconds\n",
			(double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
			(double) (tv2.tv_sec - tv1.tv_sec));


	plot_print * pp_tmp;


	attach_plots_to_command(plot_print_list,cs_info->cmnd);
	plot_print_list = NULL;



	clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);
	command_type * tmp;
	tmp = cs_info->cmnd;
	if(tmp == NULL){
		printf("Im null\n");
	}
	while(tmp != NULL){
		printf("-%d %s l:%d\n", tmp->type, command_names[tmp->type], tmp->line );
		pp_tmp = tmp->plot_list;
		while(pp_tmp != NULL){
			printf("\tPLOT: node_name: %s] [%d] l:%d\n",pp_tmp->name,pp_tmp->index, pp_tmp->line );
			pp_tmp = pp_tmp->next;
		}

		tmp = tmp->next;
	}

	//return -1;
/*
	printf("plots-------------------------------\n");
	pp_tmp = plot_print_list;
	while(1>2 && pp_tmp != NULL){
		printf("PLOT: node name: %s l:%d\n",pp_tmp->name, pp_tmp->line );
		pp_tmp = pp_tmp->next;
	}
	//exit(-1);
	printf("command_count: %d plots_count: %d\n",cs_info->commands_count, cs_info->plots_count );
	//return 1;
	printf("get_sparse: %d\n",get_sparse() );
	//exit(-1);*/

	if(init_analysis(table_counter, m2) == -1){
		fprintf(stderr, "failed to init mna_construction m2: %d n:%d\n",m2,table_counter );
		return -1;
	}
	//voltage_source 	* tmp_v = vltg_list;
	//printf("is_ac: %d\n",tmp_v->is_ac );
//	printf("0x%x transient_analysis\n",tmp_v->transient_info );

	//enables floating point exceptions
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	//fix_plot_print_list();
	printf("mna_construction\n");

	mna_construction();
//	printf_matrix();
	//free_list();
	clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);

	printf ("construction time = %10g seconds\n",
			(double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
			(double) (tv2.tv_sec - tv1.tv_sec));

	clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);

	//printf("get_freq_type: %d\n",get_freq_type() );
	//exit(-1);

	int command_id = 1;
	tmp = cs_info->cmnd;
	gsl_permutation * p = NULL;
	if(tmp == NULL){
		//DC
		//p = 
		mna_solver();
	}else{
		while(tmp != NULL){
			//reverse previously applied permutation
			if(p != NULL && tmp->type != DC_SWEEP){
				for(int i=0;i<p->size;i++){
					if(i != p->data[i] && i < p->data[i]){
						unsigned int line_1 = i;
						unsigned int line_2 = p->data[i];
						for(int j=0;j<alpha_matrix_size;j++){
							double swap_tmp;
							swap_tmp = alpha_matrix[line_1*alpha_matrix_size+j];
							alpha_matrix[line_1*alpha_matrix_size+j] 
							= alpha_matrix[line_2*alpha_matrix_size+j];
							alpha_matrix[line_2*alpha_matrix_size+j] = swap_tmp;  
						}
					}
				}
				gsl_permutation_free(p); p = NULL;
			}

			printf("RUNNING ANALYSIS %d %s l:%d\n", tmp->type, command_names[tmp->type], tmp->line );
			if(tmp->type == DC_SWEEP){
				printf("- DC_SWEEP\n");
				//printf("============================================================\n");
				//printf_matrix();
					
				if(p == NULL)
					p = mna_solver();
	
				//printf_matrix();
				//for(int i=0;i<p->size;i++){
				//	printf("%d %d\n",i,p->data[i] );
				//}
				//printf_matrix();
				//printf("============================================================\n");
				//size_t size;
  				//size_t * data;
				//tmp->plot_list = plot_print_list;
				dc_sweep_analysis(p,tmp,command_id);
			}else if(tmp->type == TRANSIENT){
				printf("TRANSIENT\n");
				p = mna_solver();
				trnsient_analysis(command_id,tmp);

			}else if(tmp->type == AC){
				printf("AC\n");
				run_freq_analysis(command_id,tmp);
			}
			command_id++;
			tmp = tmp->next;
		}
	}
/*
	printf_matrix();
	if(!get_freq_type()){
		p = mna_solver();
	}
	if(get_tran())
		trnsient_analysis();

	if(get_freq_type())
		run_freq_analysis();
	//XXX
	if(1>2)
		dc_sweep_analysis(p);

*/	
	clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);
	
	printf ("solve time = %10g seconds\n",
			(double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
			(double) (tv2.tv_sec - tv1.tv_sec));

	printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf_matrix();

	/*
	printf("===============================================\n");
	print_resistor();
	print_capacitor();
	print_coil();
	print_crnt();
	print_vltg();
	print_diode();
	print_transistor_bjt();
	print_transistor_mos();
	*/
	printf("===============================================\n");
	printf("get_analysis_type: %d %s %d\n",get_analysis_type(), get_analysis_type_text(), get_sparse());
	return 0;
}
