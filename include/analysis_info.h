#ifndef ANALYSIS_INFO_H
#define ANALYSIS_INFO_H
#include "element_structs.h"
#include <time.h>

typedef struct  command_type {
	enum { DC, DC_SWEEP, TRANSIENT, AC} type;
	enum { LU, CHOLESKY, CG, BICG } solver_type;

	unsigned int is_trapezoidal 	: 1;
	unsigned int is_backeuler 		: 1;
	unsigned int output_in_graph 	: 1;
	unsigned int output_in_list 	: 1;

	union{
		//-----------------------------------
		// TRANSIENT OPTIONS
		//-----------------------------------
		struct t_info{
			double tran_time_step;
			double tran_fin_time;
		}t_info;
		//-----------------------------------
		// DC_SWEEP OPTIONS
		//-----------------------------------
		/*struct d_info{
			char * dc_swp_node_name;
			double dc_swp_from_voltage;
			double dc_swp_to_voltage;
			double dc_swp_step;

		}d_info;*/
		dc_options * d_info;

		//-----------------------------------
		// AC_ANALYSIS OPTIONS
		//-----------------------------------
		struct a_info{
			//double ac_phase;
			//double ac_mag;
			int ac_sweep;
			int ac_steps;
			double ac_start_freq;
			double ac_end_freq;
		}a_info;
	}d;
	//struct dc_options *next;
	//} dc_options;
	char * raw_cmnd;
	unsigned int line;
	plot_print * plot_list;
	struct command_type * next;
}command_type ;



typedef struct {
	time_t timestamp;
	char * sol_dir;
	command_type * cmnd;
	command_type * end;
	unsigned int commands_count;
	unsigned int plots_count;
	unsigned int info_count;
	unsigned int info_size;
	unsigned int is_sparse 			: 1;
	//char * type_names[];
	//char * solver_names[];
} crct_sim_spec;


crct_sim_spec * cs_info;
#endif