#ifndef ELEMENT_STRUCTS_H
#define ELEMENT_STRUCTS_H

#include <stdint.h>

typedef struct resistor {
	char *name;
	uint32_t positive_node;
	uint32_t negative_node;
	double value;
	struct resistor *next;
} resistor;

typedef struct capacitor {
	char *name;
	uint32_t positive_node;
	uint32_t negative_node;
	double value;
	struct capacitor *next;
} capacitor;

typedef struct coil {
	char *name;
	uint32_t positive_node;
	uint32_t negative_node;
	double value;
	struct coil *next;
} coil;

typedef struct voltage_source {
	char *name;
	uint32_t positive_node;
	uint32_t negative_node;
	double value;
	short is_ac;
	double ac_phase;
	double ac_mag;
	int b_vector_pos;
	struct voltage_source *next;
	struct transient_spec * transient_info;
} voltage_source;



/*
XXX not use yet
typedef struct  {
	double ac_phase;
	double ac_mag;
} ac_voltage_spec;
*/

typedef struct transient_spec {
	unsigned char type;
	union{
		struct exp{
			//exp
			double i1;
			double i2;
			double td1; 
			double tc1;
			double td2;
			double tc2;
		} exp_info;
		struct sin {
			//sin
			double i1;
			double ia;
			double fr;
			double td;
			double df;
			double ph;	
		} sin_info;
		struct pulse {
			//pulse
			double i1;
			double i2;
			double td;
			double tr;
			double tf;
			double pw;
			double per;
		} pulse_info;
		struct pwl{
			double *t1;
			double *i1;
			unsigned int vector_size;
		} pwl_info;
	}spec;
} transient_spec;

typedef struct current_source {
	char *name;
	uint32_t positive_node;
	uint32_t negative_node;
	double value;
	double ac_phase;
	double ac_mag;
	short is_ac;
	struct current_source *next;
	struct transient_spec * transient_info;
} current_source;

typedef struct diode {
	char *name;	
	uint32_t positive_node;
	uint32_t negative_node;
	char *model_name;
	double area;
	struct diode *next;
} diode;

typedef struct transistor_mos {
	char *name;
	uint32_t d_node;
	uint32_t g_node;
	uint32_t s_node;
	uint32_t b_node;
	char *model_name;
	double L_value;
	double W_value;
	struct transistor_mos *next;
} transistor_mos;

typedef struct transistor_bjt{
	char *name;
	uint32_t c_node;
	uint32_t b_node;
	uint32_t e_node;
	char *model_name;
	double area;
	struct transistor_bjt *next;
	//area
} transistor_bjt;


typedef struct dc_options{
	char * node_name;
	double from_voltage;
	double to_voltage;
	double step;
	struct dc_options *next;
} dc_options;

typedef struct plot_print {
	char * name;
	int index;
	unsigned int line;
	unsigned short is_db : 1;
	struct plot_print *next;
} plot_print;


double tran_time_step, tran_fin_time;
int derivative_approximation_type;
dc_options 		* dc_option_list;
plot_print 		* plot_print_list;
transient_spec 	* tran_spec_info_temp;
resistor		* resistor_list;
capacitor 		* capacitor_list;
coil 			* coil_list;
voltage_source 	* vltg_list;
current_source 	* crnt_list;
diode 			* diode_list;
transistor_mos 	* tmos_list;
transistor_bjt 	* tbjt_list;
unsigned int transient_analysis;
unsigned int m2,n,m1,m;
unsigned int analysis_type;
double itol;
unsigned int sparse_analysis;
unsigned int non_zero;
unsigned int c_bar_non_zero;
//-------------------------------------------
// FINAL STEP OPTIONS
//-------------------------------------------
double ac_phase;
double ac_mag;
int ac_sweep;
int ac_steps;
double ac_start_freq;
double ac_end_freq;
unsigned int freq_analysis;
#endif
