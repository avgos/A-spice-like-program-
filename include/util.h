#ifndef UTIL_H
#define UTIL_H
#include "element_structs.h"

void init_lists();
int insert(char *stack[], int qlen, int ignore_line, int line_num);
void print_resistor();
void print_capacitor();
void print_coil();
void print_crnt();
void print_vltg();
void print_diode();
void print_transistor_bjt();
void print_transistor_mos();
void set_analysis_type(int type);
void set_itol(double given_itol);
char * get_analysis_type_text(void);
unsigned int get_analysis_type(void);
//plot_print * insert_plot_print(char *yytext, int line);
plot_print * insert_plot_print(char *yytext, int line, int is_db);
void fix_plot_print_list(void);
void set_sparse(void);
int get_sparse(void);
void free_list(void);
//void set_tran(void);
void set_tran(char * yytext, int line);
int get_tran(void);
int get_approx_type(void);
void set_approx_type(int type);
void parse_ac_mag_phase(char * yytext);
void parse_ac_options(char * yytext, int line);


void set_freq_type(void);
int get_freq_type(void);

void parse_transient_spec(char * yytext, unsigned char type);
#endif
