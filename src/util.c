#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "element_structs.h"
#include "node_to_int.h"
#include "analysis_info.h"


#define ALLOC_CHECK(a) if(a == NULL) { \
		fprintf(stderr, "[ERR] MALLOC FAILED FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__ );\
       	exit(-1);\
    }


void enqueue(command_type *p){
	if(p == NULL) return;
	if(cs_info == NULL) return;
	cs_info->commands_count++;
	p->plot_list = NULL;
	if(cs_info->end == NULL){
		p->next = NULL;
		cs_info->cmnd = p;
		cs_info->end  = p;
	}else{
		//if(p->type == DC_SWEEP){
		//	p->next = cs_info->cmnd;
		//	cs_info->cmnd = p;
		//}else{
		cs_info->end->next = p;
		p->next = NULL;
		cs_info->end = p;
		//}
	}
}

void init_lists(){
	resistor_list 	= NULL;
	capacitor_list 	= NULL;
	coil_list 		= NULL;
	vltg_list 		= NULL;
	crnt_list		= NULL;
	diode_list 		= NULL;
	tmos_list 		= NULL;
	tbjt_list 		= NULL;
	plot_print_list = NULL;
	dc_option_list  = NULL;
	tran_spec_info_temp = NULL;
	m2 = 0;
	m1 = 0; 
	m  = 0;
	n  = 0; 
	analysis_type = 0;
	itol=0.001;
	sparse_analysis = 0;
	non_zero = 0;
	c_bar_non_zero = 0;
	transient_analysis = 0;
	derivative_approximation_type = 1; //default trapezoidal
	freq_analysis = 0;
}

int get_sparse(void){
	return sparse_analysis;
}
void set_sparse(void){
	sparse_analysis = 1;
}
void set_itol(double given_itol){
	itol=given_itol;
}
void set_approx_type(int type){
	derivative_approximation_type = type;
}
int get_approx_type(void){
	return derivative_approximation_type;
}
void set_analysis_type(int type){
	analysis_type = type;
}
void set_freq_type(void){
	freq_analysis = 1;
}
int get_freq_type(void){
	return freq_analysis;
}



void set_tran(char * yytext, int line){

	char *t1;
	t1 = &yytext[5];
	while(*t1 == ' ' && *t1 != '\0') t1++;
	tran_time_step = atof(t1);
	//printf("A: %f\n", tran_time_step);
	while(*t1 != ' ' && *t1 != '\0') t1++;
	while(*t1 == ' ' && *t1 != '\0') t1++;
	tran_fin_time = atof(t1);
	//printf("B: %f\n", tran_fin_time); 

	command_type * cmnd_type;
	cmnd_type = malloc(sizeof(command_type));
	if(!cmnd_type){
		fprintf(stderr, "Memory allocation failed %s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
		exit(-1);
	}
	cmnd_type->line = line;
	char * raw_cmnd;
	int res;
	res = asprintf(&raw_cmnd,".TRANSIENT %.5e %.5e",tran_time_step, tran_fin_time);
	if(res == -1 || raw_cmnd == NULL){
		fprintf(stderr, "asprintf failed\n");
		exit(-1);
	}
	cmnd_type->raw_cmnd = raw_cmnd; 
	cmnd_type->type = TRANSIENT;
	cmnd_type->d.t_info.tran_time_step = tran_time_step;
	cmnd_type->d.t_info.tran_fin_time  = tran_fin_time;
	//enqueue
	enqueue(cmnd_type);
	//cmnd_type->next = cs_info->cmnd;
	//cs_info->cmnd = cmnd_type;
	transient_analysis = 1;
}


int get_tran(void){
	return transient_analysis;
}

char * get_analysis_type_text(void){
	char * names[] = {"LU","Cholesky","CG","BICG", "trapezoidal", "backward-euler"};
	return names[analysis_type];
}

unsigned int get_analysis_type(void){
	return analysis_type;
}

plot_print * insert_plot_print(char *yytext, int line, int is_db){
	static long int prev_line = -1;
	static plot_print * plot_end = NULL;
	printf("[insert_plot_print]: %s 0x%x\n",yytext, plot_end );
	if(prev_line == -1){
		cs_info->plots_count = 1;
		prev_line = line;
	}else{
		if(line != prev_line){
			prev_line = line;
			cs_info->plots_count += 1;
		}
	}
	//printf("insert_plot_print: %s\n", yytext);
	//pos == 0 -> .DC
	//pos >  0 -> val
	plot_print * temp;
	//printf("[insert_plot_print] %s| ",yytext );
	temp = malloc(sizeof(plot_print));
	temp->is_db = is_db;
	char * t;
	t = yytext;
	while((*t == ' ' || *t == '\n') && *t != '\0') t++;
	if(*t != '\0'){
		temp->name = strdup(t);
		t = temp->name;
		int i = strlen(temp->name) - 1;
		t = t + i;
		char * t_initial = t;
		while(*t == '\n' || *t == '\r' || *t == ' ') t--;
		if(t != t_initial)
			*(t+1) = '\0';
		if(is_db){
			int name_len = strlen(temp->name);
			if(temp->name[name_len-1] == ')'){
				temp->name[name_len-1-3] = ')';
				temp->name[name_len-1-2] = '\0';
			}else{
				temp->name[name_len-1-2] = '\0';
			}
		}
	}
	else
		temp->name = NULL;
	temp->line = line;
	//temp->next = plot_print_list;
	//plot_print_list = temp;
	
	if(plot_end == NULL){
		temp->next = NULL;
		plot_print_list = temp;
		plot_end = temp;
	}else{
		plot_end->next = temp;
		temp->next = NULL;
		plot_end = temp;
	}

	return NULL;
}



void parse_ac_options(char * yytext, int line){
	//Solver options
	set_freq_type();
	printf("parse_ac_options: %s\n",yytext);
	char * t1;
	t1 = &yytext[3];
	printf("parse_ac_options: %s\n",t1);
	while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
	int cmp = strncasecmp("lin", t1, 3);
	if(cmp == 0){
		//exoume lin
		ac_sweep = 0;
	}else{
		//exoume log , de mporei na einai alliws -> o parser kobei ta panta 
		//ektos apo lin/log ;)
		ac_sweep = 1;
	}
	printf("sweep: %d\n", ac_sweep);
	printf("cmp : %d\n", cmp);
	while(*t1 != ' ' && *t1 != '\t' && *t1 != '\0')	t1++;
	ac_steps = atoi(t1);
	printf("steps : %d\n", ac_steps);	
	t1++;
	while(*t1 != ' ' && *t1 != '\t' && *t1 != '\0')	t1++;
	while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;	
	ac_start_freq = atof(t1);
	printf("start_freq: %f\n", ac_start_freq);
	while(*t1 != ' ' && *t1 != '\t' && *t1 != '\0')	t1++;
	while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;	
	if(*t1 == '\0'){
		fprintf(stderr, "parsing ac options failed\nFILE: %s LINE: %d\n",__FILE__, __LINE__);
		exit(-1);
	}

	ac_end_freq = atof(t1);
	printf("end_freq: %f\n", ac_end_freq);


	command_type * cmnd_type;
	cmnd_type = malloc(sizeof(command_type));
	if(!cmnd_type){
		fprintf(stderr, "Memory allocation failed %s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
		exit(-1);
	}

	cmnd_type->line = line;
	cmnd_type->type = AC;
	char * raw_cmnd;
	asprintf(&raw_cmnd,".AC %s %d %lf %lf", ac_sweep==0?"lin":"log", ac_steps, ac_start_freq, ac_end_freq);
	cmnd_type->raw_cmnd = raw_cmnd; //trim_stdup(yytext);
	cmnd_type->d.a_info.ac_sweep = ac_sweep;
	cmnd_type->d.a_info.ac_steps = ac_steps;
	cmnd_type->d.a_info.ac_start_freq = ac_start_freq;
	cmnd_type->d.a_info.ac_end_freq = ac_end_freq;

	//enqueue
	//cmnd_type->next = cs_info->cmnd;
	//cs_info->cmnd = cmnd_type;
	enqueue(cmnd_type);
	//cmnd_type = NULL;
}


void parse_ac_mag_phase(char * yytext){
	// Voltage / Current Source Option
	set_freq_type();
	printf("ac_mag: %s\n",yytext );
	char *t1;
	t1 = &yytext[2];
	printf("|%s|\n", t1);
	while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
	ac_mag = atof(t1);
	printf("mag:%f\n", ac_mag);
	while(*t1 != ' ' && *t1 != '\0')t1++;
	while((*t1 == ' '  || *t1 == '\t') && *t1 != '\0')t1++;
	if(*t1 == '\0'){
		fprintf(stderr, "parsing ac mag phase failed\nFILE: %s LINE: %d\n",__FILE__, __LINE__);
		exit(-1);
	}
	ac_phase = atof(t1);
	printf("phase: %f\n",ac_phase);
}

void parse_transient_spec(char * yytext, unsigned char type){
	transient_spec * tran_spec_info;
	tran_spec_info = (transient_spec *)malloc(sizeof(transient_spec));
	tran_spec_info->type = type;
	char *t1;
	t1 = yytext;
	while(*t1 == ' ' || *t1 == '(') t1++;
	if(type == 'e'){
		//exp
		tran_spec_info->spec.exp_info.i1 = atof(t1);
		//printf("yytext: %f - %s\n",atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.exp_info.i2 = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.exp_info.td1 = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.exp_info.tc1 = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.exp_info.td2 = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.exp_info.tc2 = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);

		printf("i1:%f i2:%f td1:%f tc1:%f td2:%f tc2:%f\n", tran_spec_info->spec.exp_info.i1,
			tran_spec_info->spec.exp_info.i2, tran_spec_info->spec.exp_info.td1, tran_spec_info->spec.exp_info.tc1, 
			tran_spec_info->spec.exp_info.td2, tran_spec_info->spec.exp_info.tc2);
	}else if(type == 's'){
		//sin
		tran_spec_info->spec.sin_info.i1 = atof(t1);
		//printf("yytext: %f - %s\n",atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.sin_info.ia = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.sin_info.fr = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.sin_info.td = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.sin_info.df = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.sin_info.ph = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
	}else if(type == 'p'){
		//pulse
		tran_spec_info->spec.pulse_info.i1 = atof(t1);
		////printf("yytext: %f - %s\n",atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.i2 = atof(t1);
		////printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.td = atof(t1);
		////printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.tr = atof(t1);
		////printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.tf = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.pw = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
		while(*t1 != ' ' && *t1 != '\0') t1++;
		while((*t1 == ' ' || *t1 == '\t') && *t1 != '\0') t1++;
		tran_spec_info->spec.pulse_info.per = atof(t1);
		//printf("yytext: %f - %s\n", atof(t1),t1);
	}else if(type == 'w'){
		//pwl
		printf("pwl: %s\n", t1);
		char * t2;
		t2 = t1;
		int vector_size = 0;
		while(*t2 != '\0'){
			vector_size += (*t2 == ')');
			t2++;
		}
		printf("vector_size: %d\n", vector_size);
		if(vector_size == 0){
			fprintf(stderr, "[ERROR] pwl source without arguments \n");
			exit(-1);
		}
		tran_spec_info->spec.pwl_info.vector_size = vector_size;
		tran_spec_info->spec.pwl_info.t1 = malloc(sizeof(double)*vector_size);
		tran_spec_info->spec.pwl_info.i1 = malloc(sizeof(double)*vector_size);
		t2 = t1;
		int i;
		for(i=0;i<vector_size;i++){
			while((*t2 == ' ' || *t2 == '\t' || *t2 == '(') && *t2 == '\0') t2++;
			//printf("A:t2: %s\n",t2);
			tran_spec_info->spec.pwl_info.t1[i] = atof(t2);
			while(*t2 != ' ' && *t2 != '\0') t2++;
			while((*t2 == ' ' || *t2 == '\t') && *t2 != '\0') t2++;
			//printf("B:t2: %s\n",t2);
			tran_spec_info->spec.pwl_info.i1[i] = atof(t2);
			while(*t2 != ')' && *t2 != '\0') t2++;
			while((*t2 == ')' || *t2 == ' ' || *t2 == '(') && *t2 != '\0') t2++;
		}
		for(i=0;i<vector_size;i++){
			printf("i:%d t1:%f i1:%f\n", i, tran_spec_info->spec.pwl_info.t1[i],tran_spec_info->spec.pwl_info.i1[i]);
		}
	}
	printf("------------------------------------------------\n");
	//vltg_list->transient_info = tran_spec_info;
	tran_spec_info_temp  = tran_spec_info;
	//free();
}

void fix_plot_print_list(void){
	plot_print * temp;
	temp = plot_print_list;
	while(temp != NULL){
		temp->name[strlen(temp->name)-1] = '\0';
		//printf("EDW: |%s|->[%d]\n", &temp->name[2],get_key(&(temp->name[2])) );
		temp->index = get_key(&(temp->name[2]));
		temp = temp->next;
	}
}



//dc_options * insert_dc_option_node(dc_options * head, char *stack[], int qlen){
void insert_dc_option_node(char *stack[], int qlen,int line_no){
	static unsigned int line_num = -1;
	static command_type * last_node = NULL;
	int i;
	for(i=0;i<qlen;i++){
		printf("DC_OPTION: %s\n",stack[i]);
	}
	if(qlen != 5){
		fprintf(stderr, "qlen: %d\n",qlen);
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4 insert_dc_option_node\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return;
	}

	dc_options * temp = malloc(sizeof(dc_options));
	ALLOC_CHECK(temp)
	temp->node_name = strdup(stack[1]);  stack[1] = NULL;
	temp->from_voltage = atof(stack[2]); stack[2] = NULL;
	temp->to_voltage = atof(stack[3]);   stack[3] = NULL;
	temp->step = atof(stack[4]);  		 stack[4] = NULL;
	

	if(line_num != line_no){
		line_num = line_no;
		temp->next = NULL;
		command_type * cmnd_type;
		cmnd_type = malloc(sizeof(command_type));
		if(!cmnd_type){
			fprintf(stderr, "Memory allocation failed %s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
			exit(-1);
		}
		char * raw_cmnd;
		asprintf(&raw_cmnd,".DC %s %lf %lf %lf", temp->node_name, temp->from_voltage, temp->to_voltage, temp->step);
		cmnd_type->raw_cmnd = raw_cmnd; //strdup(yytext);
		cmnd_type->line = line_no;
		cmnd_type->type = DC_SWEEP;
		cmnd_type->plot_list = NULL;
		cmnd_type->d.d_info = temp;
		last_node = cmnd_type;
		enqueue(cmnd_type);
	}else{
		temp->next = last_node->d.d_info;
		last_node->d.d_info = temp;
	}
}


resistor * insert_resistor(resistor * head, char *stack[], int qlen){
	if(qlen != 4){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	m1++;
	resistor *temp = malloc(sizeof(resistor));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	//printf("RESISTOR: qlen=%d key1=[%s/%d] key2=[%s/%d]\n",qlen,stack[1],add_node(stack[1]),stack[2],add_node(stack[2]) );
	temp->value = atof(stack[3]); stack[3] = NULL;
	temp->next = head;
	if(temp->positive_node == 0 || temp->negative_node == 0){
		// + 1 nz
		non_zero += 1;
	}else{
		// +4 nz
		non_zero += 4;
	}
	return temp;
}

capacitor * insert_capacitor(capacitor * head, char *stack[], int qlen){
	if(qlen != 4){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	m1++;
	capacitor *temp = malloc(sizeof(capacitor));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	temp->value = atof(stack[3]); stack[3] = NULL;

	if(temp->positive_node == 0 || temp->negative_node == 0){
		// + 1 nz
		c_bar_non_zero += 1;
	}else{
		// +4 nz
		c_bar_non_zero += 4;
	}
	temp->next = head;
	return temp;
}

coil * insert_coil(coil * head, char *stack[], int qlen){
	if(qlen != 4){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	m2++;
	coil *temp = malloc(sizeof(coil));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	temp->value = atof(stack[3]); stack[3] = NULL;
	temp->next = head;
	if(temp->positive_node != temp->negative_node){
		c_bar_non_zero += 1;
	}
	//XXX ?????
	if(temp->positive_node == 0 || temp->negative_node == 0){
		// +2 nz
		non_zero += 2;
	}else{
		// +4 nz
		non_zero += 4;
	}	
	return temp;
}


diode * insert_diode(diode * head, char *stack[], int qlen){
	if(qlen != 5 && qlen != 4){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;		
	}
	diode *temp = malloc(sizeof(diode));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	temp->model_name = stack[3]; stack[3] = NULL;
	if(qlen == 5){
		temp->area = (stack[4][0] == '[') ? atof(&stack[4][1]) : atof(stack[4]); stack[4] = NULL;
	}
	temp->next = head;
	return temp;
}


transistor_bjt * insert_bjt(transistor_bjt * head, char * stack[], int qlen){
	if(qlen != 6 && qlen != 5){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;		
	}
	transistor_bjt *temp = malloc(sizeof(transistor_bjt));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->c_node = add_node(stack[1]); stack[1] = NULL;
	temp->b_node = add_node(stack[2]); stack[2] = NULL;
	temp->e_node = add_node(stack[3]); stack[3] = NULL;
	temp->model_name = stack[4]; stack[4] = NULL;
	if(qlen == 6){
		temp->area = (stack[5][0] == '[') ? atof(&stack[5][1]) : atof(stack[5]); stack[5] = NULL;
	}
	temp->next = head;
	return temp;
}


voltage_source * insert_voltage_source(voltage_source * head, char *stack[], int qlen){
	if(qlen != 4 && (get_freq_type() != 1 && qlen == 7)){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	m2++;
	voltage_source *temp = malloc(sizeof(voltage_source));
	ALLOC_CHECK(temp)
	if(qlen > 0){
		temp->name = strdup(stack[0]); stack[0] = NULL;
	}
	if(qlen > 1){
		temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	}else{
		temp->positive_node = 0;
	}
	
	if(qlen > 2){
		temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	}else{
		temp->negative_node = 0;
	}
	if(qlen > 3){
		temp->value = atof(stack[3]); stack[3] = NULL;
	}else{
		temp->value = 0.0;
	}
	temp->next = head;
	temp->is_ac = 0;
	temp->transient_info = tran_spec_info_temp;
	tran_spec_info_temp = NULL;
	if(get_freq_type() == 1){
		//printf("HELLO\n");
		temp->ac_phase = ac_phase;
		temp->ac_mag = ac_mag;
		temp->is_ac = 1;
	}else{
		//printf("NOT HELLO\n");
	}
	if(temp->positive_node == 0 || temp->negative_node == 0){
		// +2 nz
		non_zero += 2;
	}else{
		// +4 nz
		non_zero += 4;
	}
	return temp;
}

current_source * insert_current_source(current_source * head, char *stack[], int qlen){
	if(qlen != 4){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	m1++;
	current_source *temp = malloc(sizeof(current_source));
	ALLOC_CHECK(temp)
	temp->name = stack[0]; stack[0] = NULL;
	temp->positive_node = add_node(stack[1]); stack[1] = NULL;
	temp->negative_node = add_node(stack[2]); stack[2] = NULL;
	temp->value = atof(stack[3]); stack[3] = NULL;
	temp->next = head;
	temp->is_ac = 0;
	temp->transient_info = tran_spec_info_temp;
	tran_spec_info_temp = NULL;
	if(get_freq_type() == 1){
		temp->ac_phase = ac_phase;
		temp->ac_mag = ac_mag;
		temp->is_ac = 1;
	}
	return temp;
}

transistor_mos * insert_transistor_mos(transistor_mos * head, char *stack[], int qlen){
	if(qlen != 8){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d qlen(%d) != 4\n",__FILE__,__FUNCTION__,__LINE__,qlen);
		return head;
	}
	transistor_mos * temp = malloc(sizeof(transistor_mos));
	ALLOC_CHECK(temp)
	temp->name   = stack[0]; 			stack[0] = NULL;
	temp->d_node = add_node(stack[1]); 	stack[1] = NULL;
	temp->s_node = add_node(stack[2]); 	stack[2] = NULL;
	temp->g_node = add_node(stack[3]); 	stack[3] = NULL;
	temp->b_node = add_node(stack[4]); 	stack[4] = NULL;
	temp->model_name = stack[5]; 		stack[5] = NULL;
	temp->L_value = atof(&stack[6][2]); 	stack[6] = NULL;
	temp->W_value = atof(&stack[7][2]); 	stack[7] = NULL;
	temp->next = head;
	return temp;
}

void print_resistor(){

	resistor * temp;
	temp = resistor_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, value:%lf, name:%s\n",temp->positive_node, temp->negative_node,temp->value,temp->name);
		temp = temp->next;
	}
}

void print_capacitor(){

	capacitor * temp;
	temp = capacitor_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, value:%lf, name:%s\n",temp->positive_node, temp->negative_node,temp->value,temp->name);
		temp = temp->next;
	}
}

void print_coil(){

	coil * temp;
	temp = coil_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, value:%lf, name:%s\n",temp->positive_node, temp->negative_node,temp->value,temp->name);
		temp = temp->next;
	}
}

void print_vltg(){

	voltage_source * temp;
	temp = vltg_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, value:%lf, name:%s\n",temp->positive_node, temp->negative_node,temp->value,temp->name);
		temp = temp->next;
	}
}

void print_crnt(){

	current_source * temp;
	temp = crnt_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, value:%lf, name:%s\n",temp->positive_node, temp->negative_node,temp->value,temp->name);
		temp = temp->next;
	}
}

void print_diode(){

	diode * temp;
	temp = diode_list;
	while(temp != NULL){
		//printf("pnode: %d, nnode:%d, model name:%s, name:%s area:%lf\n",temp->positive_node, temp->negative_node,temp->model_name,temp->name,temp->area);
		temp = temp->next;
	}
}

void print_transistor_mos(){

	transistor_mos * temp;
	temp = tmos_list;
	while(temp != NULL){
		//printf("dnode: %d, gnode:%d, snode:%d, bnode:%d, model name:%s, name:%s, Lvalue:%lf, Wvalue:%lf\n",temp->d_node, temp->g_node,temp->s_node,temp->b_node,temp->model_name,temp->name,temp->L_value,temp->W_value);
		temp = temp->next;
	}
}

//model_name??
void print_transistor_bjt(){

	transistor_bjt * temp;
	temp = tbjt_list;
	while(temp != NULL){
		//printf("cnode: %d, bnode:%d, enode:%d, model name:%s, name:%s area:%lf\n",temp->c_node, temp->b_node,temp->e_node,temp->model_name,temp->name, temp->area);
		temp = temp->next;
	}
}


void free_list(void){
	void * temp, * curr;
	//temp = (dc_options *)dc_option_list;
	//while(temp != NULL){curr = temp;temp = ((dc_options *)temp)->next;free(curr);}
	//temp = (plot_print *) plot_print_list;
	//while(temp != NULL){curr = temp;temp = ((plot_print *)temp)->next;free(curr);}
	 temp = (resistor *) resistor_list; 
	while(temp != NULL){curr = temp; temp = ((resistor *)temp)->next; free(((resistor *)curr)->name); free(curr);}
	 temp = (capacitor *) capacitor_list; 
	while(temp != NULL){curr = temp;temp = ((capacitor *)temp)->next; free(((capacitor *)curr)->name); free(curr);}
	temp = (coil *) coil_list; 
	while(temp != NULL){curr = temp;temp = ((coil *)temp)->next; free(((coil *)curr)->name); free(curr);}
	temp = (voltage_source *) vltg_list; 
	while(temp != NULL){curr = temp;temp = ((voltage_source *)temp)->next; free(((voltage_source *)curr)->name); free(curr);}
	temp = (current_source *) crnt_list; 
	while(temp != NULL){curr = temp;temp = ((current_source *)temp)->next; free(((current_source *)curr)->name); free(curr);}
	temp = (diode *) diode_list; 
	while(temp != NULL){curr = temp;temp = ((diode *)temp)->next; free(((diode *)curr)->name); free(((diode *)curr)->model_name); free(curr);}
	temp = (transistor_mos *) tmos_list; 
	while(temp != NULL){curr = temp;temp = ((transistor_mos *)temp)->next; free(((transistor_mos *)curr)->name); free(((transistor_mos *)curr)->model_name); free(curr);}
	temp = (transistor_bjt *) tbjt_list; 
	while(temp != NULL){curr = temp;temp = ((transistor_bjt *)temp)->next; free(((transistor_bjt *)curr)->name); free(((transistor_bjt *)curr)->model_name); free(curr);}

}

int insert(char *stack[], int stack_pos, int ignore_line, int line_num){
	
	//printf("ignore_line:%d stack_pos:%d \n",ignore_line,stack_pos );
	//int i;
	//for(i=0;i<stack_pos;i++)
	//	printf("%d: %s\t",i,stack[i]);
	//printf("\n");
	if(ignore_line == 1){
		fprintf(stderr, "[ERR] FILE:%s FUNCTION: %s LINE: %d %d line is ignored (ignore_line = 1)\n",__FILE__,__FUNCTION__,__LINE__,line_num);
		return 1;
	}
	switch(stack[0][0]){
		case 'r':
		case 'R': resistor_list = insert_resistor(resistor_list,stack,stack_pos); 
				break;
		case 'c':
		case 'C': capacitor_list = insert_capacitor(capacitor_list,stack,stack_pos); 
				break;
		case 'l':
		case 'L': coil_list = insert_coil(coil_list,stack,stack_pos); 
				break;
		case 'm':
		case 'M': tmos_list = insert_transistor_mos(tmos_list,stack,stack_pos); 
				break;
		case 'v':
		case 'V': vltg_list = insert_voltage_source(vltg_list,stack,stack_pos); 
				break;
		case 'i':
		case 'I': crnt_list = insert_current_source(crnt_list,stack,stack_pos); 
				break;
		case 'd':
		case 'D': diode_list = insert_diode(diode_list, stack, stack_pos); 
				break;
		case 'q':
		case 'Q': tbjt_list = insert_bjt(tbjt_list, stack, stack_pos); 
				break;
		//case 'K': dc_option_list = insert_dc_option_node(dc_option_list, stack, stack_pos); break;
		case 'K': insert_dc_option_node(stack, stack_pos,line_num); 
				break;
		default: break;
	}
	/*printf("------------------------------------\n");
	printf("List pointers\n");
	printf("0x%x\n", (unsigned int)resistor_list);
	printf("0x%x\n", (unsigned int)capacitor_list);
	printf("0x%x\n", (unsigned int)coil_list);
	printf("0x%x\n", (unsigned int)tmos_list);
	printf("0x%x\n", (unsigned int)vltg_list);
	printf("0x%x\n", (unsigned int)crnt_list);
	printf("------------------------------------\n");
	*/
	return 1;
}

