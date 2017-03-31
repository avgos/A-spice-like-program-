%option noyywrap default yylineno case-insensitive
%option nounput
%option noinput
%x RLC_elem
%x cur_vol_source
%x tr_mos
%x tr_bjt
%x diode
%x analysis_option
%x dc_option
%x plot_option
%x itol_option
%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node_to_int.h"
#include "util.h"
int pos = 0;
char * stack[8];
int stack_pos = 0;
int ignore_line = 0;
int call = 0;
int exist_tran_spec = 0;
#define ALLOC_CHECK(a) if(a == NULL) { \
		fprintf(stderr, "[ERR] MALLOC FAILED FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__ );\
       	exit(-1);\
    }

#define SYNTAX_CHECK {\
		fprintf(stderr,"[ERR] WRONG SYNTAX LINE: %d parser_line: %d\n",yylineno,__LINE__);\
	}
%}


%%
 

".TRAN"(([ \t]*((([0-9]*".")?[0-9]*e("-"|"+")?[1-9]+[0-9]*)|((0|[1-9][0-9]*)("."[0-9]*)?))))+[ \t]*(\n|\r\n)? { 
	set_tran(yytext,yylineno);
}

".AC"([ \t]+)("lin"|"log"|"dec")([ \t]+)([0-9]+)(([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)){2}) {
	parse_ac_options(yytext,yylineno);
}

".DC" { BEGIN(dc_option); stack[stack_pos++] = strdup("K"); }
<dc_option>[\t ]*
<dc_option>[VI][0-9_A-Z]+ {  stack[stack_pos++] = strdup(yytext);}
<dc_option>((0|[1-9][0-9]*)"."[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*)|[0-9]* {  stack[stack_pos++] = strdup(yytext); }
<dc_option>\n|\r\n { BEGIN(INITIAL);  insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<dc_option><<EOF>> { BEGIN(INITIAL);  insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<dc_option>. {}



".PLOT"|".PRINT" { BEGIN(plot_option); call=0;}
<plot_option>[\t ]*
<plot_option>"tran" { }
<plot_option>[V]"("([0-9_A-Za-z]+)")" { insert_plot_print(yytext,yylineno,0);call++; }
<plot_option>[V]([0-9_A-Za-z]+)(" "|\n|\r\n) { insert_plot_print(yytext,yylineno,0); call++;}

<plot_option>[V]"("([0-9_A-Za-z]+)"#db)" { insert_plot_print(yytext,yylineno,1);call++; }
<plot_option>[V]([0-9_A-Za-z]+)"#db" { insert_plot_print(yytext,yylineno,1); call++;}

<plot_option>\n|\r\n { if(call == 0) insert_plot_print("",yylineno,0); BEGIN(INITIAL);}
<plot_option><<EOF>> { if(call == 0) insert_plot_print("",yylineno,0);  BEGIN(INITIAL);}
<plot_option>. { }




".opti"([(.)*]) { printf("DIE hello .opti\n");}
".width"([(.)*]) { printf("DIE hello .width\n");}


".OPTIONS"[ \t]+"SPARSE"[ \t]*(\n|\r\n)? {  set_sparse(); }
".OPTIONS"(([ \t]+"SPD"[ \t]+"iter")|([ \t]+"iter"[ \t]+"SPD"))[ \t]*(\n|\r\n)? { set_analysis_type(2); }
".OPTIONS"[ \t]+"iter"[ \t]*(\n|\r\n)? {  set_analysis_type(3); }
".OPTIONS"[ \t]+"SPD"[ \t]*(\n|\r\n)? {  set_analysis_type(1); }
".OPTIONS"[ \t]+"METHOD=TR"(\n|\r\n)? {  set_approx_type(1); /*printf("TR: %s\n", yytext);trapezoidal*/}
".OPTIONS"[ \t]+"METHOD=BE"(\n|\r\n)? {  set_approx_type(2); /*printf("BE: %s\n", yytext); backward euler*/}
".OPTIONS"[ \t]*(\n|\r\n)? {  }




"ITOL=" {BEGIN(itol_option);}
<itol_option>((0|[1-9][0-9]*)"."[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*) {set_itol(atof(yytext)); printf("************ITOL=%lf\n",atof(yytext)); BEGIN(INITIAL);}
<itol_option>[\t ]* {}


[RLC][0-9_A-Z]+ { BEGIN(RLC_elem); ignore_line = 0; stack_pos = 0; stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); }
<RLC_elem>[\t ]*
<RLC_elem>0|[1-9]+[0-9]* { stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); }
<RLC_elem>([_A-Z]*[0-9]*[_A-Z]*)*|[_A-Z]* 	{ if(stack_pos >= 1 && stack_pos <= 3){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<RLC_elem>((0|[1-9][0-9]*)"."[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*) {  if(stack_pos == 3){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<RLC_elem>\n|\r\n { BEGIN(INITIAL); if(stack_pos != 4){ ignore_line = 1; } insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<RLC_elem><<EOF>> { BEGIN(INITIAL); insert(stack, stack_pos,ignore_line,yylineno); if(stack_pos != 4){ exit(1);} stack_pos = 0;}      
<RLC_elem>. { fprintf(stderr, "unrecognized character ERROR1 [%s]\n",yytext); ignore_line = 1;}



[VI][0-9_A-Z]+ { BEGIN(cur_vol_source); exist_tran_spec = 0; ignore_line = 0; stack_pos = 0; stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); }
<cur_vol_source>[\t ]* 
<cur_vol_source>"transient" { }
<cur_vol_source>"AC"((([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)[,]?[ \t]*)){2})$ {
	parse_ac_mag_phase(yytext);
}
<cur_vol_source>"exp"([ \t]*)"("(([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)[,]?[ \t]*)+)")" {
	exist_tran_spec = 1;
	printf("voltage exp: %s\n", &yytext[3]);
	parse_transient_spec(&yytext[3], 'e');
}
<cur_vol_source>"sin"([ \t]*)"("(([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)[,]?[ \t]*)+)")" {
	exist_tran_spec = 1;
	printf("voltage sin: %s\n", &yytext[3]);
	parse_transient_spec(&yytext[3], 's');
}
<cur_vol_source>"pulse"([ \t]*)"("(([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)[,]?[ \t]*)+)")" {
	exist_tran_spec = 1;
	printf("voltage pulse: %s\n", &yytext[5]);
	parse_transient_spec(&yytext[5], 'p');
}
<cur_vol_source>"pwl"(([ \t]*)"("(([ \t]*(([0-9]*"."?)[0-9]*|[0-9]*"."?[0-9]*"e"[-+]?[0-9]*)[,]?[ \t]*)+)")")+ {
	exist_tran_spec = 1;
	printf("voltage pwl: %s\n", &yytext[3]);
	parse_transient_spec(&yytext[3], 'w');
}
<cur_vol_source>0|[1-9]+[0-9]* 	{ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); }
<cur_vol_source>([_A-Z]*[0-9]*[_A-Z]*)+|[_A-Z]* { if(stack_pos != 3 && exist_tran_spec != 1){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else {SYNTAX_CHECK} }
<cur_vol_source>((([0-9]*".")?[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*))[m]? { if(stack_pos == 3 && exist_tran_spec != 1){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); }else SYNTAX_CHECK }
<cur_vol_source>\n|\r\n 	{ BEGIN(INITIAL); if(stack_pos != 4 && exist_tran_spec != 1){ ignore_line = 1; } insert(stack, stack_pos, ignore_line,yylineno); stack_pos = 0;}
<cur_vol_source><<EOF>> { BEGIN(INITIAL); insert(stack, stack_pos,ignore_line,yylineno); if(stack_pos != 4){ exit(1);} stack_pos = 0;}
<cur_vol_source>. { fprintf(stderr, "unrecognized character ERROR %s\n", yytext); ignore_line = 1;}

[D][0-9_A-Z]+ { BEGIN(diode); ignore_line = 0; stack_pos = 0; stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}
<diode>[\t ]* 
<diode>0|[1-9]+[0-9]* { if((stack_pos >= 1 && stack_pos <= 3) || stack_pos == 4){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); } }
<diode>[_A-Z]*[0-9]*[_A-Z]*|[_A-Z]* 	{ if(stack_pos >= 1 && stack_pos <= 3){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<diode>("["?)(((0|[1-9][0-9]*)"."[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*))("]"?) { if(stack_pos==4){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<diode>\n|\r\n { BEGIN(INITIAL); if(stack_pos != 5 && stack_pos != 4){ ignore_line = 1; } insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<diode><<EOF>> { BEGIN(INITIAL); insert(stack, stack_pos,ignore_line,yylineno); if(stack_pos != 5){ exit(1); } stack_pos = 0;  }
<diode>. { fprintf(stderr, "[DIODE] unrecognized character ERROR |%s|\n",yytext); ignore_line = 1;}


[M][0-9_A-Z]+ { BEGIN(tr_mos); ignore_line = 0; stack_pos = 0; stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}
<tr_mos>[\t ]* 
<tr_mos>"L="((([0-9]*".")?[0-9]*e("-"|"+")?[1-9]+[0-9]*)|((0|[1-9][0-9]*)("."[0-9]*)?)) { if(stack_pos==6){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<tr_mos>"W="((([0-9]*".")?[0-9]*e("-"|"+")?[1-9]+[0-9]*)|((0|[1-9][0-9]*)("."[0-9]*)?)) { if(stack_pos==7){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<tr_mos>0|[1-9]+[0-9]*|[_A-Z_]*[1-9]+[0-9]*[_A-Z_]*|[_A-Z_]* 	{ if(stack_pos >= 1 && stack_pos <= 5){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<tr_mos>\n|\r\n { BEGIN(INITIAL); if(stack_pos != 8){ ignore_line = 1; } insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<tr_mos><<EOF>> { BEGIN(INITIAL); insert(stack, stack_pos,ignore_line,yylineno); if(stack_pos != 8){ exit(1); } stack_pos = 0;  } 
<tr_mos>. { fprintf(stderr, "unrecognized character ERROR\n"); ignore_line = 1;}


[Q][0-9_A-Z]+ { BEGIN(tr_bjt); ignore_line = 0; stack_pos = 0; stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}
<tr_bjt>[\t ]* 
<tr_bjt>0|[1-9]+[0-9]* { if((stack_pos >= 1 && stack_pos <= 3) || stack_pos == 5){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]); } }
<tr_bjt>[_A-Z]*[0-9]*[_A-Z_]*|[_A-Z_]* 	{ if(stack_pos >= 1 && stack_pos <= 4){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<tr_bjt>("["?)(((0|[1-9][0-9]*)"."[0-9]*)|(([0-9]*".")?[0-9]*e("-"|"+")?[0-9]*))("]"?) { if(stack_pos==5){ stack[stack_pos++] = strdup(yytext); ALLOC_CHECK(stack[stack_pos-1]);}else SYNTAX_CHECK }
<tr_bjt>\n|\r\n { BEGIN(INITIAL); if(stack_pos != 6 && stack_pos != 5){ ignore_line = 1; } insert(stack, stack_pos,ignore_line,yylineno); stack_pos = 0;}
<tr_bjt><<EOF>> { BEGIN(INITIAL); insert(stack, stack_pos,ignore_line,yylineno); if(stack_pos != 6){ exit(1);} stack_pos = 0; }
<tr_bjt>. { fprintf(stderr, "BJT unrecognized character ERROR |%s|\n",yytext); ignore_line = 1;}


"*".* 		{}
[ \t\n\r]   {}
<<EOF>> 	{printf("End Of File\n"); return 2;}
".END" {}
".op" {}
.       	{ printf("Unknown character [%c] - %s @ %d\n",yytext[0],yytext,yylineno);  return -1;    }
%%

/*
int main(int argc, char **argv)
{
	".TRAN"(([ \t]*(([0-9]*)"."[0-9]*)|([0-9]*e("-"|"+")?[0-9]*)|(([0-9]*"."?)[0-9]*e("-"|"+")?[0-9]*))+)+[ \t]* (\n|\r\n)? { set_tran(&yytext[5]);}
 yylex();
 return 0;
}
*/
