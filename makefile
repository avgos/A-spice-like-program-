IDIR =./include/
CC=gcc
CFLAGS=-I$(IDIR) -Wall -O0 -g -std=gnu99
BDIR=./src
ODIR=./obj
LIBS=-lfl -lgsl -lgslcblas -lm ./CXSparse/Lib/libcxsparse.a

_OBJ = util.o dc_analysis.o transient_analysis.o main.o node_to_int.o lex.yy.o csparse.o freq_analysis.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: cxsparse parser sol_parser


cxsparse:
	$(MAKE) -C CXSparse library

parser: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

sol_parser: obj/sol_parser.o
	$(CC) -o $@ $^ $(CFLGAS) -lm

obj/freq_analysis.o: src/freq_analysis.c
	$(CC) -c $< -o $@ $(CFLAGS)

obj/transient_analysis.o: src/transient_analysis.c
	$(CC) -c $< -o $@ $(CFLAGS)

obj/sol_parser.o: src/sol_parser.c
	$(CC) -c $< -o $@ $(CFLAGS)
obj/csparse.o: src/csparse.c
	$(CC) -c $< -o $@ $(CFLAGS)
obj/dc_analysis.o: src/dc_analysis.c
	$(CC) -c $< -o $@  $(CFLAGS)
obj/util.o: src/util.c
	$(CC) -c $< -o $@  $(CFLAGS)
obj/main.o: src/main.c
	$(CC) -c $< -o $@  $(CFLAGS)
obj/node_to_int.o: src/node_to_int.c
	$(CC) -c $< -o $@  $(CFLAGS)
obj/lex.yy.o: src/lex.yy.c
	$(CC) -c $< -o $@  $(CFLAGS)

src/lex.yy.c: src/parser.l
	cd src; flex parser.l ; cd .. ; 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o parser
	rm -f sol_parser
	rm -f LU_sol
	rm -f cholesky_sol
	rm -f cg_sol
	rm -f bicg_sol
	rm -f BiCG_sol_nodes
	rm -f CG_sol_nodes
	rm -f LU_sol_nodes
	rm -f cholesky_sol_nodes
	#$(MAKE) -C CXSparse clean

