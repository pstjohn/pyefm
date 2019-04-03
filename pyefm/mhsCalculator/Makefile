################################################################################
################################################################################
# Author:  Christian Jungreuthmayer
# Date:    Tue Oct 15 11:05:43 CEST 2013
# Company: Austrian Centre of Industrial Biotechnology (ACIB)
################################################################################

CC         = gcc
TARGET     = mhsCalculator_linprog_glpk

# CC_OPTS = -g -O0 -g -O2
CC_OPTS = -O3 -Wall -g -O2
LD_OPTS = 


GURBOI_LIB = -lgurobi55 -lm -lpthread
GUROBI_LOCAL_OBJ = mhsCalculator_linprog_gurobi.o prepro.o read_files.o printing.o qsort_routines.o

GLPK_LIB   = -lglpk -lm
GLPK_LOCAL_OBJ = mhsCalculator_linprog_glpk.o prepro.o read_files.o printing.o qsort_routines.o

CPLEX_LIB  = -lcplex -lilocplex -pthread -lm
CPLEX_LOCAL_OBJ = mhsCalculator_linprog_cplex.o prepro.o read_files.o printing.o qsort_routines.o
CPLEX_SOLPOOL_LOCAL_OBJ = mhsCalculator_linprog_cplex_solpool.o prepro.o read_files.o printing.o qsort_routines.o

BERGE_LIB  = -lpthread -lm
BERGE_LOCAL_OBJ = mhsCalculator_berge.o prepro.o read_files.o printing.o qsort_routines.o berge.o tree.o

BRUTEFORCE_LIB  = -lpthread -lm
BRUTEFORCE_LOCAL_OBJ = mhsCalculator_bruteforce.o prepro.o read_files.o printing.o qsort_routines.o bruteforce.o tree.o

all: mhsCalculator_berge mhsCalculator_bruteforce convert_text_to_bin $(TARGET)

mhsCalculator_berge: $(BERGE_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_berge $(BERGE_LOCAL_OBJ) $(BERGE_LIB) $(LD_OPTS)

mhsCalculator_bruteforce: $(BRUTEFORCE_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_bruteforce $(BRUTEFORCE_LOCAL_OBJ) $(BRUTEFORCE_LIB) $(LD_OPTS)

mhsCalculator_linprog_gurobi: $(GUROBI_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_linprog_gurobi $(GUROBI_LOCAL_OBJ) $(GURBOI_LIB) $(LD_OPTS)

mhsCalculator_linprog_glpk: $(GLPK_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_linprog_glpk $(GLPK_LOCAL_OBJ) $(GLPK_LIB) $(LD_OPTS)

mhsCalculator_linprog_cplex: $(CPLEX_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_linprog_cplex $(CPLEX_LOCAL_OBJ) $(CPLEX_LIB) $(LD_OPTS)

mhsCalculator_linprog_cplex_solpool: $(CPLEX_SOLPOOL_LOCAL_OBJ)
	gcc $(CC_OPTS) -o mhsCalculator_linprog_cplex_solpool $(CPLEX_SOLPOOL_LOCAL_OBJ) $(CPLEX_LIB) $(LD_OPTS)

mhsCalculator_berge.o: mhsCalculator_berge.c mhsCalculator_berge.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -c  mhsCalculator_berge.c

mhsCalculator_bruteforce.o: mhsCalculator_bruteforce.c mhsCalculator_bruteforce.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -c  mhsCalculator_bruteforce.c

mhsCalculator_linprog_gurobi.o: mhsCalculator_linprog.c mhsCalculator_linprog.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -DGUROBI -c  mhsCalculator_linprog.c -o mhsCalculator_linprog_gurobi.o $(GUROBI_INC)

mhsCalculator_linprog_glpk.o: mhsCalculator_linprog.c mhsCalculator_linprog.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -DGLPK -c  mhsCalculator_linprog.c -o mhsCalculator_linprog_glpk.o $(GLPK_INC)

mhsCalculator_linprog_cplex.o: mhsCalculator_linprog.c mhsCalculator_linprog.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -DCPLEX -c  mhsCalculator_linprog.c -o mhsCalculator_linprog_cplex.o $(CPLEX_INC)

mhsCalculator_linprog_cplex_solpool.o: mhsCalculator_linprog.c mhsCalculator_linprog.h type_defs.h prepro.h read_files.h printing.h qsort_routines.h
	gcc $(CC_OPTS) -DCPLEX_SOLPOOL -c  mhsCalculator_linprog.c -o mhsCalculator_linprog_cplex_solpool.o $(CPLEX_INC)

prepro.o: prepro.c prepro.h type_defs.h qsort_routines.h
	gcc $(CC_OPTS) -c  prepro.c

read_files.o: read_files.c read_files.h type_defs.h
	gcc $(CC_OPTS) -c  read_files.c

printing.o: printing.c printing.h type_defs.h
	gcc $(CC_OPTS) -c  printing.c

qsort_routines.o: qsort_routines.c qsort_routines.h type_defs.h
	gcc $(CC_OPTS) -c  qsort_routines.c

berge.o: berge.c berge.h type_defs.h berge_type_defs.h tree_type_defs.h
	gcc $(CC_OPTS) -c  berge.c

bruteforce.o: bruteforce.c bruteforce.h type_defs.h berge_type_defs.h tree_type_defs.h
	gcc $(CC_OPTS) -c  bruteforce.c

tree.o: tree.c tree.h type_defs.h berge_type_defs.h tree_type_defs.h
	gcc $(CC_OPTS) -c  tree.c

convert_text_to_bin: convert_text_to_bin.c
	gcc $(CC_OPTS) convert_text_to_bin.c -o convert_text_to_bin -lm

clean_all: clean
	rm -f Makefile config.log config.status
	rm -fr autom4te.cache

clean:
	rm -f mhsCalculator_berge mhsCalculator_bruteforce mhsCalculator_linprog_gurobi mhsCalculator_linprog_glpk mhsCalculator_linprog_cplex mhsCalculator_linprog_cplex_solpool convert_text_to_bin mhsCalculator_berge.o prepro.o read_files.o printing.o qsort_routines.o berge.o bruteforce.o tree.o mhsCalculator_linprog_gurobi.o mhsCalculator_linprog_glpk.o mhsCalculator_linprog_cplex.o mhsCalculator_linprog_cplex_solpool.o mhsCalculator_bruteforce.o
