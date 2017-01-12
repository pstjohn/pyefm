///////////////////////////////////////////////////////////////////////////////
// Author: Christian Jungreuthmayer
// Email: christian.jungreuthmayer@acib.at
// Company: Austrian Centre of Industrial Biotechnology (ACIB)
// Web: http://www.acib.at
// Copyright (C) 2012, 2013
// Published unter GNU Public License V3
////////////////////////////////////////////////////////////////////////////////
// Basic Permissions.
//
// All rights granted under this License are granted for the term of copyright
// on the Program, and are irrevocable provided the stated conditions are met.
// This License explicitly affirms your unlimited permission to run the
// unmodified Program. The output from running a covered work is covered by
// this License only if the output, given its content, constitutes a covered
// work. This License acknowledges your rights of fair use or other equivalent,
// as provided by copyright law.
//
// You may make, run and propagate covered works that you do not convey, without
// conditions so long as your license otherwise remains in force. You may convey
// covered works to others for the sole purpose of having them make modifications
// exclusively for you, or provide you with facilities for running those works,
// provided that you comply with the terms of this License in conveying all
// material for which you do not control copyright. Those thus making or running
// the covered works for you must do so exclusively on your behalf, under your
// direction and control, on terms that prohibit them from making any copies of
// your copyrighted material outside their relationship with you.
//
// Disclaimer of Warranty.
//
// THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE
// LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR
// OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND,
// EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
// SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY
// SERVICING, REPAIR OR CORRECTION.
//
// Limitation of Liability.
//
// IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL
// ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE
// PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL,
// SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
// INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR
// DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR
// A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH
// HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#define USE_GLPK           1
#define USE_GUROBI         2
#define USE_CPLEX          3
#define USE_CPLEX_SOLPOOL  4
#define YES          0
#define NO           1

#define FIND_MAXIMUM -1
#define FIND_MINIMUM  1

#define LINPROG_MIP_GAP_REL 1e-07
#define LINPROG_MIP_GAP_ABS 1e-12

#ifdef GLPK
#define LP_TOOLKIT USE_GLPK
#endif

#if GUROBI
#define LP_TOOLKIT USE_GUROBI
#endif

#if CPLEX
#define LP_TOOLKIT USE_CPLEX
#endif

#if CPLEX_SOLPOOL
#define LP_TOOLKIT USE_CPLEX_SOLPOOL
#endif

#define DO_SHORTCUT YES
// #define DO_SHORTCUT NO

#define WRITE_LINPROG_FILE NO
// #define WRITE_LINPROG_FILE YES

#define INIT_MAX_EXPECTED_KNOCKOUT_COMBINATIONS 50000

#if LP_TOOLKIT == USE_GLPK
#include <glpk.h>
#endif
#if LP_TOOLKIT == USE_GUROBI
#include "gurobi_c.h"
#endif
#if LP_TOOLKIT == USE_CPLEX
#include "ilcplex/cplex.h"
#endif
#if LP_TOOLKIT == USE_CPLEX_SOLPOOL
#include "ilcplex/cplex.h"

#define MAX_SOLUTION_ENUM 1000000
#define MAX_SOLUTION_CAPACITY 1000000
#define SOLUTION_POOL_INTENSITY 4
#define SOLUTION_POOL_GAP 0.00001
#endif

#include "type_defs.h"
#include "printing.h"
#include "read_files.h"
#include "prepro.h"
////////////////////////////////////////////////////////////////////////////////

int g_n_out_of_D;
unsigned long long int g_size_linprog_mem;
int *g_pnt_norms;
unsigned long long *g_pnt_tmp_cutset;
double *g_lp_max;
double *g_r;
unsigned int g_solutions_cnt = 0;
unsigned long g_num_efm_writes = 0;
int g_shortcut_possible = 0;
FILE *g_fh_out;
unsigned long long int g_max_sol_enhancer = 1;
unsigned int g_uncomp_solutions_cnt = 0;
unsigned int g_max_expected_knockout_combinations;

////////////////////////////////////////////////////////////////////////////////
// glpk variables
////////////////////////////////////////////////////////////////////////////////
#if LP_TOOLKIT == USE_GLPK
glp_prob *g_lp;
glp_iocp g_parm;
int *g_pnt_ia;
int *g_pnt_ja;
double *g_pnt_ar;
#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Gurobi variables
////////////////////////////////////////////////////////////////////////////////
#if LP_TOOLKIT == USE_GUROBI
GRBenv   *g_gurobi_masterenv = NULL;
GRBmodel *g_gurobi_model     = NULL;
GRBenv   *g_gurobi_modelenv  = NULL;
char *g_gurobi_col_is_integer;
double *g_gurobi_obj_cf;
int *g_gurobi_start;
int *g_gurobi_len;
int *g_gurobi_index;
double *g_gurobi_value;
char *g_gurobi_row_constr_type;
double *g_gurobi_row_righthandside;
unsigned int g_solutions_cnt_last = 0;
char **g_gurobi_varnames;
#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Cplex variables
////////////////////////////////////////////////////////////////////////////////
#if LP_TOOLKIT == USE_CPLEX
CPXENVptr g_cplex_env = NULL;
CPXLPptr g_cplex_lp = NULL;
char *g_cplex_col_is_integer;
double *g_cplex_obj_cf;
int *g_cplex_start;
int *g_cplex_len;
int *g_cplex_index;
double *g_cplex_value;
char *g_cplex_row_constr_type;
double *g_cplex_row_righthandside;
unsigned int g_solutions_cnt_last = 0;
char **g_cplex_varnames;
double *g_cplex_lb;
double *g_cplex_ub;
// transposed
int *g_cplex_start_trans;
int *g_cplex_index_trans;
double *g_cplex_value_trans;
#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Cplex variables
////////////////////////////////////////////////////////////////////////////////
#if LP_TOOLKIT == USE_CPLEX_SOLPOOL
CPXENVptr g_cplex_env = NULL;
CPXLPptr g_cplex_lp = NULL;
char *g_cplex_col_is_integer;
double *g_cplex_obj_cf;
int *g_cplex_start;
int *g_cplex_len;
int *g_cplex_index;
double *g_cplex_value;
char *g_cplex_row_constr_type;
double *g_cplex_row_righthandside;
unsigned int g_solutions_cnt_last = 0;
char **g_cplex_varnames;
double *g_cplex_lb;
double *g_cplex_ub;
// transposed
int *g_cplex_start_trans;
int *g_cplex_index_trans;
double *g_cplex_value_trans;
#endif
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void usage(char *message);
void handle_arguments(int largc, char **largv, struct_cmd_options *cmd_options);
void do_frees(struct_reac_info *reac_info, struct_mode_info *mode_info, struct_map_info *map_info);
void uncompress_and_write_solution(int sol_num, struct_map_info *map_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info);
void print_uncomp_cutset_to_file(double *new_l_r, unsigned int new_l_r_cnt, unsigned int l_unc_sol_cnt, struct_reac_info *reac_info, struct_map_info *map_info);

// glpk specific methods
#if LP_TOOLKIT == USE_GLPK
void init_glpk(struct_mode_info *mode_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_reac_info *reac_info);
void resize_arrays_glpk(struct_mode_info *mode_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_reac_info *reac_info);
void set_glpk_prob(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info *reac_info);
void fill_linprog_mem_glpk(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info);
int do_linprog_mip_glpk(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info);
void cleanup_glpk();
#endif

// gurobi specific methods
#if LP_TOOLKIT == USE_GUROBI
void init_gurobi(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void resize_arrays_gurobi(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void fill_linprog_mem_gurobi(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info);
int do_linprog_mip_gurobi(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info);
void cleanup_gurobi(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info);
#endif

// cplex specific methods
#if LP_TOOLKIT == USE_CPLEX
void init_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void resize_arrays_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void fill_linprog_mem_cplex(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info);
int do_linprog_mip_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info);
void cleanup_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info);
void print_g_cplex_arrays(struct_map_info *map_info, struct_mode_info *mode_info, struct_cmd_options *cmd_options, unsigned long long int write_cnt);
void cplex_transpose_arrays(unsigned long long int num_idx, unsigned long long int num_rows, unsigned long long int num_cols);
void cplex_clean_transpose();
#endif

// cplex specific methods
#if LP_TOOLKIT == USE_CPLEX_SOLPOOL
void init_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void resize_arrays_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info);
void fill_linprog_mem_cplex(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info);
int do_linprog_mip_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info);
void cleanup_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info);
void print_g_cplex_arrays(struct_map_info *map_info, struct_mode_info *mode_info, struct_cmd_options *cmd_options, unsigned long long int write_cnt);
void cplex_transpose_arrays(unsigned long long int num_idx, unsigned long long int num_rows, unsigned long long int num_cols);
void cplex_clean_transpose();
#endif
////////////////////////////////////////////////////////////////////////////////
