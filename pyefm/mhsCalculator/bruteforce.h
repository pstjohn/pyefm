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
#ifndef BERGE_HEADER
#define BERGE_HEADER 1
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>

#include "type_defs.h"
#include "berge_type_defs.h"
#include "tree_type_defs.h"
#include "tree.h"

#define INIT_MCS_ELEMS     100000
#define INIT_PCS_ELEMS     MAX_THREADS
#define INIT_PCS_ELEMS_THR  5000

#define VERBOSITY_CNT 100
FILE *g_fh_out;
long long int g_uncomp_solutions_cnt = 0;
unsigned long long int g_max_sol_enhancer = 1;
long int g_num_bruteforce_combinations = 0;
long int g_bruteforce_depth = 0;
long int g_invoke_cnt = 0;
long int g_found_cutsets = 0;
unsigned *g_bruteforce_cutset;
unsigned long long int *g_pnt_tmp_cutset2;
unsigned long long int *g_n_over_k;
time_t g_start_time;
struct structInThread
{
   int thread_id;
   struct_mode_info *mode_info;
   struct_map_info *map_info;
   struct_reac_info *reac_info;
   struct_cutset_info *cutset_info;
   struct_tree_info   *tree_info;
   struct_cmd_options *cmd_options;
};
struct structInThread g_strInThread[MAX_THREADS];
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

void do_mincutset_bruteforce(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_cutset_info *cutset_info);
void set_max_bruteforce_combinations(int nr_cancels, int num_reacs);
unsigned long long int n_over_k(int n, int k);
void *worker_bruteforce(void *str_ptr);
void evaluate_cutset_bruteforce(int thread_id, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info,
                                struct_mode_info *mode_info, struct_cutset_info *cutset_info, struct_tree_info *tree_info);
void fill_n_over_k(int cur_num_cancellations, int max_cancellation, int num_dupset_keepers, int num_always_zero_reacs );
void allocate_bruteforce_mem(struct_cmd_options *cmd_options, struct_map_info *map_info, struct_mode_info *mode_info);
void free_bruteforce_mem();
void uncompress_and_write_solution(int thread_id, int sol_num, struct_map_info *map_info, struct_reac_info *reac_info, struct_mode_info *mode_info, char *filename, int max_cancellations);
void print_uncomp_cutset_to_file(double *new_l_r, unsigned int new_l_r_cnt, unsigned int l_unc_sol_cnt, struct_reac_info *reac_info, struct_map_info *map_info);


void init_cs_arr_bruteforce(struct_cutset_info *cutset_info, struct_tree_info *tree_info, struct_mode_info *mode_info, struct_cmd_options cmd_options,
                       struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_reac_info reac_info);
void free_mem_cutset_info(struct_cutset_info *cutset_info, struct_cmd_options cmd_options);
void reorder_reaction_of_modes(struct_mode_info *mode_info, struct_map_info map_info, struct_bitorder_info bitorder_info, struct_cmd_options cmd_options);
int compute_precutsets(long long int e, unsigned long long int * t_cutset, unsigned long long int *m_cutset,
                       struct_mode_info *mode_info, struct_map_info *map_info, struct_cutset_info *cutset_info,    
                       struct_tree_info *tree_info, struct_cmd_options cmd_options);
int calc_num_ones(unsigned long long *cs, unsigned int num_unit_size);
void check_and_resize_arr_pcs(struct_cutset_info *cutset_info, struct_mode_info *mode_info);
void check_and_resize_arr_mcs(struct_cutset_info *cutset_info, struct_mode_info *mode_info);
int enough_good_modes_survive(unsigned long long *cs, struct_cmd_options cmd_options, struct_mode_info *mode_info);
int is_mode_killed(unsigned long long *cs, unsigned long long *t_pnt_flux, int idx, struct_mode_info *mode_info);
void fill_precutsets_to_cutsets(struct_tree_info *tree_info, struct_cutset_info *cutset_info, struct_mode_info *mode_info,
                                struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_cmd_options cmd_options);
int find_superset_in_tree(long int node_idx, long int i_pcs, int thr_id, struct_tree_info *tree_info, struct_cutset_info *cutset_inf, struct_mode_info *mode_info);
void reorder_reaction_of_cutsets(struct_cutset_info *cutset_info, struct_map_info *map_info, struct_mode_info *mode_info, struct_bitorder_info *bitorder_info);
#endif
