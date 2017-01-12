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
#ifndef PREPRO_HEADER
#define PREPRO_HEADER 1
#include <stdlib.h>
#include <math.h>

#include "type_defs.h"
#include "berge_type_defs.h"
#include "qsort_routines.h"
#include "printing.h"

#define PREPRO_REMOVE_ESSENTIAL_REACTIONS 1
#define PREPRO_REMOVE_DUPLICATE_MODES     1
#define PREPRO_REMOVE_SUPERSET_MODES      1
#define PREPRO_REMOVE_DUPLICATE_COLS      1

#define INIT_NORM_CONSTRAINTS 20000

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void find_essential_reacs(struct_mode_info *mode_info, struct_cmd_options cmd_options, struct_reac_info reac_info);
void update_essential_reacs(struct_reac_info *reac_info, struct_mode_info *mode_info);
void initialize_map_remove_essentials(struct_map_info *map_info, struct_reac_info reac_info, struct_mode_info *mode_info);
void free_map_memory(struct_map_info *map_info);
void do_preprocessing(struct_map_info *map_info, struct_reac_info *reac_info, struct_mode_info *mode_info, struct_cmd_options cmd_options);
void free_mem_map_info(struct_map_info *map_info);
void create_map_remove_essentials(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info reac_info);
void do_mem_allocation_bin_second(struct_mode_info *mode_info, struct_map_info map_info);
void restructure_flux_array_essentials(struct_reac_info reac_info, struct_mode_info *mode_info, struct_map_info map_info);
void remove_duplicate_modes(struct_mode_info *mode_info);
void init_always_zero_variables(struct_map_info *map_info, struct_mode_info mode_info, struct_reac_info reac_info);
void remove_superset_modes(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info reac_info, struct_cmd_options cmd_options);
void init_duplicate_col_is_keeper_array(struct_map_info *map_info, struct_reac_info reac_info);
void find_duplicate_columns(struct_map_info *map_info, struct_mode_info mode_info, struct_reac_info reac_info);
struct_mode* allocate_transpose_memory(int t_unit_size, struct_map_info map_info);
void create_tranpose_mode_array(struct_mode **l_mode_tranpose, int t_unit_size, struct_mode_info mode_info, struct_map_info map_info);
void identify_duplicate_cols(struct_mode *l_mode_tranpose, int t_unit_size, struct_map_info *map_info, struct_reac_info reac_info);
void restructure_flux_array_duplicate_columns(struct_mode_info *mode_info, struct_reac_info reac_info, struct_map_info map_info);
void uncompress_cutsets(struct_mode_info *mode_info, struct_cutset_info *cutset_info, struct_reac_info *reac_info, struct_map_info *map_info);
////////////////////////////////////////////////////////////////////////////////
#endif
