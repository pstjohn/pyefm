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

#ifndef TYPE_DEFS
#define TYPE_DEFS 1

#define VERSION "1.1"
#define YEARS   "2012, 2013"
#define AUTHORS "Christian Jungreuthmayer"
#define COMPANY "Austrian Centre of Industrial Biotechnology (ACIB)"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
typedef struct
{
   char **reactions;
   int num_reactions;
   int max_len_reac_name;
   char **readin_essential_reactions;
   unsigned int *readin_essential_idx;
   unsigned int num_readin_essential_reacs;
} struct_reac_info;

typedef struct
{
   unsigned long long int num_efms;
   unsigned long long int bad_efms;
   unsigned long long int *pnt_good_first;
   unsigned long long int *pnt_flux_first;
   unsigned long long int *pnt_reac_occ;
   unsigned long long int *pnt_reac_occ_bad;
   unsigned int *essential_idx;
   unsigned int *is_essential;
   unsigned int num_unit_size_first;
   unsigned int num_essential_reacs;
   unsigned long long int *pnt_good;
   unsigned long long int *pnt_flux;
   unsigned int num_unit_size;
   unsigned long long int *pnt_tmp_cutset;
   unsigned int max_norm;
   unsigned long long int max_norm_constraints;
   unsigned long long int *norm_bitmap;
   unsigned long long int num_norm_constraints;
   unsigned int num_unit_size_orig;
   unsigned long long int *pnt_good_duplicates;
   unsigned long long int good_emfs;
   unsigned long long int good_emfs_orig;
   unsigned long long int num_good_removed_by_always_zero;
   unsigned long long int num_good_removed_because_all_zero;
} struct_mode_info;

typedef struct
{
   char *o_filename;
   char *r_filename;
   char *m_filename;
   char *e_filename;
   unsigned long long int good_efms;
   unsigned long long int good_efms_wanted;
   int lin_superset_test;
   int apply_heuristics;
   int num_threads;
   int max_cancellations;
   int rand_seed_integer;
   int solution_range;
   int output_bitvector;
} struct_cmd_options;

typedef struct
{
   unsigned int map_len;
   unsigned int *map;
   unsigned int num_always_zero_reacs;
   unsigned long long int num_always_zero_modes;
   unsigned long long int *always_zero_bitmap;
   unsigned int *always_zero_map;
   unsigned int *is_always_zero;
   int *always_zero_mover;
   unsigned int *dupsets_keep_reac;
   unsigned int *dupsets_is_keeper;
   unsigned int *dupsets_map;
   int num_dupset_keepers;
   unsigned int num_dupsets;
   unsigned int *num_reac_dupset;
   unsigned int *dupsets;
   unsigned long long int unit_col_length;
} struct_map_info;

typedef struct
{
  unsigned long long int *flux;
  int norm;
  int unit_size;
} struct_flux;

typedef struct
{
  unsigned long long int *flux;
  int orig_reac;
  unsigned long long int unit_col_length;
} struct_mode;
////////////////////////////////////////////////////////////////////////////////
#endif
