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
#include "mhsCalculator_linprog.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int main(int argc, char **argv)
{
   struct timeval time_start_calc;
   struct timeval time_before_prepro;
   struct timeval time_before_mincutset;
   struct timeval time_stop_calc;

   struct_cmd_options cmd_options = {NULL, // o_filename
                                     NULL, // r_filename
                                     NULL, // m_filename
                                     NULL, // e_filename
                                     0,    // good_emfs
                                     0,    // wanted_emfs
                                     0,    // lin_superset_test
                                     0,    // apply_heuristics
                                     1,    // num_threads
                                     0,    // max_cancellations
                                     0};   // rand_seed_integer
   struct_reac_info reac_info = {NULL, // array of reaction names
                                 0,    // num_reactions
                                 0,    // max_len_reac_name
                                 NULL, // readin_essential_reactions
                                 NULL, // readin_essential_idx
                                 0};   // num_readin_essential_reacs

   struct_mode_info mode_info = {0,     // num_emfs
                                 0,     // num_bad_emfs
                                 NULL,  // pnt_good_first
                                 NULL,  // pnt_flux_first
                                 NULL,  // pnt_reac_occ
                                 NULL,  // pnt_reac_occ_bad
                                 NULL,  // essential_idx
                                 NULL,  // is_essential
                                 0,     // num_unit_size_first
                                 0,     // num_essential_reacs
                                 NULL,  // pnt_good
                                 NULL,  // pnt_flux
                                 0,     // num_unit_size
                                 NULL,  // pnt_tmp_cutset
                                 0,     // max_norm
                                 0,     // max_norm_constraints
                                 NULL,  // norm_bitmap
                                 0,     // num_norm_constraints
                                 0};    // num_unit_size_orig

   struct_map_info map_info = {0,     // map_len
                               NULL,  // map
                               0,     // num_always_zero_reacs
                               0,     // num_always_zero_modes
                               NULL,  // always_zero_bitmap
                               NULL,  // always_zero_map
                               NULL,  // is_always_zero
                               NULL,  // always_zero_mover
                               NULL,  // dupsets_keep_reac
                               NULL,  // dupsets_is_keeper
                               NULL,  // dupsets_map
                               0,     // num_dupset_keepers
                               0,     // num_dupsets
                               NULL,  // num_reac_dupset
                               NULL,  // dupsets
                               0};    // unit_col_length

   gettimeofday(&time_start_calc,NULL);

   printf("INFO: program %s started.\n",argv[0]);

   handle_arguments(argc, argv, &cmd_options);


   readin_reactions_file(cmd_options.r_filename, &reac_info);
   print_reactions(reac_info);

   readin_efm_file_bin(cmd_options, reac_info, &mode_info);
   print_reac_occ_arr(mode_info, reac_info, cmd_options);

   // find essential reactions by inspection of modes
   find_essential_reacs(&mode_info, cmd_options, reac_info);
   print_essential_reacs(mode_info, reac_info);

   // read in essential reactions provided by user
   if( cmd_options.e_filename != NULL )
   {
      readin_essential_file(cmd_options.e_filename, &reac_info);
      update_essential_reacs(&reac_info, &mode_info);
   }
   print_final_essential_reacs(mode_info, reac_info);

   gettimeofday(&time_before_prepro,NULL);

   // do pre-processing
   do_preprocessing(&map_info, &reac_info, &mode_info, cmd_options);

   gettimeofday(&time_before_mincutset,NULL);

   ////////////////////////////////////////////////////////////////////////////
   // do linear programming
   ////////////////////////////////////////////////////////////////////////////
#if LP_TOOLKIT == USE_GLPK
   init_glpk(&mode_info, &map_info, &cmd_options, &reac_info);

   do
   {
      gettimeofday(&time_stop_calc,NULL);
      printf("calculation time without reading in and pre-processing modes ");
      display_execution_time(time_stop_calc, time_before_mincutset);

      set_glpk_prob(&mode_info, &map_info, &reac_info);
      fill_linprog_mem_glpk(&mode_info, &cmd_options, &map_info, &reac_info);
   }while(do_linprog_mip_glpk(&reac_info, &map_info, &cmd_options, &mode_info));

   cleanup_glpk();
#endif

#if LP_TOOLKIT == USE_GUROBI
   init_gurobi(&mode_info, &reac_info, &cmd_options, &map_info);

   do
   {
      gettimeofday(&time_stop_calc,NULL);
      printf("calculation time without reading in and pre-processing modes ");
      display_execution_time(time_stop_calc,time_before_mincutset);

      fill_linprog_mem_gurobi(&mode_info, &cmd_options, &map_info, &reac_info);
   }while(do_linprog_mip_gurobi(&reac_info, &map_info, &cmd_options, &mode_info));

   cleanup_gurobi(&reac_info, &map_info, &mode_info);
#endif

#if LP_TOOLKIT == USE_CPLEX
   init_cplex(&mode_info, &reac_info, &cmd_options, &map_info);

   do
   {
      gettimeofday(&time_stop_calc,NULL);
      printf("calculation time without reading in and pre-processing modes ");
      display_execution_time(time_stop_calc,time_before_mincutset);

      fill_linprog_mem_cplex(&mode_info, &cmd_options, &map_info, &reac_info);
   }while(do_linprog_mip_cplex(&reac_info, &map_info, &cmd_options, &mode_info));


   cleanup_cplex(&reac_info, &map_info, &mode_info);
#endif

#if LP_TOOLKIT == USE_CPLEX_SOLPOOL
   init_cplex(&mode_info, &reac_info, &cmd_options, &map_info);

   do
   {
      gettimeofday(&time_stop_calc,NULL);
      printf("calculation time without reading in and pre-processing modes ");
      display_execution_time(time_stop_calc,time_before_mincutset);

      fill_linprog_mem_cplex(&mode_info, &cmd_options, &map_info, &reac_info);
   }while(do_linprog_mip_cplex(&reac_info, &map_info, &cmd_options, &mode_info));


   cleanup_cplex(&reac_info, &map_info, &mode_info);
#endif
   ////////////////////////////////////////////////////////////////////////////


   // print_min_max_knockouts(&mode_info, &map_info, &cutset_info);
   // printf("Number of found cutsets: %llu\n",cutset_info.mcs_idx);
   if( cmd_options.o_filename != NULL )
   {
      fclose(g_fh_out);
   }
   printf("Total number of (uncompressed) solutions: %u\n",g_uncomp_solutions_cnt);

   gettimeofday(&time_stop_calc,NULL);

   printf("Time reading files:    "); display_execution_time(time_before_prepro,time_start_calc);
   printf("Preprocessing time:    "); display_execution_time(time_before_mincutset,time_before_prepro);
   printf("Mincutset time:        "); display_execution_time(time_stop_calc,time_before_mincutset);
   printf("Total execution time:  "); display_execution_time(time_stop_calc,time_start_calc);

   // clean up
   do_frees(&reac_info, &mode_info, &map_info);
   printf("INFO: program %s stopped.\n",argv[0]);
   return(0);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void do_frees(struct_reac_info *reac_info, struct_mode_info *mode_info, struct_map_info *map_info)
{
   free_mem_reac_info(reac_info);
   free_mem_mode_info(mode_info);
   free_mem_map_info(map_info);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void handle_arguments(int largc, char **largv, struct_cmd_options *cmd_options)
{
   int r_opt;
   int i;
   char *char_num_good_modes    = NULL;
   char *char_num_wanted_modes  = NULL;
   char *char_num_threads       = NULL;
   char *char_max_cancellations = NULL;
   char *char_solution_range    = NULL;
   char *char_endpointer;

   printf("mhsCalculator_linprog Version: %s\n",VERSION);
   printf("Copyright (C) %s %s %s\n",YEARS, AUTHORS, COMPANY);

   printf("Executed command: %s\n",largv[0]);
   printf("Options:");
   for(i = 1; i < largc; i++ )
   {
      printf(" %s",largv[i]);
   }
   printf("\n");

   while(( r_opt = getopt(largc, largv, "hlkm:r:e:o:n:w:t:s:")) != -1 )
   {
      switch(r_opt)
      {
         case 'h':
            usage("");
            break;
         case 'k':
            cmd_options->apply_heuristics = 1;
            break;
         case 'm':
            cmd_options->m_filename = optarg;
            break;
         case 'r':
            cmd_options->r_filename = optarg;
            break;
         case 'e':
            cmd_options->e_filename = optarg;
            break;
         case 'o':
            cmd_options->o_filename = optarg;
            break;
         case 'n':
            char_num_good_modes = optarg;
            break;
         case 'w':
            char_num_wanted_modes = optarg;
            break;
         case 't':
            char_num_threads = optarg;
            break;
         case 's':
            char_solution_range = optarg;
            break;
         case '?':
            usage("ERROR: invalid arguments\n");
            break;
         default:
            usage("ERROR: invalid arguments\n");
      }
   }
   ////////////////////////////////////////////////////////////////////////////
   // file names
   ////////////////////////////////////////////////////////////////////////////
   if( cmd_options->m_filename == NULL )
   {
      usage("ERROR: name of file containing modes (in binary form) not provided\n");
   }

   if( cmd_options->r_filename == NULL )
   {
      usage("ERROR: name of file containing reaction names not provided\n");
   }

   if( char_num_good_modes == NULL )
   {
      usage("ERROR: number of good modes not defined (Note: good modes must be at top of mode file)\n");
   }

   if( cmd_options->o_filename != NULL )
   {
      g_fh_out = fopen(cmd_options->o_filename, "w");

      if( g_fh_out <= 0 )
      {
         printf("FATAL ERROR: open file '%s' for writing failed: %s\n",cmd_options->o_filename,strerror(errno));
         exit(EXIT_FAILURE);
      }
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number of maximum cancellations
   ////////////////////////////////////////////////////////////////////////////
   if( char_solution_range == NULL )
   {
      cmd_options->solution_range = 0;
   }
   else
   {
      cmd_options->solution_range = strtol(char_solution_range, &char_endpointer, 10);
      if( char_endpointer == char_solution_range || errno == EINVAL || errno == ERANGE )
      {
         printf("FATAL ERROR: error while converting solution range (-s %s): %s\n",char_solution_range,strerror(errno));
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
   }
   if( cmd_options->solution_range < 0 )
   {
      printf("FATAL ERROR: solution range (-s) is less than 0 (%d)\n",cmd_options->solution_range);
      exit(EXIT_FAILURE);
   }
   else if( cmd_options->solution_range == 0 )
   {
      printf("INFO: unlimited solution range\n");
   }
   else
   {
      printf("INFO: solution range: %d\n",cmd_options->solution_range);
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number of 'good' efms
   ////////////////////////////////////////////////////////////////////////////
   cmd_options->good_efms = strtol(char_num_good_modes, &char_endpointer, 10);
   if( char_endpointer == char_num_good_modes || errno == EINVAL || errno == ERANGE )
   {
      fprintf(stderr, "FATAL ERROR: error while converting number of 'good' modes (-n %s): %s\n",char_num_good_modes,strerror(errno));
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }
   printf("INFO: number of keeper modes (-n %s): %llu\n",char_num_good_modes,cmd_options->good_efms);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number of threads
   ////////////////////////////////////////////////////////////////////////////
   if( char_num_threads == NULL )
   {
      cmd_options->num_threads = 1;
   }
   else
   {
      cmd_options->num_threads = strtol(char_num_threads, &char_endpointer, 10);
      if( char_endpointer == char_num_threads || errno == EINVAL || errno == ERANGE )
      {
         fprintf(stderr, "FATAL ERROR: error while converting number of threads (-t %s): %s\n",char_num_threads,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
   }
   printf("INFO: number of threads: %d\n",cmd_options->num_threads);

   if( cmd_options->num_threads <= 0 )
   {
      fprintf(stderr, "FATAL ERROR: number of threads provided as argument (%d) is less or equal 0\n",cmd_options->num_threads);
      exit(EXIT_FAILURE);
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number of maximum cancellations
   ////////////////////////////////////////////////////////////////////////////
   if( char_max_cancellations == NULL )
   {
      cmd_options->max_cancellations = 0;
   }
   else
   {
      cmd_options->max_cancellations = strtol(char_max_cancellations, &char_endpointer, 10);
      if( char_endpointer == char_max_cancellations || errno == EINVAL || errno == ERANGE )
      {
         fprintf(stderr, "FATAL ERROR: error while converting number of cancellations (-c %s): %s\n",char_max_cancellations,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
   }
   if( cmd_options->max_cancellations < 0 )
   {
      fprintf(stderr, "FATAL ERROR: number of maximum number of cancellations is less than 0 (%d)\n",cmd_options->max_cancellations);
      exit(EXIT_FAILURE);
   }
   else if( cmd_options->max_cancellations == 0 )
   {
      printf("INFO: unlimited number of cancellations are allowed\n");
   }
   else
   {
      printf("INFO: number of maximum cancellations: %d\n",cmd_options->max_cancellations);
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number wanted 'good' modes
   ////////////////////////////////////////////////////////////////////////////
   if( char_num_wanted_modes != NULL )
   {
      cmd_options->good_efms_wanted = strtol(char_num_wanted_modes, &char_endpointer, 10);
      if( char_endpointer == char_num_wanted_modes || errno == EINVAL || errno == ERANGE )
      {
         fprintf(stderr, "FATAL ERROR: error while converting number of wanted 'good' modes (-w %s): %s\n",char_num_wanted_modes,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      printf("INFO: Number of modes that have to be kept was specified\n");
      printf("INFO: number of keeper modes that have to be kept: %llu\n",cmd_options->good_efms_wanted);
      if( cmd_options->good_efms < cmd_options->good_efms_wanted )
      {
         fprintf(stderr, "ERROR: number of keeper modes (%llu) is less than number of keeper modes that must be kept (%llu)\n",
                 cmd_options->good_efms,cmd_options->good_efms_wanted);
         fprintf(stderr, "       execution aborted.\n");
         exit(EXIT_FAILURE);
      }
   }
   else
   {
      cmd_options->good_efms_wanted = cmd_options->good_efms;
   }

   if( cmd_options->good_efms < 0 )
   {
         fprintf(stderr, "ERROR: number of keeper modes (%llu) must be 0 or greater\n",cmd_options->good_efms);
         fprintf(stderr, "       execution aborted.\n");
         exit(EXIT_FAILURE);
   }

   if( cmd_options->good_efms_wanted < 0 )
   {
         fprintf(stderr, "ERROR: number of keeper modes that must be kept (%llu) must not be less than 0\n",cmd_options->good_efms_wanted);
         fprintf(stderr, "       execution aborted.\n");
         exit(EXIT_FAILURE);
   }

   if( cmd_options->good_efms_wanted > cmd_options->good_efms )
   {
      fprintf(stderr, "ERROR: number of good modes (%llu) must be larger than or equal to wanted good modes (%llu)\n",
              cmd_options->good_efms,cmd_options->good_efms_wanted);
      fprintf(stderr, "       execution aborted.\n");
      exit(EXIT_FAILURE);
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   if( cmd_options->good_efms > 0 && cmd_options->good_efms_wanted < cmd_options->good_efms )
   {
      g_n_out_of_D = 1;
   }
   ////////////////////////////////////////////////////////////////////////////

   return;
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void usage(char *message)
{
   printf("%s",message);
#if LP_TOOLKIT == USE_GLPK
   printf("usage: mhsCalculator_linprog_glpk -m efms.sorted.bin -r rfile [-e essential_reacs.txt] [-o output_cutsets] ");
#elif LP_TOOLKIT == USE_GUROBI
   printf("usage: mhsCalculator_linprog_gurobi -m efms.sorted.bin -r rfile [-e essential_reacs.txt] [-o output_cutsets] ");
#elif LP_TOOLKIT == USE_CPLEX
   printf("usage: mhsCalculator_linprog_cplex -m efms.sorted.bin -r rfile [-e essential_reacs.txt] [-o output_cutsets] ");
#elif LP_TOOLKIT == USE_CPLEX_SOLPOOL
   printf("usage: mhsCalculator_linprog_cplex_solpool -m efms.sorted.bin -r rfile [-e essential_reacs.txt] [-o output_cutsets] ");
#endif
   printf("-n num_good_modes [-k] [-w num_wanted_good_modes] [-t num_threads] [-s solution_range]\n");
   printf("\n");
   printf("-m ..... filename containing elementary flux modes (in binary form!)\n");
   printf("-r ..... filename containing names of reactions\n");
   printf("-e ..... filename containing essential reactions\n");
   printf("-o ..... filename (output) of file containing computed minimal cutsets\n");
   printf("-n ..... number of good/keeper modes\n");
   printf("-w ..... number of wanted good modes (number of modes that must survive if cutset is applied)\n");
   printf("-t ..... number of parallel threads\n");
   printf("-k ..... bail out of duplicate mode check if it seems to be ineffective\n");
   printf("-s ..... solution range defines how many knockouts a solution may have,\n");
   printf("         a solution range of 1 means that only the solutions with\n");
   printf("         the lowest number of knockouts are calculated, e.g. 10 knockout\n");
   printf("         a solution range of 2 means that solutions with the minimum number\n");
   printf("         knockouts are calculated and the solutions with minimum plus 1\n");
   printf("         e.g. 10 knockouts and 11 knockout\n");
   printf("-h ..... print this help message\n");

   exit(EXIT_FAILURE);
}
//////////////////////////////////////////////////////

#if LP_TOOLKIT == USE_GLPK
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_glpk(struct_mode_info *mode_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_reac_info *reac_info)
{
   int num_cols;
   g_max_expected_knockout_combinations = INIT_MAX_EXPECTED_KNOCKOUT_COMBINATIONS;

   g_lp_max = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*num_cols;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_pnt_ia = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_ia == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate memory for glpk-memory (row)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_ja = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_ja == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate memory for glpk-memory (column)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_ar = (double *) malloc( (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_pnt_ar == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate memory for glpk-memory (value)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_tmp_cutset = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*sizeof(unsigned long long)));
   if( g_pnt_tmp_cutset == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate memory for temporary cutset\n");
      exit(EXIT_FAILURE);
   }

   g_r = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: init_glpk(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void resize_arrays_glpk(struct_mode_info *mode_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_reac_info *reac_info)
{
   int num_cols;
   g_max_expected_knockout_combinations *= 2;

   g_lp_max = (double *) realloc(g_lp_max,(size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: resize_arrays_glpk(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*num_cols;
   }

   g_pnt_ia = (int *) realloc(g_pnt_ia, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_ia == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for glpk-memory (row)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_ja = (int *) realloc(g_pnt_ja, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_ja == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for glpk-memory (column)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_ar = (double *) realloc(g_pnt_ar, (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_pnt_ar == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for glpk-memory (value)\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) realloc(g_pnt_norms, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_r = (double *) realloc(g_r, (size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void cleanup_glpk()
{
   free(g_lp_max);
   free(g_pnt_ia);
   free(g_pnt_ja);
   free(g_pnt_ar);
   free(g_pnt_norms);
   free(g_pnt_tmp_cutset);
   free(g_r);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void set_glpk_prob(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info *reac_info)
{
   unsigned int i;
   char col_name[100];

   g_lp = glp_create_prob();
   glp_set_prob_name(g_lp, "knockout");
   glp_set_obj_dir(g_lp, GLP_MAX);

   if( g_n_out_of_D == 1 )
   {
      glp_add_cols(g_lp,map_info->num_dupset_keepers + mode_info->good_emfs);

      for( i = 1; i <= map_info->num_dupset_keepers; i++ )
      {
         glp_set_col_name(g_lp, i, reac_info->reactions[map_info->map[map_info->dupsets_map[i-1]]]);
         glp_set_col_kind(g_lp, i, GLP_BV);
         glp_set_col_bnds(g_lp, i, GLP_DB, 0, 1);
         glp_set_obj_coef(g_lp, i, 1);
      }

      for( i = 1; i <= mode_info->good_emfs; i++ )
      {
         sprintf(col_name,"Y%06u",i);
         glp_set_col_name(g_lp, i+map_info->num_dupset_keepers, col_name);
         glp_set_col_kind(g_lp, i+map_info->num_dupset_keepers, GLP_BV);
         glp_set_col_bnds(g_lp, i+map_info->num_dupset_keepers, GLP_DB, 0, 1);
      }
   }
   else
   {
      glp_add_cols(g_lp,map_info->num_dupset_keepers);

      for( i = 1; i <= map_info->num_dupset_keepers; i++ )
      {
         glp_set_col_name(g_lp, i, reac_info->reactions[map_info->map[map_info->dupsets_map[i-1]]]);
         glp_set_col_kind(g_lp, i, GLP_BV);
         glp_set_col_bnds(g_lp, i, GLP_DB, 0, 1);
         glp_set_obj_coef(g_lp, i, 1);
      }
   }

   glp_init_iocp(&g_parm);
   g_parm.presolve = GLP_ON;
   g_parm.binarize = GLP_ON;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void fill_linprog_mem_glpk(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info)
{
   unsigned long i;
   unsigned long lp_row = 0;
   unsigned int j,k;
   int idx;
   unsigned long write_cnt;
   int unit_cnt;
   int bitmover;
   unsigned long long tmp_unit;
   int do_set;
   int norm_mode;
   char mode_name[20];
   int num_cols;
#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

   if( g_n_out_of_D == 1 )
   {
#if DO_SHORTCUT == NO
      glp_add_rows(g_lp,mode_info->bad_efms + 2*mode_info->good_emfs + 1 + 2*g_solutions_cnt);
#else
      glp_add_rows(g_lp,mode_info->bad_efms + 2*mode_info->good_emfs + 1 + g_solutions_cnt + 1);
#endif

      if( g_solutions_cnt == 0 )
      {
         write_cnt = 0;

         //////////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////////
         for( i = 0; i < mode_info->bad_efms; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  // g_pnt_tmp_cutset[mode_info->num_unit_size*sizeof(unsigned long long) + k] = 0;
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               // idx = map_info->map[j];
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               // printf("j=%d, idx=%d: ",j,idx);
               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     // printf("do_set=1: j=%d idx=%d reaction=%s\n",j,idx,reac_info->reactions[idx]);
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
                  write_cnt++;
                  if( write_cnt > g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: linear program storage not big enough!\n");
                     exit(EXIT_FAILURE);
                  }
                  g_pnt_ia[write_cnt] = lp_row;
                  g_pnt_ja[write_cnt] = j+1;
                  g_pnt_ar[write_cnt] = 1.0;
               }
            }

            if( norm_mode == 0 )
            {
               printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
               exit(EXIT_FAILURE);
            }

            g_pnt_norms[i] = norm_mode;
         }
         //////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////
         // and no the >= ||bd||
         //////////////////////////////////////////////////////////////////////
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  // g_pnt_tmp_cutset[mode_info->num_unit_size*sizeof(unsigned long long) + k] = 0;
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               // idx = map_info->map[j];
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               // printf("j=%d, idx=%d: ",j,idx);

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     // printf("do_set=1: j=%d idx=%d reaction=%s\n",j,idx,reac_info->reactions[idx]);
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
                  write_cnt++;
                  if( write_cnt > g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: linear program storage not big enough!\n");
                     exit(EXIT_FAILURE);
                  }
                  g_pnt_ia[write_cnt] = lp_row;
                  g_pnt_ja[write_cnt] = j+1;
                  g_pnt_ar[write_cnt] = 1.0;
               }
            }

            // if( norm_mode == 0 )
            // {
            //    printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
            //    exit(EXIT_FAILURE);
            // }
            g_pnt_norms[i+mode_info->bad_efms]             = norm_mode;

            write_cnt++;
            g_pnt_ia[write_cnt] = lp_row;
            g_pnt_ja[write_cnt] = map_info->num_dupset_keepers + i + 1;
            g_pnt_ar[write_cnt] = norm_mode;

         }
         //////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////
         // and now the <= 2 * ||bd|| - 1
         //////////////////////////////////////////////////////////////////////
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  // g_pnt_tmp_cutset[mode_info->num_unit_size*sizeof(unsigned long long) + k] = 0;
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               // idx = map_info->map[j];
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               // printf("j=%d, idx=%d: ",j,idx);

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     // printf("do_set=1: j=%d idx=%d reaction=%s\n",j,idx,reac_info->reactions[idx]);
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
                  write_cnt++;
                  if( write_cnt > g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: linear program storage not big enough!\n");
                     exit(EXIT_FAILURE);
                  }
                  g_pnt_ia[write_cnt] = lp_row;
                  g_pnt_ja[write_cnt] = j+1;
                  g_pnt_ar[write_cnt] = 1.0;
               }
            }

            printf("norm_mode %lu: %d\n",i,norm_mode);
            // if( norm_mode == 0 )
            // {
            //    printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
            //    exit(EXIT_FAILURE);
            // }
            g_pnt_norms[i+mode_info->bad_efms + mode_info->good_emfs]             = norm_mode;

            write_cnt++;
            g_pnt_ia[write_cnt] = lp_row;
            g_pnt_ja[write_cnt] = map_info->num_dupset_keepers + i + 1;
            g_pnt_ar[write_cnt] = norm_mode;
         }
         //////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////
         // ||y||  <= good_emfs - good_emfs_wanted
         //////////////////////////////////////////////////////////////////////
         g_pnt_norms[mode_info->bad_efms+2*mode_info->good_emfs] = mode_info->good_emfs_orig - mode_info->num_good_removed_by_always_zero - cmd_options->good_efms_wanted;
         lp_row++;
         for( j = 0; j < mode_info->good_emfs; j++ )
         {
            write_cnt++;
            g_pnt_ia[write_cnt] = lp_row;
            g_pnt_ja[write_cnt] = map_info->num_dupset_keepers + j + 1;
            // g_pnt_ar[write_cnt] = 1.0;
            g_pnt_ar[write_cnt] = mode_info->pnt_good_duplicates[j];
         }

         g_num_efm_writes = write_cnt;
      }


      lp_row = 0;

      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         lp_row++;

         sprintf(mode_name,"EFM%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(g_pnt_norms[i] - 1) );
      }

      for( i = 0; i < mode_info->good_emfs; i++ )
      {
         lp_row++;

         sprintf(mode_name,"NOOD%09lu",i);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_LO, (double)(g_pnt_norms[i+mode_info->bad_efms]), 0.0 );
      }

      for( i = 0; i < mode_info->good_emfs; i++ )
      {
         lp_row++;

         sprintf(mode_name,"NOOD%09llu",i+mode_info->good_emfs);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(2*g_pnt_norms[i+mode_info->bad_efms+mode_info->good_emfs] - 1) );
      }

      lp_row++;
      glp_set_row_name(g_lp, lp_row, "Y");
      glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(g_pnt_norms[mode_info->bad_efms+2*mode_info->good_emfs]));

      if( (mode_info->bad_efms + 2*mode_info->good_emfs + 1) != lp_row )
      {
         printf("FATAL ERROR: number of expected lin-program rows (%llu) does is not equal to set rows (%lu)\n",
                mode_info->bad_efms + 2*mode_info->good_emfs + 1,lp_row);
         exit(EXIT_FAILURE);
      }



      write_cnt = g_num_efm_writes;
#if DO_SHORTCUT == NO
      for( i = 0; i < g_solutions_cnt; i++ )
      {
         lp_row++;
         norm_mode = 0;
         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            if( g_r[i*num_cols + j] > 0 )
            {
               norm_mode++;
               write_cnt++;
               if( write_cnt > g_size_linprog_mem )
               {
                  printf("FATAL ERROR: linear program storage not big enough!\n");
                  exit(EXIT_FAILURE);
               }
               g_pnt_ia[write_cnt] = lp_row;
               g_pnt_ja[write_cnt] = j+1;
               g_pnt_ar[write_cnt] = 1.0;
            }
         }
         sprintf(mode_name,"SOLP%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(norm_mode - 1) );
      }
#else
      lp_row++;
      for( j = 0; j < map_info->num_dupset_keepers; j++ )
      {
         write_cnt++;
         if( write_cnt > g_size_linprog_mem )
         {
            printf("FATAL ERROR: linear program storage not big enough!\n");
            exit(EXIT_FAILURE);
         }
         g_pnt_ia[write_cnt] = lp_row;
         g_pnt_ja[write_cnt] = j+1;
         g_pnt_ar[write_cnt] = 1.0;
      }
      sprintf(mode_name,"SHORTCUT%09lu",lp_row);
      glp_set_row_name(g_lp, lp_row, mode_name);
      if( g_solutions_cnt == 0 )
      {
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(map_info->num_dupset_keepers) );
      }
      else
      {
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(g_lp_max[g_solutions_cnt-1]) );
      }
#endif

      for( i = 0; i < g_solutions_cnt; i++ )
      {
         lp_row++;
         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            if( g_r[i*num_cols + j] <= 0 )
            {
               write_cnt++;
               if( write_cnt > g_size_linprog_mem )
               {
                  printf("FATAL ERROR: linear program storage not big enough!\n");
                  exit(EXIT_FAILURE);
               }
               g_pnt_ia[write_cnt] = lp_row;
               g_pnt_ja[write_cnt] = j+1;
               g_pnt_ar[write_cnt] = 1.0;
            }
         }
         sprintf(mode_name,"SOLN%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_LO, 1.0, 0.0 );
      }

      printf("INFO: stored %lu rows in glpk-arrays\n",lp_row);
      printf("INFO: wrote %lu elements to glpk-arrays\n",write_cnt);
      glp_load_matrix(g_lp, write_cnt, g_pnt_ia, g_pnt_ja, g_pnt_ar);
   }
   else  // if n-out-of-D
   {
#if DO_SHORTCUT == NO
      glp_add_rows(g_lp,mode_info->bad_efms + 2*g_solutions_cnt);
#else
      glp_add_rows(g_lp,mode_info->bad_efms + g_solutions_cnt + 1);
#endif

      if( g_solutions_cnt == 0 )
      {
         write_cnt = 0;
         // for( i = 0; i < g_num_efms; i++ )
         for( i = 0; i < mode_info->bad_efms; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  // g_pnt_tmp_cutset[mode_info->num_unit_size*sizeof(unsigned long long) + k] = 0;
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               // idx = map_info->map[j];
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               // printf("j=%d, idx=%d: ",j,idx);

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     // printf("do_set=1: j=%d idx=%d reaction=%s\n",j,idx,reac_info->reactions[idx]);
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
                  write_cnt++;
                  if( write_cnt > g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: linear program storage not big enough!\n");
                     exit(EXIT_FAILURE);
                  }
                  g_pnt_ia[write_cnt] = lp_row;
                  g_pnt_ja[write_cnt] = j+1;
                  g_pnt_ar[write_cnt] = 1.0;
               }
            }

            if( norm_mode == 0 )
            {
               printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
               exit(EXIT_FAILURE);
            }

            g_pnt_norms[i] = norm_mode;
         }
         g_num_efm_writes = write_cnt;
      }


      lp_row = 0;
      // for( i = 0; i < g_num_efms; i++ )
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         lp_row++;

         sprintf(mode_name,"EFM%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(g_pnt_norms[i] - 1) );
      }

      if( mode_info->bad_efms != lp_row )
      {
         printf("FATAL ERROR: number of expected lin-program rows does is not equal to set rows\n");
         exit(EXIT_FAILURE);
      }

      write_cnt = g_num_efm_writes;
#if DO_SHORTCUT == NO
      for( i = 0; i < g_solutions_cnt; i++ )
      {
         lp_row++;
         norm_mode = 0;
         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            if( g_r[i*num_cols + j] > 0 )
            {
               norm_mode++;
               write_cnt++;
               if( write_cnt > g_size_linprog_mem )
               {
                  printf("FATAL ERROR: linear program storage not big enough!\n");
                  exit(EXIT_FAILURE);
               }
               g_pnt_ia[write_cnt] = lp_row;
               g_pnt_ja[write_cnt] = j+1;
               g_pnt_ar[write_cnt] = 1.0;
            }
         }
         sprintf(mode_name,"SOLP%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(norm_mode - 1) );
      }
#else
      lp_row++;
      for( j = 0; j < map_info->num_dupset_keepers; j++ )
      {
         write_cnt++;
         if( write_cnt > g_size_linprog_mem )
         {
            printf("FATAL ERROR: linear program storage not big enough!\n");
            exit(EXIT_FAILURE);
         }
         g_pnt_ia[write_cnt] = lp_row;
         g_pnt_ja[write_cnt] = j+1;
         g_pnt_ar[write_cnt] = 1.0;
      } 
      sprintf(mode_name,"SHORTCUT%09lu",lp_row);
      glp_set_row_name(g_lp, lp_row, mode_name);
      if( g_solutions_cnt == 0 )
      {
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(map_info->num_dupset_keepers) );
      }
      else
      {
         glp_set_row_bnds(g_lp, lp_row, GLP_UP, 0.0, (double)(g_lp_max[g_solutions_cnt-1]) );
      }
#endif

      for( i = 0; i < g_solutions_cnt; i++ )
      {
         lp_row++;
         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            if( g_r[i*num_cols + j] <= 0 )
            {
               write_cnt++;
               if( write_cnt > g_size_linprog_mem )
               {
                  printf("FATAL ERROR: linear program storage not big enough!\n");
                  exit(EXIT_FAILURE);
               }
               g_pnt_ia[write_cnt] = lp_row;
               g_pnt_ja[write_cnt] = j+1;
               g_pnt_ar[write_cnt] = 1.0;
            }
         }
         sprintf(mode_name,"SOLN%09lu",lp_row);
         glp_set_row_name(g_lp, lp_row, mode_name);
         glp_set_row_bnds(g_lp, lp_row, GLP_LO, 1.0, 0.0 );
      }

      printf("INFO: stored %lu rows in glpk-arrays\n",lp_row);
      printf("INFO: wrote %lu elements to glpk-arrays\n",write_cnt);
      glp_load_matrix(g_lp, write_cnt, g_pnt_ia, g_pnt_ja, g_pnt_ar);
   } // if n-out-of-D

#if WRITE_LINPROG_FILE == YES
   sprintf(file_name,"glpk_linprog.%05d.lp",g_solutions_cnt);
   glp_write_lp(g_lp, NULL, file_name);
#endif
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int do_linprog_mip_glpk(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info)
{
   unsigned int i;
   int ret_val;
   int num_cols;
#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

   ret_val = glp_intopt(g_lp, &g_parm);

   if( g_solutions_cnt >= g_max_expected_knockout_combinations )
   {
      resize_arrays_glpk(mode_info, map_info, cmd_options, reac_info);
      // printf("FATAL ERROR: number of solutions (%u) is larger than expected maximum (%u)\n",g_solutions_cnt,g_max_expected_knockout_combinations);
      // printf("             execution aborted.\n");
      // exit(EXIT_FAILURE);
   }

   if( ret_val == 0 )
   {
      if( glp_mip_obj_val(g_lp) == 0 )
      {
         // no optimal solution was found!
         glp_delete_prob(g_lp);
         glp_free_env();
         return(0);
      }

      g_lp_max[g_solutions_cnt] = glp_mip_obj_val(g_lp);
      printf("GLPK: g_solutions_cnt=%d maximum norm=%f, minimum number of knockouts = %f ret_val=%d\n",g_solutions_cnt,g_lp_max[g_solutions_cnt],map_info->num_dupset_keepers - g_lp_max[g_solutions_cnt] + map_info->num_always_zero_reacs,ret_val);
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         g_r[g_solutions_cnt*num_cols + i] = glp_mip_col_val(g_lp, i+1);
         // printf("%s=%3.1f ",reac_info->reactions[map_info->map[map_info->dupsets_map[i]]],g_r[g_solutions_cnt*num_cols + i]);
      }
      // printf("\n");

      // printf("valid unfolded cutset: ");
      // for( i = 0; i < map_info->num_dupset_keepers; i++ )
      // {
      //    if( g_r[g_solutions_cnt*num_cols + i] < 0.5)
      //    {
      //       printf("\"%s\" ",reac_info->reactions[map_info->map[map_info->dupsets_map[i]]]);
      //    }
      // }
      // for( i = 0; i < map_info->num_always_zero_reacs; i++ )
      // {
      //    printf("\"%s\" ",reac_info->reactions[map_info->always_zero_map[i]]);
      // }
      // printf("\n");

   if( (g_solutions_cnt > 1) && (fabs(g_lp_max[g_solutions_cnt - 1]) != fabs(g_lp_max[g_solutions_cnt - 2])) )
   {
      // if old solution uses less knockouts then this solution
      // then we can use 'short' cut
      g_shortcut_possible = 1;
   }  
   else  
   {     
      g_shortcut_possible = 0;
   }

      if( cmd_options->solution_range > 0 && g_solutions_cnt > 0 && !(g_lp_max[g_solutions_cnt] > 0 && (fabs(g_lp_max[g_solutions_cnt] - g_lp_max[0] ) < cmd_options->solution_range) ) )
      {
         printf("INFO: stopping programing: g_lp_max[g_solutions_cnt]=%f g_lp_max[0]=%f solution range=%d\n",g_lp_max[g_solutions_cnt],g_lp_max[0],cmd_options->solution_range);
         glp_delete_prob(g_lp);
         glp_free_env();
         return(0);
      }
      uncompress_and_write_solution(g_solutions_cnt, map_info, reac_info, cmd_options, mode_info);
   }
   else
   {
      printf("FATAL ERROR: glp_intopt() returned with an error!\n");
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

#if WRITE_LINPROG_FILE == YES
   sprintf(file_name,"glpk_mip_efm_outfile.%05d.txt",g_solutions_cnt);
   glp_print_mip(g_lp,file_name);
#endif

   glp_delete_prob(g_lp);
   glp_free_env();
   g_solutions_cnt++;
   return(1);
}
//////////////////////////////////////////////////////
#endif

#if LP_TOOLKIT == USE_GUROBI
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void cleanup_gurobi(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info)
{
   free(g_lp_max);
   free(g_pnt_norms);
   free(g_pnt_tmp_cutset);
   free(g_r);

   // GRBfreemodel(g_gurobi_model);
   // GRBfreeenv(g_gurobi_masterenv);
   free(g_gurobi_col_is_integer);
   free(g_gurobi_obj_cf);
   free(g_gurobi_start);
   free(g_gurobi_len);
   free(g_gurobi_index);
   free(g_gurobi_value);
   free(g_gurobi_row_constr_type);
   free(g_gurobi_row_righthandside);

   int i;
   unsigned long int num_cols;
   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }
   for( i = 0; i < num_cols; i++ )
   {
      free(g_gurobi_varnames[i]);
   }
   free(g_gurobi_varnames);

   GRBfreemodel(g_gurobi_model);
   GRBfreeenv(g_gurobi_masterenv);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_gurobi(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations = INIT_MAX_EXPECTED_KNOCKOUT_COMBINATIONS;

   g_lp_max = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: init_gurobi(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   g_pnt_tmp_cutset = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*sizeof(unsigned long long)));
   if( g_pnt_tmp_cutset == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for temporary cutset\n");
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf("map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu num_cols=%lu\n",map_info->num_dupset_keepers,mode_info->good_emfs,num_cols);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf("mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf("g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: init_gurobi(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_gurobi_varnames = (char **) malloc((size_t) sizeof(char *)*num_cols);

   if( g_gurobi_varnames == NULL )
   {
      printf("FATAL ERROR: init_gurobi(): couldn't allocate %lu bytes for *g_gurobi_varnames\n",sizeof(char *)*num_cols);
      exit(EXIT_FAILURE);
   }

   g_gurobi_col_is_integer = (char *) malloc( (size_t) (num_cols*sizeof(char)));
   if( g_gurobi_col_is_integer == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_col_is_integer array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_obj_cf = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_gurobi_obj_cf == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_obj_cf array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_start = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_gurobi_start == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_start array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_len = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_gurobi_len == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_len array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_index = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_gurobi_index == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_gurobi_value = (double *) malloc( (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_gurobi_value == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_row_constr_type = (char *) malloc( (size_t) (num_rows*sizeof(char)));
   if( g_gurobi_row_constr_type == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_row_righthandside = (double *) malloc( (size_t) (num_rows*sizeof(double)));
   if( g_gurobi_row_righthandside == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_gurobi_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }

   int i;
   for( i = 0; i < num_cols; i++ )
   {
      printf("allocate memory for g_gurobi_varnames[%d]: %lu bytes\n",i,(reac_info->max_len_reac_name+1)*sizeof(char));
      g_gurobi_varnames[i] = (char *) malloc( (size_t) ((reac_info->max_len_reac_name+1)*sizeof(char)));
      if( g_gurobi_varnames[i] == NULL )
      {
         printf("FATAL ERROR: couldn't allocate memory for g_gurobi_varnames[%d] array\n",i);
         exit(EXIT_FAILURE);
      }
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void resize_arrays_gurobi(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations *= 2;

   g_lp_max = (double *) realloc(g_lp_max, (size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf(" resize_arrays_gurobi():map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu\n",map_info->num_dupset_keepers,mode_info->good_emfs);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf(" resize_arrays_gurobi():mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf(" resize_arrays_gurobi():g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) realloc(g_r, (size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_gurobi_index = (int *) realloc( g_gurobi_index, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_gurobi_index == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate memory for g_gurobi_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_gurobi_value = (double *) realloc( g_gurobi_value, (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_gurobi_value == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate memory for g_gurobi_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) realloc(g_pnt_norms, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_row_constr_type = (char *) realloc(g_gurobi_row_constr_type, (size_t) (num_rows*sizeof(char)));
   if( g_gurobi_row_constr_type == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate memory for g_gurobi_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_gurobi_row_righthandside = (double *) realloc(g_gurobi_row_righthandside, (size_t) (num_rows*sizeof(double)));
   if( g_gurobi_row_righthandside == NULL )
   {
      printf("FATAL ERROR: resize_arrays_gurobi(): couldn't allocate memory for g_gurobi_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void fill_linprog_mem_gurobi(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info)
{
   unsigned long i;
   unsigned long lp_row = 0;
   unsigned int j,k;
   int idx;
   unsigned long write_cnt;
   int unit_cnt;
   int bitmover;
   unsigned long long tmp_unit;
   int do_set;
   int norm_mode;
   int gurobi_ret;
   char name_neg[300];
#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif



   printf("DEBUG: fill_linprog_mem_gurobi(): entered. g_shortcut_possible=%d\n",g_shortcut_possible);
   if( (g_solutions_cnt == 0) || ((DO_SHORTCUT == YES) && g_shortcut_possible == 1) )
   {
      if( g_solutions_cnt > 0 )
      {
         // delete old context
         GRBfreemodel(g_gurobi_model);
         GRBfreeenv(g_gurobi_masterenv);
      }

      gurobi_ret = GRBloadenv( &g_gurobi_masterenv, NULL);
      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: couldn't initiate gurobi environment. Returned error code = %d\n",gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      /////////////////////////////////////////////////////////////////////////
      // define all columns
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         g_gurobi_obj_cf[i] = -1.0;
         g_gurobi_col_is_integer[i] = GRB_BINARY;
         // printf("DEBUG: filling data to g_gurobi_varnames[%lu]: '%s'\n",i,reac_info->reactions[map_info->map[map_info->dupsets_map[i]]]);
         strcpy(g_gurobi_varnames[i],reac_info->reactions[map_info->map[map_info->dupsets_map[i]]]);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // define all extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            g_gurobi_obj_cf[i+map_info->num_dupset_keepers] = 0.0;
            g_gurobi_col_is_integer[i+map_info->num_dupset_keepers] = GRB_BINARY;
            sprintf(name_neg,"Y%06lu",i);
            strcpy(g_gurobi_varnames[i+map_info->num_dupset_keepers],name_neg);
         }
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // set contraint type
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         g_gurobi_row_constr_type[i]  = GRB_LESS_EQUAL;
      }

      /////////////////////////////////////////////////////////////////////////
      // set all extra constraint types because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            g_gurobi_row_constr_type[i + mode_info->bad_efms]  = GRB_GREATER_EQUAL;
         }
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            g_gurobi_row_constr_type[i + mode_info->bad_efms + mode_info->good_emfs]  = GRB_LESS_EQUAL;
         }
         // extra constraint for ||y|| <= |D| - n
         // which is equal to ||y|| <= g_good_emfs - g_good_emfs_wanted
         g_gurobi_row_constr_type[2*mode_info->good_emfs + mode_info->bad_efms]  = GRB_LESS_EQUAL;
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // calculate norm of all bad modes
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         norm_mode = 0;
         lp_row++;

         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            // reset temporary cutset array
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               g_pnt_tmp_cutset[k] = 0;
            }

            tmp_unit = 1;
            idx = j;

            bitmover = idx;

            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);

            tmp_unit <<= bitmover;

            // fill temporary cutset array with data
            g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               norm_mode++;
            }
         }

         if( norm_mode == 0 )
         {
            printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
            exit(EXIT_FAILURE);
         }

         g_pnt_norms[i] = norm_mode;
         g_gurobi_row_righthandside[i] = norm_mode - 1;
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // get the norm of the all good modes because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               // reset temporary cutset array
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;

               // fill temporary cutset array with data
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
               }
            }

            g_pnt_norms[i+mode_info->bad_efms]             = norm_mode;
            g_pnt_norms[i+mode_info->bad_efms+mode_info->good_emfs] = norm_mode;
            g_gurobi_row_righthandside[i+mode_info->bad_efms]             = norm_mode;
            g_gurobi_row_righthandside[i+mode_info->bad_efms+mode_info->good_emfs] = (2*norm_mode) - 1;
         }
         // extra constraint ||y|| <= g_good_emfs - g_good_emfs_wanted
         // g_gurobi_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = g_good_emfs_orig - g_good_emfs_wanted;
         g_gurobi_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = mode_info->good_emfs_orig - mode_info->num_good_removed_by_always_zero - cmd_options->good_efms_wanted;
      }
      /////////////////////////////////////////////////////////////////////////
      printf("number of constraints: %llu\n mode_info->bad_efms=%llu mode_info->good_emfs=%llu\n",
              mode_info->bad_efms+mode_info->good_emfs*2,mode_info->bad_efms,mode_info->good_emfs);

      printf("DEBUG: fill_linprog_mem_gurobi(): after norm calculation\n");

      /////////////////////////////////////////////////////////////////////////
      // write linear programming problem column by column
      /////////////////////////////////////////////////////////////////////////
      write_cnt = 0;
      for( j = 0; j < map_info->num_dupset_keepers; j++ )
      {
         g_gurobi_start[j] = write_cnt;
         for( k = 0; k < mode_info->num_unit_size; k++ )
         {
            g_pnt_tmp_cutset[k] = 0;
         }

         tmp_unit = 1;
         idx = j;

         bitmover = idx;

         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);

         tmp_unit <<= bitmover;
         g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

         for( i = 0; i < mode_info->bad_efms; i++ )
         {

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               if( write_cnt >= g_size_linprog_mem )
               {
                  printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                  printf("             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               g_gurobi_index[write_cnt] = i;
               g_gurobi_value[write_cnt] = 1.0;
               write_cnt++;
            }
         }

         //////////////////////////////////////////////////////////////////////
         // write column value for extra constraint because of n-out-of-D 
         //////////////////////////////////////////////////////////////////////
         if( g_n_out_of_D == 1 )
         {
            ///////////////////////////////////////////////////////////////////
            // do it for the >= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_gurobi_index[write_cnt] = i + mode_info->bad_efms;
                  g_gurobi_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
            ///////////////////////////////////////////////////////////////////
            // do it for the <= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_gurobi_index[write_cnt] = i + mode_info->bad_efms + mode_info->good_emfs;
                  g_gurobi_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
         }
         g_gurobi_len[j] = write_cnt - g_gurobi_start[j];
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // write extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( j = 0; j < mode_info->good_emfs; j++ )
         {
            g_gurobi_start[j+map_info->num_dupset_keepers] = write_cnt;

            g_gurobi_index[write_cnt] = j + mode_info->bad_efms;
            g_gurobi_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_gurobi_index[write_cnt] = j + mode_info->bad_efms + mode_info->good_emfs;
            g_gurobi_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_gurobi_index[write_cnt] = mode_info->bad_efms + 2*mode_info->good_emfs;
            // g_gurobi_value[write_cnt] = 1;
            g_gurobi_value[write_cnt] = (float)mode_info->pnt_good_duplicates[j];
            write_cnt++;

            g_gurobi_len[j+map_info->num_dupset_keepers] = write_cnt - g_gurobi_start[j+map_info->num_dupset_keepers];
         }

      }
      /////////////////////////////////////////////////////////////////////////


      g_num_efm_writes = write_cnt;
      printf("DEBUG: wrote %lu index and value entries to arrays\n",write_cnt);

      printf("DEBUG: fill_linprog_mem_gurobi(): after filling index/value values\n");

      if( g_n_out_of_D == 1 )
      {
         gurobi_ret = GRBloadmodel(g_gurobi_masterenv, &g_gurobi_model, "knockout", map_info->num_dupset_keepers + mode_info->good_emfs,
                                   mode_info->bad_efms + 2*mode_info->good_emfs + 1, FIND_MINIMUM, 0.0,
                                   g_gurobi_obj_cf, g_gurobi_row_constr_type, g_gurobi_row_righthandside, g_gurobi_start, g_gurobi_len,
                                   g_gurobi_index, g_gurobi_value, NULL, NULL, g_gurobi_col_is_integer, g_gurobi_varnames, NULL);
      }
      else
      {
         gurobi_ret = GRBloadmodel(g_gurobi_masterenv, &g_gurobi_model, "knockout", map_info->num_dupset_keepers, mode_info->bad_efms, FIND_MINIMUM, 0.0,
                                   g_gurobi_obj_cf, g_gurobi_row_constr_type, g_gurobi_row_righthandside, g_gurobi_start, g_gurobi_len,
                                   g_gurobi_index, g_gurobi_value, NULL, NULL, g_gurobi_col_is_integer, g_gurobi_varnames, NULL);
      }

      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: GRBloadmodel returned abnormally: return code = %d\n",gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      printf("DEBUG: after GRBloadmodel()\n");

      g_gurobi_modelenv = GRBgetenv(g_gurobi_model);
      if( g_gurobi_modelenv == NULL )
      {
         printf("FATAL ERROR: couldn't get model environment.\n");
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      gurobi_ret = GRBsetintparam(g_gurobi_modelenv, "THREADS",  cmd_options->num_threads);
      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: couldn't set number of gurobi threads to %d. Returned error code = %d\n",cmd_options->num_threads,gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("INFO: going to set (absolute) MIPGapAbs to %e\n", LINPROG_MIP_GAP_ABS);
      gurobi_ret = GRBsetdblparam(g_gurobi_modelenv, "MIPGapAbs", LINPROG_MIP_GAP_ABS);
      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: couldn't set (absolute) MIPGapAbs to '%e'. Returned error code = %d\n", LINPROG_MIP_GAP_ABS,gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("INFO: going to set (relative) MIPGap to %e\n", LINPROG_MIP_GAP_REL);
      gurobi_ret = GRBsetdblparam(g_gurobi_modelenv, "MIPGap", LINPROG_MIP_GAP_REL);
      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: couldn't set (relative) MIPGap to '%e'. Returned error code = %d\n", LINPROG_MIP_GAP_REL,gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      // gurobi_ret = GRBsetintparam(g_gurobi_modelenv, "PRESOLVE", 0);  // switch off presolving
      // gurobi_ret = GRBsetintparam(g_gurobi_modelenv, "MIPFOCUS", 0);  // set solving strategy
      // gurobi_ret = GRBsetintparam(g_gurobi_modelenv, "METHOD", 4);  // set solving strategy
      // if( gurobi_ret != 0 )
      // {
      //    printf("FATAL ERROR: couldn't set solve method. Returned error code = %d\n",gurobi_ret);
      //    printf("             execution aborted.\n");
      //    exit(EXIT_FAILURE);
      // }

      // write linear programming problem to file
#if WRITE_LINPROG_FILE == YES
      printf("INFO: going to write 'gurobi_linprog_refill.%05d.lp'",g_solutions_cnt);
      sprintf(file_name,"gurobi_linprog_refill.%05d.lp",g_solutions_cnt);
      GRBwrite( g_gurobi_model, file_name);
#endif

   }
}
//////////////////////////////////////////////////////



//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int do_linprog_mip_gurobi(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info)
{
   int gurobi_ret;
   unsigned int i;
   int *indices_pos;
   double *elements_pos;
   unsigned cnt_pos = 0;
   int *indices_neg;
   double *elements_neg;
   unsigned cnt_neg = 0;
   int optimstatus;
   int solcount;
   int s;
   int num_cols;

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif

#if DO_SHORTCUT == NO
   char name_pos[300];
#else
   int *indices_shortcut;
   double *elements_shortcut;
   int cnt_shortcut = 0;

   if( (indices_shortcut = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for indices_shortcut\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_shortcut = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for elements_shortcut\n");
      exit(EXIT_FAILURE);
   }
#endif

   char name_neg[300];

   if( (indices_pos = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for indices_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (indices_neg = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for indices_neg\n");
      exit(EXIT_FAILURE);
   }
   if( (elements_pos = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for elements_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_neg = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_gurobi(): couldn't locate memory for elements_neg\n");
      exit(EXIT_FAILURE);
   }

   if( (g_solutions_cnt > 0) )
   {
#if DO_SHORTCUT == YES
      if( g_shortcut_possible == 1 )
      {
         for( s = 0; s < g_solutions_cnt_last; s++ )
         {
            cnt_neg = 0;
            for( i = 0; i < map_info->num_dupset_keepers; i++ )
            {
               if( g_r[s*num_cols + i] <= 0.5 )
               {
                  indices_neg[cnt_neg] = i;
                  elements_neg[cnt_neg] = 1.0;
                  cnt_neg++;
               }
            }

            sprintf(name_neg,"NEWNEG_RESET%07d",s);
            gurobi_ret = GRBaddconstr(g_gurobi_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);

            if( gurobi_ret != 0 )
            {
               printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",gurobi_ret);
               printf("             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
         }

         cnt_shortcut = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            indices_shortcut[cnt_shortcut] = i;
            elements_shortcut[cnt_shortcut] = 1.0;
            cnt_shortcut++;
         }
         gurobi_ret = GRBaddconstr(g_gurobi_model, cnt_shortcut, indices_shortcut, elements_shortcut, GRB_LESS_EQUAL, g_lp_max[g_solutions_cnt-1],"SHORTCUT");

         if( gurobi_ret != 0 )
         {
            printf("FATAL ERROR: Adding shortcut constrain to model failed. Return code = %d\n",gurobi_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
#endif

      for( s = g_solutions_cnt_last; s < g_solutions_cnt; s++ )
      {
         cnt_pos = 0;
         cnt_neg = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            if( g_r[s*num_cols + i] > 0.5 )
            {
               indices_pos[cnt_pos] = i;
               elements_pos[cnt_pos] = 1.0;
               cnt_pos++;
            }
            else
            {
               indices_neg[cnt_neg] = i;
               elements_neg[cnt_neg] = 1.0;
               cnt_neg++;
            }
         }

         if( cnt_pos != (unsigned int)g_lp_max[s] )
         {
            printf("FATAL ERROR: number of elements of new row (%u) are not identical to norm of row (%u) for s=%d\n",cnt_pos,(unsigned int)g_lp_max[s],s);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

#if DO_SHORTCUT == NO
         sprintf(name_pos,"NEWPOS%07d",g_solutions_cnt);
         gurobi_ret = GRBaddconstr(g_gurobi_model, cnt_pos, indices_pos, elements_pos, GRB_LESS_EQUAL, (g_lp_max[s] - 1.0),name_pos);

         if( gurobi_ret != 0 )
         {
            printf("FATAL ERROR: Adding positive constrain to model failed. Return code = %d\n",gurobi_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
#endif
         sprintf(name_neg,"NEWNEG%07d",g_solutions_cnt);
         gurobi_ret = GRBaddconstr(g_gurobi_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);

         if( gurobi_ret != 0 )
         {
            printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",gurobi_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
   }

   GRBupdatemodel( g_gurobi_model );

#if WRITE_LINPROG_FILE == YES
   printf("INFO: going to write 'gurobi_linprog_addconst.%05d.lp'",g_solutions_cnt);
   sprintf(file_name,"gurobi_linprog_addconst.%05d.lp",g_solutions_cnt);
   GRBwrite( g_gurobi_model, file_name);
   GRBwriteparams( g_gurobi_modelenv, "gurobi_params");
#endif

   g_solutions_cnt_last = g_solutions_cnt;

   gurobi_ret = GRBoptimize(g_gurobi_model);

   if( gurobi_ret != 0 )
   {
      printf("FATAL ERROR: GRBoptimize failed. Return code = %d\n",gurobi_ret);
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   if( g_solutions_cnt >= g_max_expected_knockout_combinations )
   {
      resize_arrays_gurobi(mode_info, reac_info, cmd_options, map_info);
      // printf("FATAL ERROR: number of solutions (%u) is larger than expected maximum (%u)\n",g_solutions_cnt,g_max_expected_knockout_combinations);
      // printf("             execution aborted.\n");
      // exit(EXIT_FAILURE);
   }

   gurobi_ret = GRBgetintattr(g_gurobi_model, GRB_INT_ATTR_STATUS, &optimstatus);

   if( gurobi_ret != 0 )
   {
      printf("FATAL ERROR: could retrieve status of optimization. Return code = %d\n",gurobi_ret);
      printf("             execution aborted.\n");
      GRBfreemodel(g_gurobi_model);
      exit(EXIT_FAILURE);
   }

#if DO_SHORTCUT != NO
   free(indices_shortcut);
   free(elements_shortcut);
#endif
   free(indices_pos);
   free(elements_pos);
   free(indices_neg);
   free(elements_neg);

   if (optimstatus != GRB_OPTIMAL)
   {
      printf("INFO: Gurobi could not find optimal solution. Optimization status = %d\n",optimstatus);
      return(0);
   }


   gurobi_ret = GRBgetdblattr(g_gurobi_model, GRB_DBL_ATTR_OBJVAL, &g_lp_max[g_solutions_cnt]);

   if( gurobi_ret != 0 )
   {
      printf("FATAL ERROR: could retrieve status of optimization. Return code = %d\n",gurobi_ret);
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }
   g_lp_max[g_solutions_cnt] *= -1;

   printf("INFO: solution calculated by gurobi = %e\n",g_lp_max[g_solutions_cnt]);

   g_lp_max[g_solutions_cnt] = round(g_lp_max[g_solutions_cnt]);

   gurobi_ret = GRBgetintattr(g_gurobi_model, "SolCount", &solcount);
   if( gurobi_ret != 0 )
   {
      printf("FATAL ERROR: could retrieve solution cnt. Return code = %d\n",gurobi_ret);
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }
   printf("INFO: number of solutions found by gurobi = %d\n",solcount);

   // Attentention: we only grab 1 solution even if there would be more
   for( s = 0; s < 1; s++ )
   {
      gurobi_ret = GRBsetintparam(g_gurobi_modelenv, "SolutionNumber", s);
      if( gurobi_ret != 0 )
      {
         printf("FATAL ERROR: could set parameter \"SolutionNumber\" to %d. Return code = %d\n",s,gurobi_ret);
         printf("             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         gurobi_ret = GRBgetdblattrelement(g_gurobi_model, "Xn", i, &g_r[g_solutions_cnt*num_cols + i]);
         if( gurobi_ret != 0 )
         {
            printf("FATAL ERROR: Failed to obtain solution for reaction %d. Return code = %d\n",i,gurobi_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }


      printf("GUROBI: g_solutions_cnt=%d maximum norm=%f, minimum number of knockouts = %f\n",
             g_solutions_cnt,g_lp_max[g_solutions_cnt],map_info->num_dupset_keepers - g_lp_max[g_solutions_cnt] + map_info->num_always_zero_reacs);
      int norm = 0;
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         // printf("%s=%3.1f ",reac_info->reactions[map_info->map[map_info->dupsets_map[i]]],g_r[g_solutions_cnt*num_cols + i]);
         if( g_r[g_solutions_cnt*num_cols + i] > 0.5 )
         {
            norm++;
         }
      }
      printf("\n");

      if( fabs(g_lp_max[g_solutions_cnt] - (float)norm) > 1e-7 )
      {
         printf("FATAL ERROR: gurobi's maximum %d is not equal to the norm of the result vector (%d)\n",(int) g_lp_max[g_solutions_cnt],norm);
         printf("             ");
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            printf("col %d: val %f ",i,g_r[g_solutions_cnt*num_cols + i]);
         }
         printf("\n");
         exit(EXIT_FAILURE);
      }

      // uncompress last solution and write them all to file if required
      if( cmd_options->solution_range > 0 && g_solutions_cnt > 0 && !(g_lp_max[g_solutions_cnt] > 0 && (fabs(g_lp_max[g_solutions_cnt] - g_lp_max[0] ) < cmd_options->solution_range) ) )
      {
         printf("INFO: stopping programing: g_lp_max[0]=%f g_lp_max[g_solutions_cnt]=%f solution_range=%d\n",g_lp_max[0],g_lp_max[g_solutions_cnt],cmd_options->solution_range);
         return(0);
      }
      uncompress_and_write_solution(g_solutions_cnt, map_info, reac_info, cmd_options, mode_info);

      g_solutions_cnt++;
   }

   if( (g_solutions_cnt > 1) && (g_lp_max[g_solutions_cnt - 1] != g_lp_max[g_solutions_cnt - 2]) )
   {
      // if old solution uses less knockouts then this solution
      // then we can use 'short' cut
      g_shortcut_possible = 1;
   }
   else
   {
      g_shortcut_possible = 0;
   }

   return(1);
}
//////////////////////////////////////////////////////
#endif

#if LP_TOOLKIT == USE_CPLEX
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_g_cplex_arrays(struct_map_info *map_info, struct_mode_info *mode_info, struct_cmd_options *cmd_options, unsigned long long int write_cnt)
{
   unsigned int c;
   unsigned long long int r;

   for( r = 0; r < mode_info->bad_efms; r++ )
   {
      print_mode(&mode_info->pnt_flux[r], mode_info->num_unit_size);
      printf(" norm of mode %llu: %d\n",r,g_pnt_norms[r]);
   }

   printf("pointer to g_cplex_env=%p\n",g_cplex_env);
   printf("pointer to g_cplex_lp=%p\n",g_cplex_lp);

   if( g_n_out_of_D == 1 )
   {
      printf("g_solutions_cnt=%d\n",g_solutions_cnt);
      for( c = 0; c < map_info->num_dupset_keepers + mode_info->good_emfs; c++ )
      {
         printf("cplex variable name %u: %s ",c,g_cplex_varnames[c]);
         if( g_cplex_col_is_integer[c] == CPX_BINARY )
         {
            printf("BINARY\n");
         }
         else if( g_cplex_col_is_integer[c] == CPX_CONTINUOUS )
         {
            printf("GRB_CONTINUOUS\n");
         }
         else
         {
            printf("UNKNOWN!!!!!!\n");
         }
         printf("objective coeff %u: %g\n",c,g_cplex_obj_cf[c]);
         printf("lower bound %u: %g\n",c,g_cplex_lb[c]);
         printf("upper bound %u: %g\n",c,g_cplex_ub[c]);
      }
      for( r = 0; r < mode_info->bad_efms + g_solutions_cnt + 2*mode_info->good_emfs + 1; r++ )
      {
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("constraint type %llu: EQUAL ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("constraint type %llu: LESS ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf("constraint type %llu: GREATER ",r);
         }
         else
         {
            printf("constraint type %llu: UNkNOWN ",r);
         }
         // printf("right_hand_side %llu: %g\n",r,g_cplex_row_righthandside[r]);
         printf(" %g\n",g_cplex_row_righthandside[r]);
      }
      for( c = 0; c < map_info->num_dupset_keepers + mode_info->good_emfs; c++ )
      {
         printf("g_cplex_start[%d]=%d g_cplex_len[%d]=%d\n",c,g_cplex_start[c],c,g_cplex_len[c]);
      }
   }
   else
   {
      printf("g_solutions_cnt=%d\n",g_solutions_cnt);
      for( c = 0; c < map_info->num_dupset_keepers; c++ )
      {
         printf("cplex variable name %u: %s ",c,g_cplex_varnames[c]);
         if( g_cplex_col_is_integer[c] == CPX_BINARY )
         {
            printf("BINARY\n");
         }
         else if( g_cplex_col_is_integer[c] == CPX_CONTINUOUS )
         {
            printf("GRB_CONTINUOUS\n");
         }
         else
         {
            printf("UNKNOWN!!!!!!\n");
         }
         printf("objective coeff %u: %g\n",c,g_cplex_obj_cf[c]);
         printf("lower bound %u: %g\n",c,g_cplex_lb[c]);
         printf("upper bound %u: %g\n",c,g_cplex_ub[c]);
         printf("variablename %u: %s\n",c,g_cplex_varnames[c]);
      }
      for( r = 0; r < mode_info->bad_efms + g_solutions_cnt; r++ )
      {
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("constraint type %llu: EQUAL ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("constraint type %llu: LESS ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf("constraint type %llu: GREATER ",r);
         }
         else
         {
            printf("constraint type %llu: UNkNOWN ",r);
         }
         // printf("right_hand_side %llu: %g\n",r,g_cplex_row_righthandside[r]);
         printf(" %g\n",g_cplex_row_righthandside[r]);
      }
      for( c = 0; c < map_info->num_dupset_keepers; c++ )
      {
         printf("g_cplex_start[%d]=%d g_cplex_len[%d]=%d g_cplex_index[%d]=%d g_cplex_value[%d]=%g\n",
                 c,g_cplex_start[c],c,g_cplex_len[c],c,g_cplex_index[c],c,g_cplex_value[c]);
      }

      /////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////
      for( r = 0; r < mode_info->bad_efms; r++ )
      {
         int start = g_cplex_start_trans[r];
         int stop;
         if( r == mode_info->bad_efms - 1 )
         {
            stop = write_cnt;
         }
         else
         {
            stop = g_cplex_start_trans[r+1];
         }
         printf("r=%llu start=%x stop=%d\n",r,start,stop);
         for( c = start; c < stop; c++ )
         {
            printf("%g*%d:%s ",g_cplex_value_trans[c],g_cplex_index_trans[c],g_cplex_varnames[g_cplex_index_trans[c]]);
         }
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("== %g\n",g_cplex_row_righthandside[r]);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("<= %g\n",g_cplex_row_righthandside[r]);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf(">= %g\n",g_cplex_row_righthandside[r]);
         }
         else
         {
            printf("UNKNOWN %g\n",g_cplex_row_righthandside[r]);
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////
   }

   // print transpose information
   printf("DEBUG: leaving print_g_cplex_arrays().\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void cplex_transpose_arrays(unsigned long long int num_idx, unsigned long long int num_rows, unsigned long long int num_cols)
{
   unsigned long long int r, c, x;
   unsigned long long int new_idx = 0;

   printf("DEBUG: cplex_transpose_arrays(): entered.\n");

   g_cplex_start_trans = (int *) malloc( (size_t) sizeof(int)*num_rows );

   if( g_cplex_start_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_start_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_index_trans = (int *) malloc( (size_t) sizeof(int)*num_idx );

   if( g_cplex_index_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_index_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_value_trans = (double *) malloc( (size_t) sizeof(double)*num_idx );

   if( g_cplex_value_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_value_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   for( r = 0; r < num_rows; r++ )
   {
      g_cplex_start_trans[r] = new_idx;

      for( c = 0; c < num_cols; c++ )
      {
         for( x = g_cplex_start[c]; x < g_cplex_start[c] + g_cplex_len[c]; x++ )
         {
            if( r == g_cplex_index[x] )
            {
               g_cplex_value_trans[new_idx] = g_cplex_value[x];
               g_cplex_index_trans[new_idx] = c;
               new_idx++;
            }
         }
      }
   }

   if( new_idx != num_idx )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): transposing arrays failed: new_idx=%llu old_idx=%llu\n",new_idx, num_idx );
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   printf("DEBUG: cplex_transpose_arrays(): leaving. transpose %llu values\n", new_idx);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void cplex_clean_transpose()
{
   free(g_cplex_start_trans);
   free(g_cplex_index_trans);
   free(g_cplex_value_trans);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void cleanup_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info)
{
   free(g_lp_max);
   free(g_pnt_norms);
   free(g_pnt_tmp_cutset);
   free(g_r);

   free(g_cplex_col_is_integer);
   free(g_cplex_obj_cf);
   free(g_cplex_lb);
   free(g_cplex_ub);
   free(g_cplex_start);
   free(g_cplex_len);
   free(g_cplex_index);
   free(g_cplex_value);
   free(g_cplex_row_constr_type);
   free(g_cplex_row_righthandside);

   int i;
   unsigned long int num_cols;
   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }
   for( i = 0; i < num_cols; i++ )
   {
      free(g_cplex_varnames[i]);
   }
   free(g_cplex_varnames);

   printf("DEBUG: cleanup_cplex() deleting old cplex prob\n");
   CPXfreeprob (g_cplex_env, &g_cplex_lp);
   printf("DEBUG: cleanup_cplex(): deleting old cplex env\n");
   CPXcloseCPLEX (&g_cplex_env);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations = INIT_MAX_EXPECTED_KNOCKOUT_COMBINATIONS;

   g_lp_max = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   g_pnt_tmp_cutset = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*sizeof(unsigned long long)));
   if( g_pnt_tmp_cutset == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for temporary cutset\n");
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf("map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu num_cols=%lu\n",map_info->num_dupset_keepers,mode_info->good_emfs,num_cols);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf("mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf("g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_varnames = (char **) malloc((size_t) sizeof(char *)*num_cols);

   if( g_cplex_varnames == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for *g_cplex_varnames\n",sizeof(char *)*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_col_is_integer = (char *) malloc( (size_t) (num_cols*sizeof(char)));
   if( g_cplex_col_is_integer == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_col_is_integer array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_obj_cf = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_obj_cf == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_obj_cf array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_lb = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_lb == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_lb array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_ub = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_ub == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_ub array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_start = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_cplex_start == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_start array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_len = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_cplex_len == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_len array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_index = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_cplex_index == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_cplex_value = (double *) malloc( (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_cplex_value == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_constr_type = (char *) malloc( (size_t) (num_rows*sizeof(char)));
   if( g_cplex_row_constr_type == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_righthandside = (double *) malloc( (size_t) (num_rows*sizeof(double)));
   if( g_cplex_row_righthandside == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }

   int i;
   for( i = 0; i < num_cols; i++ )
   {
      printf("allocate memory for g_cplex_varnames[%d]: %lu bytes\n",i,(reac_info->max_len_reac_name+1)*sizeof(char));
      g_cplex_varnames[i] = (char *) malloc( (size_t) ((reac_info->max_len_reac_name+1)*sizeof(char)));
      if( g_cplex_varnames[i] == NULL )
      {
         printf("FATAL ERROR: couldn't allocate memory for g_cplex_varnames[%d] array\n",i);
         exit(EXIT_FAILURE);
      }
   }
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void resize_arrays_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations *= 2;

   g_lp_max = (double *) realloc(g_lp_max, (size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf(" resize_arrays_cplex():map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu\n",map_info->num_dupset_keepers,mode_info->good_emfs);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf(" resize_arrays_cplex():mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf(" resize_arrays_cplex():g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) realloc(g_r, (size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_index = (int *) realloc( g_cplex_index, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_cplex_index == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_cplex_value = (double *) realloc( g_cplex_value, (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_cplex_value == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) realloc(g_pnt_norms, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_constr_type = (char *) realloc(g_cplex_row_constr_type, (size_t) (num_rows*sizeof(char)));
   if( g_cplex_row_constr_type == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_righthandside = (double *) realloc(g_cplex_row_righthandside, (size_t) (num_rows*sizeof(double)));
   if( g_cplex_row_righthandside == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void fill_linprog_mem_cplex(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info)
{
   unsigned long i;
   unsigned long lp_row = 0;
   unsigned int j,k;
   int idx;
   unsigned long write_cnt;
   int unit_cnt;
   int bitmover;
   unsigned long long tmp_unit;
   int do_set;
   int norm_mode;
   int cplex_ret;
   char name_neg[300];
   int status;
#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif



   printf("DEBUG: fill_linprog_mem_cplex(): entered. g_shortcut_possible=%d\n",g_shortcut_possible);
   if( (g_solutions_cnt == 0) || ((DO_SHORTCUT == YES) && g_shortcut_possible == 1) )
   {
      if( g_solutions_cnt > 0 )
      {
         printf("DEBUG: deleting old cplex prob\n");
         // delete old context
         status = CPXfreeprob (g_cplex_env, &g_cplex_lp);
         if( status )
         {
            fprintf(stderr,"FATAL ERROR: couldn't free cplex environment. Returned error code = %d\n",status);
            fprintf(stderr,"             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
         printf("DEBUG: deleting old cplex environment\n");
         status = CPXcloseCPLEX (&g_cplex_env);
         if ( status )
         {
            char errmsg[1024];
            fprintf (stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring (g_cplex_env, status, errmsg);
            fprintf (stderr, "%s", errmsg);
            fprintf(stderr,"             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }

      /////////////////////////////////////////////////////////////////////////
      // get new CPLEX environment
      /////////////////////////////////////////////////////////////////////////
      printf("DEBUG: creating cplex environment\n");
      g_cplex_env = CPXopenCPLEX (&status);

      if ( g_cplex_env == NULL )
      {
         char errmsg[1024];
         fprintf (stderr, "Could not open CPLEX environment.\n");
         CPXgeterrorstring (g_cplex_env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("DEBUG: creating cplex prob\n");
      g_cplex_lp = CPXcreateprob (g_cplex_env, &status, "lpex1");

      if ( g_cplex_lp == NULL )
      {
         fprintf (stderr, "Failed to create CPLEX LP.\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_OPPORTUNISTIC);
      // status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_AUTO);
      // status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_DETERMINISTIC);
      if ( status )
      {
         fprintf (stderr, "Failed to parallel mode type\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_AUXROOTTHREADS, cmd_options->num_threads);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to use %d threads for root.\n", cmd_options->num_threads);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_THREADS, cmd_options->num_threads);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to use %d threads.\n", cmd_options->num_threads);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("INFO: going to set absolute MIP gap to %g\n", LINPROG_MIP_GAP_ABS);
      status = CPXsetdblparam( g_cplex_env, CPX_PARAM_EPAGAP, LINPROG_MIP_GAP_ABS);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to absolute MIP gap tolerance of %e.\n", LINPROG_MIP_GAP_ABS);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("INFO: going to set relative MIP gap to %g\n", LINPROG_MIP_GAP_REL);
      status = CPXsetdblparam( g_cplex_env, CPX_PARAM_EPGAP, LINPROG_MIP_GAP_REL);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to relative MIP gap tolerance of %e.\n",LINPROG_MIP_GAP_REL);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXchgobjsen (g_cplex_env, g_cplex_lp, CPX_MIN);
      if ( status )
      {
         fprintf (stderr, "Failed to tell CPLEX to minimize\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // define all columns
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         g_cplex_obj_cf[i] = -1.0;
         // g_cplex_col_is_integer[i] = GRB_BINARY;
         g_cplex_col_is_integer[i] = CPX_BINARY;
         g_cplex_lb[i] = 0.0;
         g_cplex_ub[i] = 1.0;
         strcpy(g_cplex_varnames[i],reac_info->reactions[map_info->map[map_info->dupsets_map[i]]]);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // define all extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            g_cplex_obj_cf[i+map_info->num_dupset_keepers] = 0.0;
            // g_cplex_col_is_integer[i+map_info->num_dupset_keepers] = GRB_BINARY;
            g_cplex_col_is_integer[i+map_info->num_dupset_keepers] = CPX_BINARY;
            g_cplex_lb[i+map_info->num_dupset_keepers] = 0.0;
            g_cplex_ub[i+map_info->num_dupset_keepers] = 1.0;
            sprintf(name_neg,"Y%06lu",i);
            strcpy(g_cplex_varnames[i+map_info->num_dupset_keepers],name_neg);
         }
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // set contraint type
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         // g_cplex_row_constr_type[i]  = GRB_LESS_EQUAL;
         g_cplex_row_constr_type[i]  = 'L';
      }

      /////////////////////////////////////////////////////////////////////////
      // set all extra constraint types because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            // g_cplex_row_constr_type[i + mode_info->bad_efms]  = GRB_GREATER_EQUAL;
            g_cplex_row_constr_type[i + mode_info->bad_efms]  = 'G';
         }
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            // g_cplex_row_constr_type[i + mode_info->bad_efms + mode_info->good_emfs]  = GRB_LESS_EQUAL;
            g_cplex_row_constr_type[i + mode_info->bad_efms + mode_info->good_emfs]  = 'L';
         }
         // extra constraint for ||y|| <= |D| - n
         // which is equal to ||y|| <= g_good_emfs - g_good_emfs_wanted
         g_cplex_row_constr_type[2*mode_info->good_emfs + mode_info->bad_efms]  = 'L';
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // calculate norm of all bad modes
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         norm_mode = 0;
         lp_row++;

         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            // reset temporary cutset array
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               g_pnt_tmp_cutset[k] = 0;
            }

            tmp_unit = 1;
            idx = j;

            bitmover = idx;

            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);

            tmp_unit <<= bitmover;

            // fill temporary cutset array with data
            g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               norm_mode++;
            }
         }

         if( norm_mode == 0 )
         {
            printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
            exit(EXIT_FAILURE);
         }

         g_pnt_norms[i] = norm_mode;
         g_cplex_row_righthandside[i] = norm_mode - 1;
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // get the norm of the all good modes because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               // reset temporary cutset array
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;

               // fill temporary cutset array with data
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
               }
            }

            g_pnt_norms[i+mode_info->bad_efms]             = norm_mode;
            g_pnt_norms[i+mode_info->bad_efms+mode_info->good_emfs] = norm_mode;
            g_cplex_row_righthandside[i+mode_info->bad_efms]             = norm_mode;
            g_cplex_row_righthandside[i+mode_info->bad_efms+mode_info->good_emfs] = (2*norm_mode) - 1;
         }
         // extra constraint ||y|| <= g_good_emfs - g_good_emfs_wanted
         // g_cplex_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = g_good_emfs_orig - g_good_emfs_wanted;
         g_cplex_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = mode_info->good_emfs_orig - mode_info->num_good_removed_by_always_zero - cmd_options->good_efms_wanted;
      }
      /////////////////////////////////////////////////////////////////////////
      printf("number of constraints: %llu\n mode_info->bad_efms=%llu mode_info->good_emfs=%llu\n",
              mode_info->bad_efms+mode_info->good_emfs*2,mode_info->bad_efms,mode_info->good_emfs);

      printf("DEBUG: fill_linprog_mem_cplex(): after norm calculation\n");

      /////////////////////////////////////////////////////////////////////////
      // write linear programming problem column by column
      /////////////////////////////////////////////////////////////////////////
      write_cnt = 0;
      for( j = 0; j < map_info->num_dupset_keepers; j++ )
      {
         g_cplex_start[j] = write_cnt;
         for( k = 0; k < mode_info->num_unit_size; k++ )
         {
            g_pnt_tmp_cutset[k] = 0;
         }

         tmp_unit = 1;
         idx = j;

         bitmover = idx;

         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);

         tmp_unit <<= bitmover;
         g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

         for( i = 0; i < mode_info->bad_efms; i++ )
         {

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               if( write_cnt >= g_size_linprog_mem )
               {
                  printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                  printf("             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               g_cplex_index[write_cnt] = i;
               g_cplex_value[write_cnt] = 1.0;
               write_cnt++;
            }
         }

         //////////////////////////////////////////////////////////////////////
         // write column value for extra constraint because of n-out-of-D 
         //////////////////////////////////////////////////////////////////////
         if( g_n_out_of_D == 1 )
         {
            ///////////////////////////////////////////////////////////////////
            // do it for the >= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_cplex_index[write_cnt] = i + mode_info->bad_efms;
                  g_cplex_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
            ///////////////////////////////////////////////////////////////////
            // do it for the <= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_cplex_index[write_cnt] = i + mode_info->bad_efms + mode_info->good_emfs;
                  g_cplex_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
         }
         g_cplex_len[j] = write_cnt - g_cplex_start[j];
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // write extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( j = 0; j < mode_info->good_emfs; j++ )
         {
            g_cplex_start[j+map_info->num_dupset_keepers] = write_cnt;

            g_cplex_index[write_cnt] = j + mode_info->bad_efms;
            g_cplex_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_cplex_index[write_cnt] = j + mode_info->bad_efms + mode_info->good_emfs;
            g_cplex_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_cplex_index[write_cnt] = mode_info->bad_efms + 2*mode_info->good_emfs;
            // g_cplex_value[write_cnt] = 1;
            g_cplex_value[write_cnt] = (float)mode_info->pnt_good_duplicates[j];
            write_cnt++;

            g_cplex_len[j+map_info->num_dupset_keepers] = write_cnt - g_cplex_start[j+map_info->num_dupset_keepers];
         }

      }
      /////////////////////////////////////////////////////////////////////////


      g_num_efm_writes = write_cnt;

      printf("DEBUG: fill_linprog_mem_cplex(): after filling index/value values\n");


      if( g_n_out_of_D == 1 )
      {
         cplex_transpose_arrays(write_cnt, mode_info->bad_efms + 2*mode_info->good_emfs + 1, map_info->num_dupset_keepers + mode_info->good_emfs);
         // print_g_cplex_arrays(map_info, mode_info, cmd_options, write_cnt);

         cplex_ret = CPXnewcols(g_cplex_env, g_cplex_lp, map_info->num_dupset_keepers + mode_info->good_emfs, g_cplex_obj_cf, g_cplex_lb, g_cplex_ub, g_cplex_col_is_integer, g_cplex_varnames);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXnewcols returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         cplex_ret = CPXaddrows (g_cplex_env, g_cplex_lp, 0, mode_info->bad_efms + 2*mode_info->good_emfs + 1,
                                 write_cnt, g_cplex_row_righthandside,
                                 g_cplex_row_constr_type, g_cplex_start_trans, g_cplex_index_trans, g_cplex_value_trans, NULL, NULL);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXaddrows returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
      else
      {
         printf("DEBUG: number of rows going to be added %llu\n",mode_info->bad_efms);
         printf("DEBUG: wrote %lu index and value entries to arrays\n",write_cnt);
         cplex_transpose_arrays(write_cnt, mode_info->bad_efms, map_info->num_dupset_keepers);
         // print_g_cplex_arrays(map_info, mode_info, cmd_options, write_cnt);

         cplex_ret = CPXnewcols(g_cplex_env, g_cplex_lp, map_info->num_dupset_keepers, g_cplex_obj_cf, g_cplex_lb, g_cplex_ub, g_cplex_col_is_integer, g_cplex_varnames);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXnewcols returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, mode_info->bad_efms, write_cnt, g_cplex_row_righthandside,
                                g_cplex_row_constr_type, g_cplex_start_trans, g_cplex_index_trans, g_cplex_value_trans, NULL, NULL);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXaddrows returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }

      cplex_clean_transpose();


      // write linear programming problem to file
#if WRITE_LINPROG_FILE == YES
      sprintf(file_name,"cplex_linprog_refill.%05d.lp",g_solutions_cnt);
      CPXwriteprob( g_cplex_env, g_cplex_lp, file_name, "LP");
#endif

   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int do_linprog_mip_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info)
{
   int cplex_ret;
   unsigned int i;
   int *indices_pos;
   double *elements_pos;
   unsigned cnt_pos = 0;
   int *indices_neg;
   double *elements_neg;
   unsigned cnt_neg = 0;
   int s;
   int indices_start[1] = {0};
   double rhs_value[1] = {1.0};
   int num_cols;
   int solstat;
   char name_neg[300];

#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif

#if DO_SHORTCUT == NO
   char name_pos[300];
#else
   int *indices_shortcut;
   double *elements_shortcut;
   int cnt_shortcut = 0;

   if( (indices_shortcut = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_shortcut\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_shortcut = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_shortcut\n");
      exit(EXIT_FAILURE);
   }
#endif

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

   if( (indices_pos = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (indices_neg = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_neg\n");
      exit(EXIT_FAILURE);
   }
   if( (elements_pos = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_neg = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_neg\n");
      exit(EXIT_FAILURE);
   }

   if( (g_solutions_cnt > 0) )
   {
#if DO_SHORTCUT == YES
      if( g_shortcut_possible == 1 )
      {
         for( s = 0; s < g_solutions_cnt_last; s++ )
         {
            cnt_neg = 0;
            for( i = 0; i < map_info->num_dupset_keepers; i++ )
            {
               if( g_r[s*num_cols + i] <= 0.5 )
               {
                  indices_neg[cnt_neg] = i;
                  elements_neg[cnt_neg] = 1.0;
                  cnt_neg++;
               }
            }

            sprintf(name_neg,"NEWNEG_RESET%07d",s);
            // cplex_ret = GRBaddconstr(g_cplex_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);
            rhs_value[0] = 1.0;
            cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_neg, rhs_value, "G", indices_start, indices_neg, elements_neg, NULL, NULL);

            if( cplex_ret != 0 )
            {
               printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",cplex_ret);
               printf("             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
         }

         cnt_shortcut = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            indices_shortcut[cnt_shortcut] = i;
            elements_shortcut[cnt_shortcut] = 1.0;
            cnt_shortcut++;
         }
         // cplex_ret = GRBaddconstr(g_cplex_model, cnt_shortcut, indices_shortcut, elements_shortcut, GRB_LESS_EQUAL, g_lp_max[g_solutions_cnt-1],"SHORTCUT");
         rhs_value[0] = g_lp_max[g_solutions_cnt-1];
         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_shortcut, rhs_value, "L", indices_start, indices_shortcut, elements_shortcut, NULL, NULL);

         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: Adding shortcut constrain to model failed. Return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
#endif

      for( s = g_solutions_cnt_last; s < g_solutions_cnt; s++ )
      {
         cnt_pos = 0;
         cnt_neg = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            if( g_r[s*num_cols + i] > 0.5 )
            {
               indices_pos[cnt_pos] = i;
               elements_pos[cnt_pos] = 1.0;
               cnt_pos++;
            }
            else
            {
               indices_neg[cnt_neg] = i;
               elements_neg[cnt_neg] = 1.0;
               cnt_neg++;
            }
         }

         if( cnt_pos != (unsigned int)g_lp_max[s] )
         {
            printf("FATAL ERROR: number of elements of new row (%u) are not identical to norm of row (%u) for s=%d\n",cnt_pos,(unsigned int)g_lp_max[s],s);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

#if DO_SHORTCUT == NO
            sprintf(name_pos,"NEWPOS%07d",g_solutions_cnt);
            // cplex_ret = GRBaddconstr(g_cplex_model, cnt_pos, indices_pos, elements_pos, GRB_LESS_EQUAL, (g_lp_max[s] - 1.0),name_pos);
            rhs_value[0] = g_lp_max[s] - 1.0;
            cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_pos, rhs_value, "L", indices_start, indices_pos, elements_pos, NULL, NULL);

            if( cplex_ret != 0 )
            {
               printf("FATAL ERROR: Adding positive constrain to model failed. Return code = %d\n",cplex_ret);
               printf("             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
#endif
         sprintf(name_neg,"NEWNEG%07d",g_solutions_cnt);
         // cplex_ret = GRBaddconstr(g_cplex_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);
         rhs_value[0] = 1.0;
         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_neg, rhs_value, "G", indices_start, indices_neg, elements_neg, NULL, NULL);

         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
   }


#if WRITE_LINPROG_FILE == YES
   sprintf(file_name,"cplex_linprog_addconst.%05d.lp",g_solutions_cnt);
   printf("INFO: writing lp file '%s'\n",file_name);
   CPXwriteprob( g_cplex_env, g_cplex_lp, file_name, "LP");;
#endif

   g_solutions_cnt_last = g_solutions_cnt;

   cplex_ret = CPXmipopt (g_cplex_env, g_cplex_lp);

   if( cplex_ret != 0 )
   {
      printf("INFO: CPLEX could not find optimal solution. Return code = %d\n",cplex_ret);
      return(0);
   }

   if( g_solutions_cnt >= g_max_expected_knockout_combinations )
   {
      resize_arrays_cplex(mode_info, reac_info, cmd_options, map_info);
   }

#if DO_SHORTCUT != NO
   free(indices_shortcut);
   free(elements_shortcut);
#endif
   free(indices_pos);
   free(elements_pos);
   free(indices_neg);
   free(elements_neg);

   cplex_ret = CPXsolution (g_cplex_env, g_cplex_lp, &solstat, &g_lp_max[g_solutions_cnt],&g_r[g_solutions_cnt*num_cols] , NULL, NULL, NULL);
   if ( cplex_ret )
   {
      fprintf(stderr, "WARNING: failed to obtain solution.\n");
      return(0);;
   }
   printf ("\nSolution status = %d\n", solstat);
   printf ("Solution value = %f\n", g_lp_max[g_solutions_cnt]);


   g_lp_max[g_solutions_cnt] *= -1;

   printf("INFO: solution calculated by cplex = %e\n",g_lp_max[g_solutions_cnt]);

   g_lp_max[g_solutions_cnt] = round(g_lp_max[g_solutions_cnt]);

   printf("CPLEX: g_solutions_cnt=%d maximum norm=%f, minimum number of knockouts = %f\n",
          g_solutions_cnt,g_lp_max[g_solutions_cnt],map_info->num_dupset_keepers - g_lp_max[g_solutions_cnt] + map_info->num_always_zero_reacs);

   int norm = 0;
   for( i = 0; i < map_info->num_dupset_keepers; i++ )
   {
      // printf("%s=%3.1f ",reac_info->reactions[map_info->map[map_info->dupsets_map[i]]],g_r[g_solutions_cnt*num_cols + i]);
      if( g_r[g_solutions_cnt*num_cols + i] > 0.5 )
      {
         norm++;
      }
   }
   printf("\n");

   if( fabs(g_lp_max[g_solutions_cnt] - (float)norm) > 1e-7 )
   {
      printf("FATAL ERROR: cplex's maximum %d is not equal to the norm of the result vector (%d)\n",(int) g_lp_max[g_solutions_cnt],norm);
      printf("             ");
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         printf("col %d: val %f ",i,g_r[g_solutions_cnt*num_cols + i]);
      }
      printf("\n");
      exit(EXIT_FAILURE);
   }

   // uncompress last solution and write them all to file if required
   if( cmd_options->solution_range > 0 && g_solutions_cnt > 0 && !(g_lp_max[g_solutions_cnt] > 0 && (fabs(g_lp_max[g_solutions_cnt] - g_lp_max[0] ) < cmd_options->solution_range) ) )
   {
      printf("INFO: stopping programing: g_lp_max[0]=%f g_lp_max[g_solutions_cnt]=%f solution_range=%d\n",g_lp_max[0],g_lp_max[g_solutions_cnt],cmd_options->solution_range);
      return(0);
   }
   uncompress_and_write_solution(g_solutions_cnt, map_info, reac_info, cmd_options, mode_info);

   g_solutions_cnt++;

   if( (g_solutions_cnt > 1) && (g_lp_max[g_solutions_cnt - 1] != g_lp_max[g_solutions_cnt - 2]) )
   {
      // if old solution uses less knockouts then this solution
      // then we can use 'short' cut
      g_shortcut_possible = 1;
   }
   else
   {
      g_shortcut_possible = 0;
   }

   return(1);
}
//////////////////////////////////////////////////////
#endif

#if LP_TOOLKIT == USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void print_g_cplex_arrays(struct_map_info *map_info, struct_mode_info *mode_info, struct_cmd_options *cmd_options, unsigned long long int write_cnt)
{
   unsigned int c;
   unsigned long long int r;

   for( r = 0; r < mode_info->bad_efms; r++ )
   {
      print_mode(&mode_info->pnt_flux[r], mode_info->num_unit_size);
      printf(" norm of mode %llu: %d\n",r,g_pnt_norms[r]);
   }

   printf("pointer to g_cplex_env=%p\n",g_cplex_env);
   printf("pointer to g_cplex_lp=%p\n",g_cplex_lp);

   if( g_n_out_of_D == 1 )
   {
      printf("g_solutions_cnt=%d\n",g_solutions_cnt);
      for( c = 0; c < map_info->num_dupset_keepers + mode_info->good_emfs; c++ )
      {
         printf("cplex variable name %u: %s ",c,g_cplex_varnames[c]);
         if( g_cplex_col_is_integer[c] == CPX_BINARY )
         {
            printf("BINARY\n");
         }
         else if( g_cplex_col_is_integer[c] == CPX_CONTINUOUS )
         {
            printf("GRB_CONTINUOUS\n");
         }
         else
         {
            printf("UNKNOWN!!!!!!\n");
         }
         printf("objective coeff %u: %g\n",c,g_cplex_obj_cf[c]);
         printf("lower bound %u: %g\n",c,g_cplex_lb[c]);
         printf("upper bound %u: %g\n",c,g_cplex_ub[c]);
      }
      for( r = 0; r < mode_info->bad_efms + g_solutions_cnt + 2*mode_info->good_emfs + 1; r++ )
      {
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("constraint type %llu: EQUAL ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("constraint type %llu: LESS ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf("constraint type %llu: GREATER ",r);
         }
         else
         {
            printf("constraint type %llu: UNkNOWN ",r);
         }
         // printf("right_hand_side %llu: %g\n",r,g_cplex_row_righthandside[r]);
         printf(" %g\n",g_cplex_row_righthandside[r]);
      }
      for( c = 0; c < map_info->num_dupset_keepers + mode_info->good_emfs; c++ )
      {
         printf("g_cplex_start[%d]=%d g_cplex_len[%d]=%d\n",c,g_cplex_start[c],c,g_cplex_len[c]);
      }
   }
   else
   {
      printf("g_solutions_cnt=%d\n",g_solutions_cnt);
      for( c = 0; c < map_info->num_dupset_keepers; c++ )
      {
         printf("cplex variable name %u: %s ",c,g_cplex_varnames[c]);
         if( g_cplex_col_is_integer[c] == CPX_BINARY )
         {
            printf("BINARY\n");
         }
         else if( g_cplex_col_is_integer[c] == CPX_CONTINUOUS )
         {
            printf("GRB_CONTINUOUS\n");
         }
         else
         {
            printf("UNKNOWN!!!!!!\n");
         }
         printf("objective coeff %u: %g\n",c,g_cplex_obj_cf[c]);
         printf("lower bound %u: %g\n",c,g_cplex_lb[c]);
         printf("upper bound %u: %g\n",c,g_cplex_ub[c]);
         printf("variablename %u: %s\n",c,g_cplex_varnames[c]);
      }
      for( r = 0; r < mode_info->bad_efms + g_solutions_cnt; r++ )
      {
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("constraint type %llu: EQUAL ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("constraint type %llu: LESS ",r);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf("constraint type %llu: GREATER ",r);
         }
         else
         {
            printf("constraint type %llu: UNkNOWN ",r);
         }
         // printf("right_hand_side %llu: %g\n",r,g_cplex_row_righthandside[r]);
         printf(" %g\n",g_cplex_row_righthandside[r]);
      }
      for( c = 0; c < map_info->num_dupset_keepers; c++ )
      {
         printf("g_cplex_start[%d]=%d g_cplex_len[%d]=%d g_cplex_index[%d]=%d g_cplex_value[%d]=%g\n",
                 c,g_cplex_start[c],c,g_cplex_len[c],c,g_cplex_index[c],c,g_cplex_value[c]);
      }

      /////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////
      for( r = 0; r < mode_info->bad_efms; r++ )
      {
         int start = g_cplex_start_trans[r];
         int stop;
         if( r == mode_info->bad_efms - 1 )
         {
            stop = write_cnt;
         }
         else
         {
            stop = g_cplex_start_trans[r+1];
         }
         printf("r=%llu start=%x stop=%d\n",r,start,stop);
         for( c = start; c < stop; c++ )
         {
            printf("%g*%d:%s ",g_cplex_value_trans[c],g_cplex_index_trans[c],g_cplex_varnames[g_cplex_index_trans[c]]);
         }
         if( g_cplex_row_constr_type[r] == 'E' )
         {
            printf("== %g\n",g_cplex_row_righthandside[r]);
         }
         else if( g_cplex_row_constr_type[r] == 'L' )
         {
            printf("<= %g\n",g_cplex_row_righthandside[r]);
         }
         else if( g_cplex_row_constr_type[r] == 'G' )
         {
            printf(">= %g\n",g_cplex_row_righthandside[r]);
         }
         else
         {
            printf("UNKNOWN %g\n",g_cplex_row_righthandside[r]);
         }
      }
      /////////////////////////////////////////////////////////////////////////////////////////
   }

   // print transpose information
   printf("DEBUG: leaving print_g_cplex_arrays().\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void cplex_transpose_arrays(unsigned long long int num_idx, unsigned long long int num_rows, unsigned long long int num_cols)
{
   unsigned long long int r, c, x;
   unsigned long long int new_idx = 0;

   printf("DEBUG: cplex_transpose_arrays(): entered.\n");

   g_cplex_start_trans = (int *) malloc( (size_t) sizeof(int)*num_rows );

   if( g_cplex_start_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_start_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_index_trans = (int *) malloc( (size_t) sizeof(int)*num_idx );

   if( g_cplex_index_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_index_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_value_trans = (double *) malloc( (size_t) sizeof(double)*num_idx );

   if( g_cplex_value_trans == NULL )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): failed to allocate memory for g_cplex_value_trans\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   for( r = 0; r < num_rows; r++ )
   {
      g_cplex_start_trans[r] = new_idx;

      for( c = 0; c < num_cols; c++ )
      {
         for( x = g_cplex_start[c]; x < g_cplex_start[c] + g_cplex_len[c]; x++ )
         {
            if( r == g_cplex_index[x] )
            {
               g_cplex_value_trans[new_idx] = g_cplex_value[x];
               g_cplex_index_trans[new_idx] = c;
               new_idx++;
            }
         }
      }
   }

   if( new_idx != num_idx )
   {
      fprintf(stderr, "FATAL ERROR: cplex_transpose_arrays(): transposing arrays failed: new_idx=%llu old_idx=%llu\n",new_idx, num_idx );
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   printf("DEBUG: cplex_transpose_arrays(): leaving. transpose %llu values\n", new_idx);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void cplex_clean_transpose()
{
   free(g_cplex_start_trans);
   free(g_cplex_index_trans);
   free(g_cplex_value_trans);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void cleanup_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_mode_info *mode_info)
{
   free(g_lp_max);
   free(g_pnt_norms);
   free(g_pnt_tmp_cutset);
   free(g_r);

   free(g_cplex_col_is_integer);
   free(g_cplex_obj_cf);
   free(g_cplex_lb);
   free(g_cplex_ub);
   free(g_cplex_start);
   free(g_cplex_len);
   free(g_cplex_index);
   free(g_cplex_value);
   free(g_cplex_row_constr_type);
   free(g_cplex_row_righthandside);

   int i;
   unsigned long int num_cols;
   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }
   for( i = 0; i < num_cols; i++ )
   {
      free(g_cplex_varnames[i]);
   }
   free(g_cplex_varnames);

   printf("DEBUG: cleanup_cplex() deleting old cplex prob\n");
   CPXfreeprob (g_cplex_env, &g_cplex_lp);
   printf("DEBUG: cleanup_cplex(): deleting old cplex env\n");
   CPXcloseCPLEX (&g_cplex_env);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void init_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations = INIT_MAX_EXPECTED_KNOCKOUT_COMBINATIONS;

   g_lp_max = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   g_pnt_tmp_cutset = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*sizeof(unsigned long long)));
   if( g_pnt_tmp_cutset == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for temporary cutset\n");
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf("map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu num_cols=%lu\n",map_info->num_dupset_keepers,mode_info->good_emfs,num_cols);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf("mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf("g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) malloc((size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_varnames = (char **) malloc((size_t) sizeof(char *)*num_cols);

   if( g_cplex_varnames == NULL )
   {
      printf("FATAL ERROR: init_cplex(): couldn't allocate %lu bytes for *g_cplex_varnames\n",sizeof(char *)*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_col_is_integer = (char *) malloc( (size_t) (num_cols*sizeof(char)));
   if( g_cplex_col_is_integer == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_col_is_integer array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_obj_cf = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_obj_cf == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_obj_cf array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_lb = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_lb == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_lb array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_ub = (double *) malloc( (size_t) (num_cols*sizeof(double)));
   if( g_cplex_ub == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_ub array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_start = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_cplex_start == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_start array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_len = (int *) malloc( (size_t) ((num_cols)*sizeof(int)));
   if( g_cplex_len == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_len array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_index = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_cplex_index == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_cplex_value = (double *) malloc( (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_cplex_value == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) malloc( (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_constr_type = (char *) malloc( (size_t) (num_rows*sizeof(char)));
   if( g_cplex_row_constr_type == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_righthandside = (double *) malloc( (size_t) (num_rows*sizeof(double)));
   if( g_cplex_row_righthandside == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for g_cplex_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }

   int i;
   for( i = 0; i < num_cols; i++ )
   {
      printf("allocate memory for g_cplex_varnames[%d]: %lu bytes\n",i,(reac_info->max_len_reac_name+1)*sizeof(char));
      g_cplex_varnames[i] = (char *) malloc( (size_t) ((reac_info->max_len_reac_name+1)*sizeof(char)));
      if( g_cplex_varnames[i] == NULL )
      {
         printf("FATAL ERROR: couldn't allocate memory for g_cplex_varnames[%d] array\n",i);
         exit(EXIT_FAILURE);
      }
   }
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void resize_arrays_cplex(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_map_info *map_info)
{
   unsigned long int num_cols;
   unsigned long int num_rows;

   g_max_expected_knockout_combinations *= 2;

   g_lp_max = (double *) realloc(g_lp_max, (size_t) sizeof(double)*g_max_expected_knockout_combinations);
   if( g_lp_max == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations);
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
      printf(" resize_arrays_cplex():map_info->num_dupset_keepers=%u mode_info->good_emfs=%llu\n",map_info->num_dupset_keepers,mode_info->good_emfs);

      num_rows = mode_info->good_emfs*2 + 1 + mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      printf(" resize_arrays_cplex():mode_info->bad_efms=%llu g_max_expected_knockout_combinations=%d num_rows=%lu\n",mode_info->bad_efms,g_max_expected_knockout_combinations,num_rows);

      g_size_linprog_mem = (mode_info->good_emfs*2 + 1 + mode_info->bad_efms + g_max_expected_knockout_combinations)*(map_info->num_dupset_keepers + 1) + 2*mode_info->good_emfs;
      printf(" resize_arrays_cplex():g_size_linprog_mem=%llu\n",g_size_linprog_mem);
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
      num_rows = mode_info->bad_efms + map_info->num_dupset_keepers + g_max_expected_knockout_combinations;
      g_size_linprog_mem = (mode_info->bad_efms + g_max_expected_knockout_combinations)*map_info->num_dupset_keepers;
   }

   g_r = (double *) realloc(g_r, (size_t) sizeof(double)*g_max_expected_knockout_combinations*num_cols);
   if( g_r == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate %lu bytes for g_r\n",sizeof(double)*g_max_expected_knockout_combinations*num_cols);
      exit(EXIT_FAILURE);
   }

   g_cplex_index = (int *) realloc( g_cplex_index, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_cplex_index == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_index array: %s\n",strerror(errno));
      printf("             g_size_linprog_mem=%llu\n",g_size_linprog_mem);
      exit(EXIT_FAILURE);
   }

   g_cplex_value = (double *) realloc( g_cplex_value, (size_t) (g_size_linprog_mem*sizeof(double)));
   if( g_cplex_value == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_value array\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_norms = (int *) realloc(g_pnt_norms, (size_t) (g_size_linprog_mem*sizeof(int)));
   if( g_pnt_norms == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for norm array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_constr_type = (char *) realloc(g_cplex_row_constr_type, (size_t) (num_rows*sizeof(char)));
   if( g_cplex_row_constr_type == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_row_constr_type array\n");
      exit(EXIT_FAILURE);
   }

   g_cplex_row_righthandside = (double *) realloc(g_cplex_row_righthandside, (size_t) (num_rows*sizeof(double)));
   if( g_cplex_row_righthandside == NULL )
   {
      printf("FATAL ERROR: resize_arrays_cplex(): couldn't allocate memory for g_cplex_row_righthandside array\n");
      exit(EXIT_FAILURE);
   }
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
void fill_linprog_mem_cplex(struct_mode_info *mode_info, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info)
{
   unsigned long i;
   unsigned long lp_row = 0;
   unsigned int j,k;
   int idx;
   unsigned long write_cnt;
   int unit_cnt;
   int bitmover;
   unsigned long long tmp_unit;
   int do_set;
   int norm_mode;
   int cplex_ret;
   char name_neg[300];
   int status;
#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif



   printf("DEBUG: fill_linprog_mem_cplex(): entered. g_shortcut_possible=%d\n",g_shortcut_possible);
   if( (g_solutions_cnt == 0) || ((DO_SHORTCUT == YES) && g_shortcut_possible == 1) )
   {
      if( g_solutions_cnt > 0 )
      {
         printf("DEBUG: deleting old cplex prob\n");
         // delete old context
         status = CPXfreeprob (g_cplex_env, &g_cplex_lp);
         if( status )
         {
            fprintf(stderr,"FATAL ERROR: couldn't free cplex environment. Returned error code = %d\n",status);
            fprintf(stderr,"             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
         printf("DEBUG: deleting old cplex env\n");
         status = CPXcloseCPLEX (&g_cplex_env);
         if ( status )
         {
            char errmsg[1024];
            fprintf (stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring (g_cplex_env, status, errmsg);
            fprintf (stderr, "%s", errmsg);
            fprintf(stderr,"             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }

      /////////////////////////////////////////////////////////////////////////
      // get new CPLEX environment
      /////////////////////////////////////////////////////////////////////////
      printf("DEBUG: creating cplex environment\n");
      g_cplex_env = CPXopenCPLEX (&status);

      if ( g_cplex_env == NULL )
      {
         char errmsg[1024];
         fprintf (stderr, "Could not open CPLEX environment.\n");
         CPXgeterrorstring (g_cplex_env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      printf("DEBUG: creating cplex prob\n");
      g_cplex_lp = CPXcreateprob (g_cplex_env, &status, "lpex1");

      if ( g_cplex_lp == NULL )
      {
         fprintf (stderr, "Failed to create CPLEX LP.\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_OPPORTUNISTIC);
      // status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_AUTO);
      // status = CPXsetintparam( g_cplex_env, CPX_PARAM_PARALLELMODE, CPX_PARALLEL_DETERMINISTIC);
      if ( status )
      {
         fprintf (stderr, "Failed to parallel mode type\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_AUXROOTTHREADS, cmd_options->num_threads);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to use %d threads for root.\n", cmd_options->num_threads);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_THREADS, cmd_options->num_threads);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to use %d threads.\n", cmd_options->num_threads);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_SCRIND, CPX_ON);
      if ( status )
      {
         fprintf (stderr, "Failed to turn on CPLEX screen indicator.\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_POPULATELIM, MAX_SOLUTION_ENUM);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to enumarate %d solutions.\n", MAX_SOLUTION_ENUM);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_SOLNPOOLCAPACITY, MAX_SOLUTION_CAPACITY);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to a pool capacity of %d.\n", MAX_SOLUTION_CAPACITY);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetintparam( g_cplex_env, CPX_PARAM_SOLNPOOLINTENSITY, SOLUTION_POOL_INTENSITY);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to solution pool intensity of %d.\n", SOLUTION_POOL_INTENSITY);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetdblparam( g_cplex_env, CPX_PARAM_SOLNPOOLGAP, SOLUTION_POOL_GAP);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to solution pool gap of %f.\n", SOLUTION_POOL_GAP);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetdblparam( g_cplex_env, CPX_PARAM_EPAGAP, LINPROG_MIP_GAP_ABS);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to absolute MIP gap tolerance of %f.\n", LINPROG_MIP_GAP_ABS);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXsetdblparam( g_cplex_env, CPX_PARAM_EPGAP, LINPROG_MIP_GAP_REL);
      if ( status )
      {
         fprintf (stderr, "Failed to configure CPLEX to relative MIP gap tolerance of %f.\n", LINPROG_MIP_GAP_REL);
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      status = CPXchgobjsen (g_cplex_env, g_cplex_lp, CPX_MIN);
      if ( status )
      {
         fprintf (stderr, "Failed to tell CPLEX to minimize\n");
         fprintf(stderr,"             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // define all columns
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         g_cplex_obj_cf[i] = -1.0;
         // g_cplex_col_is_integer[i] = GRB_BINARY;
         g_cplex_col_is_integer[i] = CPX_BINARY;
         g_cplex_lb[i] = 0.0;
         g_cplex_ub[i] = 1.0;
         strcpy(g_cplex_varnames[i],reac_info->reactions[map_info->map[map_info->dupsets_map[i]]]);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // define all extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            g_cplex_obj_cf[i+map_info->num_dupset_keepers] = 0.0;
            // g_cplex_col_is_integer[i+map_info->num_dupset_keepers] = GRB_BINARY;
            g_cplex_col_is_integer[i+map_info->num_dupset_keepers] = CPX_BINARY;
            g_cplex_lb[i+map_info->num_dupset_keepers] = 0.0;
            g_cplex_ub[i+map_info->num_dupset_keepers] = 1.0;
            sprintf(name_neg,"Y%06lu",i);
            strcpy(g_cplex_varnames[i+map_info->num_dupset_keepers],name_neg);
         }
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // set contraint type
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         // g_cplex_row_constr_type[i]  = GRB_LESS_EQUAL;
         g_cplex_row_constr_type[i]  = 'L';
      }

      /////////////////////////////////////////////////////////////////////////
      // set all extra constraint types because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            // g_cplex_row_constr_type[i + mode_info->bad_efms]  = GRB_GREATER_EQUAL;
            g_cplex_row_constr_type[i + mode_info->bad_efms]  = 'G';
         }
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            // g_cplex_row_constr_type[i + mode_info->bad_efms + mode_info->good_emfs]  = GRB_LESS_EQUAL;
            g_cplex_row_constr_type[i + mode_info->bad_efms + mode_info->good_emfs]  = 'L';
         }
         // extra constraint for ||y|| <= |D| - n
         // which is equal to ||y|| <= g_good_emfs - g_good_emfs_wanted
         g_cplex_row_constr_type[2*mode_info->good_emfs + mode_info->bad_efms]  = 'L';
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // calculate norm of all bad modes
      /////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         norm_mode = 0;
         lp_row++;

         for( j = 0; j < map_info->num_dupset_keepers; j++ )
         {
            // reset temporary cutset array
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               g_pnt_tmp_cutset[k] = 0;
            }

            tmp_unit = 1;
            idx = j;

            bitmover = idx;

            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);

            tmp_unit <<= bitmover;

            // fill temporary cutset array with data
            g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               norm_mode++;
            }
         }

         if( norm_mode == 0 )
         {
            printf("FATAL ERROR: norm of mode %d is equal to 0\n",norm_mode);
            exit(EXIT_FAILURE);
         }

         g_pnt_norms[i] = norm_mode;
         g_cplex_row_righthandside[i] = norm_mode - 1;
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // get the norm of the all good modes because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( i = 0; i < mode_info->good_emfs; i++ )
         {
            norm_mode = 0;
            lp_row++;

            for( j = 0; j < map_info->num_dupset_keepers; j++ )
            {
               // reset temporary cutset array
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  g_pnt_tmp_cutset[k] = 0;
               }

               tmp_unit = 1;
               idx = j;

               bitmover = idx;

               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);

               tmp_unit <<= bitmover;

               // fill temporary cutset array with data
               g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  norm_mode++;
               }
            }

            g_pnt_norms[i+mode_info->bad_efms]             = norm_mode;
            g_pnt_norms[i+mode_info->bad_efms+mode_info->good_emfs] = norm_mode;
            g_cplex_row_righthandside[i+mode_info->bad_efms]             = norm_mode;
            g_cplex_row_righthandside[i+mode_info->bad_efms+mode_info->good_emfs] = (2*norm_mode) - 1;
         }
         // extra constraint ||y|| <= g_good_emfs - g_good_emfs_wanted
         // g_cplex_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = g_good_emfs_orig - g_good_emfs_wanted;
         g_cplex_row_righthandside[mode_info->bad_efms+mode_info->good_emfs*2] = mode_info->good_emfs_orig - mode_info->num_good_removed_by_always_zero - cmd_options->good_efms_wanted;
      }
      /////////////////////////////////////////////////////////////////////////
      printf("number of constraints: %llu\n mode_info->bad_efms=%llu mode_info->good_emfs=%llu\n",
              mode_info->bad_efms+mode_info->good_emfs*2,mode_info->bad_efms,mode_info->good_emfs);

      printf("DEBUG: fill_linprog_mem_cplex(): after norm calculation\n");

      /////////////////////////////////////////////////////////////////////////
      // write linear programming problem column by column
      /////////////////////////////////////////////////////////////////////////
      write_cnt = 0;
      for( j = 0; j < map_info->num_dupset_keepers; j++ )
      {
         g_cplex_start[j] = write_cnt;
         for( k = 0; k < mode_info->num_unit_size; k++ )
         {
            g_pnt_tmp_cutset[k] = 0;
         }

         tmp_unit = 1;
         idx = j;

         bitmover = idx;

         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);

         tmp_unit <<= bitmover;
         g_pnt_tmp_cutset[unit_cnt] |= tmp_unit;

         for( i = 0; i < mode_info->bad_efms; i++ )
         {

            do_set = 0;
            for( k = 0; k < mode_info->num_unit_size; k++ )
            {
               if( (mode_info->pnt_flux[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
               {
                  do_set = 1;
                  break;
               }
            }

            if( do_set == 1 )
            {
               if( write_cnt >= g_size_linprog_mem )
               {
                  printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                  printf("             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               g_cplex_index[write_cnt] = i;
               g_cplex_value[write_cnt] = 1.0;
               write_cnt++;
            }
         }

         //////////////////////////////////////////////////////////////////////
         // write column value for extra constraint because of n-out-of-D 
         //////////////////////////////////////////////////////////////////////
         if( g_n_out_of_D == 1 )
         {
            ///////////////////////////////////////////////////////////////////
            // do it for the >= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_cplex_index[write_cnt] = i + mode_info->bad_efms;
                  g_cplex_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
            ///////////////////////////////////////////////////////////////////
            // do it for the <= constraints of the n-out-of-D feature
            ///////////////////////////////////////////////////////////////////
            for( i = 0; i < mode_info->good_emfs; i++ )
            {

               do_set = 0;
               for( k = 0; k < mode_info->num_unit_size; k++ )
               {
                  if( (mode_info->pnt_good[i*mode_info->num_unit_size + k] & g_pnt_tmp_cutset[k]) != 0 )
                  {
                     do_set = 1;
                     break;
                  }
               }

               if( do_set == 1 )
               {
                  if( write_cnt >= g_size_linprog_mem )
                  {
                     printf("FATAL ERROR: index and value array not large enough (write_cnt=%lu,g_size_linprog_mem=%llu)\n",write_cnt,g_size_linprog_mem);
                     printf("             execution aborted.\n");
                     exit(EXIT_FAILURE);
                  }
                  g_cplex_index[write_cnt] = i + mode_info->bad_efms + mode_info->good_emfs;
                  g_cplex_value[write_cnt] = 1.0;
                  write_cnt++;
               }
            }
         }
         g_cplex_len[j] = write_cnt - g_cplex_start[j];
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // write extra columns because of n-out-of-D feature
      /////////////////////////////////////////////////////////////////////////
      if( g_n_out_of_D == 1 )
      {
         for( j = 0; j < mode_info->good_emfs; j++ )
         {
            g_cplex_start[j+map_info->num_dupset_keepers] = write_cnt;

            g_cplex_index[write_cnt] = j + mode_info->bad_efms;
            g_cplex_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_cplex_index[write_cnt] = j + mode_info->bad_efms + mode_info->good_emfs;
            g_cplex_value[write_cnt] = g_pnt_norms[j+mode_info->bad_efms];
            write_cnt++;

            g_cplex_index[write_cnt] = mode_info->bad_efms + 2*mode_info->good_emfs;
            // g_cplex_value[write_cnt] = 1;
            g_cplex_value[write_cnt] = (float)mode_info->pnt_good_duplicates[j];
            write_cnt++;

            g_cplex_len[j+map_info->num_dupset_keepers] = write_cnt - g_cplex_start[j+map_info->num_dupset_keepers];
         }

      }
      /////////////////////////////////////////////////////////////////////////


      g_num_efm_writes = write_cnt;

      printf("DEBUG: fill_linprog_mem_cplex(): after filling index/value values\n");


      if( g_n_out_of_D == 1 )
      {
         cplex_transpose_arrays(write_cnt, mode_info->bad_efms + 2*mode_info->good_emfs + 1, map_info->num_dupset_keepers + mode_info->good_emfs);
         // print_g_cplex_arrays(map_info, mode_info, cmd_options, write_cnt);

         cplex_ret = CPXnewcols(g_cplex_env, g_cplex_lp, map_info->num_dupset_keepers + mode_info->good_emfs, g_cplex_obj_cf, g_cplex_lb, g_cplex_ub, g_cplex_col_is_integer, g_cplex_varnames);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXnewcols returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         cplex_ret = CPXaddrows (g_cplex_env, g_cplex_lp, 0, mode_info->bad_efms + 2*mode_info->good_emfs + 1,
                                 write_cnt, g_cplex_row_righthandside,
                                 g_cplex_row_constr_type, g_cplex_start_trans, g_cplex_index_trans, g_cplex_value_trans, NULL, NULL);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXaddrows returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
      else
      {
         printf("DEBUG: number of rows going to be added %llu\n",mode_info->bad_efms);
         printf("DEBUG: wrote %lu index and value entries to arrays\n",write_cnt);
         cplex_transpose_arrays(write_cnt, mode_info->bad_efms, map_info->num_dupset_keepers);
         // print_g_cplex_arrays(map_info, mode_info, cmd_options, write_cnt);

         cplex_ret = CPXnewcols(g_cplex_env, g_cplex_lp, map_info->num_dupset_keepers, g_cplex_obj_cf, g_cplex_lb, g_cplex_ub, g_cplex_col_is_integer, g_cplex_varnames);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXnewcols returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, mode_info->bad_efms, write_cnt, g_cplex_row_righthandside,
                                g_cplex_row_constr_type, g_cplex_start_trans, g_cplex_index_trans, g_cplex_value_trans, NULL, NULL);
         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: CPXaddrows returned abnormally: return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }

      cplex_clean_transpose();


      // write linear programming problem to file
#if WRITE_LINPROG_FILE == YES
      sprintf(file_name,"cplex_linprog_refill.%05d.lp",g_solutions_cnt);
      CPXwriteprob( g_cplex_env, g_cplex_lp, file_name, "LP");
#endif

   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// USE_CPLEX_SOLPOOL
//////////////////////////////////////////////////////
int do_linprog_mip_cplex(struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info)
{
   int cplex_ret;
   unsigned int i;
   int *indices_pos;
   double *elements_pos;
   unsigned cnt_pos = 0;
   int *indices_neg;
   double *elements_neg;
   unsigned cnt_neg = 0;
   int s;
   int indices_start[1] = {0};
   double rhs_value[1] = {1.0};
   int num_cols;
   int solstat;
   char name_neg[300];

#if WRITE_LINPROG_FILE == YES
   char file_name[300];
#endif

   int *indices_shortcut;
   double *elements_shortcut;
   int cnt_shortcut = 0;

   if( (indices_shortcut = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_shortcut\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_shortcut = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_shortcut\n");
      exit(EXIT_FAILURE);
   }

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

   if( (indices_pos = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (indices_neg = (int *) malloc((size_t)sizeof(int)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for indices_neg\n");
      exit(EXIT_FAILURE);
   }
   if( (elements_pos = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_pos\n");
      exit(EXIT_FAILURE);
   }

   if( (elements_neg = (double *) malloc((size_t)sizeof(double)*reac_info->num_reactions)) == NULL )
   {
      printf("FATAL ERROR: do_linprog_mip_cplex(): couldn't locate memory for elements_neg\n");
      exit(EXIT_FAILURE);
   }

   if( (g_solutions_cnt > 0) )
   {
      if( g_shortcut_possible == 1 )
      {
         for( s = 0; s < g_solutions_cnt_last; s++ )
         {
            cnt_neg = 0;
            for( i = 0; i < map_info->num_dupset_keepers; i++ )
            {
               if( g_r[s*num_cols + i] <= 0.5 )
               {
                  indices_neg[cnt_neg] = i;
                  elements_neg[cnt_neg] = 1.0;
                  cnt_neg++;
               }
            }

            sprintf(name_neg,"NEWNEG_RESET%07d",s);
            // cplex_ret = GRBaddconstr(g_cplex_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);
            rhs_value[0] = 1.0;
            cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_neg, rhs_value, "G", indices_start, indices_neg, elements_neg, NULL, NULL);

            if( cplex_ret != 0 )
            {
               printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",cplex_ret);
               printf("             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
         }

         cnt_shortcut = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            indices_shortcut[cnt_shortcut] = i;
            elements_shortcut[cnt_shortcut] = 1.0;
            cnt_shortcut++;
         }
         // cplex_ret = GRBaddconstr(g_cplex_model, cnt_shortcut, indices_shortcut, elements_shortcut, GRB_LESS_EQUAL, g_lp_max[g_solutions_cnt-1],"SHORTCUT");
         rhs_value[0] = g_lp_max[g_solutions_cnt-1];
         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_shortcut, rhs_value, "L", indices_start, indices_shortcut, elements_shortcut, NULL, NULL);

         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: Adding shortcut constrain to model failed. Return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }

      for( s = g_solutions_cnt_last; s < g_solutions_cnt; s++ )
      {
         cnt_pos = 0;
         cnt_neg = 0;
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            if( g_r[s*num_cols + i] > 0.5 )
            {
               indices_pos[cnt_pos] = i;
               elements_pos[cnt_pos] = 1.0;
               cnt_pos++;
            }
            else
            {
               indices_neg[cnt_neg] = i;
               elements_neg[cnt_neg] = 1.0;
               cnt_neg++;
            }
         }

         if( cnt_pos != (unsigned int)g_lp_max[s] )
         {
            printf("FATAL ERROR: number of elements of new row (%u) are not identical to norm of row (%u) for s=%d\n",cnt_pos,(unsigned int)g_lp_max[s],s);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         sprintf(name_neg,"NEWNEG%07d",g_solutions_cnt);
         // cplex_ret = GRBaddconstr(g_cplex_model, cnt_neg, indices_neg, elements_neg, GRB_GREATER_EQUAL, 1.0,name_neg);
         rhs_value[0] = 1.0;
         cplex_ret = CPXaddrows(g_cplex_env, g_cplex_lp, 0, 1, cnt_neg, rhs_value, "G", indices_start, indices_neg, elements_neg, NULL, NULL);

         if( cplex_ret != 0 )
         {
            printf("FATAL ERROR: Adding negative constrain to model failed. Return code = %d\n",cplex_ret);
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }
      }
   }


#if WRITE_LINPROG_FILE == YES
   sprintf(file_name,"cplex_linprog_addconst.%05d.lp",g_solutions_cnt);
   printf("INFO: writing lp file '%s'\n",file_name);
   CPXwriteprob( g_cplex_env, g_cplex_lp, file_name, "LP");;
#endif

   g_solutions_cnt_last = g_solutions_cnt;

   cplex_ret = CPXpopulate(g_cplex_env, g_cplex_lp);

   if( cplex_ret != 0 )
   {
      printf("INFO: CPLEX could not find optimal solution. Return code = %d\n",cplex_ret);
      free(indices_shortcut);
      free(elements_shortcut);
      free(indices_pos);
      free(elements_pos);
      free(indices_neg);
      free(elements_neg);
      return(0);
   }

   solstat = CPXgetstat(g_cplex_env, g_cplex_lp);
   printf("INFO: solution status: %d\n",solstat);

   if( solstat == CPXMIP_INFEASIBLE )
   {
      printf("INFO: CPLEX prolbem is infeasible. solstat = %d\n",solstat);
      free(indices_shortcut);
      free(elements_shortcut);
      free(indices_pos);
      free(elements_pos);
      free(indices_neg);
      free(elements_neg);
      return(0);
   }

   free(indices_shortcut);
   free(elements_shortcut);
   free(indices_pos);
   free(elements_pos);
   free(indices_neg);
   free(elements_neg);


   int numsolns = CPXgetsolnpoolnumsolns(g_cplex_env, g_cplex_lp);
   printf("INFO: CPXpopulate computed %d solutions\n",numsolns);

   int numreplaced = CPXgetsolnpoolnumreplaced(g_cplex_env, g_cplex_lp);
   printf("INFO: Solutions replace from pool: %d\n",numreplaced);

   int scn;
   double sol_obj_val;
   for( scn = 0; scn < numsolns; scn++ )
   {
      if( g_solutions_cnt >= g_max_expected_knockout_combinations )
      {
         resize_arrays_cplex(mode_info, reac_info, cmd_options, map_info);
      }

      cplex_ret = CPXgetsolnpoolobjval(g_cplex_env, g_cplex_lp,scn, &sol_obj_val);
      if( cplex_ret )
      {
         fprintf(stderr, "WARNING: failed to obtain objective value form solution pool for solution %d of %d solutions.\n",scn,numsolns);
         exit(EXIT_FAILURE);
      }
      printf("INFO: objective value form solution pool for solution %d of %d solutions: %f\n",scn,numsolns,sol_obj_val);

      g_lp_max[g_solutions_cnt] = sol_obj_val;
      g_lp_max[g_solutions_cnt] *= -1;
      g_lp_max[g_solutions_cnt] = round(g_lp_max[g_solutions_cnt]);
      printf("CPLEX: g_solutions_cnt=%d maximum norm=%f, minimum number of knockouts = %f\n",
             g_solutions_cnt,g_lp_max[g_solutions_cnt],map_info->num_dupset_keepers - g_lp_max[g_solutions_cnt] + map_info->num_always_zero_reacs);

      cplex_ret = CPXgetsolnpoolx(g_cplex_env, g_cplex_lp,scn,&g_r[g_solutions_cnt*num_cols], 0, num_cols - 1);
      if( cplex_ret )
      {
         fprintf(stderr, "WARNING: failed to obtain solution form solution pool for solution %d of %d solutions.\n",scn,numsolns);
         exit(EXIT_FAILURE);
      }

      int norm = 0;
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {
         printf("%s=%3.1f ",reac_info->reactions[map_info->map[map_info->dupsets_map[i]]],g_r[g_solutions_cnt*num_cols + i]);
         if( g_r[g_solutions_cnt*num_cols + i] > 0.5 )
         {
            norm++;
         }
      }
      printf("\n");

      if( fabs(g_lp_max[g_solutions_cnt] - (float)norm) > 1e-7 )
      {
         printf("FATAL ERROR: cplex's maximum %d is not equal to the norm of the result vector (%d)\n",(int) g_lp_max[g_solutions_cnt],norm);
         printf("             ");
         for( i = 0; i < map_info->num_dupset_keepers; i++ )
         {
            printf("col %d: val %f ",i,g_r[g_solutions_cnt*num_cols + i]);
         }
         printf("\n");
         exit(EXIT_FAILURE);
      }

      // uncompress last solution and write them all to file if required
      if( cmd_options->solution_range > 0 && g_solutions_cnt > 0 && !(g_lp_max[g_solutions_cnt] > 0 && (fabs(g_lp_max[g_solutions_cnt] - g_lp_max[0] ) < cmd_options->solution_range) ) )
      {
         printf("INFO: stopping programing: g_lp_max[0]=%f g_lp_max[g_solutions_cnt]=%f solution_range=%d\n",g_lp_max[0],g_lp_max[g_solutions_cnt],cmd_options->solution_range);
         return(0);
      }
      uncompress_and_write_solution(g_solutions_cnt, map_info, reac_info, cmd_options, mode_info);

      g_solutions_cnt++;
   }

   g_shortcut_possible = 1;

   return(1);
}
//////////////////////////////////////////////////////
#endif

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void uncompress_and_write_solution(int sol_num, struct_map_info *map_info, struct_reac_info *reac_info, struct_cmd_options *cmd_options, struct_mode_info *mode_info)
{
   int i;
   int d;
   int m;
   int x;
   int z;
   int new_l_r_cnt = 0;
   double *new_cs;
   double *new_cs_dup;
   double *new_l_r;
   unsigned long long int t_max_sol_enhancer = g_max_sol_enhancer;
   int num_cols;

   printf("DEBUG: entered uncompress_and_write_solution(): sol_num=%d\n",sol_num);

   if( g_n_out_of_D == 1 )
   {
      num_cols = map_info->num_dupset_keepers + mode_info->good_emfs;
   }
   else
   {
      num_cols = map_info->num_dupset_keepers;
   }

   new_l_r = (double *) malloc((size_t) t_max_sol_enhancer*map_info->map_len*sizeof(double));
   if( new_l_r == NULL )
   {
      printf("FATAL ERROR: uncompress_and_write_solution(): couldn't allocate memory for local new solution array\n");
      printf("             number of bytes tried to allocate: %llu\n",g_max_sol_enhancer*map_info->map_len*sizeof(double));
      printf("             g_max_sol_enhancer: %llu\n",g_max_sol_enhancer);
      printf("             map_info->map_len: %d\n",map_info->map_len);
      printf("             sizeof(double): %lu\n",sizeof(double));
      exit(EXIT_FAILURE);
   }

   if( (new_cs = (double *) malloc( (size_t)sizeof(double)*reac_info->num_reactions) ) == NULL )
   {
      printf("FATAL ERROR: uncompress_and_write_solution(): couldn't allocate memory for new_cs\n");
      exit(EXIT_FAILURE);
   }

   if( (new_cs_dup = (double *) malloc( (size_t)sizeof(double)*reac_info->num_reactions) ) == NULL )
   {
      printf("FATAL ERROR: uncompress_and_write_solution(): couldn't allocate memory for new_cs\n");
      exit(EXIT_FAILURE);
   }

   // int duplicate_start = new_g_r_cnt;
   int duplicate_start = new_l_r_cnt;

   // initialize unfolded cutset
   for( i = 0; i < map_info->map_len; i++ )
   {
      new_cs[i] = 1.0;
      new_l_r[new_l_r_cnt*map_info->map_len + i] = 1.0;
   }

   // uncompress cutset
   for( i = 0; i < map_info->num_dupset_keepers; i++ )
   {
      if( g_r[sol_num*num_cols + i] < 0.5)
      {
         new_cs[map_info->dupsets_map[i]] = 0.0;
         new_l_r[new_l_r_cnt*map_info->map_len + map_info->dupsets_map[i]] = 0.0;
      }
   }

   if( cmd_options->o_filename != NULL )
   {
      print_uncomp_cutset_to_file(new_l_r,new_l_r_cnt,g_uncomp_solutions_cnt, reac_info, map_info);
   }

   g_uncomp_solutions_cnt++;
   new_l_r_cnt++;
   if( new_l_r_cnt >= t_max_sol_enhancer )
   {
      t_max_sol_enhancer *= 2;
      new_l_r = (double *) realloc(new_l_r, (size_t) t_max_sol_enhancer*map_info->map_len*sizeof(double));
      if( new_l_r == NULL )
      {
         printf("FATAL ERROR: uncompress_and_write_solution(): couldn't reallocate memory for new_l_r\n");
         exit(EXIT_FAILURE);
      }
   }

   int duplicate_stop = new_l_r_cnt;

   for( d = 0; d < map_info->num_dupsets; d++ )
   {

      int duplicated_col = map_info->dupsets[d*reac_info->num_reactions + 0];

      if( new_cs[duplicated_col] == 0 )
      {
         for( m = duplicate_start; m < duplicate_stop; m++ )
         {
            // printf("   m=%d duplicate_start=%d duplicate_stop=%d\n",m,duplicate_start,duplicate_stop);
            for( x = 1; x < map_info->num_reac_dupset[d]; x++ )
            {
               for( z = 0; z < map_info->map_len; z++ )
               {
                  new_cs_dup[z] = new_l_r[m*map_info->map_len + z];
               }

               if( new_cs_dup[map_info->dupsets[d*reac_info->num_reactions + 0]] == 1 )
               {
                  printf("FATAL ERROR: uncompress_cutsets()\n");
                  exit(EXIT_FAILURE);
               }
               new_cs_dup[map_info->dupsets[d*reac_info->num_reactions + 0]] = 1;
               new_cs_dup[map_info->dupsets[d*reac_info->num_reactions + x]] = 0;

               for( z = 0; z < map_info->map_len; z++ )
               {
                  new_l_r[new_l_r_cnt*map_info->map_len + z] = new_cs_dup[z];
               }

               if( cmd_options->o_filename != NULL )
               {
                  print_uncomp_cutset_to_file(new_l_r,new_l_r_cnt,g_uncomp_solutions_cnt, reac_info, map_info);
               }
               g_uncomp_solutions_cnt++;
               new_l_r_cnt++;
               if( new_l_r_cnt >= t_max_sol_enhancer )
               {
                  t_max_sol_enhancer *= 2;
                  new_l_r = (double *) realloc(new_l_r, (size_t) t_max_sol_enhancer*map_info->map_len*sizeof(double));
                  if( new_l_r == NULL )
                  {
                     printf("FATAL ERROR: uncompress_and_write_solution(): couldn't reallocate memory for new_l_r\n");
                     exit(EXIT_FAILURE);
                  }
               }
            }
         }
      }
      duplicate_stop = new_l_r_cnt;
   }

   printf("INFO: uncompress_and_write_solution(): Total number of minimal cutsets after uncompressing: %d\n",g_uncomp_solutions_cnt);
   free(new_l_r);
   free(new_cs);
   free(new_cs_dup);

   fflush(g_fh_out);
}
/////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_uncomp_cutset_to_file(double *new_l_r, unsigned int new_l_r_cnt, unsigned int l_unc_sol_cnt, struct_reac_info *reac_info, struct_map_info *map_info)
{
   int i;
   int num_knockouts = 0;

   fprintf(g_fh_out,"%06u: ",l_unc_sol_cnt+1);

   for( i = 0; i < map_info->map_len; i++ )
   {
      if( new_l_r[new_l_r_cnt*map_info->map_len + i] < 0.5)
      {
         fprintf(g_fh_out,"\"%s\" ",reac_info->reactions[map_info->map[i]]);
         num_knockouts++;
      }
   }
   for( i = 0; i < map_info->num_always_zero_reacs; i++ )
   {
      fprintf(g_fh_out,"\"%s\" ",reac_info->reactions[map_info->always_zero_map[i]]);
      num_knockouts++;
   }
   fprintf(g_fh_out,"  num_knockouts=%d\n",num_knockouts);
}
//////////////////////////////////////////////////////
