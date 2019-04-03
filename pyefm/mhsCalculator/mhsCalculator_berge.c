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
#include "mhsCalculator_berge.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int main(int argc, char **argv)
{
   struct timeval time_start_calc;
   struct timeval time_before_prepro;
   struct timeval time_before_mincutset;
   struct timeval time_before_postpro;
   struct timeval time_before_printing;
   struct timeval time_stop_calc;
   struct_cmd_options cmd_options = {NULL, // o_filename
                                     NULL, // r_filename
                                     NULL, // m_filename
                                     NULL, // e_filename
                                     0,    // good_efms
                                     0,    // wanted_emfs
                                     0,    // lin_superset_test
                                     0,    // apply_heuristics
                                     1,    // num_threads
                                     0,    // max_cancellations
                                     0,    // rand_seed_integer
                                     0,    // solution range
                                     0};   // output_bitvector
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

   struct_cutset_info cutset_info = { NULL,     // arr_pcs
                                      INIT_PCS_ELEMS,        // max_pcs_idx
                                      0,        // pcs_idx
                                      NULL,     // arr_mcs
                                      INIT_MCS_ELEMS,        // max_mcs_idx
                                      0,        // mcs_idx
                                      NULL,     // arr_mcs2node
                                      {0},      // pcs_idx_thr[MAX_THREADS];
                                      {0},      // max_pcs_idx_thr[MAX_THREADS]
                                      {NULL}}; // arr_pcs_thr[MAX_THREADS

   gettimeofday(&time_start_calc,NULL);

   printf("INFO: program %s started.\n",argv[0]);

   handle_arguments(argc, argv, &cmd_options);


   readin_reactions_file(cmd_options.r_filename, &reac_info);
   print_reactions(reac_info);

   readin_efm_file_bin(cmd_options, reac_info, &mode_info);

#if PRINT_DEBUG_OUTPUT != 0
#endif
   print_reac_occ_arr(mode_info, reac_info, cmd_options);

   // find essential reactions by inspection of modes
   find_essential_reacs(&mode_info, cmd_options, reac_info);

#if PRINT_DEBUG_OUTPUT != 0
#endif
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

   // perform Berge algorithm
   do_mincutset_berge(&mode_info, reac_info, &map_info, cmd_options, &cutset_info);

   gettimeofday(&time_before_postpro,NULL);

   uncompress_cutsets(&mode_info, &cutset_info, &reac_info, &map_info);

   gettimeofday(&time_before_printing,NULL);

   if( cmd_options.o_filename != NULL )
   {
      print_all_cutsets_text(&cmd_options, &cutset_info, &mode_info, &reac_info, &map_info);
   }
   print_min_max_knockouts(&mode_info, &map_info, &cutset_info);
   printf("Number of found cutsets: %llu\n",cutset_info.mcs_idx);
   // printf("Number of comparison events: %lld\n",g_num_cmp_events_total);

   gettimeofday(&time_stop_calc,NULL);
#if PRINT_DEBUG_OUTPUT != 0
   // printf("Time tree management:      %lld s %lld us\n",g_usec_tree_handling/1000000,g_usec_tree_handling%1000000);
   // printf("Time candidate generation: %lld s %lld us\n",g_usec_candidate_generation/1000000,g_usec_candidate_generation%1000000);
   // printf("Time filtering:            %lld s %lld us\n",g_usec_filtering/1000000,g_usec_filtering%1000000);
#endif
   printf("Time reading files:    "); display_execution_time(time_before_prepro,time_start_calc);
   printf("Preprocessing time:    "); display_execution_time(time_before_mincutset,time_before_prepro);
   printf("Mincutset time:        "); display_execution_time(time_before_postpro,time_before_mincutset);
   printf("Postprocessing time:   "); display_execution_time(time_before_printing,time_before_postpro);
   printf("Time printing cutsets: "); display_execution_time(time_stop_calc,time_before_printing);
   printf("Total execution time:  "); display_execution_time(time_stop_calc,time_start_calc);

   // clean up
   do_frees(&reac_info, &mode_info, &map_info, &cutset_info, cmd_options);
   printf("INFO: program %s stopped.\n",argv[0]);
   return(0);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void do_frees(struct_reac_info *reac_info, struct_mode_info *mode_info, struct_map_info *map_info,
             struct_cutset_info *cutset_info, struct_cmd_options cmd_options)
{
   free_mem_reac_info(reac_info);
   free_mem_mode_info(mode_info);
   free_mem_map_info(map_info);
   free_mem_cutset_info(cutset_info, cmd_options);
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
   char *char_rand_seed         = NULL;
   char *char_endpointer;

   printf("mhsCalculator_berge Version: %s\n",VERSION);
   printf("Copyright (C) %s %s %s\n",YEARS, AUTHORS, COMPANY);

   printf("Executed command: %s\n",largv[0]);
   printf("Options:");
   for(i = 1; i < largc; i++ )
   {
      printf(" %s",largv[i]);
   }
   printf("\n");

   while(( r_opt = getopt(largc, largv, "hblkm:r:e:o:n:w:t:c:s:")) != -1 )
   {
      switch(r_opt)
      {
         case 'h':
            usage("");
            break;
         case 'b':
            cmd_options->output_bitvector = 1;
            break;
         case 'l':
            cmd_options->lin_superset_test = 1;
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
         case 'c':
            char_max_cancellations = optarg;
            break;
         case 's':
            char_rand_seed = optarg;
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
   if( cmd_options->num_threads > MAX_THREADS )
   {
      fprintf(stderr, "FATAL ERROR: number of threads provided as argument (%d) is larger than MAX_THREADS (%d)\n",cmd_options->num_threads,MAX_THREADS);
      fprintf(stderr, "             Increase MAX_THREADS and recompile\n");
      exit(EXIT_FAILURE);
   }
   else if( cmd_options->num_threads <= 0 )
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
   // inter for seeding randum number generation
   ////////////////////////////////////////////////////////////////////////////
   if( char_rand_seed != NULL )
   {
      cmd_options->rand_seed_integer = strtol(char_rand_seed, &char_endpointer, 10);
      if( char_endpointer == char_rand_seed || errno == EINVAL || errno == ERANGE )
      {
         fprintf(stderr, "FATAL ERROR: error while converting random seed number (-s %s): %s\n",char_rand_seed,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }
      printf("INFO: integer value for seeding random number generation: %d\n",cmd_options->rand_seed_integer);
      srandom(cmd_options->rand_seed_integer);
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
   return;
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void usage(char *message)
{
   printf("%s",message);
   printf("usage: mhsCalculator_berge -m efms.sorted.bin -r rfile [-e essential_reacs.txt] [-o output_cutsets] ");
   printf("-n num_good_modes [-l] [-k] [-w num_wanted_good_modes] [-t num_threads] [-c max cancellations] [-s integer] [-b]\n");
   printf("\n");
   printf("-m ..... filename containing elementary flux modes (in binary form!)\n");
   printf("-r ..... filename containing names of reactions\n");
   printf("-e ..... filename containing essential reactions\n");
   printf("-o ..... filename (output) of file containing computed minimal cutsets\n");
   printf("-n ..... number of good/keeper modes\n");
   printf("-w ..... number of wanted good modes (number of modes that must survive if cutset is applied)\n");
   printf("-t ..... number of parallel threads\n");
   printf("-k ..... bail out of duplicate mode check if it seems to be ineffective\n");
   printf("-b ..... print cutsets as bitvector string\n");
   printf("-c ..... maximum number of cancelations\n");
   printf("-l ..... use linear approach to find subsets, by default a tree search approach is used\n");
   printf("-s ..... seed value for random number generation\n");

   printf("-h ..... print this help message\n");

   exit(EXIT_FAILURE);
}
//////////////////////////////////////////////////////
