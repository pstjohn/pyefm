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

#include "bruteforce.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void do_mincutset_bruteforce(struct_mode_info *mode_info, struct_reac_info *reac_info, struct_map_info *map_info, struct_cmd_options *cmd_options, struct_cutset_info *cutset_info)
{
   int c,i;
   int iret[MAX_THREADS];
   pthread_t thread[MAX_THREADS];
   unsigned long long int *t_cutset;
   unsigned long long int *m_cutset;
   struct_tree_info tree_info = {NULL,    // nodes
                                 0,       // node_max_idx
                                 0,       // node_idx
                                 NULL };  // bitorder

   struct_bitorder_info bitorder_info = {NULL,  // struct_bitorder
                                         NULL}; // reac_occ

   if( cmd_options->o_filename != NULL )
   {
      g_fh_out = fopen(cmd_options->o_filename, "w");

      if( g_fh_out <= 0 )
      {
         printf("FATAL ERROR: open file '%s' for writing failed: %s\n",cmd_options->o_filename,strerror(errno));
         exit(EXIT_FAILURE);
      }
   }

   allocate_bruteforce_mem(cmd_options, map_info, mode_info);

   init_cs_arr_bruteforce(cutset_info, &tree_info, mode_info, *cmd_options, *map_info, &bitorder_info, *reac_info);

   // change order of reactions of modes
   // reorder_reaction_of_modes(mode_info, *map_info, bitorder_info, cmd_options);

   t_cutset = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int) );
   m_cutset = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int) );
   if( t_cutset == NULL || m_cutset == NULL )
   {
      fprintf(stderr, "FATAL ERROR: do_mincutset_berge(): couldn't allocate memory for either t_cutset pr m_cutset\n");
      exit(EXIT_FAILURE);
   }
   ///////////////////////////////////////////////////


   for( c = 1; c <= cmd_options->max_cancellations - map_info->num_always_zero_reacs; c++ )
   {
      printf("INFO: trying to find cutsets with %d knockouts\n", c + map_info->num_always_zero_reacs);

      // fill global n-over-k array for speed up
      fill_n_over_k(c, cmd_options->max_cancellations, map_info->num_dupset_keepers, map_info->num_always_zero_reacs );

      set_max_bruteforce_combinations(c,map_info->num_dupset_keepers);
      g_bruteforce_depth = c;
      g_invoke_cnt = 0;

      for( i = 0; i < cmd_options->num_threads; i++ )
      {
         g_strInThread[i].thread_id   = i;
         g_strInThread[i].mode_info   = mode_info;
         g_strInThread[i].map_info    = map_info;
         g_strInThread[i].reac_info   = reac_info;
         g_strInThread[i].cutset_info = cutset_info;
         g_strInThread[i].tree_info   = &tree_info;
         g_strInThread[i].cmd_options = cmd_options;
      }

      // create threads
      for( i = 0; i < cmd_options->num_threads; i++ )
      {
         iret[i] = pthread_create(&thread[i], NULL, worker_bruteforce, (void *) &g_strInThread[i]);
      }


      // wait for all threads
      for( i = 0; i < cmd_options->num_threads; i++ )
      {
         pthread_join( thread[i], NULL);
      }

      printf("found %ld compressed sets of knockouts resulting in %lld cutsets\n",g_found_cutsets,g_uncomp_solutions_cnt);
      fill_precutsets_to_cutsets(&tree_info, cutset_info, mode_info, *map_info, &bitorder_info, *cmd_options);
   }

   free_bruteforce_mem();

   if( cmd_options->o_filename != NULL )
   {
      fclose(g_fh_out);
   }

   // reorder_reaction_of_cutsets(cutset_info, map_info, mode_info, &bitorder_info);

   free(m_cutset);
   free(t_cutset);
   free_mem_tree_info(&tree_info);
   free_mem_bitorder_info(&bitorder_info);
   return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void allocate_bruteforce_mem(struct_cmd_options *cmd_options, struct_map_info *map_info, struct_mode_info *mode_info)
{
   g_bruteforce_cutset = (unsigned *) malloc(sizeof(unsigned)*cmd_options->num_threads*cmd_options->max_cancellations);
   if( g_bruteforce_cutset == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_bruteforce_mem(): couldn't allocate memory for g_bruteforce_cutset\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_n_over_k = (unsigned long long int *) malloc(sizeof(unsigned long long int)*(cmd_options->max_cancellations - map_info->num_always_zero_reacs + 1)*(map_info->num_dupset_keepers + 1));
   if( g_n_over_k == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_bruteforce_mem(): couldn't allocate memory for g_n_over_k\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   g_pnt_tmp_cutset2 = (unsigned long long int *) malloc(sizeof(unsigned long long int)*cmd_options->num_threads*mode_info->num_unit_size);
   if( g_pnt_tmp_cutset2 == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_bruteforce_mem(): couldn't allocate memory for g_pnt_tmp_cutset2\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void free_bruteforce_mem()
{
   free(g_bruteforce_cutset);
   free(g_n_over_k);
   free(g_pnt_tmp_cutset2);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
unsigned long long int n_over_k(int n, int k)
{
   int i;
   unsigned long long int num_cmbs = 1;
   int n_help = n;

   for( i = 0; i < k; i++ )
   {
      num_cmbs *= n_help;
      n_help--;
      num_cmbs /= (i+1);
   }

   return(num_cmbs);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void *worker_bruteforce(void *str_ptr)
{
   unsigned long long int c;
   int thr_id;
   struct_mode_info   *mode_info;
   struct_map_info    *map_info;
   struct_reac_info   *reac_info;
   struct_cutset_info *cutset_info;
   struct_tree_info   *tree_info;
   struct_cmd_options *cmd_options;
   int depth;
   struct structInThread *str;
   str = (struct structInThread *) str_ptr;
   thr_id       = (*str).thread_id;
   mode_info    = (*str).mode_info;
   map_info     = (*str).map_info;
   reac_info    = (*str).reac_info;
   cutset_info  = (*str).cutset_info;
   tree_info    = (*str).tree_info;
   cmd_options  = (*str).cmd_options;


   // printf("INFO: in worker_bruteforce() subroutine for thread %d\n",thr_id);

   int n,k,j;
   unsigned long long int limit;
   unsigned long long int z;
   time_t t_cur_time;
   double diff_time;
   double expect_runtime;
   g_start_time = time(NULL);

   for( c = 0 + thr_id; c < g_num_bruteforce_combinations; c += cmd_options->num_threads )
   {
      // printf("thread %d: c=%llu g_num_bruteforce_combinations=%lu g_bruteforce_depth=%d\n",thr_id,c,g_num_bruteforce_combinations,g_bruteforce_depth);
      n = map_info->num_dupset_keepers - 1;

      k = g_bruteforce_depth - map_info->num_always_zero_reacs;
      z = c;
      depth = 0;
      for( j = 0; j < map_info->num_dupset_keepers && depth < g_bruteforce_depth; j++ )
      {
         // limit = n_over_k(n,k);
         limit = g_n_over_k[n*(cmd_options->max_cancellations - map_info->num_always_zero_reacs + 1) + k];
         if( z < limit )
         {
            n--;
         }
         else
         {
            // g_bruteforce_cutset[thr_id*mode_info->num_unit_size + depth] = j;
            g_bruteforce_cutset[thr_id*cmd_options->max_cancellations + depth] = j;
            depth++;
            n--;
            k--;
            z -= limit;
         }
         // printf("INFO: thread=%d c=%llu j=%d g_map_len=%d limit=%llu n=%d k=%d z=%llu depth=%d\n",thr_id,c,j,g_map_len,limit,n,k,z,depth);
      }

      evaluate_cutset_bruteforce(thr_id, cmd_options, map_info, reac_info, mode_info, cutset_info, tree_info);

      if( c > 0 && c%500000 == 0 )
      {
         t_cur_time = time(NULL);
         diff_time = difftime( t_cur_time, g_start_time );
         expect_runtime = diff_time*g_num_bruteforce_combinations/(float)c/3600;

         printf("thr_id=%d processed %llu of %lu combinations (%5.2f%%) spent time=%es expected runtime=%e h\n",thr_id,c,g_num_bruteforce_combinations,100*(float)(c)/(float)g_num_bruteforce_combinations,diff_time,expect_runtime);
      }
   }

   return((void *)NULL);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void evaluate_cutset_bruteforce(int thread_id, struct_cmd_options *cmd_options, struct_map_info *map_info, struct_reac_info *reac_info,
                                struct_mode_info *mode_info, struct_cutset_info *cutset_info, struct_tree_info *tree_info)
{
   int i,d;
   int unit_cnt;
   int idx;
   int bitmover;
   int knockedout;
   unsigned long long tmp_unit;
   unsigned long m;
   unsigned long long bad_kept  = 0;
   unsigned long long good_kept = 1;
   int i_pcs = thread_id;
   int max_cancellations = cmd_options->max_cancellations;
   char *filename = cmd_options->o_filename;

   // reset global temporary bit-stored cutset
   for( i = 0; i < mode_info->num_unit_size; i++ )
   {
      // g_pnt_tmp_cutset[thread_id*UNIT_LENGTH*sizeof(unsigned long long) + i] = 0;
      g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size + i] = 0;
   }

   // create bitmap of current cutset
   // printf("candidate: ");
   for( d = 0; d < g_bruteforce_depth; d++ )
   {
      tmp_unit = 1;
      // idx = g_bruteforce_cutset[thread_id*mode_info->num_unit_size + d];
      idx = g_bruteforce_cutset[thread_id*max_cancellations + d];
      // printf("%d ",idx);
      bitmover = idx;

      unit_cnt = bitmover/(8*sizeof(unsigned long long));
      bitmover -= unit_cnt*8*sizeof(unsigned long long);

      tmp_unit <<= bitmover;
      g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size + unit_cnt] |= tmp_unit;
   }
   // printf(": %lld\n",g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size]);

   // iterate through all modes
   for( m = 0; m < mode_info->bad_efms; m++ )
   {
      knockedout = 0;
      for( i = 0; i < mode_info->num_unit_size; i++ )
      {
         if( (mode_info->pnt_flux[m*mode_info->num_unit_size + i] & g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size + i]) != 0 )
         {
            knockedout = 1;
            break;
         }
      }

      if( knockedout == 0 )
      {
         bad_kept = 1;
         break;
      }
   }

   if( bad_kept == 0 )
   {
      // bad mode check was ok -> not a single bad mode survived
      // check if enough good modes will survive cutset
      if( cmd_options->good_efms_wanted > 0 && cmd_options->good_efms_wanted < cmd_options->good_efms )
      {
         // note that i_pcs is set to thread_id, as we always check only one mode
         // which is store at the position 'thread_id'
         if( enough_good_modes_survive(&g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size], *cmd_options, mode_info) )
         {
            good_kept = 1;
         }
         else
         {
            good_kept = 0;
         }
      }
   }

   if( bad_kept == 0 && good_kept == 1 )
   {
      // printf("INFO: we found a potential cutset -> superset test needs to be done\n");

      // copy currently check cutset to the arr_pcs[] at position i_pcs = thr_id
      // the do a superset search
      for( i = 0; i < mode_info->num_unit_size; i++ )
      {
         // cutset_info->arr_pcs[i_pcs*mode_info->num_unit_size + i] = g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size + i];
         cutset_info->arr_pcs[i_pcs*mode_info->num_unit_size + i] = g_pnt_tmp_cutset2[thread_id*mode_info->num_unit_size + i];
      }

      // we found a proper candidate cutset -> do superset test
      int is_not_superset = 1;

      if( cutset_info->mcs_idx > 0 )
      {
         // printf("INFO: going to call find_superset_in_tree(): cutset_info->mcs_idx=%lld\n",cutset_info->mcs_idx);
         // only do superset check if there where any already computed cutsets to compare with
         if( find_superset_in_tree(0, i_pcs, thread_id, tree_info, cutset_info, mode_info) )
         {
            is_not_superset = 0;
         }
      }

      if( is_not_superset == 1)
      {
         //////////////////////////////////////////////////////////
         // store new cutset in arr_pcs_thr[][] of this thread
         //////////////////////////////////////////////////////////
         for( i = 0; i < mode_info->num_unit_size; i++ )
         {
            cutset_info->arr_pcs_thr[thread_id][cutset_info->pcs_idx_thr[thread_id]*mode_info->num_unit_size + i] = cutset_info->arr_pcs[i_pcs*mode_info->num_unit_size + i];
         }
         cutset_info->pcs_idx_thr[thread_id]++;
         //////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////
         // if array is not big enough anymore -> resize it
         //////////////////////////////////////////////////////////
         if( cutset_info->pcs_idx_thr[thread_id] >= cutset_info->max_pcs_idx_thr[thread_id] )
         {
            cutset_info->arr_pcs_thr[thread_id] = (unsigned long long int *) realloc(cutset_info->arr_pcs_thr[thread_id],
                                     (size_t) mode_info->num_unit_size*(cutset_info->max_pcs_idx_thr[thread_id]*2)*sizeof(unsigned long long int) );
            if( cutset_info->arr_pcs_thr[thread_id] == NULL )
            {
               fprintf(stderr, "FATAL ERROR: cutset_info->arr_pcs_thr[%d] is not big enough\n",thread_id);
               fprintf(stderr, "             cutset_info->max_pcs_idx_thr[%d]=%llu\n",thread_id,cutset_info->max_pcs_idx_thr[thread_id]);
               fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
               exit(EXIT_FAILURE);
            }
            cutset_info->max_pcs_idx_thr[thread_id] *= 2;
         }
         //////////////////////////////////////////////////////////

         // new cutset is not a superset of any existing cutset -> add it
         pthread_mutex_lock( &mutex1 );
         g_found_cutsets++;
         // printf("found a working set of reactions thr_id=%d: %ld:\n",thread_id,g_found_cutsets);
         // int k;
         // for( k = 0; k < g_bruteforce_depth; k++ )
         // {
         //     printf("\"%s\" ",reac_info->reactions[map_info->dupsets_map[g_bruteforce_cutset[thread_id*max_cancellations + k]]]);
         // }
         // for( k = 0; k < map_info->num_always_zero_reacs; k++ )
         // {
         //    printf("\"%s\" ",reac_info->reactions[map_info->always_zero_map[k]]);
         // }
         // printf("\n");

         uncompress_and_write_solution(thread_id, g_found_cutsets, map_info, reac_info, mode_info, filename, max_cancellations);

         pthread_mutex_unlock( &mutex1 );
      }
   }

   // printf("   good_kept=%llu good_rmvd=%llu bad_kept=%llu bad_rmvd=%llu\n",good_kept,good_rmvd,bad_kept,bad_rmvd);
   return;
}
//////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void fill_n_over_k(int cur_num_cancellations, int max_cancellation, int num_dupset_keepers, int num_always_zero_reacs )
{
   int n;
   int k;

   for( n = 0; n <= num_dupset_keepers; n++ )
   {
      // for( k = 0; k <= max_cancellation - num_always_zero_reacs; k++ )
      for( k = 0; k <= cur_num_cancellations - num_always_zero_reacs; k++ )
      {
         // g_n_over_k[n][k] = n_over_k(n,k);
         g_n_over_k[n*(max_cancellation - num_always_zero_reacs + 1) + k] = n_over_k(n,k);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void set_max_bruteforce_combinations(int nr_cancels, int num_reacs)
{
   int i;
   int t_num_reacs;

   t_num_reacs = num_reacs;
   // initialize global combination number
   g_num_bruteforce_combinations = 1;

   for( i = 0; i < nr_cancels; i++ )
   {
      g_num_bruteforce_combinations *= num_reacs;
      num_reacs--;
      g_num_bruteforce_combinations /= (i+1);
   }
   printf("INFO: %d cancellations of %d relevant reaction results in %lu (brute force) combinations\n",nr_cancels,t_num_reacs,g_num_bruteforce_combinations);
}
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_cs_arr_bruteforce(struct_cutset_info *cutset_info, struct_tree_info *tree_info, struct_mode_info *mode_info, struct_cmd_options cmd_options,
     struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_reac_info reac_info)
{
   int j;
   int u;

   printf("DEBUG: init_cs_arr_bruteforce(): entered\n");

   cutset_info->arr_pcs = (unsigned long long int *) malloc( (size_t) cutset_info->max_pcs_idx*mode_info->num_unit_size*sizeof(unsigned long long int));
   if( cutset_info->arr_pcs == NULL )
   {
      fprintf(stderr, "FATAL ERROR: could not allocate memory for cutset_info->arr_pcs\n");
      exit(EXIT_FAILURE);
   }

   cutset_info->arr_mcs = (unsigned long long int *) malloc( (size_t) cutset_info->max_mcs_idx*mode_info->num_unit_size*sizeof(unsigned long long int));
   if( cutset_info->arr_mcs == NULL )
   {
      fprintf(stderr, "FATAL ERROR: could not allocate memory for cutset_info->arr_mcs\n");
      exit(EXIT_FAILURE);
   }

   cutset_info->arr_mcs2node = (unsigned long long int *) malloc( (size_t) cutset_info->max_mcs_idx*sizeof(unsigned long long int));
   if( cutset_info->arr_mcs2node == NULL )
   {
      fprintf(stderr, "FATAL ERROR: could not allocate memory for cutset_info->arr_mcs2node\n");
      exit(EXIT_FAILURE);
   }

   cutset_info->xor_list = (struct_xor_list *) malloc((size_t) sizeof(struct_xor_list)*reac_info.num_reactions);
   if( cutset_info->xor_list == NULL )
   {
      fprintf(stderr, "FATAL ERROR: init_cs_arr_bruteforce(): couldn't allocate %lu bytes for cutset_info->xor_list\n",sizeof(struct_xor_list)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }


   determine_tree_bit_order(map_info, *mode_info, bitorder_info);

   ////////////////////////////////////////////////////////////////////////////
   // alocate memory for tree
   ////////////////////////////////////////////////////////////////////////////
   tree_info->node_max_idx = INIT_MAX_NODES;
   tree_info->nodes = (struct_tree_node *) malloc( (size_t) tree_info->node_max_idx*sizeof(struct_tree_node));
   if( tree_info->nodes == NULL )
   {
      fprintf(stderr, "FATAL ERROR: init_cs_arr_bruteforce(): couldn't allocate memory for tree\n");
      exit(EXIT_FAILURE);
   }
   tree_info->node_idx = 0;
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // create root node
   ////////////////////////////////////////////////////////////////////////////
   tree_info->nodes[tree_info->node_idx].parent   = ROOT;
   tree_info->nodes[tree_info->node_idx].left     = EMPTY;
   tree_info->nodes[tree_info->node_idx].right    = EMPTY;
   tree_info->nodes[tree_info->node_idx].type     = NON_LEAF;
   tree_info->nodes[tree_info->node_idx].depth    = 0;
   tree_info->nodes[tree_info->node_idx].bitpos   = 0;
   tree_info->nodes[tree_info->node_idx].unit_cnt = 0;
   tree_info->nodes[tree_info->node_idx].bit_map  = 1;
   tree_info->nodes[tree_info->node_idx].i_mcs    = -1; // the root node is never assigned to a cutset
   tree_info->nodes[tree_info->node_idx].cs     = (unsigned long long int*) calloc((size_t) mode_info->num_unit_size, (size_t) sizeof(unsigned long long int));
   if( tree_info->nodes[tree_info->node_idx].cs == NULL )
   {
      fprintf(stderr, "FATAL ERROR: init_cs_arr_bruteforce(): couldn't allocate memory for cutset (cs) of tree node %lu\n",tree_info->node_idx);
      exit(EXIT_FAILURE);
   }
   for( u = 0; u < mode_info->num_unit_size; u++ )
   {
      tree_info->nodes[tree_info->node_idx].cs[u] = 0xFFFFFFFFFFFFFFFFULL;
   }

   tree_info->node_idx++;

   CHECK_AND_RESIZE_G_NODE
   ////////////////////////////////////////////////////////////////////////////


   for( j = 0; j < cmd_options.num_threads; j++ )
   {
      cutset_info->max_pcs_idx_thr[j] = INIT_PCS_ELEMS_THR;
      printf("DEBUG: init_cs_arr_bruteforce(): cutset_info->max_pcs_idx_thr[%d]=%llu\n",j,cutset_info->max_pcs_idx_thr[j]);
      cutset_info->arr_pcs_thr[j] = (unsigned long long int *) malloc( (size_t) cutset_info->max_pcs_idx_thr[j]*mode_info->num_unit_size*sizeof(unsigned long long int) );
      if( cutset_info->arr_pcs_thr[j] == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocating memory for arr_pcs_thr[%d] failed\n",j);
         exit(EXIT_FAILURE);
      }
   }
   printf("DEBUG: init_cs_arr_bruteforce(): leaving\n");
   fflush(stdout);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_cutset_info(struct_cutset_info *cutset_info, struct_cmd_options cmd_options)
{
   int i;
   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      free(cutset_info->arr_pcs_thr[i]);
   }
   free(cutset_info->arr_pcs);
   free(cutset_info->arr_mcs);
   free(cutset_info->arr_mcs2node);
   free(cutset_info->xor_list);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void reorder_reaction_of_modes(struct_mode_info *mode_info, struct_map_info map_info, struct_bitorder_info bitorder_info, struct_cmd_options cmd_options)
{
   long int e;
   int u;
   int r_in, r_out;
   int unit_cnt_in, unit_cnt_out;
   int bitmover_in, bitmover_out;
   unsigned long long int bitmask_in, bitmask_out;
   unsigned long long int *store;


   store = (unsigned long long int *) malloc(mode_info->num_unit_size*sizeof(unsigned long long int));

   if( store == NULL )
   {
      fprintf(stderr,"FATAL ERROR: reorder_reaction_of_modes(): couldn't allocate memory for tmp\n");
      exit(EXIT_FAILURE);
   }

   for( e = 0; e < mode_info->bad_efms; e++ )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         store[u] = mode_info->pnt_flux[e*mode_info->num_unit_size + u];
         mode_info->pnt_flux[e*mode_info->num_unit_size + u] = 0;
      }
      for( r_out = 0; r_out < map_info.num_dupset_keepers; r_out++ )
      {
         r_in = bitorder_info.bitorder[r_out].bitpos;
         bitmask_in = 1;
         bitmover_in = r_in;
         unit_cnt_in = bitmover_in/(8*sizeof(unsigned long long));
         bitmover_in -= unit_cnt_in*8*sizeof(unsigned long long);
         bitmask_in <<= bitmover_in;

         if( store[unit_cnt_in] & bitmask_in )
         {
            // bit was set -> move it
            bitmask_out = 1;
            bitmover_out = r_out;
            unit_cnt_out = bitmover_out/(8*sizeof(unsigned long long));
            bitmover_out -= unit_cnt_out*8*sizeof(unsigned long long);
            bitmask_out <<= bitmover_out;
            mode_info->pnt_flux[e*mode_info->num_unit_size + unit_cnt_out] |= bitmask_out;
         }
      }
   }

   for( e = 0; e < mode_info->good_emfs; e++ )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         store[u] = mode_info->pnt_good[e*mode_info->num_unit_size + u];
         mode_info->pnt_good[e*mode_info->num_unit_size + u] = 0;
      }
      for( r_out = 0; r_out < map_info.num_dupset_keepers; r_out++ )
      {
         r_in = bitorder_info.bitorder[r_out].bitpos;
         bitmask_in = 1;
         bitmover_in = r_in;
         unit_cnt_in = bitmover_in/(8*sizeof(unsigned long long));
         bitmover_in -= unit_cnt_in*8*sizeof(unsigned long long);
         bitmask_in <<= bitmover_in;

         if( store[unit_cnt_in] & bitmask_in )
         {
            // bit was set -> move it
            bitmask_out = 1;
            bitmover_out = r_out;
            unit_cnt_out = bitmover_out/(8*sizeof(unsigned long long));
            bitmover_out -= unit_cnt_out*8*sizeof(unsigned long long);
            bitmask_out <<= bitmover_out;
            mode_info->pnt_good[e*mode_info->num_unit_size + unit_cnt_out] |= bitmask_out;
         }
      }
   }

   free(store);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int compute_precutsets(long long int e, unsigned long long int * t_cutset, unsigned long long int *m_cutset,
                       struct_mode_info *mode_info, struct_map_info *map_info, struct_cutset_info *cutset_info, 
                       struct_tree_info *tree_info, struct_cmd_options cmd_options)
{
   int hit;
   int u,v,b;
   int t_array_changed = 0;
   long int i_mcs;
   int num_ones;
   long unsigned int eliminated_mcs    = 0;
   long unsigned int new_potential_pcs = 0;

   for( u = 0; u < mode_info->num_unit_size; u++ )
   {
      t_cutset[u] = 0;
      m_cutset[u] = 0;
   }

   for(i_mcs = 0; i_mcs < cutset_info->mcs_idx; i_mcs++ )
   {
      hit = 0;
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         if( ((mode_info->pnt_flux[e*mode_info->num_unit_size + u] & cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + u]) != 0 ) )
         {
            hit = 1;
            break;
         }
      }

      if( hit == 0 )
      {
         t_array_changed = 1;
         eliminated_mcs++;

         if( cmd_options.lin_superset_test == 0 )
         {
            delete_from_tree(cutset_info->arr_mcs2node[i_mcs],i_mcs, tree_info, mode_info, cutset_info);
         }

         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            m_cutset[u] = cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + u];
         }

         if( i_mcs != cutset_info->mcs_idx - 1 )
         {
            if( cmd_options.lin_superset_test == 0 )
            {
               cutset_info->arr_mcs2node[i_mcs] = cutset_info->arr_mcs2node[(cutset_info->mcs_idx - 1)];
               tree_info->nodes[cutset_info->arr_mcs2node[i_mcs]].i_mcs = i_mcs;
            }
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + u] = cutset_info->arr_mcs[(cutset_info->mcs_idx - 1)*mode_info->num_unit_size + u];
            }
            i_mcs--;
         }
         cutset_info->mcs_idx--;
         num_ones = calc_num_ones(m_cutset, mode_info->num_unit_size);

         if( (cmd_options.max_cancellations == 0) || (num_ones < cmd_options.max_cancellations - map_info->num_always_zero_reacs) )
         {
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               if( mode_info->pnt_flux[e*mode_info->num_unit_size + u] != 0 )
               {
                  unsigned long long tmp_unit = 0x8000000000000000ULL;
                  for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
                  {
                     if( mode_info->pnt_flux[e*mode_info->num_unit_size + u] & tmp_unit )
                     {
                        new_potential_pcs++;
                        for( v = 0; v < mode_info->num_unit_size; v++ )
                        {
                           t_cutset[v] = m_cutset[v];
                        }
                        t_cutset[u] = t_cutset[u] | tmp_unit;
                        for( v = 0; v < mode_info->num_unit_size; v++ )
                        {
                           cutset_info->arr_pcs[cutset_info->pcs_idx*mode_info->num_unit_size + v] = t_cutset[v];
                        }
                        cutset_info->pcs_idx++;
                        check_and_resize_arr_pcs(cutset_info, mode_info);
                     }
                     tmp_unit >>= 1;
                  }
               }
            }
         }
      }
   }

   return t_array_changed;

}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int calc_num_ones(unsigned long long *cs, unsigned int num_unit_size)
{
   int u;
   int b;
   int num_ones = 0;

   if( num_unit_size == 1 )
   {
      unsigned long long tmp_unit = 0x8000000000000000ULL;

      for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
      {
         if( cs[0] & tmp_unit )
         {
            num_ones++;
         }
         tmp_unit >>= 1;
      }
   }
   else
   {
      for( u = num_unit_size - 1; u >= 0; u-- )
      {
         unsigned long long tmp_unit = 0x8000000000000000ULL;
         for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
         {
            if( cs[u] & tmp_unit )
            {
               num_ones++;
            }
            tmp_unit >>= 1;
         }
      }
   }

   return(num_ones);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void check_and_resize_arr_pcs(struct_cutset_info *cutset_info, struct_mode_info *mode_info)
{
   // printf("DEBUG: entered check_and_resize_arr_pcs()\n");
   if( cutset_info->pcs_idx >= cutset_info->max_pcs_idx )
   {
      cutset_info->max_pcs_idx *= 2;
      cutset_info->arr_pcs = (unsigned long long int *) realloc(cutset_info->arr_pcs, sizeof(unsigned long long int)*(cutset_info->max_pcs_idx)*mode_info->num_unit_size);
      if( cutset_info->arr_pcs == NULL )
      {
         fprintf(stderr, "FATAL ERROR: cutset_info->arr_pcs is not big enough!\n");
         fprintf(stderr, "             reallocating memory failed\n");
         exit(EXIT_FAILURE);
      }
   }
   // printf("DEBUG: leaving check_and_resize_arr_pcs()\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void check_and_resize_arr_mcs(struct_cutset_info *cutset_info, struct_mode_info *mode_info)
{
   // printf("DEBUG: entered check_and_resize_arr_mcs()\n");
   if( cutset_info->mcs_idx >= cutset_info->max_mcs_idx )
   {
      cutset_info->max_mcs_idx *= 2;
      cutset_info->arr_mcs2node = (unsigned long long int *) realloc(cutset_info->arr_mcs2node, (size_t) (cutset_info->max_mcs_idx)*sizeof(unsigned long long int));
      cutset_info->arr_mcs = (unsigned long long int *) realloc(cutset_info->arr_mcs, (size_t) (cutset_info->max_mcs_idx)*mode_info->num_unit_size_orig*sizeof(unsigned long long int));
      if( cutset_info->arr_mcs == NULL || cutset_info->arr_mcs2node == NULL )
      {
         fprintf(stderr, "FATAL ERROR: could not re-allocate memory for cutset_info->arr_mcs\n");
         exit(EXIT_FAILURE);
      }
   }

   // printf("DEBUG: leaving check_and_resize_arr_mcs()\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int enough_good_modes_survive(unsigned long long *cs, struct_cmd_options cmd_options, struct_mode_info *mode_info)
{
   unsigned long long int e;
   unsigned long long int survived_modes = mode_info->num_good_removed_because_all_zero;
   unsigned long long int killed_modes   = mode_info->num_good_removed_by_always_zero;

   if( mode_info->good_emfs == 0 )
   {
      // because of pre-processing all good modes might be gone
      // this means that there will always be enough good modes
      // whatever 'potential' cut set is applied
      return(1);
   }

   // loop over all good modes
   for( e = 0; e < mode_info->good_emfs; e++ )
   {

      if( is_mode_killed(cs,mode_info->pnt_good,e, mode_info) )
      {
         killed_modes += mode_info->pnt_good_duplicates[e];
      }
      else
      {
         // survived_modes++;
         survived_modes += mode_info->pnt_good_duplicates[e];
      }

      if( survived_modes >= cmd_options.good_efms_wanted )
      {
         return(1);
      }

      if( killed_modes > (cmd_options.good_efms - cmd_options.good_efms_wanted) )
      {
         return(0);
      }
   }

   // we should never end up here
   fprintf(stderr, "FATAL INTERNAL ERROR: enough_good_modes_survive()\n");
   fprintf(stderr, "                      e=%llu, survived_modes=%llu, killed_modes=%llu\n",e,survived_modes,killed_modes);
   fprintf(stderr, "                      num_good_removed_because_all_zero=%llu\n",mode_info->num_good_removed_because_all_zero);
   fprintf(stderr, "                      num_good_removed_by_always_zero=%llu\n",mode_info->num_good_removed_by_always_zero);
   fprintf(stderr, "                      mode_info->good_emfs=%llu\n",mode_info->good_emfs);
   fprintf(stderr, "                      cmd_options.good_emfs=%llu, cmd_options.good_emfs_wanted=%llu\n",cmd_options.good_efms,cmd_options.good_efms_wanted);
   fprintf(stderr, "                      execution aborted.\n");
   exit(EXIT_FAILURE);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int is_mode_killed(unsigned long long *cs, unsigned long long *t_pnt_flux, int idx, struct_mode_info *mode_info)
{
   int u;

   if( mode_info->num_unit_size == 1 )
   {
      if( (cs[0] & t_pnt_flux[idx]) != 0 )
      {
         // mode is killed
         return(1);
      }
      else
      {
         // mode is not killed
         return(0);
      }
   }
   else
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         if( (cs[u] & t_pnt_flux[idx*mode_info->num_unit_size + u]) != 0 )
         {
            // mode is killed
            return(1);
         }
      }
      // all AND-operation gave 0 ->
      // no cutset-reaction coincides with a mode reaction ->
      // mode is not killed
      return(0);
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// us1 stands for unit size is equal to 1
//////////////////////////////////////////////////////
void fill_precutsets_to_cutsets_us1(struct_tree_info *tree_info, struct_cutset_info *cutset_info, struct_mode_info *mode_info,
                                    struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_cmd_options cmd_options)
{
   int i_thr;
   long int i_pcs;
   unsigned long long int filled_in_cutsets = 0;

   // printf("DEBUG: fill_precutsets_to_cutsets_us1(): entered\n");
   // fflush(stdout);
   for( i_thr = 0; i_thr < cmd_options.num_threads; i_thr++ )
   {
      for( i_pcs = 0; i_pcs < cutset_info->pcs_idx_thr[i_thr]; i_pcs++ )
      {
         cutset_info->arr_mcs[cutset_info->mcs_idx] = cutset_info->arr_pcs_thr[i_thr][i_pcs];
         filled_in_cutsets++;

         // add new mcs to tree
         if( cmd_options.lin_superset_test == 0 )
         {
            add_to_tree(0,cutset_info->mcs_idx, tree_info, cutset_info, mode_info, map_info, bitorder_info);
         }

         cutset_info->mcs_idx++;

         check_and_resize_arr_mcs(cutset_info, mode_info);
         // CHECK_AND_RESIZE_ARR_MCS
      }
      cutset_info->pcs_idx_thr[i_thr] = 0;
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void fill_precutsets_to_cutsets(struct_tree_info *tree_info, struct_cutset_info *cutset_info, struct_mode_info *mode_info,
                                struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_cmd_options cmd_options)
{
   int u;
   int i_thr;
   long int i_pcs;
   unsigned long long int filled_in_cutsets = 0;

   // printf("DEBUG: fill_precutsets_to_cutsets(): entered\n");
   // fflush(stdout);

   for( i_thr = 0; i_thr < cmd_options.num_threads; i_thr++ )
   {
      // printf("DEBUG: i_thr=%d cutset_info->pcs_idx_thr[i_thr]=%d\n",i_thr,cutset_info->pcs_idx_thr[i_thr]);
      for( i_pcs = 0; i_pcs < cutset_info->pcs_idx_thr[i_thr]; i_pcs++ )
      {
         //if(i_pcs%10000 == 0 )
         //{
         //  printf("thr_id:%d/%d %llu/%llu=%04.1f%% --> %04.1f%%\n",
         //            i_thr,cmd_options.num_threads,i_pcs,cutset_info->pcs_idx_thr[i_thr],100.0*i_pcs/cutset_info->pcs_idx_thr[i_thr],100.0*i_pcs/cutset_info->pcs_idx_thr[i_thr]/cmd_options.num_threads + 100.0*i_thr/cmd_options.num_threads);
         //}
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            cutset_info->arr_mcs[cutset_info->mcs_idx*mode_info->num_unit_size + u] = cutset_info->arr_pcs_thr[i_thr][i_pcs*mode_info->num_unit_size + u];
         }
         filled_in_cutsets++;

         // add new mcs to tree
         add_to_tree(0,cutset_info->mcs_idx, tree_info, cutset_info, mode_info, map_info, bitorder_info);

         cutset_info->mcs_idx++;
         check_and_resize_arr_mcs(cutset_info, mode_info);
         // CHECK_AND_RESIZE_ARR_MCS
      }
      cutset_info->pcs_idx_thr[i_thr] = 0;
   }

}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int find_superset_in_tree(long int node_idx, long int i_pcs, int thr_id, struct_tree_info *tree_info,
                          struct_cutset_info *cutset_info, struct_mode_info *mode_info)
{
   int u;

   int hit = 0;

   if( tree_info->nodes[node_idx].type == LEAF )
   {
      // g_num_cmp_events[thr_id]++;
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         if( tree_info->nodes[node_idx].cs[u] & ~cutset_info->arr_pcs[i_pcs*mode_info->num_unit_size + u] )
         {
            hit = 1;
            break;
         }
      }

      if( hit == 0 )
      {
         // g_num_filter_events[thr_id]++;
         // leaf was not a subset
         return(1);
      }
      else
      {
         return(0);
      }
   }
   else
   {
      // g_num_cmp_events[thr_id]++;
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         if( tree_info->nodes[node_idx].cs[u] & ~cutset_info->arr_pcs[i_pcs*mode_info->num_unit_size + u] )
         {
            hit = 1;
            break;
         }
      }

      if( hit == 0 )
      {
         // there might be a subset down there
         if( tree_info->nodes[node_idx].left != EMPTY  && tree_info->nodes[node_idx].right != EMPTY )
         {
            if( find_superset_in_tree( tree_info->nodes[node_idx].left, i_pcs, thr_id, tree_info, cutset_info, mode_info) != 0 )
            {
               // we found a subset -> go back -> no further search required
               return(1);
            }
            else
            {
               return find_superset_in_tree( tree_info->nodes[node_idx].right, i_pcs, thr_id, tree_info, cutset_info, mode_info);
            }
         }
         else if( tree_info->nodes[node_idx].left != EMPTY )
         {
            return find_superset_in_tree( tree_info->nodes[node_idx].left, i_pcs, thr_id, tree_info, cutset_info, mode_info);
         }
         else if( tree_info->nodes[node_idx].right != EMPTY )
         {
            return find_superset_in_tree( tree_info->nodes[node_idx].right, i_pcs, thr_id, tree_info, cutset_info, mode_info);
         }
         else
         {
            fprintf(stderr, "FATAL ERROR: we are in a non-leaf node and both kids (left and right) are empty\n");
            fprintf(stderr, "             node_idx=%ld\n",node_idx);
            exit(EXIT_FAILURE);
         }
      }
      else
      {
         // everything down this branch is not a subset
         // do not go down there! :-)
         return(0);
      }
   }

   // should never end up here
   return(0);
}
//////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void reorder_reaction_of_cutsets(struct_cutset_info *cutset_info, struct_map_info *map_info, struct_mode_info *mode_info, struct_bitorder_info *bitorder_info)
{
   long int e;
   int u;
   int r_in, r_out;
   int unit_cnt_in, unit_cnt_out;
   int bitmover_in, bitmover_out;
   unsigned long long int bitmask_in, bitmask_out;
   unsigned long long int *store;

   for( r_in = 0; r_in < map_info->num_dupset_keepers; r_in++ )
   {
      printf("reorder_reaction_of_cutsets(): r=%d bitorder[r].occ=%lu bitorder[r].bitpos=%d\n",r_in,bitorder_info->bitorder[r_in].occ,bitorder_info->bitorder[r_in].bitpos);
   }
   fflush(stdout);

   store = (unsigned long long int *) malloc(mode_info->num_unit_size*sizeof(unsigned long long int));

   if( store == NULL )
   {
      fprintf(stderr,"FATAL ERROR: reorder_reaction_of_modes(): couldn't allocate memory for tmp\n");
      exit(EXIT_FAILURE);
   }

   for( e = 0; e < cutset_info->mcs_idx; e++ )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         store[u] = cutset_info->arr_mcs[e*mode_info->num_unit_size + u];
         cutset_info->arr_mcs[e*mode_info->num_unit_size + u] = 0;
      }
      for( r_in = 0; r_in < map_info->num_dupset_keepers; r_in++ )
      {
         bitmask_in = 1;
         bitmover_in = r_in;
         unit_cnt_in = bitmover_in/(8*sizeof(unsigned long long));
         bitmover_in -= unit_cnt_in*8*sizeof(unsigned long long);
         bitmask_in <<= bitmover_in;
         r_out = bitorder_info->bitorder[r_in].bitpos;

         if( store[unit_cnt_in] & bitmask_in )
         {
            // bit was set -> move it
            bitmask_out = 1;
            bitmover_out = r_out;
            unit_cnt_out = bitmover_out/(8*sizeof(unsigned long long));
            bitmover_out -= unit_cnt_out*8*sizeof(unsigned long long);
            bitmask_out <<= bitmover_out;
            cutset_info->arr_mcs[e*mode_info->num_unit_size + unit_cnt_out] |= bitmask_out;
         }
      }
   }

   free(store);
}
/////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void uncompress_and_write_solution(int thread_id, int sol_num, struct_map_info *map_info, struct_reac_info *reac_info, struct_mode_info *mode_info, char *filename, int max_cancellations)
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
   // int num_cols;
   unsigned long long int t_max_sol_enhancer = g_max_sol_enhancer;

    // printf("DEBUG: entered uncompress_and_write_solution(): sol_num=%d\n",sol_num);

   // num_cols = map_info->num_dupset_keepers;

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

   int duplicate_start = new_l_r_cnt;

   for( i = 0; i < map_info->map_len; i++ )
   {
      new_cs[i] = 1.0;
      new_l_r[new_l_r_cnt*map_info->map_len + i] = 1.0;
   }

   for( i = 0; i < g_bruteforce_depth; i++ )
   {
      new_cs[map_info->dupsets_map[g_bruteforce_cutset[thread_id*max_cancellations + i]]] = 0.0;
      new_l_r[new_l_r_cnt*map_info->map_len + map_info->dupsets_map[g_bruteforce_cutset[thread_id*max_cancellations + i]]] = 0.0;
   }

   if( filename != NULL )
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
      // printf("d=%d duplicated_col=%d new_cs[duplicated_col]=%f\n",d,duplicated_col,new_cs[duplicated_col]);

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

               if( filename != NULL )
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

   // printf("INFO: uncompress_and_write_solution(): Total number of minimal cutsets after uncompressing: %lld\n",g_uncomp_solutions_cnt);
   free(new_l_r);
   free(new_cs);
   free(new_cs_dup);
}
/////////////////////////////////////////////////////////


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
