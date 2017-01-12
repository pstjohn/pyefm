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

#include "berge.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void do_mincutset_berge(struct_mode_info *mode_info, struct_reac_info reac_info, struct_map_info *map_info, struct_cmd_options cmd_options, struct_cutset_info *cutset_info)
{
   unsigned long long int e;
   struct timeval t_time_start;
   struct timeval t_time_stop;

   int t_array_changed;
   unsigned long long int *t_cutset;
   unsigned long long int *m_cutset;

   struct_tree_info tree_info = {NULL,    // nodes
                                 0,       // node_max_idx
                                 0,       // node_idx
                                 NULL };  // bitorder

   struct_bitorder_info bitorder_info = {NULL,  // struct_bitorder
                                         NULL}; // reac_occ

   init_cs_arr_berge(cutset_info, &tree_info, mode_info, cmd_options, *map_info, &bitorder_info, reac_info);

   // change order of reactions of modes
   reorder_reaction_of_modes(mode_info, *map_info, bitorder_info, cmd_options);

   t_cutset = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int) );
   m_cutset = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int) );
   if( t_cutset == NULL || m_cutset == NULL )
   {
      fprintf(stderr, "FATAL ERROR: do_mincutset_berge(): couldn't allocate memory for either t_cutset pr m_cutset\n");
      exit(EXIT_FAILURE);
   }


   gettimeofday(&t_time_start,NULL);
   for( e = 0; e < mode_info->bad_efms; e++ )
   {
      cutset_info->pcs_idx = 0;

      if( mode_info->num_unit_size == 1 )
      {
         t_array_changed = compute_precutsets_us1(e, t_cutset, m_cutset, mode_info, map_info, cutset_info, &tree_info, cmd_options);
      }
      else
      {
         t_array_changed = compute_precutsets(e, t_cutset, m_cutset, mode_info, map_info, cutset_info, &tree_info, cmd_options);
      }

      if( t_array_changed == 1 )
      {
         if( cmd_options.lin_superset_test == 1 )
         {
            do_superset_check_linear(e, cmd_options, mode_info, cutset_info);
         }
         else
         {
            do_superset_check_tree(e, cmd_options, mode_info, cutset_info, &tree_info);
         }
      } // if( t_array_changed != 0 )

      if( mode_info->num_unit_size == 1 )
      {
         fill_precutsets_to_cutsets_us1(&tree_info, cutset_info, mode_info, *map_info, &bitorder_info, cmd_options);
      }
      else
      {
         fill_precutsets_to_cutsets(&tree_info, cutset_info, mode_info, *map_info, &bitorder_info, cmd_options);
      }

      if( e%VERBOSITY_CNT == 0 )
      {
         gettimeofday(&t_time_stop,NULL);
         printf("e=%llu/%llu cutset_info->mcs_idx=%llu cutset_info->pcs_idx=%llu ",e,mode_info->bad_efms,cutset_info->mcs_idx,cutset_info->pcs_idx);
         display_execution_time(t_time_stop,t_time_start);
         fflush(stdout);
      }
   }

   printf("number of cutsets before expanding: %llu\n",cutset_info->mcs_idx);

   reorder_reaction_of_cutsets(cutset_info, map_info, mode_info, &bitorder_info);

   free(m_cutset);
   free(t_cutset);
   free_mem_tree_info(&tree_info);
   free_mem_bitorder_info(&bitorder_info);
   return;
}
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_cs_arr_berge(struct_cutset_info *cutset_info, struct_tree_info *tree_info, struct_mode_info *mode_info, struct_cmd_options cmd_options,
     struct_map_info map_info, struct_bitorder_info *bitorder_info, struct_reac_info reac_info)
{
   int j;
   int u;

   printf("DEBUG: init_cs_arr_berge(): entered\n");

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
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for cutset_info->xor_list\n",sizeof(struct_xor_list)*reac_info.num_reactions);
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
      fprintf(stderr, "FATAL ERROR: init_cs_arr_berge(): couldn't allocate memory for tree\n");
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
      fprintf(stderr, "FATAL ERROR: init_cs_arr_berge(): couldn't allocate memory for cutset (cs) of tree node %lu\n",tree_info->node_idx);
      exit(EXIT_FAILURE);
   }
   for( u = 0; u < mode_info->num_unit_size; u++ )
   {
      tree_info->nodes[tree_info->node_idx].cs[u] = 0xFFFFFFFFFFFFFFFFULL;
   }

   tree_info->node_idx++;

   CHECK_AND_RESIZE_G_NODE
   ////////////////////////////////////////////////////////////////////////////

   // there very first cutset is an 'empty' cutset
   for( j = 0; j < mode_info->num_unit_size; j++ )
   {
      cutset_info->arr_mcs[ mode_info->num_unit_size*cutset_info->mcs_idx + j] = 0;
   }
   // add first node to tree
   add_to_tree(0,cutset_info->mcs_idx, tree_info, cutset_info, mode_info, map_info, bitorder_info);

   cutset_info->mcs_idx++;
   ////////////////////////////////////////////////////////////////////////////

   for( j = 0; j < cmd_options.num_threads; j++ )
   {
      cutset_info->max_pcs_idx_thr[j] = INIT_PCS_ELEMS_THR;
      printf("DEBUG: init_cs_arr_berge(): cutset_info->max_pcs_idx_thr[%d]=%llu\n",j,cutset_info->max_pcs_idx_thr[j]);
      cutset_info->arr_pcs_thr[j] = (unsigned long long int *) malloc( (size_t) cutset_info->max_pcs_idx_thr[j]*mode_info->num_unit_size*sizeof(unsigned long long int) );
      if( cutset_info->arr_pcs_thr[j] == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocating memory for arr_pcs_thr[%d] failed\n",j);
         exit(EXIT_FAILURE);
      }
   }
   printf("DEBUG: init_cs_arr_berge(): leaving\n");
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
// us1 stands for unit size is equal to 1
//////////////////////////////////////////////////////
int compute_precutsets_us1(long long int e, unsigned long long int * t_cutset, unsigned long long int *m_cutset,
                           struct_mode_info *mode_info, struct_map_info *map_info, struct_cutset_info *cutset_info, 
                           struct_tree_info *tree_info, struct_cmd_options cmd_options)
{
   int b;
   int t_array_changed = 0;
   long int i_mcs;
   int num_ones;
   long unsigned int eliminated_mcs    = 0;
   long unsigned int new_potential_pcs = 0;

   t_cutset[0] = 0;
   m_cutset[0] = 0;


   for(i_mcs = 0; i_mcs < cutset_info->mcs_idx; i_mcs++ )
   {
      if( ((mode_info->pnt_flux[e] & cutset_info->arr_mcs[i_mcs]) == 0 ) )
      {
         t_array_changed = 1;
         m_cutset[0] = cutset_info->arr_mcs[i_mcs];

         eliminated_mcs++;

         if( cmd_options.lin_superset_test == 0 )
         {
            delete_from_tree(cutset_info->arr_mcs2node[i_mcs],i_mcs, tree_info, mode_info, cutset_info);
         }

         if( i_mcs != cutset_info->mcs_idx - 1 )
         {
            if( cmd_options.lin_superset_test == 0 )
            {
               cutset_info->arr_mcs2node[i_mcs] = cutset_info->arr_mcs2node[cutset_info->mcs_idx - 1];
               tree_info->nodes[cutset_info->arr_mcs2node[i_mcs]].i_mcs = i_mcs;
            }
            cutset_info->arr_mcs[i_mcs] = cutset_info->arr_mcs[cutset_info->mcs_idx - 1];
            i_mcs--;
            if( cutset_info->mcs_idx == 0 )
            {
               printf("FATAL ERROR: invalid number for cutset_info->mcs_idx: %llu\n",cutset_info->mcs_idx);
            }
         }
         cutset_info->mcs_idx--;

         num_ones = calc_num_ones(m_cutset, mode_info->num_unit_size);

         if( (cmd_options.max_cancellations == 0) || (num_ones < (cmd_options.max_cancellations - map_info->num_always_zero_reacs)) )
         {
            unsigned long long tmp_unit = 1;

            for( b = map_info->num_dupset_keepers - 1; b >= 0; b-- )
            {
               if( mode_info->pnt_flux[e] & tmp_unit )
               {
                  t_cutset[0] = m_cutset[0] | tmp_unit;

                  new_potential_pcs++;
                  cutset_info->arr_pcs[cutset_info->pcs_idx] = t_cutset[0];
                  cutset_info->pcs_idx++;
                  check_and_resize_arr_pcs(cutset_info, mode_info);
               }
               tmp_unit <<= 1;
            }
         }
      }
   }

   return t_array_changed;
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
void do_superset_check_linear(long long int e, struct_cmd_options cmd_options, struct_mode_info *mode_info, struct_cutset_info *cutset_info)
{
   int i;
   int iret[MAX_THREADS];
   int thread_id[MAX_THREADS];
   struct timeval t_time_start;
   struct timeval t_time_stop;
   struct_worker_relay worker_relay[MAX_THREADS];
   pthread_t thread[MAX_THREADS];

   // unsigned long long int num_filter_events = 0;
   // unsigned long long int num_cmp_events    = 0;


   gettimeofday(&t_time_start,NULL);

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      thread_id[i] = i;
      worker_relay[i].thread_id   = i;
      worker_relay[i].cmd_options = cmd_options;
      worker_relay[i].mode_info   = mode_info;
      worker_relay[i].cutset_info = cutset_info;
   }

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      // iret[i] = pthread_create(&thread[i], NULL, do_superset_check_linear_worker, (void *) &thread_id[i]);
      iret[i] = pthread_create(&thread[i], NULL, do_superset_check_linear_worker, (void *) &worker_relay[i]);
      // printf("Thread creation %i returns: %d\n",i,iret[i]);
   }

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      pthread_join( thread[i], NULL);
   }

   // for( i = 0; i < cmd_options.num_threads; i++ )
   // {
   //    num_filter_events += g_num_filter_events[i];
   //    num_cmp_events    += g_num_cmp_events[i];
   // }
   // g_num_cmp_events_total += num_cmp_events;

   gettimeofday(&t_time_stop,NULL);

   // if( e%VERBOSITY_CNT == 0 )
   // {
   //    printf("INFO: Superset test of candidates - linear:  "); display_execution_time(t_time_stop,t_time_start);
   //    printf("INFO: do_superset_check_linear(): number of comparison events: %llu/%llu\n",num_cmp_events,g_num_cmp_events_total);
   //    printf("INFO: do_superset_check_linear(): number of found subset: %llu/%lu (new subsets: %llu)\n",num_filter_events,g_pcs_idx,(g_pcs_idx-num_filter_events));
   // }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void *do_superset_check_linear_worker(void *ptr_worker_relay)
{
   struct_worker_relay *worker_relay = (struct_worker_relay *)ptr_worker_relay;

   int u;
   long int i_pcs, i_mcs;
   int is_subset;

   int thr_id                        = worker_relay->thread_id;
   struct_cmd_options  t_cmd_options = worker_relay->cmd_options;
   struct_mode_info   *t_mode_info   = worker_relay->mode_info;
   struct_cutset_info *t_cutset_info = worker_relay->cutset_info;

   t_cutset_info->pcs_idx_thr[thr_id] = 0;

   // g_num_cmp_events[thr_id] = 0;
   // g_num_filter_events[thr_id] = 0;

   if( t_mode_info->num_unit_size == 1 )
   {
      for( i_pcs = thr_id; i_pcs < t_cutset_info->pcs_idx; i_pcs+= t_cmd_options.num_threads )
      {
         is_subset = 0;
         for( i_mcs = 0; i_mcs < t_cutset_info->mcs_idx; i_mcs++ )
         {
            // g_num_cmp_events[thr_id]++;

            if( ((t_cutset_info->arr_mcs[i_mcs] & t_cutset_info->arr_pcs[i_pcs]) == t_cutset_info->arr_mcs[i_mcs] ) )
            {
               // g_num_filter_events[thr_id]++;
               is_subset = 1;
               break;
            }
         }

         if( is_subset == 0 )
         {
            // (pre-)cutset is ok -> check if we must test 
            // (pre-)cutset against the good modes
            if( t_cmd_options.good_efms_wanted > 0 && t_cmd_options.good_efms_wanted < t_cmd_options.good_efms )
            {
               // max number of wanted modes was specified
               // test if enough of the good modes surive
               if( enough_good_modes_survive(&t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size], t_cmd_options, t_mode_info) )
               {
                  // this is a precutset we keep
                  t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]] = t_cutset_info->arr_pcs[i_pcs];
                  t_cutset_info->pcs_idx_thr[thr_id]++;
                  if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
                  {
                     t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                              (size_t) t_mode_info->num_unit_size*(t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                     if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                     {
                        printf("FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                        printf("             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                        printf("             reallocating memory failed: %s\n",strerror(errno));
                     }
                     t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
                  }
               }
            }
            else
            {
               // this is a precutset we keep
               t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]] = t_cutset_info->arr_pcs[i_pcs];
               t_cutset_info->pcs_idx_thr[thr_id]++;
               if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
               {
                  t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                           (size_t) t_mode_info->num_unit_size*(t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                  if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                  {
                     fprintf(stderr, "FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                     fprintf(stderr, "             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                     fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
                     exit(EXIT_FAILURE);
                  }
                  t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
               }
            }
         }  // is subset != 1
      }
   }
   else  // t_mode_info->num_unit_size > 1
   {
      for( i_pcs = thr_id; i_pcs < t_cutset_info->pcs_idx; i_pcs+= t_cmd_options.num_threads )
      {
         is_subset = 0;
         for( i_mcs = 0; i_mcs < t_cutset_info->mcs_idx; i_mcs++ )
         {
            // printf("testing super set: i_mcs=%ld t_cutset_info->mcs_idx=%lu i_pcs=%ld, t_cutset_info->pcs_idx=%lu\n",i_mcs,t_cutset_info->mcs_idx,i_pcs,t_cutset_info->pcs_idx);
            // g_num_cmp_events[thr_id]++;
            int equ = 1;
            for( u = 0; u < t_mode_info->num_unit_size; u++ )
            {
               // if( (t_elem_mcs->cutset[u] != 0) && ((t_elem_mcs->cutset[u] & t_elem_pcs1->cutset[u]) != t_elem_mcs->cutset[u] ) )
               if( ((t_cutset_info->arr_mcs[i_mcs*t_mode_info->num_unit_size + u] & t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size + u]) != t_cutset_info->arr_mcs[i_mcs*t_mode_info->num_unit_size + u] ) )
               {
                  equ = 0;
                  break;
               }
            }

            if( equ == 1 )
            {
               // g_num_filter_events[thr_id]++;
               is_subset = 1;
               // t_elem_mcs = NULL;
               break;
            }
         }

         if( is_subset == 0 )
         {
            // (pre-)cutset is ok -> check if we must test 
            // (pre-)cutset against the good modes
            if( t_cmd_options.good_efms_wanted > 0 && t_cmd_options.good_efms_wanted < t_cmd_options.good_efms )
            {
               // max number of wanted modes was specified
               // test if enough of the good modes surive
               // if( enough_good_modes_survive(t_elem_pcs1->cutset) )
               if( enough_good_modes_survive(&t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size], t_cmd_options, t_mode_info) )
               {
                  for( u = 0; u < t_mode_info->num_unit_size; u++ )
                  {
                     t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]*t_mode_info->num_unit_size + u] = t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size + u];
                  }
                  t_cutset_info->pcs_idx_thr[thr_id]++;
                  if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
                  {
                     t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                              (size_t) t_mode_info->num_unit_size*(t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                     if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                     {
                        fprintf(stderr, "FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                        fprintf(stderr, "             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                        fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
                        exit(EXIT_FAILURE);
                     }
                     t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
                  }
               }
            }
            else
            {
               for( u = 0; u < t_mode_info->num_unit_size; u++ )
               {
                  t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]*t_mode_info->num_unit_size + u] = t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size + u];
               }
               t_cutset_info->pcs_idx_thr[thr_id]++;
               if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
               {
                  t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                           (size_t) t_mode_info->num_unit_size*(t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                  if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                  {
                     fprintf(stderr, "FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                     fprintf(stderr, "             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                     fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
                     exit(EXIT_FAILURE);
                  }
                  t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
               }
            }
         }  // is subset != 1
      }
   }

   return((void *)NULL);
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
      for( i_pcs = 0; i_pcs < cutset_info->pcs_idx_thr[i_thr]; i_pcs++ )
      {
         //if(i_pcs%10000 == 0 )
         //{
         //   printf("thr_id:%d/%lu %lu/%lu=%04.1f%% --> %04.1f%%\n",
         //           i_thr,g_num_threads,i_pcs,g_pcs_idx_thr[i_thr],100.0*i_pcs/g_pcs_idx_thr[i_thr],100.0*i_pcs/g_pcs_idx_thr[i_thr]/g_num_threads + 100.0*i_thr/g_num_threads);
         //}
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            cutset_info->arr_mcs[cutset_info->mcs_idx*mode_info->num_unit_size + u] = cutset_info->arr_pcs_thr[i_thr][i_pcs*mode_info->num_unit_size + u];
         }
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
void do_superset_check_tree(long long int e, struct_cmd_options cmd_options, struct_mode_info *mode_info, struct_cutset_info *cutset_info, struct_tree_info *tree_info)
{
   pthread_t thread[MAX_THREADS];
   int iret[MAX_THREADS];
   int thread_id[MAX_THREADS];
   struct_worker_relay_tree worker_relay[MAX_THREADS];
   int i;
   // struct timeval t_time_start;
   // struct timeval t_time_stop;

   // unsigned long long int num_filter_events = 0;
   // unsigned long long int num_cmp_events    = 0;


   // gettimeofday(&t_time_start,NULL);

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      thread_id[i] = i;
      worker_relay[i].thread_id   = i;
      worker_relay[i].cmd_options = cmd_options;
      worker_relay[i].mode_info   = mode_info;
      worker_relay[i].cutset_info = cutset_info;
      worker_relay[i].tree_info   = tree_info;
   }

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      // iret[i] = pthread_create(&thread[i], NULL, do_superset_check_linear_worker, (void *) &thread_id[i]);
      iret[i] = pthread_create(&thread[i], NULL, do_superset_check_tree_worker, (void *) &worker_relay[i]);
      // printf("Thread creation %i returns: %d\n",i,iret[i]);
   }

   for( i = 0; i < cmd_options.num_threads; i++ )
   {
      pthread_join( thread[i], NULL);
   }

   // for( i = 0; i < cmd_options.num_threads; i++ )
   // {
   //    num_filter_events += g_num_filter_events[i];
   //    num_cmp_events    += g_num_cmp_events[i];
   // }
   // g_num_cmp_events_total += num_cmp_events;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void *do_superset_check_tree_worker(void *ptr_worker_relay)
{
   struct_worker_relay_tree *worker_relay = (struct_worker_relay_tree *)ptr_worker_relay;

   int u;
   long int i_pcs;
   int is_subset;
   int add_pcs;

   int thr_id                        = worker_relay->thread_id;
   struct_cmd_options  t_cmd_options = worker_relay->cmd_options;
   struct_mode_info   *t_mode_info   = worker_relay->mode_info;
   struct_cutset_info *t_cutset_info = worker_relay->cutset_info;
   struct_tree_info   *t_tree_info   = worker_relay->tree_info;

   t_cutset_info->pcs_idx_thr[thr_id] = 0;

   // g_num_cmp_events[thr_id] = 0;
   // g_num_filter_events[thr_id] = 0;

   // printf("DEBUG: entered do_superset_check_tree_worker(): thr_id=%d\n",thr_id);

   if(  t_mode_info->num_unit_size == 1 )
   {
      for( i_pcs = thr_id; i_pcs < t_cutset_info->pcs_idx; i_pcs+= t_cmd_options.num_threads )
      {
         is_subset = 0;

         // subset search makes only sense
         // if there are already cutsets available
         if( t_cutset_info->mcs_idx > 0 )
         {
            if( find_superset_in_tree_us1(0,i_pcs,thr_id, t_tree_info, t_cutset_info) )
            {
               is_subset = 1;
            }
         }

         if( is_subset == 0 )
         {
            // (pre-)cutset is ok -> check if we must test 
            // (pre-)cutset against the good modes
            add_pcs = 1;
            if( t_cmd_options.good_efms_wanted > 0 && t_cmd_options.good_efms_wanted < t_cmd_options.good_efms )
            {
               // max number of wanted modes was specified
               // test if enough of the good modes surive
               // if( enough_good_modes_survive(t_elem_pcs1->cutset) )
               if( ! enough_good_modes_survive(&t_cutset_info->arr_pcs[i_pcs], t_cmd_options, t_mode_info) )
               {
                  add_pcs = 0;
               }
            }

            if( add_pcs )
            {
               t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]] = t_cutset_info->arr_pcs[i_pcs];
               t_cutset_info->pcs_idx_thr[thr_id]++;

               if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
               {
                  t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                           (size_t) (t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                  if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                  {
                     fprintf(stderr, "FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                     fprintf(stderr, "             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                     fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
                     exit(EXIT_FAILURE);
                  }
                  t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
               }
            }
         }  // is subset != 1
      } // loop over all precutsets
   }
   else
   {
      for( i_pcs = thr_id; i_pcs < t_cutset_info->pcs_idx; i_pcs+= t_cmd_options.num_threads )
      {
         is_subset = 0;

         // subset search makes only sense
         // if there are already cutsets available
         if( t_cutset_info->mcs_idx > 0 )
         {
            if( find_superset_in_tree(0,i_pcs,thr_id, t_tree_info, t_cutset_info, t_mode_info) )
            {
               is_subset = 1;
            }
         }


         if( is_subset == 0 )
         {
            // (pre-)cutset is ok -> check if we must test 
            // (pre-)cutset against the good modes
            add_pcs = 1;
            if( t_cmd_options.good_efms_wanted > 0 && t_cmd_options.good_efms_wanted < t_cmd_options.good_efms )
            {
               // max number of wanted modes was specified
               // test if enough of the good modes surive
               if( ! enough_good_modes_survive(&t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size], t_cmd_options, t_mode_info) )
               {
                  add_pcs = 0;
               }
            }

            if( add_pcs )
            {
               for( u = 0; u < t_mode_info->num_unit_size; u++ )
               {
                  t_cutset_info->arr_pcs_thr[thr_id][t_cutset_info->pcs_idx_thr[thr_id]*t_mode_info->num_unit_size + u] = t_cutset_info->arr_pcs[i_pcs*t_mode_info->num_unit_size + u];
               }
               t_cutset_info->pcs_idx_thr[thr_id]++;
               if( t_cutset_info->pcs_idx_thr[thr_id] >= t_cutset_info->max_pcs_idx_thr[thr_id] )
               {
                  t_cutset_info->arr_pcs_thr[thr_id] = (unsigned long long int *) realloc(t_cutset_info->arr_pcs_thr[thr_id],
                                           (size_t) t_mode_info->num_unit_size*(t_cutset_info->max_pcs_idx_thr[thr_id]*2)*sizeof(unsigned long long int) );
                  if( t_cutset_info->arr_pcs_thr[thr_id] == NULL )
                  {
                     fprintf(stderr, "FATAL ERROR: t_cutset_info->arr_pcs_thr[%d] is not big enough\n",thr_id);
                     fprintf(stderr, "             t_cutset_info->max_pcs_idx_thr[%d]=%llu\n",thr_id,t_cutset_info->max_pcs_idx_thr[thr_id]);
                     fprintf(stderr, "             reallocating memory failed: %s\n",strerror(errno));
                     exit(EXIT_FAILURE);
                  }
                  t_cutset_info->max_pcs_idx_thr[thr_id] *= 2;
               }
            }
         }  // is subset != 1
      } // loop over all precutsets
   } // is t_mode_info->num_unit_size == 1

#if PRINT_DEBUG_OUTPUT != 0
   gettimeofday(&t_stop_time_calc,NULL);
   printf("DEBUG: leaving do_superset_check_tree_worker(): thr_id=%d: ",thr_id); display_execution_time(t_stop_time_calc,t_start_time_calc);
#endif

   return((void *)NULL);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// us1 stands for unit size is equal to 1
//////////////////////////////////////////////////////
int find_superset_in_tree_us1(long int node_idx, long int i_pcs, int thr_id, struct_tree_info *tree_info, struct_cutset_info *cutset_info)
{
   if( tree_info->nodes[node_idx].type == LEAF )
   {
      // g_num_cmp_events[thr_id]++;

      if( tree_info->nodes[node_idx].cs[0] & ~cutset_info->arr_pcs[i_pcs])
      {
         return(0);
      }
      else
      {
         // g_num_filter_events[thr_id]++;
         // leaf was not a subset
         return(1);
      }
   }
   else
   {
      // g_num_cmp_events[thr_id]++;

      if( tree_info->nodes[node_idx].cs[0] & ~cutset_info->arr_pcs[i_pcs] )
      {
         // everything down this branch is not a subset
         // do not go down there! :-)
         return(0);
      }
      else
      {
         // there might be a subset down there
         if( tree_info->nodes[node_idx].left != EMPTY  && tree_info->nodes[node_idx].right != EMPTY )
         {
            if( find_superset_in_tree_us1( tree_info->nodes[node_idx].left, i_pcs, thr_id, tree_info, cutset_info) != 0 )
            {
               // we found a subset -> go back -> no further search required
               return(1);
            }
            else
            {
               return find_superset_in_tree_us1( tree_info->nodes[node_idx].right, i_pcs, thr_id, tree_info, cutset_info);
            }
         }
         else if( tree_info->nodes[node_idx].left != EMPTY )
         {
            return find_superset_in_tree_us1( tree_info->nodes[node_idx].left, i_pcs, thr_id, tree_info, cutset_info);
         }
         else if( tree_info->nodes[node_idx].right != EMPTY )
         {
            return find_superset_in_tree_us1( tree_info->nodes[node_idx].right, i_pcs, thr_id, tree_info, cutset_info);
         }
         else
         {
            fprintf(stderr, "FATAL ERROR: we are in a non-leaf node and both kids (left and right) are empty\n");
            fprintf(stderr, "             node_idx=%ld\n",node_idx);
            exit(EXIT_FAILURE);
         }
      }
   }

   // should never end up here
   return(0);
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
