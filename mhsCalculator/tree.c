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
#include "tree.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void determine_tree_bit_order(struct_map_info map_info, struct_mode_info mode_info, struct_bitorder_info *bitorder_info)
{
   int r;
   unsigned long long int e;
   unsigned long long int tmp_unit;
   int bitmover;
   int unit_cnt;
   int *bitorder_mixer;

   bitorder_info->reac_occ = (unsigned long int *) calloc((size_t) map_info.num_dupset_keepers, (size_t) sizeof(unsigned long int));
   if( bitorder_info->reac_occ == NULL )
   {
      fprintf(stderr, "FATAL ERROR: determine_tree_bit_order(): couldn't allocate memory for bitorder_info->reac_occ\n");
      exit(EXIT_FAILURE);
   }

   bitorder_mixer = (int *) calloc((size_t) map_info.num_dupset_keepers, (size_t) sizeof(int));
   if( bitorder_mixer == NULL )
   {
      fprintf(stderr, "FATAL ERROR: determine_tree_bit_order(): couldn't allocate memory for bitorder_mixer\n");
      exit(EXIT_FAILURE);
   }
   for( r = 0; r < map_info.num_dupset_keepers; r++ )
   {
      bitorder_mixer[r] = r;
   }

   bitorder_info->bitorder = (struct_bitorder *) calloc((size_t) map_info.num_dupset_keepers, (size_t) sizeof(struct_bitorder));
   if( bitorder_info->bitorder == NULL )
   {
      fprintf(stderr, "FATAL ERROR: determine_tree_bit_order(): couldn't allocate memory for bitorder_info->bitorder\n");
      exit(EXIT_FAILURE);
   }

   // initialize bitorder array
   for( r = 0; r < map_info.num_dupset_keepers; r++ )
   {
      bitorder_info->bitorder[r].occ      = 0;
      bitorder_info->bitorder[r].bitpos   = bitorder_mixer[r];
      bitorder_info->bitorder[r].unit_cnt = bitorder_mixer[r]/(8*sizeof(unsigned long long int));
      bitorder_info->bitorder[r].bit_map  = 1;
      bitorder_info->bitorder[r].bit_map  <<= (bitorder_mixer[r] - bitorder_info->bitorder[r].unit_cnt*8*sizeof(unsigned long long int));
   }

   for( e = 0; e < mode_info.bad_efms; e++ )
   {
      for( r = 0; r < map_info.num_dupset_keepers; r++ )
      {
         tmp_unit = 1;
         bitmover = r;
         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);
         tmp_unit <<= bitmover;

         if( mode_info.pnt_flux[e*mode_info.num_unit_size + unit_cnt] & tmp_unit )
         {
            bitorder_info->bitorder[r].occ++;
         }
      }
   }

   qsort(bitorder_info->bitorder, map_info.num_dupset_keepers, sizeof(struct_bitorder), efm_cmp_bitorder);
   // qsort(bitorder_info->bitorder, map_info.num_dupset_keepers, sizeof(struct_bitorder), efm_cmp_bitorder_rev);

   for( e = 0; e < mode_info.bad_efms; e++ )
   {
      for( r = 0; r < map_info.num_dupset_keepers; r++ )
      {
         tmp_unit = 1;
         bitmover = r;
         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);
         tmp_unit <<= bitmover;

         if( mode_info.pnt_flux[e*mode_info.num_unit_size + unit_cnt] & tmp_unit )
         {
            bitorder_info->reac_occ[r]++;
         }
      }
   }

   // debug output
   for( r = 0; r < map_info.num_dupset_keepers; r++ )
   {
      printf("r=%d g_bitorder[r].occ=%lu g_bitorder[r].bitpos=%d\n",r,bitorder_info->bitorder[r].occ,bitorder_info->bitorder[r].bitpos);
   }
   for( r = 0; r < map_info.num_dupset_keepers; r++ )
   {
      printf("r=%d bitorder_info->reac_occ[r]=%lu\n",r,bitorder_info->reac_occ[r]);
   }

   free(bitorder_mixer);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_tree_info(struct_tree_info *tree_info)
{
   int i;

   for( i = 0; i < tree_info->node_idx; i++ )
   {
      free(tree_info->nodes[i].cs);
   }
   free(tree_info->nodes);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_bitorder_info(struct_bitorder_info *bitorder_info)
{
   free(bitorder_info->bitorder);
   free(bitorder_info->reac_occ);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void delete_from_tree(long int node_idx, long int i_mcs, struct_tree_info *tree_info, struct_mode_info *mode_info, struct_cutset_info *cutset_info)
{
   int u;
   long int parent = tree_info->nodes[node_idx].parent;
   long int left   = tree_info->nodes[node_idx].left;
   long int right  = tree_info->nodes[node_idx].right;
   long int parents_left;
   long int parents_right;
   long int last_idx;

   // printf("DEBUG: delete_from_tree(): entered: node_idx=%ld i_mcs=%ld left=%ld right=%ld parent=%ld\n",node_idx,i_mcs,left,right,parent);

   if( node_idx == 0 )
   {
      // we are not going to remove root node

      if( left == EMPTY && right == EMPTY )
      {
         // however we make sure that cs is set to '1'
         for( u = 0; u < mode_info->num_unit_size; u++ )
         {
            tree_info->nodes[node_idx].cs[u] = 0xFFFFFFFFFFFFFFFFULL;
         }
      }
      else if( left == EMPTY )
      {
         for( u = 0; u < mode_info->num_unit_size; u++ )
         {
            tree_info->nodes[node_idx].cs[u] = tree_info->nodes[right].cs[u];
         }

      }
      else if( right == EMPTY )
      {
         for( u = 0; u < mode_info->num_unit_size; u++ )
         {
            tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u];
         }

      }
      else
      {
         fprintf(stderr, "FATAL ERROR: delete_from_tree(): why did we end up here?\n");
         exit(EXIT_FAILURE);
      }

      // printf("DEBUG: delete_from_tree(): leaving (not going to delete root node)\n");
      return;
   }

   parents_left  = tree_info->nodes[parent].left;
   parents_right = tree_info->nodes[parent].right;

   last_idx = tree_info->node_idx - 1;

   if( tree_info->nodes[node_idx].type == LEAF )
   {
      // free memory for nodes cutset data
      free(tree_info->nodes[node_idx].cs);

      /////////////////////////////////////////////////////////////////////////
      // detach node
      /////////////////////////////////////////////////////////////////////////
      if( parents_left == node_idx )
      {
         // printf("leaf node was a left -> parents left is set to EMPTY\n");
         tree_info->nodes[parent].left = EMPTY;
      }
      else if( parents_right == node_idx )
      {
         // printf("leaf node was a right -> parents right is set to EMPTY\n");
         tree_info->nodes[parent].right = EMPTY;
      }
      else
      {
         fprintf(stderr, "FATAL ERROR: delete_from_tree(): 01: leaf node: inconsistent tree\n");
         fprintf(stderr, "             node_idx=%ld parent=%ld left=%ld right=%ld\n",node_idx,parent,left,right);
         fprintf(stderr, "             parents_left=%ld parents_right=%ld\n",parents_left,parents_right);
         exit(EXIT_FAILURE);
      }
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // free occupied space in node-array
      /////////////////////////////////////////////////////////////////////////
      if( node_idx != last_idx )
      {
         move_last_element(last_idx, node_idx, tree_info, cutset_info);
      }
      tree_info->node_idx--;
      /////////////////////////////////////////////////////////////////////////

      // check if we have moved parent of current node
      if( parent == last_idx )
      {
         delete_from_tree(node_idx,-1, tree_info, mode_info, cutset_info);
      }
      else
      {
         delete_from_tree(parent,-1, tree_info, mode_info, cutset_info);
      }
   }
   else // if( tree_info->nodes[node_idx].type == LEAF )
   {
      // this is a non-leaf
      if( left == EMPTY && right == EMPTY )
      {
         // both branches are empty -> we can delete
         // free memory for nodes cutset data
         free(tree_info->nodes[node_idx].cs);

         //////////////////////////////////////////////////////////////////////
         // detach node
         //////////////////////////////////////////////////////////////////////
         if( parents_left == node_idx )
         {
            tree_info->nodes[parent].left = EMPTY;
         }
         else if( parents_right == node_idx )
         {
            tree_info->nodes[parent].right = EMPTY;
         }
         else
         {
            fprintf(stderr, "FATAL ERROR: delete_from_tree(): 03: non-leaf node: inconsistent tree\n");
            fprintf(stderr, "             node_idx=%ld parent=%ld left=%ld right=%ld\n",node_idx,parent,left,right);
            fprintf(stderr, "             parents_left=%ld parents_right=%ld\n",parents_left,parents_right);
            exit(EXIT_FAILURE);
         }

         //////////////////////////////////////////////////////////////////////
         // free occupied space in node-array
         //////////////////////////////////////////////////////////////////////
         if( node_idx != last_idx )
         {
            move_last_element(last_idx, node_idx, tree_info, cutset_info);
         }
         tree_info->node_idx--;
         //////////////////////////////////////////////////////////////////////


         // check if we have moved parent of current node
         if( parent == last_idx )
         {
            delete_from_tree(node_idx,-1, tree_info, mode_info, cutset_info);
         }
         else
         {
            delete_from_tree(parent,-1, tree_info, mode_info, cutset_info);
         }
      }
      else
      {
         // one of the branches is not empty ->
         // reset cs-values
         update_cs_values(node_idx, tree_info, mode_info);
      }
   }
   // printf("DEBUG: delete_from_tree(): leaving\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void add_to_tree(long int cur_parent, long int i_mcs, struct_tree_info *tree_info, struct_cutset_info *cutset_info, struct_mode_info *mode_info,
                 struct_map_info map_info, struct_bitorder_info *bitorder_info)
{
   int i;
   long int left  = tree_info->nodes[cur_parent].left;
   long int right = tree_info->nodes[cur_parent].right;
   int new_bitpos;
   int new_unit_cnt;
   unsigned long long int new_bit_map;

   // printf("DEBUG: add_to_tree(): cur_parent=%ld i_mcs=%ld cutset_info->mcs_idx=%ld\n",cur_parent,i_mcs,cutset_info->mcs_idx);
   // printf("       left=%ld right=%ld unit_cnt=%d bitpos=%d\n",left,right,tree_info->nodes[cur_parent].unit_cnt,tree_info->nodes[cur_parent].bitpos);
   // printf("       tree_info->nodes[cur_parent].bit_map=     ");print_cutset(&tree_info->nodes[cur_parent].bit_map);
   // printf("       tree_info->nodes[cur_parent].cs=          ");print_cutset(tree_info->nodes[cur_parent].cs);
   // printf("       tree_info->nodes[cur_parent].cs=          ");print_cutset_reordered(tree_info->nodes[cur_parent].cs);
   // printf("       cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size]=");print_cutset(&cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size]);
   // printf("       cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size]=");print_cutset_reordered(&cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size]);
   // fflush(stdout);

   if( tree_info->nodes[cur_parent].bit_map & cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + tree_info->nodes[cur_parent].unit_cnt] )
   {
      // printf("       => set! => go to left branch\n");
      // bit is set -> we need to go to the left
      if( left == EMPTY )
      {
         // printf("       left branch is empty -> new idx: %ld\n",tree_info->node_idx);
         // there is no node -> add this mcs
         tree_info->nodes[tree_info->node_idx].parent   = cur_parent;
         tree_info->nodes[tree_info->node_idx].type     = LEAF;
         tree_info->nodes[tree_info->node_idx].left     = EMPTY;
         tree_info->nodes[tree_info->node_idx].right    = EMPTY;
         tree_info->nodes[tree_info->node_idx].i_mcs    = i_mcs;
         cutset_info->arr_mcs2node[i_mcs]        = tree_info->node_idx;

         tree_info->nodes[tree_info->node_idx].cs       = (unsigned long long int*) calloc((size_t) mode_info->num_unit_size, (size_t) sizeof(unsigned long long int));

         if( tree_info->nodes[tree_info->node_idx].cs == NULL )
         {
            fprintf(stderr, "FATAL ERROR: add_to_tree(): couldn't allocate memory for cutset (cs) of tree node %lu\n",tree_info->node_idx);
            exit(EXIT_FAILURE);
         }

         for( i = 0; i < mode_info->num_unit_size; i++ )
         {
            tree_info->nodes[tree_info->node_idx].cs[i] = cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + i];
         }

         // update current root
         tree_info->nodes[cur_parent].left  = tree_info->node_idx;

         tree_info->node_idx++;
         // check_and_resize_tree_info->node();
         CHECK_AND_RESIZE_G_NODE
      }
      else
      {
         // printf("       left branch is not empty!\n");
         if( tree_info->nodes[left].type == NON_LEAF )
         {
            // printf("       left node is a NON_LEAF!\n");
            // step down the tree
            add_to_tree(left, i_mcs, tree_info, cutset_info, mode_info, map_info, bitorder_info );
         }
         else
         {
            // printf("       left node is a LEAF!\n");

            // set xor_list and g_num_xor_elemes
            set_list_XOR_of_leaf_and_new(tree_info->nodes[left].cs, &cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size],
                                        cutset_info, mode_info, map_info, bitorder_info);

            new_bitpos   = cutset_info->xor_list[0].reac_id;
            new_unit_cnt = new_bitpos/(8*sizeof(unsigned long long int));
            new_bit_map  = 1;
            new_bit_map  <<= (new_bitpos - new_unit_cnt*8*sizeof(unsigned long long int));
            // printf("new_bitpos=%d new_unit_cnt=%d new_bit_map=%llu\n",new_bitpos,new_unit_cnt,new_bit_map);

            // create new non-leaf node
            int new_depth = tree_info->nodes[cur_parent].depth + 1;
            long int new_parent = tree_info->node_idx;
            tree_info->nodes[new_parent].i_mcs    = -1;
            tree_info->nodes[new_parent].parent   = cur_parent;
            tree_info->nodes[new_parent].left     = EMPTY;
            tree_info->nodes[new_parent].right    = EMPTY;
            tree_info->nodes[new_parent].type     = NON_LEAF;
            tree_info->nodes[new_parent].depth    = new_depth;
            tree_info->nodes[new_parent].bitpos   = new_bitpos;
            tree_info->nodes[new_parent].unit_cnt = new_unit_cnt;
            tree_info->nodes[new_parent].bit_map  = new_bit_map;

            tree_info->nodes[new_parent].cs     = (unsigned long long int*) calloc((size_t) mode_info->num_unit_size, (size_t) sizeof(unsigned long long int));

            if( tree_info->nodes[new_parent].cs == NULL )
            {
               fprintf(stderr, "FATAL ERROR: add_to_tree(): couldn't allocate memory for cutset (cs) of tree node %ld\n",new_parent);
               exit(EXIT_FAILURE);
            }

            // attach new non-leaf node to 'old' parent node of leaf node
            tree_info->nodes[cur_parent].left = tree_info->node_idx;

            // re-attach existing leaf node to new parent
            if( tree_info->nodes[new_parent].bit_map & tree_info->nodes[left].cs[tree_info->nodes[new_parent].unit_cnt] )
            {
               // bit is set -> needs to be attached on left side
               tree_info->nodes[new_parent].left = left;
            }
            else
            {
               // bit is set -> needs to be attached on right side
               tree_info->nodes[new_parent].right = left;
            }

            // set the new parent of the existing leaf node
            tree_info->nodes[left].parent = new_parent;

            tree_info->node_idx++;
            // check_and_resize_g_node();
            CHECK_AND_RESIZE_G_NODE

            add_to_tree(new_parent,i_mcs, tree_info, cutset_info, mode_info, map_info, bitorder_info);
         }
      }
   }
   else
   {
      // printf("       => is NOT set! => go to right branch\n");
      // bit is not set -> we need to go to the right
      if( right == EMPTY )
      {
         // printf("       right branch is empty -> new idx: %ld\n",tree_info->node_idx);
         // there is no node -> add this mcs
         tree_info->nodes[tree_info->node_idx].parent   = cur_parent;
         tree_info->nodes[tree_info->node_idx].type     = LEAF;
         tree_info->nodes[tree_info->node_idx].i_mcs    = i_mcs;
         tree_info->nodes[tree_info->node_idx].left     = EMPTY;
         tree_info->nodes[tree_info->node_idx].right    = EMPTY;
         cutset_info->arr_mcs2node[i_mcs]        = tree_info->node_idx;

         tree_info->nodes[tree_info->node_idx].cs     = (unsigned long long int*) calloc((size_t) mode_info->num_unit_size, (size_t) sizeof(unsigned long long int));

         if( tree_info->nodes[tree_info->node_idx].cs == NULL )
         {
            fprintf(stderr, "FATAL ERROR: add_to_tree(): couldn't allocate memory for cutset (cs) of tree node %lu\n",tree_info->node_idx);
            exit(EXIT_FAILURE);
         }

         for( i = 0; i < mode_info->num_unit_size; i++ )
         {
            tree_info->nodes[tree_info->node_idx].cs[i] = cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size + i];
         }

         // update current root
         tree_info->nodes[cur_parent].right  = tree_info->node_idx;

         tree_info->node_idx++;
         // check_and_resize_g_node();
         CHECK_AND_RESIZE_G_NODE
      }
      else
      {
         // printf("       right branch is not empty!\n");
         if( tree_info->nodes[right].type == NON_LEAF )
         {
            // printf("       right node is a NON_LEAF!\n");
            // step down the tree
            add_to_tree(right,i_mcs, tree_info, cutset_info, mode_info, map_info, bitorder_info);
         }
         else
         {
            // printf("       right node is a LEAF!\n");
            // set xor_list and g_num_xor_elemes
            set_list_XOR_of_leaf_and_new(tree_info->nodes[right].cs, &cutset_info->arr_mcs[i_mcs*mode_info->num_unit_size],
                                        cutset_info, mode_info, map_info, bitorder_info);

            new_bitpos   = cutset_info->xor_list[0].reac_id;
            new_unit_cnt = new_bitpos/(8*sizeof(unsigned long long int));
            new_bit_map  = 1;
            new_bit_map  <<= (new_bitpos - new_unit_cnt*8*sizeof(unsigned long long int));
            // printf("new_bitpos=%d new_unit_cnt=%d new_bit_map=%llu\n",new_bitpos,new_unit_cnt,new_bit_map);

            int new_depth = tree_info->nodes[cur_parent].depth + 1;
            long int new_parent = tree_info->node_idx;
            tree_info->nodes[new_parent].i_mcs    = -1;
            tree_info->nodes[new_parent].parent   = cur_parent;
            tree_info->nodes[new_parent].left     = EMPTY;
            tree_info->nodes[new_parent].right    = EMPTY;
            tree_info->nodes[new_parent].type     = NON_LEAF;
            tree_info->nodes[new_parent].depth    = new_depth;
            tree_info->nodes[new_parent].bitpos   = new_bitpos;
            tree_info->nodes[new_parent].unit_cnt = new_unit_cnt;
            tree_info->nodes[new_parent].bit_map  = new_bit_map;

            tree_info->nodes[new_parent].cs     = (unsigned long long int*) calloc((size_t) mode_info->num_unit_size, (size_t) sizeof(unsigned long long int));

            if( tree_info->nodes[new_parent].cs == NULL )
            {
               fprintf(stderr, "FATAL ERROR: add_to_tree(): couldn't allocate memory for cutset (cs) of tree node %ld\n",new_parent);
               exit(EXIT_FAILURE);
            }

            // attach new non-leaf node to 'old' parent node of leaf node
            tree_info->nodes[cur_parent].right = tree_info->node_idx;

            // re-attach existing leaf node to new parent
            if( tree_info->nodes[new_parent].bit_map & tree_info->nodes[right].cs[tree_info->nodes[new_parent].unit_cnt] )
            {
               // bit is set -> needs to be attached on left side
               tree_info->nodes[new_parent].left = right;
            }
            else
            {
               // bit is set -> needs to be attached on right side
               tree_info->nodes[new_parent].right = right;
            }

            // set the new parent of the existing leaf node
            tree_info->nodes[right].parent = new_parent;

            tree_info->node_idx++;
            // check_and_resize_g_node();
            CHECK_AND_RESIZE_G_NODE

            add_to_tree(new_parent,i_mcs, tree_info, cutset_info, mode_info, map_info, bitorder_info);
         }
      }
   }

   update_cutsets_in_tree(cur_parent, tree_info, mode_info);

   // printf("DEBUG: add_to_tree(): leaving\n");
   // fflush(stdout);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void update_cutsets_in_tree(long int node_idx, struct_tree_info *tree_info, struct_mode_info *mode_info)
{
   long int left;
   long int right;
   int u;

   // printf("DEBUG: entered update_cutsets_in_tree()\n");

   left   = tree_info->nodes[node_idx].left;
   right  = tree_info->nodes[node_idx].right;


   if( left != EMPTY && right != EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u] & tree_info->nodes[right].cs[u];
      }
   }
   else if( left != EMPTY && right == EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u];
      }
   }
   else if( left == EMPTY && right != EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[right].cs[u];
      }
   }
   else
   {
      fprintf(stderr, "FATAL ERROR: both branches of node (left and right) are EMPTY\n");
      fprintf(stderr, "             node number: node_idx=%lu\n",node_idx);
      exit(EXIT_FAILURE);
   }

   // if( parent != ROOT )
   // {
   //    update_cutsets_in_tree(parent);
   // }

   // printf("DEBUG: leaving update_cutsets_in_tree()\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void move_last_element(long int last_idx, long int new_idx, struct_tree_info *tree_info, struct_cutset_info *cutset_info)
{
   // copy last to new free spot
   tree_info->nodes[new_idx] = tree_info->nodes[last_idx];

   // set parent of left branch of moved node
   if( tree_info->nodes[new_idx].type == NON_LEAF)
   {
      if( tree_info->nodes[new_idx].left != EMPTY )
      {
         tree_info->nodes[tree_info->nodes[new_idx].left].parent = new_idx;
      }

      // set parent of left branch of moved node
      if( tree_info->nodes[new_idx].right != EMPTY )
      {
         tree_info->nodes[tree_info->nodes[new_idx].right].parent = new_idx;
      }
   }

   // set kid of moved node's parent
   long int moved_parent        = tree_info->nodes[new_idx].parent;
   long int moved_parents_left  = tree_info->nodes[moved_parent].left;
   long int moved_parents_right = tree_info->nodes[moved_parent].right;
   if( moved_parents_left == last_idx )
   {
      tree_info->nodes[moved_parent].left = new_idx;
   }
   else if( moved_parents_right == last_idx )
   {
      tree_info->nodes[moved_parent].right = new_idx;
   }
   else
   {
      fprintf(stderr, "FATAL ERROR: move_last_element(): inconsistent tree\n");
      fprintf(stderr, "             new_idx=%ld\n",new_idx);
      fprintf(stderr, "             last_idx=%ld\n",last_idx);
      fprintf(stderr, "             moved_parent=%ld moved_parents_left=%ld moved_parents_right=%ld\n",moved_parent,moved_parents_left,moved_parents_right);
      exit(EXIT_FAILURE);
   }

   if( tree_info->nodes[last_idx].type == LEAF )
   {
      // we are moving a leaf node ->
      // update cutset_info->arr_mcs2node
      long int l_i_mcs = tree_info->nodes[last_idx].i_mcs;
      cutset_info->arr_mcs2node[l_i_mcs] = new_idx;
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void update_cs_values(long int node_idx, struct_tree_info *tree_info, struct_mode_info *mode_info)
{
   int u;
   long int left   = tree_info->nodes[node_idx].left;
   long int right  = tree_info->nodes[node_idx].right;
   long int parent = tree_info->nodes[node_idx].parent;

   // printf("DEBUG: update_cs_values(): entered: node_idx=%ld left=%ld right=%ld parent=%ld\n",node_idx,left,right,parent);

   if( node_idx == 0 )
   {
      fprintf(stderr, "FATAL ERROR: update_cs_values(): this routine should never be invoked for root node\n");
      exit(EXIT_FAILURE);
   }

   if( left == EMPTY && right == EMPTY )
   {
      fprintf(stderr, "FATAL ERROR: update_cs_values(): invalid state: both kids are empty\n");
      exit(EXIT_FAILURE);
   }
   else if( left != EMPTY && right != EMPTY )
   {
      fprintf(stderr, "FATAL ERROR: update_cs_values(): invalid state: both kids are non-empty\n");
      exit(EXIT_FAILURE);
   }
   else if( left != EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u];
      }
   }
   else
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[right].cs[u];
      }
   }

   propagate_cs_value(parent, tree_info, mode_info);

   // printf("DEBUG: update_cs_values(): leaving\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void propagate_cs_value(long int node_idx, struct_tree_info *tree_info, struct_mode_info *mode_info)
{
   int u;
   long int left   = tree_info->nodes[node_idx].left;
   long int right  = tree_info->nodes[node_idx].right;
   long int parent = tree_info->nodes[node_idx].parent;

   // printf("DEBUG: propagate_cs_value(): entered: node_idx=%ld left=%ld right=%ld parent=%ld\n",node_idx,left,right,parent);

   if( left == EMPTY && right == EMPTY )
   {
      fprintf(stderr, "FATAL ERROR: propagate_cs_value(): invalid state: both kids are empty\n");
      exit(EXIT_FAILURE);
   }
   else if( left != EMPTY && right != EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u] & tree_info->nodes[right].cs[u];
      }
   }
   else if( left != EMPTY )
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[left].cs[u];
      }
   }
   else
   {
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         tree_info->nodes[node_idx].cs[u] = tree_info->nodes[right].cs[u];
      }
   }

   if( node_idx != 0 )
   {
      propagate_cs_value(parent, tree_info, mode_info);
   }

   // printf("DEBUG: propagate_cs_value(): leaving\n");
}
//////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void set_list_XOR_of_leaf_and_new(unsigned long long int *cs1, unsigned long long int *cs2, struct_cutset_info *cutset_info,
                                  struct_mode_info *mode_info, struct_map_info map_info, struct_bitorder_info *bitorder_info)
{
   cutset_info->num_xor_elems = 0;
   int u;
   int b;
   int r = 0;
   unsigned long long int xor;

   // printf("DEBUG: set_global_list_XOR_of_leaf_and_new(): entered\n");
   // printf("       r=%d\n",r);
   // printf("       cs1: "); print_cutset(cs1);
   // printf("       cs2: "); print_cutset(cs2);

   for( u = 0; u < mode_info->num_unit_size; u++ )
   {
      xor = cs1[u] ^ cs2[u];

      unsigned long long tmp_unit = 1;

      for( b = 0; b < 8*sizeof(unsigned long long); b++ )
      {
         if( xor & tmp_unit )
         {
            cutset_info->xor_list[cutset_info->num_xor_elems].reac_id = r;
            cutset_info->xor_list[cutset_info->num_xor_elems].occ     = bitorder_info->reac_occ[r];
            cutset_info->num_xor_elems++;
            // leave loop as soon as we found a bit position
            // where set1 und set2 are not equal
            // this is good enough, as reactions are ordered the proper way
            // printf("Found: u=%d, b=%d, r=%d, xor=%llu\n",u,b,r,xor);
            return;
         }
         // r--;
         r++;
         // tmp_unit >>= 1;
         tmp_unit <<= 1;
         if( r > map_info.num_dupset_keepers )
         {
            break;
         }
      }
   }

   if( cutset_info->num_xor_elems == 0 )
   {
      fprintf(stdout, "FATAL ERROR: set_global_list_XOR_of_leaf_and_new(): generated xor-list is empty!\n");
      printf("              set1="); print_mode(cs1,mode_info->num_unit_size); printf("\n");
      printf("              set2="); print_mode(cs2,mode_info->num_unit_size); printf("\n");
      exit(EXIT_FAILURE);
   }
}
/////////////////////////////////////////////////////////

