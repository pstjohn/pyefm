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
#include "prepro.h"
#include <stdio.h>


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void do_preprocessing(struct_map_info *map_info, struct_reac_info *reac_info, struct_mode_info *mode_info, struct_cmd_options cmd_options)
{
   initialize_map_remove_essentials(map_info, *reac_info, mode_info);
   if( PREPRO_REMOVE_ESSENTIAL_REACTIONS )
   {
      create_map_remove_essentials(mode_info, map_info, *reac_info);
   }
   restructure_flux_array_essentials(*reac_info, mode_info, *map_info);

   if( PREPRO_REMOVE_DUPLICATE_MODES )
   {
      remove_duplicate_modes(mode_info);
   }

   init_always_zero_variables(map_info, *mode_info, *reac_info);
   if( PREPRO_REMOVE_SUPERSET_MODES )
   {
      remove_superset_modes(mode_info, map_info, *reac_info, cmd_options);
   }

   init_duplicate_col_is_keeper_array(map_info, *reac_info);
   if( PREPRO_REMOVE_DUPLICATE_COLS )
   {
      find_duplicate_columns(map_info, *mode_info, *reac_info);
   }

   restructure_flux_array_duplicate_columns(mode_info, *reac_info, *map_info);
}
//////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void find_essential_reacs(struct_mode_info *mode_info,
                          struct_cmd_options cmd_options,
                          struct_reac_info reac_info)
{
   int i = 0;
   int essential_comp_val = 0;
   mode_info->num_essential_reacs = 0;

   if( cmd_options.good_efms_wanted > 0 )
   {
      essential_comp_val = cmd_options.good_efms - cmd_options.good_efms_wanted;
   }

   for( i = 0; i < reac_info.num_reactions; i++ )
   {
      if( mode_info->pnt_reac_occ[i] > essential_comp_val || mode_info->pnt_reac_occ_bad[i] == 0 )
      {
         mode_info->essential_idx[mode_info->num_essential_reacs] = i;
         printf("INFO: find_essential_reacs(): Identified essential reaction: \"%s\" index=%d\n",
                 reac_info.reactions[i],mode_info->essential_idx[mode_info->num_essential_reacs]);
         mode_info->num_essential_reacs++;
      }
   }
   printf("INFO: find_essential_reacs(): found %u essential reactions\n",mode_info->num_essential_reacs);

   return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void update_essential_reacs(struct_reac_info *reac_info, struct_mode_info *mode_info)
{
   int i,j;
   int found;
   int added = 0;

   printf("incorporating essential reactions read from file\n");

   for( i = 0; i < reac_info->num_readin_essential_reacs; i++ )
   {
      found = 0;
      for( j = 0; j < mode_info->num_essential_reacs; j++ )
      {
         if( reac_info->readin_essential_idx[i] == mode_info->essential_idx[j] )
         {
            found = 1;
            break;
         }
      }

      if( found == 0 )
      {
         // essential reaction specified in file
         // is not in set of implicit essential reactions
         // -> add it to total essential set
         mode_info->essential_idx[mode_info->num_essential_reacs+added] = reac_info->readin_essential_idx[i];
         added++;
         printf("INFO: essential reaction \"%s\" (index=%d) read from file is added to total set of essential reactions\n",
                reac_info->readin_essential_reactions[i],reac_info->readin_essential_idx[i]);
      }
      else
      {
         printf("INFO: essential reaction \"%s\" (index=%d) read from file is already an element of set of implicit essential reactions\n",
                reac_info->readin_essential_reactions[i],reac_info->readin_essential_idx[i]);
      }
   }

   // adjust length of essential reaction array
   mode_info->num_essential_reacs += added;

   // sort index array of essential reactions
   qsort(mode_info->essential_idx, mode_info->num_essential_reacs, sizeof(unsigned int), int_cmp);

   return;
}
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void initialize_map_remove_essentials(struct_map_info *map_info, struct_reac_info reac_info, struct_mode_info *mode_info)
{
   int i;

   map_info->map_len = reac_info.num_reactions;

   map_info->map = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->map == NULL )
   {
      fprintf(stderr, "FATAL ERROR: initialize_map_remove_essentials(): couldn't allocate %lu bytes for map\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   mode_info->is_essential = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( mode_info->is_essential == NULL )
   {
      fprintf(stderr, "FATAL ERROR: initialize_map_remove_essentials(): couldn't allocate %lu bytes for is_essential\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < map_info->map_len; i++ )
   {
      map_info->map[i] = i;
      mode_info->is_essential[i] = 0;
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_map_info(struct_map_info *map_info)
{
   free(map_info->map);
   if( map_info->always_zero_bitmap != NULL )
   {
      free(map_info->always_zero_bitmap);
      free(map_info->always_zero_map);
      free(map_info->is_always_zero);
      free(map_info->always_zero_mover);
   }

   free(map_info->num_reac_dupset);
   free(map_info->dupsets_keep_reac);
   free(map_info->dupsets_is_keeper);
   free(map_info->dupsets_map);
   free(map_info->dupsets);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void create_map_remove_essentials(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info reac_info)
{
   int i,j;
   int found;
   int found_cnt = 0;

   map_info->map_len = reac_info.num_reactions - mode_info->num_essential_reacs;

   printf("INFO: creating map ...\n");

   for( i = 0; i < reac_info.num_reactions; i++ )
   {
      found = 0;
      for( j = 0; j < mode_info->num_essential_reacs; j++ )
      {
         if( i == mode_info->essential_idx[j] )
         {
            found = 1;
            break;
         }
      }

      if( found == 0 )
      {
         map_info->map[found_cnt] = i;
         found_cnt++;
         mode_info->is_essential[i] = 0;
      }
      else
      {
         mode_info->is_essential[i] = 1;
      }
   }

   if( found_cnt != map_info->map_len )
   {
      fprintf(stderr, "FATAL ERROR: inconsistency occurred while creating map\n");
      exit(EXIT_FAILURE);
   }

   // print map
   for( i = 0; i < map_info->map_len; i++ )
   {
      printf("DEBUG: map %d -> %d (%s)\n",i,map_info->map[i],reac_info.reactions[map_info->map[i]]);
   }
   printf("DEBUG: length of map: %d\n",map_info->map_len);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// restructure_flux_array_essentials() creates a new
// array of modes that does not contain columns of
// essential reactions -> the size of the entire
// 2d array is, therefore, reduced.
//////////////////////////////////////////////////////
void restructure_flux_array_essentials(struct_reac_info reac_info, struct_mode_info *mode_info, struct_map_info map_info)
{
   int j,u,r;
   unsigned long long int i;
   unsigned long long int *new_emf;
   unsigned long long int tmp_unit;
   int bitmover;
   int unit_cnt;
   int mov_cnt;
   int *t_mover;

   t_mover = (int *) malloc((size_t)sizeof(int)*reac_info.num_reactions);
   if( t_mover == NULL )
   {
      fprintf(stderr, "FATAL ERROR: restructure_flux_array_essentials(): couldn't allocate memory for t_mover\n");
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   {
      printf("restructure_flux_array_essentials(): before restructuring: bad mode %06llu: ",i);
      print_mode(&(mode_info->pnt_flux_first[i*mode_info->num_unit_size_first]),mode_info->num_unit_size_first);
      printf("\n");
   }

   printf("INFO: restructuring and moving binary flux array\n");
   do_mem_allocation_bin_second(mode_info, map_info);

   printf("INFO: num_unit_size_first=%d num_unit_size=%d\n",mode_info->num_unit_size_first,mode_info->num_unit_size);

   new_emf = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int));
   if( new_emf == NULL )
   {
      fprintf(stderr, "FATAL ERROR: restructure_flux_array_essentials(): couldn't allocate memory for emf\n");
      exit(EXIT_FAILURE);
   }

   mov_cnt = 0;
   for( u = 0; u < reac_info.num_reactions; u++ )
   {
      if( mode_info->is_essential[u] != 0 )
      {
         mov_cnt++;
      }
      t_mover[u] = mov_cnt;
      printf("DEBUG: restructure_flux_array_essentials(): t_mover[%d]=%d mode_info->is_essential[u]=%d\n",u,t_mover[u],mode_info->is_essential[u]);
   }

   /////////////////////////////////////////////////
   // restructure flux array of bad modes
   /////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      // reset temp emf
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         new_emf[u] = 0;
      }

      for( j = 0; j < reac_info.num_reactions; j++ )
      {
         if( mode_info->is_essential[j] == 0 )
         {
            tmp_unit = 1;
            bitmover = j;
            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);
            tmp_unit <<= bitmover;
            if( unit_cnt >= mode_info->num_unit_size_first )
            {
               fprintf(stderr, "FATAL ERROR: unit_cnt(%d) is greater than num_unit_size_first (%d)\n",unit_cnt,mode_info->num_unit_size_first);
               fprintf(stderr, "             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
            if( (mode_info->pnt_flux_first[i*mode_info->num_unit_size_first + unit_cnt] & tmp_unit) != 0 )
            {
               // this bit is set
               r = j - t_mover[j];
               tmp_unit = 1;
               bitmover = r;
               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);
               tmp_unit <<= bitmover;
               if( unit_cnt >= mode_info->num_unit_size )
               {
                  fprintf(stderr, "FATAL ERROR: unit_cnt(%d) is greater than num_unit_size (%d)\n",unit_cnt,mode_info->num_unit_size);
                  fprintf(stderr, "             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               new_emf[unit_cnt] |= tmp_unit;
            }

         }
      }

      // set new emf value in global array
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         mode_info->pnt_flux[i*mode_info->num_unit_size + u] = new_emf[u];
      }

      if( (i+1)%100000 == 0 )
      {
         printf("DEBUG: restructured %llu of %llu remover modes: %llu%%\n",i,mode_info->bad_efms,100*i/mode_info->bad_efms);
      }
   }
   /////////////////////////////////////////////////

   /////////////////////////////////////////////////
   // restructure flux array of good modes
   /////////////////////////////////////////////////
   for( i = 0; i < mode_info->good_emfs; i++ )
   {
      // printf("restructuring mode %ld\n",i);
      // reset temp emf
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         new_emf[u] = 0;
      }

      for( j = 0; j < reac_info.num_reactions; j++ )
      {
         if( mode_info->is_essential[j] == 0 )
         {
            tmp_unit = 1;
            bitmover = j;
            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);
            tmp_unit <<= bitmover;
            if( unit_cnt >= mode_info->num_unit_size_first )
            {
               fprintf(stderr, "FATAL ERROR: unit_cnt(%d) is greater than mode_info->num_unit_size_first (%d)\n",unit_cnt,mode_info->num_unit_size_first);
               fprintf(stderr, "             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
            if( (mode_info->pnt_good_first[i*mode_info->num_unit_size_first + unit_cnt] & tmp_unit) != 0 )
            {
               // this bit is set
               r = j - t_mover[j];
               tmp_unit = 1;
               bitmover = r;
               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);
               tmp_unit <<= bitmover;
               if( unit_cnt >= mode_info->num_unit_size )
               {
                  fprintf(stderr, "FATAL ERROR: unit_cnt(%d) is greater than mode_info->num_unit_size (%d)\n",unit_cnt,mode_info->num_unit_size);
                  fprintf(stderr, "             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               new_emf[unit_cnt] |= tmp_unit;
            }

         }
      }

      // set new emf value in global array
      for( u = 0; u < mode_info->num_unit_size; u++ )
      {
         mode_info->pnt_good[i*mode_info->num_unit_size + u] = new_emf[u];
      }

      if( (i+1)%1000 == 0 )
      {
         printf("restructured %llu of %llu good modes: %llu%%\n",i,mode_info->good_emfs,100*i/mode_info->good_emfs);
      }
   }
   /////////////////////////////////////////////////

   free(new_emf);
   free(mode_info->pnt_flux_first);
   free(mode_info->pnt_good_first);

   for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   {
      printf("restructure_flux_array_essentials(): after restructuring: bad mode %06llu: ",i);
      print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size]),mode_info->num_unit_size);
      printf("\n");
   }

   for( i = 0; i < 50 && i < mode_info->good_emfs; i++ )
   {
      printf("restructure_flux_array_essentials(): after restructuring good mode %06llu: ",i);
      print_mode(&(mode_info->pnt_good[i*mode_info->num_unit_size]),mode_info->num_unit_size);
      printf("\n");
   }

   free(t_mover);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void do_mem_allocation_bin_second(struct_mode_info *mode_info, struct_map_info map_info)
{
   unsigned long i;
   unsigned long long t_zero = 0;

   ////////////////////////////////////////////////
   // determine unit size of 'reduced' system
   ////////////////////////////////////////////////
   mode_info->num_unit_size = ceil(((float)map_info.map_len)/(8*sizeof(unsigned long long)));
   printf("DEBUG: a single mode consists of %d 'unsigned long long' elements with a size of %lu bytes\n",mode_info->num_unit_size,sizeof(unsigned long long));
      printf("DEBUG: map_info.map_len=%u\n",map_info.map_len);
   printf("DEBUG: going to allocate %llu bytes for %llu 'unsigned long long' for flux values\n",mode_info->num_unit_size*mode_info->bad_efms*sizeof(unsigned long long),mode_info->num_unit_size*mode_info->bad_efms);
   ////////////////////////////////////////////////


   ////////////////////////////////////////////////
   // memory containing value of flux of bad modes
   ////////////////////////////////////////////////
   mode_info->pnt_flux = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*mode_info->bad_efms*sizeof(unsigned long long)));

   if( mode_info->pnt_flux == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for binary flux value (mode_info->pnt_flux)\n");
      exit(EXIT_FAILURE);
   }

   // initialize new memory
   for( i = 0; i < mode_info->num_unit_size*mode_info->bad_efms; i++ )
   {
      mode_info->pnt_flux[i] = t_zero;
   }
   ////////////////////////////////////////////////

   ////////////////////////////////////////////////
   // memory containing value of flux of good modes
   ////////////////////////////////////////////////
   mode_info->pnt_good = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*mode_info->good_emfs*sizeof(unsigned long long)));

   if( mode_info->pnt_good == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for binary flux value of keeper modes (mode_info->pnt_good)\n");
      exit(EXIT_FAILURE);
   }

   // initialize new memory
   for( i = 0; i < mode_info->num_unit_size*mode_info->good_emfs; i++ )
   {
      mode_info->pnt_good[i] = t_zero;
   }
   ////////////////////////////////////////////////

   mode_info->pnt_tmp_cutset = (unsigned long long*) malloc( (size_t) (mode_info->num_unit_size*sizeof(unsigned long long)));
   if( mode_info->pnt_tmp_cutset == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for temporary cutset\n");
      exit(EXIT_FAILURE);
   }

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void remove_duplicate_modes(struct_mode_info *mode_info)
{
   unsigned long long int i;
   unsigned long long max_good_bad_emfs;
   int u;

   struct_flux *t_pnt_struct_flux;

   for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   {
      printf("remove_duplicate_modes(): before restructuring: bad mode %06llu: ",i);
      print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size_first]),mode_info->num_unit_size_first);
      printf("\n");
   }

   ////////////////////////////////////////////////////////////////////////////
   // check if there is a mode with a norm of 0
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      int norm = 0;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         int b;
         unsigned long long tmp_unit = 0x8000000000000000ULL;
         for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
         {
            if( mode_info->pnt_flux[i*mode_info->num_unit_size + u] & tmp_unit )
            {
               norm++;
            }
            tmp_unit >>= 1;
         }
      }

      if( norm == 0 )
      {
         fprintf(stderr, "FATAL ERROR: remove_duplicate_modes(): norm (%d) of mode %llu to remove is 0!\n",norm,i);
         fprintf(stderr, "             this means that the node can never be knocked out!\n");
         fprintf(stderr, "             execution aborted.");
         print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size]), mode_info->num_unit_size);
         printf("\n");
         exit(EXIT_FAILURE);
      }
   }
   ////////////////////////////////////////////////////////////////////////////

   max_good_bad_emfs = mode_info->bad_efms;

   if( max_good_bad_emfs < mode_info->good_emfs )
   {
      max_good_bad_emfs = mode_info->good_emfs;
   }
   ////////////////////////////////////////////////////////////////////////////
   // allocate memory for temporary flux array
   ////////////////////////////////////////////////////////////////////////////
   t_pnt_struct_flux = (struct_flux*) malloc( (size_t) (max_good_bad_emfs*sizeof(struct_flux)));
   if( t_pnt_struct_flux == NULL )
   {
      fprintf(stderr, "FATAL ERROR: preprocess_efm_arr(): could not allocate memory for temporary flux array t_pnt_struct_flux\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < max_good_bad_emfs; i++ )
   {
      t_pnt_struct_flux[i].flux = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int));
      if( t_pnt_struct_flux[i].flux == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocating temporary flux memory (%llu) failed\n",i);
         exit(EXIT_FAILURE);
      }
   }
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // fill temporary flux array for sorting
   // note: the norm-parameter of structure is not used in this subroutine
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      t_pnt_struct_flux[i].unit_size = mode_info->num_unit_size;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         t_pnt_struct_flux[i].flux[u] = mode_info->pnt_flux[i*mode_info->num_unit_size + u];
      }
   }
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // sort modes
   ////////////////////////////////////////////////////////////////////////////
   printf("sorting mode array (remove_duplicate_modes) ...\n");
   // sort modes
   qsort(t_pnt_struct_flux, mode_info->bad_efms, sizeof(struct_flux), efm_cmp2);
   ////////////////////////////////////////////////////////////////////////////

   for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   {
      printf("remove_duplicate_modes(): after sorting: bad mode %06llu: ",i);
      print_mode(t_pnt_struct_flux[i].flux,mode_info->num_unit_size);
      printf("\n");
   }

   ////////////////////////////////////////////////////////////////////////////
   //  find duplicate array elements
   ////////////////////////////////////////////////////////////////////////////
   unsigned long long int *last_good;
   last_good = (unsigned long long int *) calloc( (size_t)mode_info->num_unit_size, (size_t)sizeof(unsigned long long int) );
   if( last_good == NULL )
   {
      fprintf(stderr, "FATAL ERROR: preprocess_efm_arr(): could not allocate memory for 'last_good'\n");
      exit(EXIT_FAILURE);
   }
   unsigned long int good_cnt;
   unsigned long long int duplicate_cnt = 0;
   for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
   {
      last_good[u] = t_pnt_struct_flux[0].flux[u];
   }
   good_cnt = 1;

   for( i = 1; i < mode_info->bad_efms; i++ )
   {
      int identical = 1;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         if( last_good[u] != t_pnt_struct_flux[i].flux[u] )
         {
            identical = 0;
            break;
         }
      }

      if( identical == 0 )
      {
         t_pnt_struct_flux[good_cnt].norm = t_pnt_struct_flux[i].norm;
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            t_pnt_struct_flux[good_cnt].flux[u] = t_pnt_struct_flux[i].flux[u];
            last_good[u] = t_pnt_struct_flux[i].flux[u];
         }
         good_cnt++;
      }
      else
      {
         duplicate_cnt++;
      }
   }
   printf("INFO: %llu duplicate modes (which do not contain essential reactions) have been removed\n",duplicate_cnt);
   printf("INFO: old mode_info->bad_efms=%llu\n",mode_info->bad_efms);
   mode_info->bad_efms = good_cnt;
   printf("INFO: new mode_info->bad_efms=%llu\n",mode_info->bad_efms);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // copy surviving modes back to global modes array
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         mode_info->pnt_flux[i*mode_info->num_unit_size + u] = t_pnt_struct_flux[i].flux[u];
      }
   }
   ////////////////////////////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////////////////////
   //  find duplicate array elements in good modes
   /////////////////////////////////////////////////////////////////////////////
   if( mode_info->good_emfs > 0 )
   {
       /////////////////////////////////////////////////////////////////////////
       // fill temporary flux array for sorting
       // note: the norm-parameter of structure is not used in this subroutine
       /////////////////////////////////////////////////////////////////////////
       for( i = 0; i < mode_info->good_emfs; i++ )
       {
          t_pnt_struct_flux[i].unit_size = mode_info->num_unit_size;
          for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
          {
             t_pnt_struct_flux[i].flux[u] = mode_info->pnt_good[i*mode_info->num_unit_size + u];
          }
       }
       /////////////////////////////////////////////////////////////////////////

      mode_info->pnt_good_duplicates = (unsigned long long int*) calloc(mode_info->good_emfs,sizeof(unsigned long long int));
      if( mode_info->pnt_good_duplicates == NULL )
      {
         printf("FATAL ERROR: couldn't allocate memory for mode_info->pnt_good_duplicates\n");
         exit(EXIT_FAILURE);
      }
      printf("sorting good mode array (1) ...\n");
      qsort(t_pnt_struct_flux, mode_info->good_emfs, sizeof(struct_flux), efm_cmp2);

      duplicate_cnt = 0;
      good_cnt = 0;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         last_good[u] = t_pnt_struct_flux[0].flux[u];
      }
      mode_info->pnt_good_duplicates[good_cnt] = 1;
      good_cnt++;

      for( i = 1; i < mode_info->good_emfs; i++ )
      {
         int identical = 1;
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            if( last_good[u] != t_pnt_struct_flux[i].flux[u] )
            {
               identical = 0;
               break;
            }
         }

         if( identical == 0 )
         {
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               t_pnt_struct_flux[good_cnt].flux[u] = t_pnt_struct_flux[i].flux[u];
               last_good[u] = t_pnt_struct_flux[i].flux[u];
            }
            mode_info->pnt_good_duplicates[good_cnt] = 1;
            good_cnt++;
         }
         else
         {
            mode_info->pnt_good_duplicates[good_cnt - 1]++;
            duplicate_cnt++;
         }
      }
      printf("INFO: %llu duplicate good modes have been removed\n",duplicate_cnt);
      printf("INFO: old mode_info->good_emfs=%llu ",mode_info->good_emfs);
      mode_info->good_emfs_orig = mode_info->good_emfs;
      mode_info->good_emfs = good_cnt;
      printf(", new mode_info->good_emfs=%llu\n",mode_info->good_emfs);

      ////////////////////////////////////////////////////////////////////////////
      // copy surviving modes back to global modes array
      ////////////////////////////////////////////////////////////////////////////
      for( i = 0; i < mode_info->good_emfs; i++ )
      {
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            mode_info->pnt_good[i*mode_info->num_unit_size + u] = t_pnt_struct_flux[i].flux[u];
         }
      }
   }
   /////////////////////////////////////////////////////////////////////////////

   for( i = 0; i < max_good_bad_emfs; i++ )
   {
      free(t_pnt_struct_flux[i].flux);
   }

   free(t_pnt_struct_flux);
   free(last_good);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_always_zero_variables(struct_map_info *map_info, struct_mode_info mode_info, struct_reac_info reac_info)
{
   int u;

   map_info->num_always_zero_modes = 0;
   map_info->num_always_zero_reacs = 0;

   ///////////////////////////////////////////////////
   // allocate memory
   ///////////////////////////////////////////////////
   map_info->always_zero_bitmap = (unsigned long long int *) malloc((size_t) mode_info.num_unit_size*sizeof(unsigned long long int));
   if( map_info->always_zero_bitmap == NULL )
   {
      fprintf(stderr, "FATAL ERROR: could not allocate memory for always zero bitmap\n");
      exit(EXIT_FAILURE);
   }

   map_info->always_zero_map = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->always_zero_map == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->always_zero_map\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->is_always_zero = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->is_always_zero == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->is_always_zero\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->always_zero_mover = (int *) malloc((size_t) sizeof(int)*reac_info.num_reactions);
   if( map_info->always_zero_mover == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->always_zero_mover\n",sizeof(int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }
   ///////////////////////////////////////////////////



   for( u = mode_info.num_unit_size - 1; u >= 0; u-- )
   {
      map_info->always_zero_bitmap[u] = 0;
   }
   for( u = 0; u < reac_info.num_reactions; u++ )
   {
      map_info->is_always_zero[u] = 0;
      map_info->always_zero_map[u] = 0;
      map_info->always_zero_mover[u] = 0;
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void remove_superset_modes(struct_mode_info *mode_info, struct_map_info *map_info, struct_reac_info reac_info, struct_cmd_options cmd_options)
{
   unsigned long long int i;
   unsigned long long int j;
   unsigned long long int t_new_bad_emfs = 0;
   unsigned long long int num_allocated = mode_info->bad_efms;
   int u;
   int b;
   int n;
   int r = 0;

   mode_info->max_norm_constraints = INIT_NORM_CONSTRAINTS;

   struct_flux *t_pnt_struct_flux;

   ////////////////////////////////////////////////////////////////////////////
   // allocate memory for temporary flux array
   ////////////////////////////////////////////////////////////////////////////
   t_pnt_struct_flux = (struct_flux*) malloc( (size_t) (num_allocated*sizeof(struct_flux)));
   if( t_pnt_struct_flux == NULL )
   {
      fprintf(stderr, "FATAL ERROR: preprocess_efm_arr(): could not allocate memory for temporary flux array t_pnt_struct_flux\n");
      fprintf(stderr, "             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      t_pnt_struct_flux[i].flux = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size*sizeof(unsigned long long int));
      if( t_pnt_struct_flux[i].flux == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocating temporary flux memory (%llu) failed\n",i);
         exit(EXIT_FAILURE);
      }
   }
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // fill temporary flux/norm array for sorting
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      int norm = 0;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         t_pnt_struct_flux[i].flux[u] = mode_info->pnt_flux[i*mode_info->num_unit_size + u];
         unsigned long long tmp_unit = 0x8000000000000000ULL;
         for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
         {
            if( mode_info->pnt_flux[i*mode_info->num_unit_size + u] & tmp_unit )
            {
               norm++;
            }
            tmp_unit >>= 1;
         }
      }
      t_pnt_struct_flux[i].norm = norm;
      t_pnt_struct_flux[i].unit_size = mode_info->num_unit_size;

      if( norm == 0 )
      {
         fprintf(stderr, "FATAL ERROR: norm (%d) of mode %llu to remove is 0!\n",norm,i);
         fprintf(stderr, "             this means that the node can never be knocked out!\n");
         fprintf(stderr, "             execution aborted.");
         print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size]), mode_info->num_unit_size);
         printf("\n");
         exit(EXIT_FAILURE);
      }
   }
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // sort temporary modes array
   ////////////////////////////////////////////////////////////////////////////
   printf("sorting mode array (remove_superset_modes) ...\n");
   // sort modes
   qsort(t_pnt_struct_flux, mode_info->bad_efms, sizeof(struct_flux), efm_cmp2);
   // sort modes by norm of mode
   printf("sorting mode array (remove_superset_modes) by norm of mode ...\n");
   qsort(t_pnt_struct_flux, mode_info->bad_efms, sizeof(struct_flux), efm_cmp_by_norm);
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // find all constraint that have a norm of 1
   // which means that there is only on reaction involved
   // this reactions must be knocked out in all cutsets
   // hence, the number of constraints can be reduced
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      // printf("norm=%d: ",t_pnt_struct_flux[i].norm);
      // print_bit_cutset(t_pnt_struct_flux[i].flux);
      if( t_pnt_struct_flux[i].norm == 1 )
      {
         map_info->num_always_zero_modes++;
         r = 0;
         for( u = 0; u < mode_info->num_unit_size; u++ )
         {
            unsigned long long tmp_unit = 0x0000000000000001ULL;
            for( b = 0; b < 8*sizeof(unsigned long long); b++ )
            {
               if( t_pnt_struct_flux[i].flux[u] & tmp_unit )
               {
                  // we found the single reaction of this mode
                  map_info->always_zero_bitmap[u] |= tmp_unit;
                  map_info->is_always_zero[map_info->map[r]] = 1;
                  map_info->always_zero_map[map_info->num_always_zero_reacs] = map_info->map[r];
                  map_info->num_always_zero_reacs++;
               }
               r++;
               tmp_unit <<= 1;
            }
         }
      }

      if( t_pnt_struct_flux[i].norm > mode_info->max_norm )
      {
         mode_info->max_norm = t_pnt_struct_flux[i].norm;
      }
   }

   printf("INFO: map_info->num_always_zero_modes=%llu map_info->num_always_zero_reacs=%d\n",map_info->num_always_zero_modes,map_info->num_always_zero_reacs);
   for( u = 0; u < map_info->num_always_zero_reacs; u++ )
   {
      printf("INFO: map_info->always_zero_map[%d]=%d\n",u,map_info->always_zero_map[u]);
   }
   for( u = 0; u < reac_info.num_reactions; u++ )
   {
      printf("INFO: map_info->is_always_zero[%d]=%d\n",u,map_info->is_always_zero[u]);
   }
   //printf("map_info->always_zero_bitmap: ");print_cutset(map_info->always_zero_bitmap);
   ////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////////
   // remove all modes that contain an always-zero reaction
   ////////////////////////////////////////////////////////////////////////////
   unsigned long long int t_removed_bad_emfs = 0;
   t_new_bad_emfs = 0;
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      if( t_pnt_struct_flux[i].norm > 1 )
      {
         int contains_zero_reaction = 0;
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            if( (map_info->always_zero_bitmap[u] != 0) && (t_pnt_struct_flux[i].flux[u] != 0) && (t_pnt_struct_flux[i].flux[u] & map_info->always_zero_bitmap[u]) )
            {
               contains_zero_reaction = 1;
               // printf("Mode %llu contains a zero reaction\n",i);
               // printf("   map_info->always_zero_bitmap: ");print_cutset(map_info->always_zero_bitmap);
               // printf("   flux:                 ");print_cutset(t_pnt_struct_flux[i].flux);
               break;
            }
         }
         if( contains_zero_reaction == 0 )
         {
             t_pnt_struct_flux[t_new_bad_emfs].norm = t_pnt_struct_flux[i].norm;
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               // t_pnt_struct_flux[t_new_bad_emfs*mode_info->num_unit_size + u] = t_pnt_struct_flux[i*mode_info->num_unit_size + u];
               t_pnt_struct_flux[t_new_bad_emfs].flux[u] = t_pnt_struct_flux[i].flux[u];
            }
            t_new_bad_emfs++;
         }
         else
         {
            t_removed_bad_emfs++;
         }

      }
      else
      {
         t_removed_bad_emfs++;
      }
   }
   printf("INFO: old bad_efms=%llu new bad_efms=%llu t_removed_bad_emfs=%llu after removing modes containing a zero-reaction\n",
           mode_info->bad_efms,t_new_bad_emfs,t_removed_bad_emfs);
   mode_info->bad_efms = t_new_bad_emfs;
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // remove all good modes that contain a zero reaction
   ////////////////////////////////////////////////////////////////////////////
   unsigned long long int t_removed_good_emfs = 0;
   unsigned long long int t_new_good_emfs = 0;
   for( i = 0; i < mode_info->good_emfs; i++ )
   {
      int contains_zero_reaction = 0;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         if( mode_info->pnt_good[i*mode_info->num_unit_size + u] & map_info->always_zero_bitmap[u] )
         {
            contains_zero_reaction = 1;
            break;
         }
      }
      if( contains_zero_reaction == 0 )
      {
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            mode_info->pnt_good[t_new_good_emfs*mode_info->num_unit_size + u] = mode_info->pnt_good[i*mode_info->num_unit_size + u];
         }
         t_new_good_emfs++;
      }
      else
      {
         t_removed_good_emfs++;
      }
   }
   printf("INFO: old good_emfs=%llu new good_emfs=%llu removed_good_emfs=%llu (because they contained an always-zero-reaction)\n",mode_info->good_emfs,t_new_good_emfs,t_removed_good_emfs);
   mode_info->good_emfs = t_new_good_emfs;
   mode_info->num_good_removed_by_always_zero = t_removed_good_emfs;
   printf("DEBUG: g_num_good_removed_by_always_zero=%llu\n",mode_info->num_good_removed_by_always_zero);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // remove all good modes that only contain zero elements
   ////////////////////////////////////////////////////////////////////////////
   t_removed_good_emfs = 0;
   t_new_good_emfs = 0;
   for( i = 0; i < mode_info->good_emfs; i++ )
   {
      int all_zero = 1;
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         if( mode_info->pnt_good[i*mode_info->num_unit_size + u] )
         {
            all_zero = 0;
            break;
         }
      }
      if( all_zero == 0 )
      {
         for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
         {
            mode_info->pnt_good[t_new_good_emfs*mode_info->num_unit_size + u] = mode_info->pnt_good[i*mode_info->num_unit_size + u];
         }
         t_new_good_emfs++;
      }
      else
      {
         t_removed_good_emfs++;
      }
   }
   printf("INFO: old good_emfs=%llu new good_emfs=%llu removed_good_emfs=%llu (because they contained an always-zero-reaction)\n",mode_info->good_emfs,t_new_good_emfs,t_removed_good_emfs);
   mode_info->good_emfs = t_new_good_emfs;
   mode_info->num_good_removed_because_all_zero = t_removed_good_emfs;
   printf("DEBUG: num_good_removed_because_all_zero=%llu\n",mode_info->num_good_removed_because_all_zero);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // find supersets for all sub-modes with norm > 1
   ////////////////////////////////////////////////////////////////////////////
   unsigned long int orig_bad_emfs = mode_info->bad_efms;
   mode_info->norm_bitmap = (unsigned long long int *) malloc((size_t) mode_info->max_norm_constraints*mode_info->num_unit_size*sizeof(unsigned long long int) );
   if( mode_info->norm_bitmap == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocating memory for mode_info->norm_bitmap failed\n");
      exit(EXIT_FAILURE);
   }
   for( n = 2; n < mode_info->max_norm; n++ )
   {
      mode_info->num_norm_constraints = 0;
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         if( t_pnt_struct_flux[i].norm == n )
         {
            if( mode_info->num_norm_constraints >= mode_info->max_norm_constraints )
            {
               // reallocate memory -> double the size of the array
               printf("INFO: reallocating memory for mode_info->norm_bitmap: mode_info->max_norm_constraints=%llu\n",mode_info->max_norm_constraints);
               mode_info->norm_bitmap = (unsigned long long int *) realloc(mode_info->norm_bitmap,(size_t) (mode_info->max_norm_constraints*2)*mode_info->num_unit_size*sizeof(unsigned long long int) );
               if( mode_info->norm_bitmap == NULL )
               {
                  fprintf(stderr, "FATAL ERROR: re-allocating memory for mode_info->norm_bitmap failed\n");
                  exit(EXIT_FAILURE);
               }
               mode_info->max_norm_constraints *= 2;
            }
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               mode_info->norm_bitmap[mode_info->num_norm_constraints*mode_info->num_unit_size + u] = 0;
            }
            for( u = 0; u < mode_info->num_unit_size; u++ )
            {
               unsigned long long tmp_unit = 0x0000000000000001ULL;
               for( b = 0; b < 8*sizeof(unsigned long long); b++ )
               {
                  if( t_pnt_struct_flux[i].flux[u] & tmp_unit )
                  {
                     mode_info->norm_bitmap[mode_info->num_norm_constraints*mode_info->num_unit_size + u] |= tmp_unit;
                  }
                  r++;
                  tmp_unit <<= 1;
               }
            }
            mode_info->num_norm_constraints++;
         }
         else if( t_pnt_struct_flux[i].norm > n )
         {
            break;
         }
      }

      // printf("INFO: mode_info->num_norm_constraints=%lu\n",mode_info->num_norm_constraints);

      t_new_bad_emfs = 0;
      for( i = 0; i < mode_info->bad_efms; i++ )
      {
         // if( i%10000 == 0 )
         // {
         //    printf("DEBUG: remove_superset_modes(): norm: %d/%d mode: %llu/%lu %03llu%%\n",n,mode_info->max_norm,i,mode_info->bad_efms,100*i/mode_info->bad_efms);
         // }

         int no_norm_constraint = 1;
         if( t_pnt_struct_flux[i].norm > n )
         {
            for( j = 0; j < mode_info->num_norm_constraints; j++ )
            {
               int found_const = 1;
               for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
               {
                  if( (mode_info->norm_bitmap[j*mode_info->num_unit_size + u] != 0) && ((t_pnt_struct_flux[i].flux[u] & mode_info->norm_bitmap[j*mode_info->num_unit_size + u]) !=
                       mode_info->norm_bitmap[j*mode_info->num_unit_size + u]) )
                  {
                     found_const = 0;
                     break;
                  }
               }
               if( found_const == 1 )
               {
                  // printf("found norm_constraint!: \n"); print_bit_cutset(mode_info->norm_bitmap[j]); print_bit_cutset(t_pnt_struct_flux[i].flux);
                  no_norm_constraint = 0;
                  break;
               }
            }
            // printf("i=%llu no_norm_constraint=%d\n",i,no_norm_constraint);
         }

         if( no_norm_constraint == 1 )
         {
            t_pnt_struct_flux[t_new_bad_emfs].norm = t_pnt_struct_flux[i].norm;
            for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
            {
               t_pnt_struct_flux[t_new_bad_emfs].flux[u] = t_pnt_struct_flux[i].flux[u];
            }
            t_new_bad_emfs++;
         }
         else
         {
            t_removed_bad_emfs++;
         }
      }
      printf("INFO: n=%d old bad_emfs=%llu new bad_emfs=%llu\n",n,mode_info->bad_efms,t_new_bad_emfs);
      mode_info->bad_efms = t_new_bad_emfs;

      if( cmd_options.apply_heuristics == 1 )
      {
         // do some heuristic to break out from subset-superset system reduction
         if( (((float) n)/((float) mode_info->max_norm)) > 0.2 )
         {
            // we have tested 20% of the available norms
            if( ((float)(orig_bad_emfs - mode_info->bad_efms))/((float) orig_bad_emfs) < 0.3  )
            {
               // apparently subset-superset system reduction is not really working
               printf("INFO: exiting subset-superset system reduction -> as not a lot of reduction is achieved\n");
               break;
            }
         }
      }
   }
   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   // copy remaining temporary modes back to global modes array
   ////////////////////////////////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      for( u = mode_info->num_unit_size - 1; u >= 0; u-- )
      {
         mode_info->pnt_flux[i*mode_info->num_unit_size + u] = t_pnt_struct_flux[i].flux[u];
      }
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // create mover-array
   ////////////////////////////////////////////////////////////////////////////
   int mov_cnt = 0;
   for( u = 0; u < reac_info.num_reactions; u++ )
   {
      if( map_info->is_always_zero[u] == 1 )
      {
         mov_cnt++;
      }
      map_info->always_zero_mover[u] = mov_cnt;
   }
   for( u = 0; u < reac_info.num_reactions; u++ )
   {
      printf("DEBUG: map_info->always_zero_mover[%d]=%d\n",u,map_info->always_zero_mover[u]);
   }
   ////////////////////////////////////////////////////////////////////////////

   // print map for debugging purposes
   for( i = 0; i < map_info->map_len; i++ )
   {
      printf("DEBUG: map %llu -> %d (%s)\n",i,map_info->map[i],reac_info.reactions[map_info->map[i]]);
   }
   printf("DEBUG: length of map: %d\n",map_info->map_len);

   // for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   // {
   //    printf("bad modes %06llu: ",i);
   //    print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size]),mode_info->num_unit_size);
   //    printf("\n");
   // }

   for( i = 0; i < num_allocated; i++ )
   {
      free(t_pnt_struct_flux[i].flux);
   }
   free(t_pnt_struct_flux);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void init_duplicate_col_is_keeper_array(struct_map_info *map_info, struct_reac_info reac_info)
{
   int m;

   ///////////////////////////////////////////////////
   // alocate memory
   ///////////////////////////////////////////////////
   map_info->dupsets = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions*reac_info.num_reactions);
   if( map_info->dupsets == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->dupsets\n",
              sizeof(unsigned int)*reac_info.num_reactions*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->num_reac_dupset = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->num_reac_dupset == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->num_reac_dupset\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->dupsets_keep_reac = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->dupsets_keep_reac == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->dupsets_keep_reac\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->dupsets_is_keeper = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->dupsets_is_keeper == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->dupsets_is_keeper\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }

   map_info->dupsets_map = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( map_info->dupsets_map == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for map_info->dupsets_map\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }
   ///////////////////////////////////////////////////

   for( m = 0; m < map_info->map_len; m++ )
   {
      map_info->dupsets_map[m] = m;
      map_info->dupsets_is_keeper[m] = 1;
   }

   map_info->num_dupset_keepers = map_info->map_len;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void find_duplicate_columns(struct_map_info *map_info, struct_mode_info mode_info, struct_reac_info reac_info)
{
   int t_unit_size;
   int r;
   struct_mode *l_mode_transpose = NULL;

   map_info->num_dupset_keepers = 0;

   // current size of mode set
   t_unit_size = ceil(((float)(mode_info.bad_efms+mode_info.good_emfs))/(8*sizeof(unsigned long long int)));
   map_info->unit_col_length = t_unit_size;

   l_mode_transpose = allocate_transpose_memory(t_unit_size, *map_info);

   create_tranpose_mode_array(&l_mode_transpose, t_unit_size, mode_info, *map_info);

   // printf("Transpose mode matrix before sorting:\n");
   // print_transpose_mode(l_mode_transpose,t_unit_size, mode_info, *map_info);

   qsort(l_mode_transpose, map_info->map_len, sizeof(struct_mode), efm_by_cols);

   // printf("Transpose mode matrix after sorting:\n");
   // print_transpose_mode(l_mode_transpose,t_unit_size, mode_info, *map_info);

   identify_duplicate_cols(l_mode_transpose, t_unit_size, map_info, reac_info);


   // free all alocated memories
   for( r = 0; r < map_info->map_len; r++ )
   {
      free(l_mode_transpose[r].flux);
   }
   free(l_mode_transpose);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
struct_mode* allocate_transpose_memory(int t_unit_size, struct_map_info map_info)
{
   unsigned long int i;
   int r;
   struct_mode *l_mode_transpose;

   // allocate new memory where we stored modes in tranposed way
   l_mode_transpose = (struct_mode*) malloc( (size_t) (sizeof(struct_mode)*map_info.map_len));
   if( l_mode_transpose == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for transpose mode array\n");
      exit(EXIT_FAILURE);
   }
   else
   {
      printf("DEBUG: allocate %lu bytes for tranpose memory array\n",t_unit_size*sizeof(unsigned long long)*map_info.map_len);
   }

   for( r = 0; r < map_info.map_len; r++ )
   {
      l_mode_transpose[r].flux = (unsigned long long*) malloc( (size_t) (t_unit_size*sizeof(unsigned long long)));

      for( i = 0; i < t_unit_size; i++ )
      {
         l_mode_transpose[r].flux[i] = 0;
      }
      l_mode_transpose[r].orig_reac = map_info.map_len - 1 - r;
      l_mode_transpose[r].unit_col_length = map_info.unit_col_length;
   }

   printf("DEBUG: leaving allocate_transpose_memory()\n");
   return(l_mode_transpose);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void create_tranpose_mode_array(struct_mode **l_mode_tranpose, int t_unit_size, struct_mode_info mode_info, struct_map_info map_info)
{
   unsigned long int m;
   unsigned int r;
   unsigned long long int tmp_unit_src;
   int unit_cnt_src = 0;
   int bitmover_src = 0;
   // unsigned long long int *tranpose;
   struct_mode *tranpose;

   unsigned long long int tmp_unit_des;
   int unit_cnt_des  = 0;
   int bitmover_des  = 0;

   int num_ones = 0;

   tranpose = *l_mode_tranpose;

   printf("DEBUG: entered create_tranpose_mode_array()\n");

   for( m = 0; m < mode_info.good_emfs; m++ )
   {
      tmp_unit_des = 1;
      bitmover_des = m;
      unit_cnt_des = m/(8*sizeof(unsigned long long));
      bitmover_des -= unit_cnt_des*8*sizeof(unsigned long long);
      tmp_unit_des <<= bitmover_des;

      for( r = 0; r < map_info.map_len; r++ )
      {
         tmp_unit_src  = 1;
         bitmover_src  = r;
         unit_cnt_src  = r/(8*sizeof(unsigned long long));
         bitmover_src  -= unit_cnt_src*8*sizeof(unsigned long long);
         tmp_unit_src  <<= bitmover_src;

         if (mode_info.pnt_good[m*mode_info.num_unit_size + unit_cnt_src] & tmp_unit_src )
         {
            // we found a position in mode matrix containing a '1'
            (tranpose[map_info.map_len - 1 - r]).flux[unit_cnt_des] |= tmp_unit_des;

            num_ones++;
         }
      }
   }

   for( m = 0; m < mode_info.bad_efms; m++ )
   {
      tmp_unit_des = 1;
      // bitmover_des = m;
      // unit_cnt_des = m/(8*sizeof(unsigned long long));
      bitmover_des = (mode_info.good_emfs + m);
      unit_cnt_des = (mode_info.good_emfs + m)/(8*sizeof(unsigned long long));
      bitmover_des -= unit_cnt_des*8*sizeof(unsigned long long);
      tmp_unit_des <<= bitmover_des;

      for( r = 0; r < map_info.map_len; r++ )
      {
         tmp_unit_src = 1;
         bitmover_src = r;
         unit_cnt_src = r/(8*sizeof(unsigned long long));
         bitmover_src -= unit_cnt_src*8*sizeof(unsigned long long);
         tmp_unit_src <<= bitmover_src;

         if (mode_info.pnt_flux[m*mode_info.num_unit_size + unit_cnt_src] & tmp_unit_src )
         {
            // we found a position in mode matrix containing a '1'
            (tranpose[map_info.map_len - 1 - r]).flux[unit_cnt_des] |= tmp_unit_des;

            num_ones++;
         }
      }
   }

   printf("DEBUG: create_tranpose_mode_array(): Set %d elements to '1'\n",num_ones);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// unsigned int g_num_dupsets = 0;
//////////////////////////////////////////////////////
void identify_duplicate_cols(struct_mode *l_mode_tranpose, int t_unit_size, struct_map_info *map_info, struct_reac_info reac_info)
{
   int r,m;
   int last_r = 0;
   int num_of_identical = 0;
   int identical = 1;
   int tmp_num = 0;

   // we definitely keep the first column of the sorted set
   // as it cannot be a duplicate of any other column
   map_info->dupsets_keep_reac[map_info->num_dupset_keepers++] = l_mode_tranpose[0].orig_reac;

   for( r = 1; r < map_info->map_len; r++ )
   {
      identical = 1;
      for( m = 0; m < t_unit_size; m++ )
      {
         if( l_mode_tranpose[r].flux[m] != l_mode_tranpose[last_r].flux[m] )
         {
            identical = 0;
            break;
         }
      }

      if( identical == 1 ){ tmp_num++; }
      printf("DEBUG: identify_duplicate_cols: tmp_nu%d r=%d last_r=%d identical=%d t_unit_size=%d\n",
             tmp_num,r,last_r,identical,t_unit_size);

      if( identical == 1 )
      {
         if( num_of_identical == 0 )
         {
            // we are dealing with a new set of identical columns
            map_info->dupsets[map_info->num_dupsets*reac_info.num_reactions + num_of_identical++] = l_mode_tranpose[last_r].orig_reac;
            map_info->dupsets[map_info->num_dupsets*reac_info.num_reactions + num_of_identical++] = l_mode_tranpose[r].orig_reac;
         }
         else
         {
            map_info->dupsets[map_info->num_dupsets*reac_info.num_reactions + num_of_identical++] = l_mode_tranpose[r].orig_reac;
         }
      }
      else
      {
         if( num_of_identical > 0 )
         {
            map_info->num_reac_dupset[map_info->num_dupsets] = num_of_identical;
            map_info->num_dupsets++;
         }
         num_of_identical = 0;
         last_r = r;
         map_info->dupsets_keep_reac[map_info->num_dupset_keepers++] = l_mode_tranpose[r].orig_reac;
      }
   }

   // check if the last two rows were identical
   if( identical == 1 )
   {
      map_info->num_reac_dupset[map_info->num_dupsets] = num_of_identical;
      map_info->num_dupsets++;
   }

   // create is_keeper-array using array of all keeper columns
   int keepers_cnt = 0;
   for( m = 0; m < map_info->map_len; m++ )
   {
      int is_keeper = 0;
      for( r = 0; r < map_info->num_dupset_keepers; r++ )
      {
         if( map_info->dupsets_keep_reac[r] == m )
         {
            is_keeper = 1;
            break;
         }
      }
      if( is_keeper == 1 )
      {
         map_info->dupsets_map[keepers_cnt++] = m;
         map_info->dupsets_is_keeper[m] = 1;
      }
      else
      {
         map_info->dupsets_is_keeper[m] = 0;
      }
   }

   // do some debugging/info output 
   printf("Results of duplicate idendification process:\n");
   for( m = 0; m < map_info->num_dupsets; m++ )
   {
      printf("duplicate set %d:",m);
      for( r = 0; r < map_info->num_reac_dupset[m]; r++ )
      {
         // printf(" %d",map_info->dupsets[m][r]);
         printf(" %d",map_info->dupsets[m*reac_info.num_reactions + r]);
      }
      printf(": %d columns\n",map_info->num_reac_dupset[m]);
   }

   printf("We are keeping %u columns\n",map_info->num_dupset_keepers);

   for( m = 0; m < map_info->num_dupset_keepers; m++ )
   {
      printf("keeper %d: column %u\n",m,map_info->dupsets_keep_reac[m]);
   }

   printf("Content of is_keeper array\n");
   for( m = 0; m < map_info->map_len; m++ )
   {
      printf("map_info->dupsets_is_keeper[%u]=%u\n",m,map_info->dupsets_is_keeper[m]);
   }

   printf("Duplicate columns map\n");
   for( m = 0; m < map_info->num_dupset_keepers; m++ )
   {
      printf("DEBUG: duplicate map %u -> %u -> %u (%s)\n",m,map_info->dupsets_map[m],map_info->map[map_info->dupsets_map[m]],
             reac_info.reactions[map_info->map[map_info->dupsets_map[m]]]);
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void restructure_flux_array_duplicate_columns(struct_mode_info *mode_info, struct_reac_info reac_info, struct_map_info map_info)
{
   int j,u,r;
   unsigned long long int i;
   unsigned long long int *new_emf;
   unsigned long long int tmp_unit;
   int bitmover;
   int unit_cnt;
   int mov_cnt;
   unsigned int t_num_unit_size_next;
   int *t_mover;

   t_mover = (int *) malloc((size_t)sizeof(int)*reac_info.num_reactions);
   if( t_mover == NULL )
   {
      fprintf(stderr, "FATAL ERROR: restructure_flux_array_essentials(): couldn't allocate memory for t_mover\n");
      exit(EXIT_FAILURE);
   }

   printf("INFO: restructuring binary flux array in order to remove duplicate columns\n");

   t_num_unit_size_next = ceil(((float)map_info.num_dupset_keepers)/(8*sizeof(unsigned long long)));
   printf("INFO: mode_info->num_unit_size=%d\n",mode_info->num_unit_size);
   printf("INFO: next unit_size: t_num_unit_size_next=%d\n",t_num_unit_size_next);

   new_emf = (unsigned long long int *) malloc( (size_t) t_num_unit_size_next*sizeof(unsigned long long int));

   mov_cnt = 0;
   for( u = 0; u < map_info.map_len; u++ )
   {
      if( map_info.dupsets_is_keeper[u] == 0 )
      {
         mov_cnt++;
      }
      t_mover[u] = mov_cnt;
      printf("DEBUG: restructure_flux_array_duplicate_columns(): t_mover[%d]=%d map_info.dupsets_is_keeper[u]=%d\n",u,t_mover[u],map_info.dupsets_is_keeper[u]);
   }

   /////////////////////////////////////////////////
   // restructure flux array of bad modes
   /////////////////////////////////////////////////
   for( i = 0; i < mode_info->bad_efms; i++ )
   {
      // reset temp emf
      for( u = 0; u < t_num_unit_size_next; u++ )
      {
         new_emf[u] = 0;
      }

      for( j = 0; j < map_info.map_len; j++ )
      {
         if( map_info.dupsets_is_keeper[j] == 1 )
         {
            // we want to keep this column
            tmp_unit = 1;
            bitmover = j;
            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);
            tmp_unit <<= bitmover;
            if( unit_cnt >= mode_info->num_unit_size )
            {
               fprintf(stderr, "FATAL ERROR: 01: unit_cnt(%d) is greater than mode_info->num_unit_size (%d)\n",unit_cnt,mode_info->num_unit_size);
               fprintf(stderr, "             execution aborted.\n");
               exit(EXIT_FAILURE);
            }

            if( (mode_info->pnt_flux[i*mode_info->num_unit_size + unit_cnt] & tmp_unit) != 0 )
            {
               // this bit is set -> deal with it
               r = j - t_mover[j];
               tmp_unit = 1;
               bitmover = r;
               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);
               tmp_unit <<= bitmover;

               if( unit_cnt >= t_num_unit_size_next )
               {
                  fprintf(stderr, "FATAL ERROR: 02: unit_cnt(%d) is greater than t_num_unit_size_next (%d)\n",unit_cnt,t_num_unit_size_next);
                  fprintf(stderr, "             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               new_emf[unit_cnt] |= tmp_unit;
            }

         }
      }

      // set new emf value in global array
      for( u = 0; u < t_num_unit_size_next; u++ )
      {
         mode_info->pnt_flux[i*t_num_unit_size_next + u] = new_emf[u];
      }

      if( (i+1)%100000 == 0 )
      {
         printf("restructured %llu of %llu remover modes: %llu%%\n",i,mode_info->bad_efms,100*i/mode_info->bad_efms);
      }
   }
   /////////////////////////////////////////////////

   /////////////////////////////////////////////////
   // restructure flux array of good modes
   /////////////////////////////////////////////////
   for( i = 0; i < mode_info->good_emfs; i++ )
   {
      // printf("restructuring mode %ld\n",i);
      // reset temp emf
      for( u = 0; u < t_num_unit_size_next; u++ )
      {
         new_emf[u] = 0;
      }

      for( j = 0; j < map_info.map_len; j++ )
      {
         if( map_info.dupsets_is_keeper[j] == 1 )
         {
            tmp_unit = 1;
            bitmover = j;
            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);
            tmp_unit <<= bitmover;
            if( unit_cnt >= mode_info->num_unit_size )
            {
               fprintf(stderr, "FATAL ERROR: 03: unit_cnt(%d) is greater than num_unit_size (%d)\n",unit_cnt,mode_info->num_unit_size);
               fprintf(stderr, "             execution aborted.\n");
               exit(EXIT_FAILURE);
            }
            if( (mode_info->pnt_good[i*mode_info->num_unit_size + unit_cnt] & tmp_unit) != 0 )
            {
               // this bit is set
               r = j - t_mover[j];
               tmp_unit = 1;
               bitmover = r;
               unit_cnt = bitmover/(8*sizeof(unsigned long long));
               bitmover -= unit_cnt*8*sizeof(unsigned long long);
               tmp_unit <<= bitmover;

               if( unit_cnt >= t_num_unit_size_next )
               {
                  fprintf(stderr, "FATAL ERROR: 04: unit_cnt(%d) is greater than t_num_unit_size_next (%d)\n",unit_cnt,t_num_unit_size_next);
                  fprintf(stderr, "             execution aborted.\n");
                  exit(EXIT_FAILURE);
               }
               new_emf[unit_cnt] |= tmp_unit;
            }
         }
      }

      // set new emf value in global array
      for( u = 0; u < t_num_unit_size_next; u++ )
      {
         mode_info->pnt_good[i*t_num_unit_size_next + u] = new_emf[u];
      }

      if( (i+1)%1000 == 0 )
      {
         printf("restructured %llu of %llu good modes: %llu%%\n",i,mode_info->good_emfs,100*i/mode_info->good_emfs);
      }
   }
   /////////////////////////////////////////////////
   free(new_emf);

   // resize array of bad modes and good modes
   if( mode_info->num_unit_size != t_num_unit_size_next )
   {
      if( mode_info->good_emfs > 0 )
      {
         mode_info->pnt_good = (unsigned long long*) realloc(mode_info->pnt_good, (size_t) (t_num_unit_size_next*mode_info->good_emfs*sizeof(unsigned long long)));
         if( mode_info->pnt_good == NULL )
         {
            fprintf(stderr, "FATAL ERROR: couldn't rellocate memory for binary flux value of keeper modes (mode_info->pnt_good)\n");
            exit(EXIT_FAILURE);
         }
      }
      mode_info->pnt_flux = (unsigned long long*) realloc(mode_info->pnt_flux, (size_t) (t_num_unit_size_next*mode_info->bad_efms *sizeof(unsigned long long)));
      if( mode_info->pnt_flux == NULL )
      {
         fprintf(stderr, "FATAL ERROR: couldn't rellocate memory for binary flux value of bad modes (mode_info->pnt_flux)\n");
         exit(EXIT_FAILURE);
      }
   }

   // save orig num_unit_size for unfolding
   mode_info->num_unit_size_orig = mode_info->num_unit_size;
   printf("DEBUG: mode_info->num_unit_size_orig=%u\n",mode_info->num_unit_size_orig);
   // set new global unit size
   mode_info->num_unit_size = t_num_unit_size_next;

   for( i = 0; i < 50 && i < mode_info->bad_efms; i++ )
   {
      printf("restructure_flux_array_duplicate_columns(): bad mode %06llu: ",i);
      print_mode(&(mode_info->pnt_flux[i*mode_info->num_unit_size]),mode_info->num_unit_size);
      printf("\n");
   }

   for( i = 0; i < 50 && i < mode_info->good_emfs; i++ )
   {
      printf("restructure_flux_array_duplicate_columns(): good mode %06llu: ",i);
      print_mode(&(mode_info->pnt_good[i*mode_info->num_unit_size]),mode_info->num_unit_size);
      printf("\n");
   }

   free(t_mover);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void uncompress_cutsets(struct_mode_info *mode_info, struct_cutset_info *cutset_info, struct_reac_info *reac_info, struct_map_info *map_info)
{
   long int c;
   long int orig_num_mcs;
   // unsigned long long int *bit_template;
   unsigned long long int *new_mcs;
   unsigned long long int *new_mcs_dup;
   unsigned long long int tmp_unit;
   int unit_cnt;
   int bitmover;
   unsigned long long int tmp_unit_dup;
   int unit_cnt_dup;
   int bitmover_dup;
   int u;
   long int i;
   long int d;
   long int x;
   unsigned long long int *new_arr_mcs;
   long int new_mcs_cnt = 0;
   long int max_new_mcs_idx = cutset_info->max_mcs_idx;

   printf("INFO: Uncompressing mincutsets ...\n");
   new_mcs      = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size_orig*sizeof(unsigned long long int));
   new_mcs_dup  = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size_orig*sizeof(unsigned long long int));

   // if( bit_template == NULL || new_mcs == NULL || new_mcs_dup == NULL )
   if( new_mcs == NULL || new_mcs_dup == NULL )
   {
      fprintf(stderr, "FATAL ERROR: in uncompress_cutsets(): couldn't allocate memory for bit_template or new_mcs or new_mcs_dup\n");
      exit(EXIT_FAILURE);
   }

   new_arr_mcs  = (unsigned long long int *) malloc( (size_t) mode_info->num_unit_size_orig*sizeof(unsigned long long int)*max_new_mcs_idx);

   if( new_arr_mcs == NULL )
   {
      fprintf(stderr, "FATAL ERROR: in uncompress_cutsets(): couldn't allocate memory for new_arr_mcs\n");
      exit(EXIT_FAILURE);
   }

   orig_num_mcs = cutset_info->mcs_idx;
   printf("number of cutsets before uncompressing: %llu\n",cutset_info->mcs_idx);
   for( c = 0; c < orig_num_mcs; c++ )
   {
      // reset temp mcs
      for( u = 0; u < mode_info->num_unit_size_orig; u++ )
      {
         new_mcs[u] = 0;
      }
      for( i = 0; i < map_info->num_dupset_keepers; i++ )
      {

         tmp_unit = 1;
         bitmover = i;
         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);
         tmp_unit <<= bitmover;

         if( (cutset_info->arr_mcs[c*mode_info->num_unit_size + unit_cnt] & tmp_unit) != 0 )
         {
            // this bit is set
            int r = map_info->dupsets_map[i];
            tmp_unit = 1;
            bitmover = r;
            unit_cnt = bitmover/(8*sizeof(unsigned long long));
            bitmover -= unit_cnt*8*sizeof(unsigned long long);
            tmp_unit <<= bitmover;

            new_mcs[unit_cnt] |= tmp_unit;
         }
      }
      // printf("Bit template:         ");print_cutset(bit_template);
      // printf("After uncompressing:  ");print_cutset(new_mcs);

      // here new_mcs mcs contains the uncompressed cutset
      // that is not yet duplicated

      long int duplicate_start = new_mcs_cnt;
      // set new mcs value in new mcs array
      for( u = 0; u < mode_info->num_unit_size_orig; u++ )
      {
         new_arr_mcs[new_mcs_cnt*mode_info->num_unit_size_orig + u] = new_mcs[u];
      }
      new_mcs_cnt++;

      if( new_mcs_cnt >= max_new_mcs_idx )
      {
         new_arr_mcs  = (unsigned long long int *) realloc(new_arr_mcs, (size_t) mode_info->num_unit_size_orig*sizeof(unsigned long long int)*(max_new_mcs_idx*2));
         if( new_arr_mcs == NULL )
         {
            fprintf(stderr, "FATAL ERROR: size of memory used to store uncompressed duplicated cutsets not large enough\n");
            fprintf(stderr, "             reallocating memory failed!\n");
            exit(EXIT_FAILURE);
         }
         max_new_mcs_idx *= 2;
      }

      long int duplicate_stop = new_mcs_cnt;

      // walk through all duplicate sets
      for( d = 0; d < map_info->num_dupsets; d++ )
      {
         tmp_unit = 1;
         // int r = map_info->dupsets_map[i];
         // bitmover = map_info->dupsets[d][0];
         bitmover = map_info->dupsets[d*reac_info->num_reactions + 0];
         unit_cnt = bitmover/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);
         tmp_unit <<= bitmover;

         // check if a duplicated column is active
         // in unduplicated cutset
         if( new_mcs[unit_cnt] & tmp_unit )
         {
            // we found a column that needs to be duplicated
            for( i = duplicate_start; i < duplicate_stop; i++ )
            {

               for( x = 1; x < map_info->num_reac_dupset[d]; x++ )
               {
                  // make a copy of the currently duplicated mincutset
                  for( u = 0; u < mode_info->num_unit_size_orig; u++ )
                  {
                     new_mcs_dup[u] = new_arr_mcs[i*mode_info->num_unit_size_orig + u];
                  }

                  if( (new_mcs_dup[unit_cnt] & tmp_unit) == 0 )
                  {
                     fprintf(stderr, "c=%ld d=%ld i=%ld x=%ld\n",c,d,i,x);
                     fprintf(stderr, "orig: %u dup: %u len: %d\n",map_info->dupsets[d*reac_info->num_reactions + 0],map_info->dupsets[d*reac_info->num_reactions + x],map_info->num_reac_dupset[d]);
                     fprintf(stderr, "bitmover=%d\n",bitmover);
                     fprintf(stderr, "unit_cnt=%d\n",unit_cnt);
                     fprintf(stderr, "new_mcs_cnt=%ld\n",new_mcs_cnt);
                     fprintf(stderr, "tmp_unit:    ");print_mode(&tmp_unit,1);
                     fprintf(stderr, "cutset_info->arr_mcs:   ");print_mode(&cutset_info->arr_mcs[c], mode_info->num_unit_size);
                     fprintf(stderr, "new_mcs:     ");print_mode(new_mcs, mode_info->num_unit_size_orig);
                     fprintf(stderr, "new_arr_mcs: ");print_mode(&new_arr_mcs[i*mode_info->num_unit_size_orig],mode_info->num_unit_size_orig);
                     fprintf(stderr, "new_mcs_dup: ");print_mode(new_mcs_dup, mode_info->num_unit_size_orig);
                     fflush(stderr);
                     exit(EXIT_FAILURE);
                  }
                  unsigned long long int inv_tmp_unit = ~tmp_unit;
                  // delete bit of duplicate column
                  new_mcs_dup[unit_cnt] &= inv_tmp_unit;

                  tmp_unit_dup = 1;
                  bitmover_dup = map_info->dupsets[d*reac_info->num_reactions + x];
                  unit_cnt_dup = bitmover_dup/(8*sizeof(unsigned long long));
                  bitmover_dup -= unit_cnt_dup*8*sizeof(unsigned long long);
                  tmp_unit_dup <<= bitmover_dup;
                  // set duplicated bit
                  new_mcs_dup[unit_cnt_dup] |= tmp_unit_dup;
                  // store new cutset in new-cutset array
                  for( u = 0; u < mode_info->num_unit_size_orig; u++ )
                  {
                     new_arr_mcs[new_mcs_cnt*mode_info->num_unit_size_orig + u] = new_mcs_dup[u];
                  }
                  new_mcs_cnt++;

                  if( new_mcs_cnt >= max_new_mcs_idx )
                  {
                     new_arr_mcs  = (unsigned long long int *) realloc(new_arr_mcs, (size_t) mode_info->num_unit_size_orig*sizeof(unsigned long long int)*(max_new_mcs_idx*2));
                     if( new_arr_mcs == NULL )
                     {
                        fprintf(stderr, "FATAL ERROR: size of memory used to store uncompressed duplicated cutsets not large enough\n");
                        fprintf(stderr, "             reallocating memory failed!\n");
                        exit(EXIT_FAILURE);
                     }
                     max_new_mcs_idx *= 2;
                  }
               }
            }
         }
         duplicate_stop = new_mcs_cnt;
      }
      // printf("c=%u: Total number of cutsets after uncompressing: %u\n",c,new_mcs_cnt);
   }

   printf("Total number of cutsets after uncompressing: %ld\n",new_mcs_cnt);
   free(cutset_info->arr_mcs);
   cutset_info->arr_mcs = (unsigned long long int *) malloc( (size_t) (new_mcs_cnt+1)*mode_info->num_unit_size_orig*sizeof(unsigned long long int));
   if( cutset_info->arr_mcs == NULL )
   {
      fprintf(stderr, "FATAL ERROR: unfolding cutsets: couldn't allocate cutset_info->arr_mcs\n");
      exit(EXIT_FAILURE);
   }
   cutset_info->max_mcs_idx = new_mcs_cnt + 1;

   // copy mincutsets from temp-array back to global array
   cutset_info->mcs_idx = 0;
   for( c = 0; c < new_mcs_cnt; c++ )
   {
      for( u = 0; u < mode_info->num_unit_size_orig; u++ )
      {
         cutset_info->arr_mcs[c*mode_info->num_unit_size_orig + u] = new_arr_mcs[c*mode_info->num_unit_size_orig + u];
      }
      cutset_info->mcs_idx++;

      if( cutset_info->mcs_idx >= cutset_info->max_mcs_idx )
      {
         cutset_info->arr_mcs = (unsigned long long int *) realloc(cutset_info->arr_mcs, (size_t) (cutset_info->max_mcs_idx*2)*mode_info->num_unit_size_orig*sizeof(unsigned long long int));
         if( cutset_info->arr_mcs == NULL )
         {
            fprintf(stderr, "FATAL ERROR: could not re-allocate memory for cutset_info->arr_mcs\n");
            exit(EXIT_FAILURE);
         }
         cutset_info->max_mcs_idx *= 2;
      }
   }

   // free(bit_template);
   free(new_mcs);
   free(new_mcs_dup);
   free(new_arr_mcs);
}
//////////////////////////////////////////////////////
