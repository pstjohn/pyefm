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

#include "printing.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void display_execution_time(struct timeval t_stop,struct timeval t_start)
{
   unsigned long long int t_diff_time_usec = 0;
   unsigned long long int t_diff_time_msec = 0;
   unsigned long long int t_diff_time_sec  = 0;
   unsigned long long int t_diff_time_min  = 0;
   unsigned long long int t_diff_time_hour = 0;

   t_diff_time_usec = (t_stop.tv_sec - t_start.tv_sec)*1000000 + (t_stop.tv_usec - t_start.tv_usec);

   t_diff_time_hour  = t_diff_time_usec/1000000/3600;
   t_diff_time_min   = ((t_diff_time_usec/1000000)%3600)/60;
   t_diff_time_sec   = (t_diff_time_usec/1000000)%60;
   t_diff_time_msec  = (t_diff_time_usec/1000)%1000;
   
   
   
   printf(" %012llu usec = %04lluh %02llum %02llus %03llums\n",
          t_diff_time_usec, t_diff_time_hour, t_diff_time_min, t_diff_time_sec, t_diff_time_msec);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_reactions(struct_reac_info reac_info)
{
   int nr;

   printf("Reactions: ");
   for( nr = 0; nr < reac_info.num_reactions; nr++ )
   {
      printf("%s ",reac_info.reactions[nr]);
   }
   printf("\n");

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_reac_occ_arr(struct_mode_info mode_info, struct_reac_info reac_info, struct_cmd_options cmd_options)
{
   int i = 0;
   printf("Reaction occurrence: good_emfs=%llu good_emfs_wanted=%llu\n",cmd_options.good_efms,cmd_options.good_efms_wanted);
   for( i = 0; i < reac_info.num_reactions; i++ )
   {
      printf("reaction index=%d, reaction name=\"%s\", occurrence_in_good=%llu occurrence_in_bad=%llu of %llu\n",
              i,reac_info.reactions[i],mode_info.pnt_reac_occ[i],mode_info.pnt_reac_occ_bad[i],mode_info.num_efms);
   }

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_essential_reacs(struct_mode_info mode_info, struct_reac_info reac_info)
{
   int i = 0;
   unsigned int idx;
   printf("INFO: print_essential_reacs(): Essential reactions:\n"); 
   for( i = 0; i < mode_info.num_essential_reacs; i++ )
   {
      idx = mode_info.essential_idx[i];
      printf("INFO: print_essential_reacs():reaction index=%d, reaction name=\"%s\", occurrence_in_good=%llu occurrence_in_bad=%llu of %llu\n",
             idx,reac_info.reactions[idx],mode_info.pnt_reac_occ[idx],mode_info.pnt_reac_occ_bad[idx],mode_info.num_efms);
   } 
   printf("INFO: print_essential_reacs(): In total %d essential reactions found\n",mode_info.num_essential_reacs);
   
   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_final_essential_reacs(struct_mode_info mode_info, struct_reac_info reac_info)
{
   int i = 0;
   unsigned int idx;
   printf("INFO: print_final_essential_reacs(): Final set of essential reactions\n");
   for( i = 0; i < mode_info.num_essential_reacs; i++ )
   {
      idx = mode_info.essential_idx[i];
      printf("INFO: print_final_essential_reacs(): reaction index=%d, reaction name=\"%s\"\n",idx,reac_info.reactions[idx]);
   }
   printf("INFO: print_final_essential_reacs(): In total %d essential reactions found\n",mode_info.num_essential_reacs);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_mode(unsigned long long int *mode, unsigned int unit_size)
{
   int u,b;

   for( u = unit_size - 1; u >= 0; u-- )
   {
      unsigned long long tmp_unit = 0x8000000000000000ULL;
      for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
      {
         if( mode[u] & tmp_unit )
         {
            printf("1");
         }
         else
         {
            printf("0");
         }
         tmp_unit >>= 1;
         if( b%8 == 0 )
         {
            printf("|");
         }
      }
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_transpose_mode(struct_mode *l_mode_tranpose, int t_unit_size, struct_mode_info mode_info,struct_map_info map_info)
{
   unsigned long int m;
   unsigned int r;
   unsigned long long int tmp_unit_des;
   int unit_cnt_des = 0;
   int bitmover_des = 0;


   for( m = 0; m < mode_info.bad_efms + mode_info.good_emfs; m++ )
   {
      tmp_unit_des = 1;
      bitmover_des = m;
      unit_cnt_des = m/(8*sizeof(unsigned long long));
      bitmover_des -= unit_cnt_des*8*sizeof(unsigned long long);
      tmp_unit_des <<= bitmover_des;

      for( r = 0; r < map_info.map_len; r++ )
      {
         if( ((l_mode_tranpose[r]).flux[unit_cnt_des]) & tmp_unit_des )
         {
            printf("1");
         }
         else
         {
            printf("0");
         }
      }
      printf("\n");

   }
   fflush(stdout);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_all_cutsets_text(struct_cmd_options *cmd_options, struct_cutset_info *cutset_info, struct_mode_info *mode_info,
                           struct_reac_info *reac_info, struct_map_info *map_info)
{
   unsigned long long int i;
   FILE *fh;
   printf("DEBUG: printing cutsets to file '%s'",cmd_options->o_filename);

   fh = fopen(cmd_options->o_filename, "w");

   if( fh <= 0 )
   {
      fprintf(stderr, "FATAL ERROR: open file '%s' for writing failed: %s\n",cmd_options->o_filename,strerror(errno));
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < cutset_info->mcs_idx; i++ )
   {
      print_cutset_text(fh, &cutset_info->arr_mcs[i*mode_info->num_unit_size_orig], mode_info, reac_info, map_info, cmd_options->output_bitvector);
   }
   printf("\n");
   printf("Total number of cutsets: %llu\n",cutset_info->mcs_idx);

   fclose(fh);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_cutset_text(FILE *fh, unsigned long long *cs, struct_mode_info *mode_info, struct_reac_info *reac_info,
                       struct_map_info *map_info, int use_bitvector)
{
   int u,b,i,m;
   int v;
   int r = 0;
   int num_knockouts = 0;

   if( use_bitvector == 0 )
   {
      for( u = mode_info->num_unit_size_orig - 1; u >= 0; u-- )
      {
         unsigned long long tmp_unit = 0x8000000000000000ULL;
         for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
         {
            if( cs[u] & tmp_unit )
            {
               fprintf(fh,"%s ",reac_info->reactions[map_info->map[mode_info->num_unit_size_orig*8*sizeof(unsigned long long) - 1 - r]]);
               num_knockouts++;
            }
            r++;
            tmp_unit >>= 1;
         }
      }

      for( i = 0; i < map_info->num_always_zero_reacs; i++ )
      {
         fprintf(fh,"%s ",reac_info->reactions[map_info->always_zero_map[i]]);
         num_knockouts++;
      }
   }
   else
   {
      // write minimal cut set as bitvector string to file
      // WARNING: this is a very slow implementation
      //          of printing the minimal cutsets as bitvectors!

      for( m = 0; m < reac_info->num_reactions; m++ )
      {
         v = 0;
         r = 0;
         for( u = mode_info->num_unit_size_orig - 1; u >= 0; u-- )
         {
            unsigned long long tmp_unit = 0x8000000000000000ULL;
            for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
            {
               if( cs[u] & tmp_unit )
               {
                  if( map_info->map[mode_info->num_unit_size_orig*8*sizeof(unsigned long long) - 1 - r] == m )
                  {
                     // fprintf(fh,"%s-%d ",reac_info->reactions[map_info->map[mode_info->num_unit_size_orig*8*sizeof(unsigned long long) - 1 - r]], r);
                     v = 1;
                     break;
                  }
               }
               r++;
               tmp_unit >>= 1;
            }
            if( v == 1 ) break;
         }

         if( v == 0 && map_info->num_always_zero_reacs > 0 )
         {
            for( i = 0; i < map_info->num_always_zero_reacs; i++ )
            {
               if( map_info->always_zero_map[i] == m )
               {
                  v = 1;
                  break;
               }
            }
         }

         fprintf(fh,"%1d",v);
      }
   }

   fprintf(fh,"\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void print_min_max_knockouts(struct_mode_info *mode_info, struct_map_info *map_info, struct_cutset_info *cutset_info)
{
   long int i;
   int min_knockouts = 999999;
   int max_knockouts = -1;
   int num_knockouts;

   for( i = 0; i < cutset_info->mcs_idx; i++ )
   {
      num_knockouts = get_num_knockouts(&cutset_info->arr_mcs[i*mode_info->num_unit_size_orig], mode_info, map_info);

      if( num_knockouts <= 0 )
      {
         printf("WARNING: get_num_knockouts() returned (%d): i=%ld g_mcs_idx=%llu\n",num_knockouts,i,cutset_info->mcs_idx);
      }

      if( min_knockouts > num_knockouts )
      {
         min_knockouts = num_knockouts;
      }

      if( max_knockouts < num_knockouts )
      {
         max_knockouts = num_knockouts;
      }
   }

   printf("Minimum number of knockouts %3d\n",min_knockouts);
   printf("Maximum number of knockouts %3d\n",max_knockouts);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int get_num_knockouts(unsigned long long *cs, struct_mode_info *mode_info, struct_map_info *map_info)
{
   int u,b,i;

   int num_knockouts = 0;

   for( u = mode_info->num_unit_size_orig - 1; u >= 0; u-- )
   {
      unsigned long long tmp_unit = 0x8000000000000000ULL;
      for( b = 8*sizeof(unsigned long long) - 1; b >= 0; b-- )
      {
         if( cs[u] & tmp_unit )
         {
            num_knockouts++;
         }
         tmp_unit >>= 1;
      }
   }
   for( i = 0; i < map_info->num_always_zero_reacs; i++ )
   {
      num_knockouts++;
   }

   return num_knockouts;
}
//////////////////////////////////////////////////////
