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
#include "read_files.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void readin_reactions_file(char *filename, struct_reac_info *reac_info)
{
   FILE *reac_stream;
   char *line = NULL;
   char *str_start = NULL;
   char *reac;
   int num_reacs = 0;
   int max_len_reaction_name = 0;

   str_start = line;
   ///////////////////////////////////////////////////
   // open file for reading                         //
   ///////////////////////////////////////////////////
   printf("DEBUG: going to open file '%s' ...\n",filename);
   if( (reac_stream = fopen(filename,"r")) == (FILE *)0 )
   {
      fprintf(stderr, "FATAL ERROR: open file '%s' for reading failed: %s\n",filename,strerror(errno));
      exit(EXIT_FAILURE);
   }

   line = get_line(reac_stream,filename);
   str_start = line;

   if( (reac = (char *) malloc( (size_t) (sizeof(char)*(strlen(line)+1)))) == NULL )
   {
      fprintf(stderr, "FATAL ERROR: readin_reactions_file(): couldn't allocate memory for reac container\n");
      exit(EXIT_FAILURE);
   }

   // determine number of reactions
   // maximum length of reaction name
   while( get_next_reaction(&str_start,reac) >= 0 )
   {
      num_reacs++;
      // printf("DEBUG: Analyzing reaction '%s', length: %d\n",reac,(int)strlen(reac));
      if( (int)strlen(reac) > max_len_reaction_name )
      {
         max_len_reaction_name = (int)strlen(reac);
      }
   }
   reac_info->num_reactions = num_reacs;
   reac_info->max_len_reac_name = max_len_reaction_name + 1;

   printf("DEBUG: readin_reactions_file(): reac_info->num_reactions=%u\n",reac_info->num_reactions);
   printf("INFO: found %d reactions in reaction file '%s'\n",num_reacs,filename);
   printf("INFO: maximum length of reaction name: %d\n",max_len_reaction_name);

   allocate_mem_reac_info(reac_info);

   str_start = line;
   num_reacs = 0;
   while( get_next_reaction(&str_start,reac) >= 0 )
   {
      printf("DEBUG: found reaction number %d: %s\n",num_reacs,reac);
      if( num_reacs >= reac_info->num_reactions )
      {
         fprintf(stderr, "FATAL ERROR: Race condition occurred? Has reaction file '%s' changed during read process?\n",filename);
         fprintf(stderr, "             Number of reactions found in reactions is different to previous read process!\n");
         exit(EXIT_FAILURE);
      }
      strcpy(reac_info->reactions[num_reacs],reac);
      num_reacs++;
   }

   fclose(reac_stream);
   free(line);
   free(reac);

   printf("DEBUG: finished processing file '%s'.\n",filename);

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int get_next_reaction(char **str_start_inout, char *reac)
{
   char *str_dblquote_start;
   char *str_dblquote_stop;
   char *str_start;
   str_start = *str_start_inout;

   str_dblquote_start = strstr(str_start,"\"");

   if( str_dblquote_start == NULL )
   {
      // it seems there are no reactions left
      return(-1);
   }

   str_dblquote_stop = strstr(str_dblquote_start + 1,"\"");

   if( str_dblquote_stop == NULL )
   {
      fprintf(stderr, "FATAL ERROR: invalid syntax in reaction file: missing double quote\n");
      exit(EXIT_FAILURE);
   }

   if( str_dblquote_stop == str_dblquote_start + 1 )
   {
      fprintf(stderr, "FATAL ERROR: invalid syntax in reaction file: reaction name with length 0\n");
      exit(EXIT_FAILURE);
   }

   strncpy(reac, str_dblquote_start + 1, str_dblquote_stop - str_dblquote_start - 1);
   reac[str_dblquote_stop - str_dblquote_start -1] = '\0';

   *str_start_inout = str_dblquote_stop + 1;


   return(0);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
char *get_line(FILE * f, char *filename)
{
    size_t size = 0;
    size_t len  = 0;
    size_t last = 0;
    char * buf  = NULL;

    printf("DEBUG: entered get_line()\n");

    do
    {
        size += LINE_BUFFER_INCREMENT;
        buf = realloc(buf,size);
        if( buf == NULL )
        {
           fprintf(stderr, "FATAL ERROR: get_line: couldn't allocate %lu butes for buf\n",size);
           exit(EXIT_FAILURE);
        }


        if( fgets(buf+len,LINE_BUFFER_INCREMENT,f) == (char *)0 )
        {
           fprintf(stderr, "FATAL ERROR: failed to read line from file '%s': %s\n",filename,strerror(errno));
           exit(EXIT_FAILURE);
        }

        len = strlen(buf);
        last = len - 1;
    } while (!feof(f) && buf[last]!='\n');

    printf("DEBUG: leaving get_line()\n");

    return buf;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void allocate_mem_reac_info(struct_reac_info *reac_info)
{
   int i;

   // printf("DEBUG: entered allocate_mem_reac_info(): num_reactions: %d\n",reac_info->num_reactions);
   // printf("                                         max_len_reac_name: %d\n",reac_info->max_len_reac_name);

   reac_info->reactions = (char **) malloc((size_t) sizeof(char *)*reac_info->num_reactions);
   if( reac_info->reactions == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for reac_info->reactions\n",
      sizeof(char *)*reac_info->num_reactions);
      exit(EXIT_FAILURE);
   }


   for( i = 0; i < reac_info->num_reactions; i++ )
   {
      reac_info->reactions[i] = (char *) malloc((size_t) sizeof(char)*reac_info->max_len_reac_name);
      if( reac_info->reactions[i] == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for reac_info->reactions[%d]\n",
                 sizeof(char)*reac_info->max_len_reac_name,i);
         fflush(stderr);
         exit(EXIT_FAILURE);
      }
   }
   // printf("DEBUG: leaving allocate_mem_reac_info()\n");
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_reac_info(struct_reac_info *reac_info)
{
   int i;

   // printf("DEBUG: entered free_mem_reac_info(): num_reactions: %d\n",reac_info->num_reactions);

   for( i = 0; i < reac_info->num_reactions; i++ )
   {
      free(reac_info->reactions[i]);
   }
   // printf("DEBUG: leaving free_mem_reac_info()\n");

   if( reac_info->readin_essential_idx != NULL )
   {
      free(reac_info->readin_essential_idx);
   }

   if( reac_info->readin_essential_reactions != NULL )
   {
      for( i = 0; i < reac_info->num_reactions; i++ )
      {
         free(reac_info->readin_essential_reactions[i]);
      }
      free(reac_info->readin_essential_reactions);
   }

   free(reac_info->reactions);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void free_mem_mode_info(struct_mode_info *mode_info)
{
   printf("DEBUG: entered free_mem_mode_info()\n");

   free(mode_info->pnt_reac_occ);
   free(mode_info->pnt_reac_occ_bad);
   free(mode_info->essential_idx);
   free(mode_info->is_essential);
   free(mode_info->pnt_tmp_cutset);
   free(mode_info->pnt_flux);
   free(mode_info->pnt_good);
   free(mode_info->pnt_good_duplicates);
   free(mode_info->norm_bitmap);

   printf("DEBUG: leaving free_mem_mode_info()\n");
}
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void readin_efm_file_bin(struct_cmd_options cmd_options, struct_reac_info reac_info, struct_mode_info *mode_info)
{
   int fh;
   int r;
   long int m = 0;
   int bytes_read;
   int read_len;
   unsigned long long int tmp_unit;
   int unit_cnt = 0;
   int bitmover = 0;
   unsigned int t_num_reactions;

   printf("DEBUG: going to process file '%s' ...\n",cmd_options.m_filename);
   ///////////////////////////////////////////////////
   // open file for reading                         //
   ///////////////////////////////////////////////////
   printf("DEBUG: going to open file '%s' ...\n",cmd_options.m_filename);
   fh = open(cmd_options.m_filename, O_RDONLY);

   if( fh <= 0 )
   {
      fprintf(stderr, "FATAL ERROR: open file '%s' for reading failed: %s\n",cmd_options.m_filename,strerror(errno));
      exit(EXIT_FAILURE);
   }
   ///////////////////////////////////////////////////

   ///////////////////////////////////////////////////
   // read number of modes
   ///////////////////////////////////////////////////
   if( read(fh, &(mode_info->num_efms), 8) != 8 )
   {
      fprintf(stderr, "FATAL ERROR: couldn't read 8 bytes for number of emfs\n");
      exit(EXIT_FAILURE);
   }
   printf("INFO: number of efms: %llu\n",mode_info->num_efms);
   mode_info->good_emfs      = cmd_options.good_efms;
   mode_info->good_emfs_orig = cmd_options.good_efms;
   mode_info->bad_efms = mode_info->num_efms - cmd_options.good_efms;
   mode_info->num_good_removed_by_always_zero = 0;
   printf("INFO: number of good efms: %llu\n",cmd_options.good_efms);
   printf("INFO: number of bad efms: %llu\n",mode_info->bad_efms);
   ///////////////////////////////////////////////////

   ///////////////////////////////////////////////////
   // read number of reactions
   ///////////////////////////////////////////////////
   if( read(fh, &t_num_reactions, 4) != 4 )
   {
      fprintf(stderr, "FATAL ERROR: couldn't read 4 bytes for number of reactions\n");
      exit(EXIT_FAILURE);
   }
   printf("Read from bin-file number of reactions: %u\n",t_num_reactions);
   ///////////////////////////////////////////////////
   if( t_num_reactions != reac_info.num_reactions )
   {
      fprintf(stderr, "FATAL ERROR: number of reaction (%d) read from mode file '%s' ",t_num_reactions,cmd_options.m_filename);
      fprintf(stderr, " is not equal to number of reactions (%d) ",reac_info.num_reactions);
      fprintf(stderr, " read from reactions file '%s'\n",cmd_options.r_filename);
      exit(EXIT_FAILURE);
   }

   if( t_num_reactions == 0 )
   {
      fprintf(stderr, "FATAL ERROR: number of reaction (%d) read from mode file '%s' is equal 0",t_num_reactions,cmd_options.m_filename);
      fprintf(stderr, "             check your input files\n");
      exit(EXIT_FAILURE);
   }

   do_mem_allocation_bin_first(mode_info, reac_info, cmd_options);

   ///////////////////////////////////////////////////
   ///////////////////////////////////////////////////
   read_len = mode_info->num_unit_size_first*sizeof(unsigned long long int);

   for( m = 0; m < cmd_options.good_efms; m++ )
   {
      bytes_read = read( fh, &mode_info->pnt_good_first[mode_info->num_unit_size_first*m], read_len);
      if( bytes_read != read_len )
      {
         fprintf(stderr, "FATAL ERROR: could not read %d bytes at mode %lu: %d: %s\n",read_len,m,bytes_read,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      for( r = 0; r < reac_info.num_reactions; r++ )
      {
         tmp_unit = 1;
         bitmover = r;
         unit_cnt = r/(8*sizeof(unsigned long long int));
         bitmover -= unit_cnt*8*sizeof(unsigned long long int);
         tmp_unit <<= bitmover;

         if (mode_info->pnt_good_first[mode_info->num_unit_size_first*m + unit_cnt] & tmp_unit )
         {
            mode_info->pnt_reac_occ[r]++;
         }
      }
   }
   for( m = 0; m < mode_info->bad_efms; m++ )
   {
      bytes_read = read( fh, &mode_info->pnt_flux_first[mode_info->num_unit_size_first*m], read_len);
      if( bytes_read != read_len )
      {
         fprintf(stderr, "FATAL ERROR: could not read %d bytes at mode %lu: %d: %s\n",read_len,m,bytes_read,strerror(errno));
         fprintf(stderr, "             execution aborted.\n");
         exit(EXIT_FAILURE);
      }

      int norm = 0;
      for( r = 0; r < reac_info.num_reactions; r++ )
      {
         tmp_unit = 1;
         bitmover = r;
         unit_cnt = r/(8*sizeof(unsigned long long));
         bitmover -= unit_cnt*8*sizeof(unsigned long long);
         tmp_unit <<= bitmover;

         if (mode_info->pnt_flux_first[mode_info->num_unit_size_first*m + unit_cnt] & tmp_unit )
         {
            mode_info->pnt_reac_occ_bad[r] += 1;
            norm++;
         }
      }

      if( norm == 0 )
      {
         fprintf(stderr, "FATAL ERROR: norm of mode (%ld) is equal to zero!\n", m);
         fprintf(stderr, "             this bad mode can never be knocked out!\n");
         exit(EXIT_FAILURE);
      }

      if( (m+1)%100000 == 0 )
      {
         printf("INFO: processed %09lu of %09llu modes: %03llu%%\n",m,mode_info->num_efms,(m+1)*100/mode_info->num_efms);
      }
   }
   ///////////////////////////////////////////////////

   close(fh);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void do_mem_allocation_bin_first(struct_mode_info *mode_info, struct_reac_info reac_info, struct_cmd_options cmd_options)
{
   unsigned long i;
   unsigned long long t_zero = 0;

   mode_info->essential_idx = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info.num_reactions);
   if( mode_info->essential_idx == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for essential_idx\n",sizeof(unsigned int)*reac_info.num_reactions);
      exit(EXIT_FAILURE);
   }


   mode_info->pnt_reac_occ = (unsigned long long int *) malloc( (size_t) (reac_info.num_reactions*sizeof(unsigned long long int) ) );
   if( mode_info->pnt_reac_occ == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for reaction occurrence array\n");
      exit(EXIT_FAILURE);
   }
   for( i = 0; i < reac_info.num_reactions; i++ )
   {
      mode_info->pnt_reac_occ[i] = 0;
   }

   mode_info->pnt_reac_occ_bad = (unsigned long long int *) malloc( (size_t) (reac_info.num_reactions*sizeof(unsigned long long int) ) );
   if( mode_info->pnt_reac_occ_bad == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for reaction occurrence array of modes to be removed\n");
      exit(EXIT_FAILURE);
   }
   for( i = 0; i < reac_info.num_reactions; i++ )
   {
      mode_info->pnt_reac_occ_bad[i] = 0;
   }

   mode_info->num_unit_size_first = ceil(((float)reac_info.num_reactions)/(8*sizeof(unsigned long long)));
   printf("DEBUG: a single mode consists of %d 'unsigned long long' elements with a size of %lu bytes\n",
          mode_info->num_unit_size_first,sizeof(unsigned long long));
   printf("DEBUG: going to allocate %llu bytes for %llu 'unsigned long long' for flux values\n",
           mode_info->num_unit_size_first*mode_info->num_efms*sizeof(unsigned long long),mode_info->num_unit_size_first*mode_info->num_efms);

   mode_info->pnt_flux_first = (unsigned long long int*) calloc( (size_t) (mode_info->num_unit_size_first*mode_info->num_efms),
                                                                 sizeof(unsigned long long int));

   if( mode_info->pnt_flux_first == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for binary flux value (mode_info->pnt_flux_first)\n");
      exit(EXIT_FAILURE);
   }

   // initialize new memory
   for( i = 0; i < mode_info->num_unit_size_first*mode_info->num_efms; i++ )
   {
      mode_info->pnt_flux_first[i] = t_zero;
   }

   mode_info->pnt_good_first = (unsigned long long int*) calloc( (size_t) (mode_info->num_unit_size_first*cmd_options.good_efms),
                                                                 sizeof(unsigned long long int));

   if( mode_info->pnt_good_first == NULL )
   {
      fprintf(stderr, "FATAL ERROR: couldn't allocate memory for binary flux value (mode_info->pnt_flux_first)\n");
      exit(EXIT_FAILURE);
   }

   // initialize new memory
   for( i = 0; i < mode_info->num_unit_size_first*cmd_options.good_efms; i++ )
   {
      mode_info->pnt_good_first[i] = t_zero;
   }
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void readin_essential_file(char *filename, struct_reac_info *reac_info)
{           
   FILE *essential_stream;
   char *reac;
   char *line = NULL;
   char *str_start = NULL;
   int num_reacs = 0;
   int i,j;
   int found;

   allocate_readin_essential_memory(reac_info);

   str_start = line;
   ///////////////////////////////////////////////////
   // open file for reading                         //
   ///////////////////////////////////////////////////
   printf("DEBUG: going to open file '%s' ...\n",filename);
   if( (essential_stream = fopen(filename,"r")) == (FILE *)0 )
   {
      fprintf(stderr, "FATAL ERROR: open file '%s' for reading failed: %s\n",filename,strerror(errno));
      exit(EXIT_FAILURE);
   }
   
   line = get_line(essential_stream,filename);
   str_start = line;
   
   if( (reac = (char *) malloc( (size_t) (sizeof(char)*(strlen(line)+1)))) == NULL )
   {  
      fprintf(stderr, "FATAL ERROR: readin_reactions_file(): couldn't allocate memory for reac container\n");
      exit(EXIT_FAILURE);
   }

   while( get_next_reaction(&str_start,reac) >= 0 )
   {
      printf("DEBUG: found essential reaction number %d: %s\n",num_reacs,reac);
      if( num_reacs >= reac_info->num_reactions )
      {
         fprintf(stderr, "FATAL ERROR: number of reaction (%d) read in from '%s' is larger than number of reactions (%d) read from raections file\n",num_reacs,filename,reac_info->num_reactions);
         exit(EXIT_FAILURE);
      }
      strcpy(reac_info->readin_essential_reactions[num_reacs],reac);
      num_reacs++;
   }

   reac_info->num_readin_essential_reacs = num_reacs;
   printf("INFO: %d essential reactions specified in file '%s'\n",reac_info->num_readin_essential_reacs,filename);

   if( reac_info->num_readin_essential_reacs >= reac_info->num_reactions)
   {
      fprintf(stderr, "FATAL ERROR: number of essential reactions (%d) is larger than or equal to number of reactions (%d)",
              reac_info->num_readin_essential_reacs,reac_info->num_reactions);
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < reac_info->num_readin_essential_reacs; i++ )
   {
      found = 0;
      for( j = 0; j < reac_info->num_reactions; j++ )
      {
         if( strcmp(reac_info->reactions[j],reac_info->readin_essential_reactions[i]) == 0 )
         {
            found = 1;
            reac_info->readin_essential_idx[i] = j;
            printf("DEBUG: index of essential reaction \"%s\" is %d\n",reac_info->readin_essential_reactions[i],reac_info->readin_essential_idx[i]);
            break;
         }
      }
      if( found == 0 )
      {
         fprintf(stderr, "FATAL ERROR: could find essential reaction \"%s\" in array contain all reactions\n",reac_info->readin_essential_reactions[i]);
         exit(EXIT_FAILURE);
      }
   }


   fclose(essential_stream);
   free(line);
   free(reac);

   printf("DEBUG: finished processing file '%s'.\n",filename);
   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void allocate_readin_essential_memory(struct_reac_info *reac_info)
{
   int i;

   reac_info->readin_essential_reactions = (char **) malloc((size_t) sizeof(char *)*reac_info->num_reactions);
   if( reac_info->readin_essential_reactions == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_readin_essential_memory(): couldn't allocate %lu bytes for reac_info->readin_essential_reactions\n",
              sizeof(char *)*reac_info->num_reactions);
      exit(EXIT_FAILURE);
   }

   for( i = 0; i < reac_info->num_reactions; i++ )
   {
      reac_info->readin_essential_reactions[i] = (char *) malloc((size_t) sizeof(char)*reac_info->max_len_reac_name);
      if( reac_info->readin_essential_reactions[i] == NULL )
      {
         fprintf(stderr, "FATAL ERROR: allocate_readin_essential_memory(): couldn't allocate %lu bytes for reac_info->readin_essential_reactions[%d]\n",
                 sizeof(char)*reac_info->max_len_reac_name, i);
         exit(EXIT_FAILURE);
      }
   }

   reac_info->readin_essential_idx = (unsigned int *) malloc((size_t) sizeof(unsigned int)*reac_info->num_reactions);
   if( reac_info->readin_essential_idx == NULL )
   {
      fprintf(stderr, "FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for reac_info->readin_essential_idx\n",
              sizeof(unsigned int)*reac_info->num_reactions);
      exit(EXIT_FAILURE);
   }

}
//////////////////////////////////////////////////////
