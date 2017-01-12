////////////////////////////////////////////////////////////////////////////////
// Author: Christian Jungreuthmayer
// Email: christian.jungreuthmayer@acib.at
// Company: Austrian Centre of Industrial Biotechnology (ACIB)
// Web: http://www.acib.at
// Copyright (C) 2012, 2013
// Published unter GNU Public License V3 (http://www.gnu.org/licenses/gpl.html)
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
#include <stdio.h>
#include <stdlib.h>
#include <sys/file.h>
#include <sys/time.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define VERSION "1.0"
#define YEARS   "2012, 2013"
#define AUTHORS "Christian Jungreuthmayer"
#define COMPANY "Austrian Centre of Industrial Biotechnology (ACIB)"

#define CONSIDERED_ZERO 1.0e-07
#define LINE_BUFFER_INCREMENT 2048

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
char **g_reactions;
unsigned long long int *g_modes_bin;
unsigned int g_num_reactions = 0;
int g_len_elem_reac_name = 0;
FILE *g_fho;
long int g_num_modes = 0;
int g_output_format = 0;
unsigned int g_num_unit_size = 0;
char *o_filename        = NULL;
char *r_filename        = NULL;
char *m_filename        = NULL;
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void usage(char *);
int get_next_reaction(char **str_start_inout, char *reac);
FILE* open_outputfile(char *filename);
void readin_reactions_file(char *filename);
void read_efm_file(char *filename);
void print_modes();
void handle_arguments(int largc, char **largv, char **lm_filename, char **lr_filename, char **lo_filename);
void allocate_misc_mem();
char *get_line(FILE * f, char *filename);
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
   printf("program '%s' started.\n",argv[0]);

   ////////////////////////////////////////////////////////////////////////////
   // deal with input arguments
   ////////////////////////////////////////////////////////////////////////////
   handle_arguments(argc, argv, &m_filename, &r_filename, &o_filename);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // open output file handle -> this is done here in order to avoid aborting
   // the execution after reading all the modes because of an invalid output
   // file name
   ////////////////////////////////////////////////////////////////////////////
   g_fho = open_outputfile(o_filename);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   readin_reactions_file(r_filename);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // read in file containing modes
   ////////////////////////////////////////////////////////////////////////////
   read_efm_file(m_filename);
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // write sorted modes to file
   ////////////////////////////////////////////////////////////////////////////
   printf("INFO: going to write converted modes back to file '%s'\n",o_filename);
   fflush(g_fho);
   print_modes();
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // clean up
   ////////////////////////////////////////////////////////////////////////////
   fclose(g_fho);

   ////////////////////////////////////////////////////////////////////////////

   printf("program '%s' stopped.\n",argv[0]);
   return(0);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void do_frees()
{
   int i;

   free(g_modes_bin);

   for( i = 0; i < g_num_reactions; i++ )
   {
      free(g_reactions[i]);
   }
   free(g_reactions);

}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void print_modes()
{
   long int r;
   struct timeval t_start_time;
   struct timeval t_stop_time;
   unsigned long long int t_diff_time_usec = 0;

   gettimeofday(&t_start_time,NULL);

   fwrite(&g_num_modes, 1, sizeof(g_num_modes), g_fho);
   fwrite(&g_num_reactions, 1, sizeof(g_num_reactions),    g_fho);

   printf("INFO: print_modes(): going to write %lu modes\n",g_num_modes);
   for( r = 0; r < g_num_modes; r++ )
   {
      fwrite(&g_modes_bin[r*g_num_unit_size], g_num_unit_size, sizeof(unsigned long long int), g_fho);

      if( (r+1)%10000 == 0 )
      {
         gettimeofday(&t_stop_time,NULL);
         t_diff_time_usec = (t_stop_time.tv_sec - t_start_time.tv_sec)*1000000 + (t_stop_time.tv_usec - t_start_time.tv_usec);
         printf("wrote line %ld of %ld = %ld%% runtime=%llus exp. total runtime=%llus still to go=%llus\n",
                r,g_num_modes,100*r/g_num_modes,t_diff_time_usec/1000000,t_diff_time_usec/1000000*g_num_modes/r,
                t_diff_time_usec/1000000*(g_num_modes-r)/r);

      }
   }
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void read_efm_file(char *filename)
{
   int r;
   long int l_read;
   FILE *mode_stream;
   long int num_bytes;
   int ret;
   char *mode_line;
   char *str_start;
   struct timeval t_start_time;
   struct timeval t_stop_time;
   unsigned long long int t_diff_time_usec = 0;
   double val;
   unsigned long long int tmp_unit = 1;
   int bit_cnt = 0;
   unsigned long int unit_cnt = 0;


   g_num_unit_size = ceil(((float)g_num_reactions)/(8*sizeof(unsigned long long int)));
   num_bytes = g_num_modes*g_num_unit_size*sizeof(unsigned long long int);
   printf("INFO: g_num_unit_size=%u\n",g_num_unit_size);
   printf("DEBUG: going to allocate %ld bytes = %ld Megabytes = %ld Gigabytes\n",num_bytes,num_bytes/1024/1024,num_bytes/1024/1024/1024);

   g_modes_bin = (unsigned long long int *) calloc((size_t) (g_num_modes*g_num_unit_size), sizeof(unsigned long long int));

   if( g_modes_bin == NULL )
   {
      printf("FATAL ERROR: couldn't allocate memory for mode values\n");
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   printf("DEBUG: going to open file '%s'\n",filename);
   if( (mode_stream = fopen(filename,"r")) == (FILE *)0 )
   {
      printf("FATAL ERROR: open file '%s' for reading failed: %s\n",filename,strerror(errno));
      exit(EXIT_FAILURE);
   }

   gettimeofday(&t_start_time,NULL);
   l_read = 0;
   while( l_read < g_num_modes && (mode_line = get_line(mode_stream,filename)) != NULL )
   {
      // printf("DEBUG: line read: %ld\n",l_read);

      str_start = mode_line;

      l_read++;

      tmp_unit = 1;
      bit_cnt = 0;

      for( r = 0; r < g_num_reactions; r++ )
      {
         // get rid of leading space and tabulator characters
         while( str_start[0] == 9 || str_start[0] == 32 )
         {
            str_start++;
         }

         ret = sscanf(str_start,"%le",&val);
         if( ret != 1 || ret == EOF)
         {
            printf("FATAL ERROR: couldn't read carbon flux value #%d in line %ld: rest of line '%s': %s\n",r,l_read,str_start,strerror(errno));
            printf("             execution aborted.\n");
            exit(EXIT_FAILURE);
         }

         // printf("val(r=%d)=%le ",r,val);

         // move forward to next white space character
         while( str_start[0] != 9 && str_start[0] != 32 && str_start[0] != 0 )
         {
            str_start++;
         }
 
         if( fabs(val) > CONSIDERED_ZERO )
         {
            if( unit_cnt >= g_num_unit_size*g_num_modes )
            {
               printf("FATAL ERROR: number of stored flux value (%lu) exceeds number of expected flux values (%lu)\n",unit_cnt,g_num_unit_size*g_num_modes);
               printf("             unit_cnt=%lu bit_cnt=%d val=%f r=%d g_num_unit_size=%d g_num_modes=%lu\n",unit_cnt,bit_cnt,val,r,g_num_unit_size,g_num_modes);
               exit(EXIT_FAILURE);
            }
            g_modes_bin[unit_cnt] |= tmp_unit;
         }
         tmp_unit <<= 1;
         bit_cnt++;
         if( bit_cnt%(8*sizeof(unsigned long long)) == 0 )
         {
            unit_cnt++;
            bit_cnt = 0;
            tmp_unit = 1;
         }
      }

      unit_cnt = l_read*g_num_unit_size;

      // printf("\n");

      if( l_read%10000 == 0 )
      {
         gettimeofday(&t_stop_time,NULL);
         t_diff_time_usec = (t_stop_time.tv_sec - t_start_time.tv_sec)*1000000 + (t_stop_time.tv_usec - t_start_time.tv_usec);
         printf("processed line %ld of %ld = %ld%% runtime=%llus exp. total runtime=%llus still to go=%llus\n",
                l_read,g_num_modes,100*l_read/g_num_modes,t_diff_time_usec/1000000,t_diff_time_usec/1000000*g_num_modes/l_read,
                t_diff_time_usec/1000000*(g_num_modes-l_read)/l_read);
      }

      free(mode_line);
   }

   if( l_read != g_num_modes )
   {
      printf("ERROR: number of lines read (%ld) is not equal to number of modes (%ld)\n",l_read,g_num_modes);
      printf("       execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   printf("INFO: number of loaded modes: %ld\n",g_num_modes);

   fclose(mode_stream);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
FILE* open_outputfile(char *filename)
{
   FILE *fh;

   printf("INFO: open_outputfile(): going to open output file '%s'\n",filename);
   fh = fopen(filename, "we");

   if( fh == (FILE *)0 )
   {
      printf("FATAL ERROR: open_outputfile(): couldn't open file '%s' for writing: %s\n",filename,strerror(errno));
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }

   return(fh);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void usage(char *message)
{
   printf("%s",message);
   printf("usage: convert_text_to_bin -m in_modes.txt -o out_modes.bin -n number_of_modes -r rfile\n");
   exit(EXIT_FAILURE);
}
////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
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
      printf("FATAL ERROR: invalid syntax in reaction file: missing double quote\n");
      exit(EXIT_FAILURE);
   }

   if( str_dblquote_stop == str_dblquote_start + 1 )
   {
      printf("FATAL ERROR: invalid syntax in reaction file: reaction name with length 0\n");
      exit(EXIT_FAILURE);
   }

   //if( ((str_dblquote_stop - 1) - (str_dblquote_start + 1) + 1) > MAX_LEN_REACTION_NAME )
   //{
   //   printf("FATAL ERROR: invalid syntax in reaction file: reaction name is too long: ...'%s'\n",str_dblquote_start + 1);
   //   exit(EXIT_FAILURE);
   //}

   strncpy(reac, str_dblquote_start + 1, str_dblquote_stop - str_dblquote_start - 1);
   reac[str_dblquote_stop - str_dblquote_start -1] = '\0';

   *str_start_inout = str_dblquote_stop + 1;


   return(0);
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void handle_arguments(int largc, char **largv, char **lm_filename, char **lr_filename, char **lo_filename)
{
   int r_opt;
   int i;
   char *num_modes = NULL;
   char *endpointer;

   printf("convert_text_to_bin Version: %s\n",VERSION);
   printf("Copyright (C) %s %s %s\n",YEARS, AUTHORS, COMPANY);
   printf("Executed command: %s\n",largv[0]);
   printf("Options:");
   for(i = 1; i < largc; i++ )
   {
      printf(" %s",largv[i]);
   }
   printf("\n");

   while(( r_opt = getopt(largc, largv, "hm:r:o:n:")) != -1 )
   {
      switch(r_opt)
      {
         case 'h':
            usage("");
            break;
         case 'm':
            *lm_filename = optarg;
            break;
         case 'r':
            *lr_filename = optarg;
            break;
         case 'o':
            *lo_filename = optarg;
            break;
         case 'n':
            num_modes = optarg;
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
   if( *lm_filename == NULL )
   {
      usage("ERROR: name of file containing modes (in text form) not provided\n");
   }

   if( *lr_filename == NULL )
   {
      usage("ERROR: name of file containing reaction names not provided\n");
   }

   if( num_modes == NULL )
   {
      usage("ERROR: number of modes not defined\n");
   }
   ////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////
   // number of 'good' efms
   ////////////////////////////////////////////////////////////////////////////
   g_num_modes = strtol(num_modes, &endpointer, 10);
   if( endpointer == num_modes || errno == EINVAL || errno == ERANGE )
   {
      printf("FATAL ERROR: error while converting number of modes (-n %s)\n",num_modes);
      printf("             execution aborted.\n");
      exit(EXIT_FAILURE);
   }
   printf("INFO: number of expected modes: %ld\n",g_num_modes);
   ////////////////////////////////////////////////////////////////////////////

   return;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void readin_reactions_file(char *filename)
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
      printf("FATAL ERROR: open file '%s' for reading failed: %s\n",filename,strerror(errno));
      exit(EXIT_FAILURE);
   }

   line = get_line(reac_stream,filename);
   str_start = line;

   if( (reac = (char *) malloc( (size_t) (sizeof(char)*(strlen(line)+1)))) == NULL )
   {
      printf("FATAL ERROR: readin_reactions_file(): couldn't allocate memory for reac container\n");
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
   g_num_reactions = num_reacs;
   g_len_elem_reac_name = max_len_reaction_name + 1;

   printf("DEBUG: readin_reactions_file(): g_num_reactions=%u\n",g_num_reactions);
   printf("INFO: found %d reactions in reaction file '%s'\n",num_reacs,filename);
   printf("INFO: maximum length of reaction name: %d\n",max_len_reaction_name);

   allocate_misc_mem();

   str_start = line;
   num_reacs = 0;
   while( get_next_reaction(&str_start,reac) >= 0 )
   {
      printf("DEBUG: found reaction number %d: %s\n",num_reacs,reac);
      if( num_reacs >= g_num_reactions )
      {
         printf("FATAL ERROR: Race condition occurred? Has reaction file '%s' changed during read process?\n",filename);
         printf("             Number of reactions found in reactions is different to previous read process!\n");
         exit(EXIT_FAILURE);
      }
      strcpy(g_reactions[num_reacs],reac);
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
// Code from http://sucs.org/Knowledge/Help/Program%20Advisory/Reading%20an%20arbitrarily%20long%20line%20in%20C
//////////////////////////////////////////////////////
char *get_line(FILE * f, char *filename)
{
    size_t size = 0;
    size_t len  = 0;
    size_t last = 0;
    char * buf  = NULL;

    // printf("DEBUG: entered get_line()\n");

    do
    {
        size += LINE_BUFFER_INCREMENT;
        buf = realloc(buf,size);
        if( buf == NULL )
        {
           printf("FATAL ERROR: get_line: couldn't allocate %lu butes for buf\n",size);
           exit(EXIT_FAILURE);
        }


        if( fgets(buf+len,LINE_BUFFER_INCREMENT,f) == (char *)0 )
        {
           printf("FATAL ERROR: get_line failed to read line from file '%s': %s\n",filename,strerror(errno));
           exit(EXIT_FAILURE);
        }

        len = strlen(buf);
        last = len - 1;
    } while (!feof(f) && buf[last]!='\n');

    // printf("DEBUG: leaving get_line(): read '%s'\n",buf);

    return buf;
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void allocate_misc_mem()
{
   int i;
   g_reactions = (char **) malloc((size_t) sizeof(char *)*g_num_reactions);
   if( g_reactions == NULL )
   {
      printf("FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for g_reactions\n",sizeof(char *)*g_num_reactions);
      exit(EXIT_FAILURE);
   }


   for( i = 0; i < g_num_reactions; i++ )
   {
      g_reactions[i] = (char *) malloc((size_t) sizeof(char)*g_len_elem_reac_name);
      if( g_reactions[i] == NULL )
      {
         printf("FATAL ERROR: allocate_misc_mem(): couldn't allocate %lu bytes for g_reactions[%d]\n",sizeof(char)*g_len_elem_reac_name,i);
         exit(EXIT_FAILURE);
      }
   }
}
//////////////////////////////////////////////////////
