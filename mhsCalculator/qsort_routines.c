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
#include "qsort_routines.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int int_cmp(const void *a, const void *b)
{
    const unsigned int *ia = (const unsigned int *)a; // casting pointer types 
    const unsigned int *ib = (const unsigned int *)b;
    return *ia  > *ib;
        /* integer comparison: returns negative if b > a 
        and positive if a > b */
}
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// sort modes by row
//////////////////////////////////////////////////////
int efm_cmp2(const void *a, const void *b)
{
    const struct_flux *ia = (const struct_flux *)a; // casting pointer types 
    const struct_flux *ib = (const struct_flux *)b;
    int u;

    for( u = ia->unit_size - 1; u >= 0; u-- )
    {
       if( ia->flux[u] > ib->flux[u])
       {
          // return(ia->flux[u] - ib->flux[u]);
          return(1);
       }
       else if( ia->flux[u] < ib->flux[u])
       {
          return(-1);
       }
    }

    return(0);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// sort modes by norm
//////////////////////////////////////////////////////
int efm_cmp_by_norm(const void *a, const void *b)
{
    const struct_flux *ia = (const struct_flux *)a; // casting pointer types 
    const struct_flux *ib = (const struct_flux *)b;

    return(ia->norm - ib->norm);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// sort modes by columns
//////////////////////////////////////////////////////
int efm_by_cols(const void *a, const void *b)
{
    const struct_mode *ia = (const struct_mode *)a; // casting pointer types 
    const struct_mode *ib = (const struct_mode *)b;
    int u;

    for( u = ia->unit_col_length - 1; u >= 0; u-- )
    {
       if( ia->flux[u] > ib->flux[u] )
       {
          return(1);
       }
       else if( ia->flux[u] < ib->flux[u] )
       {
          return(-1);
       }
    }

    return(0);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// sort bitorder array
//////////////////////////////////////////////////////
int efm_cmp_bitorder(const void *a, const void *b)
{
    const struct_bitorder *ia = (const struct_bitorder *)a; // casting pointer types 
    const struct_bitorder *ib = (const struct_bitorder *)b;

    return(ia->occ < ib->occ);
}
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
// sort bitorder array reverse
//////////////////////////////////////////////////////
int efm_cmp_bitorder_rev(const void *a, const void *b)
{
    const struct_bitorder *ia = (const struct_bitorder *)a; // casting pointer types 
    const struct_bitorder *ib = (const struct_bitorder *)b;

    return(ia->occ > ib->occ);
}
//////////////////////////////////////////////////////
