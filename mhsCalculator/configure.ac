################################################################################
################################################################################
# Author:  Christian Jungreuthmayer
# Date:    Tue Oct 15 09:13:38 CEST 2013
# Company: Austrian Centre of Industrial Biotechnology (ACIB)
################################################################################

AC_INIT(mhsCalculator, version-0.1,christian.jungreuthmayer@boku.ac.at)
AC_MSG_NOTICE([Going to create Makefile for the toolbox 'mhsCalculator'])

BUILD_GLPK_EXEC=yes
BUILD_CPLEX_EXEC=yes
BUILD_GUROBI_EXEC=yes

echo "Testing for a C compiler"
AC_PROG_CC

echo "Testing for general header files"
AC_CHECK_HEADERS(stdio.h)
AC_CHECK_HEADERS(stdlib.h)
AC_CHECK_HEADERS(errno.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(string.h)
AC_CHECK_HEADERS(sys/time.h)


echo "Testing if math library is installed"
AC_CHECK_HEADERS(math.h)
AC_CHECK_LIB([m],[fabs],[LIBM=true], [LIBM=false], [])
if test "$LIBM" == false
then
	AC_MSG_ERROR([mathematic library 'libm' not found!]);
fi

echo "Testing if pthread library is installed"
AC_CHECK_HEADERS(pthread.h)
AC_CHECK_LIB([pthread], [pthread_create], [LIBPTHREAD=true], [LIBPTHREAD=false], [])

echo "Testing if GLPK is installed"
AC_CHECK_HEADERS(glpk.h)
AC_CHECK_LIB([glpk], [glp_create_prob], [LIBGLPK=true], [LIBGLPK=false], [-lm])

echo "Testing if CPLEX is installed"
AC_CHECK_HEADERS(ilcplex/cplex.h)
AC_CHECK_LIB([cplex],[CPXopenCPLEX], [LIBCPLEX=true], [LIBCPLEX=false], [-lpthread -lm])

echo "Testing if GUROBI is installed"
AC_CHECK_HEADERS(gurobi_c.h)
AC_CHECK_LIB([gurobi75],[GRBloadenv], [LIBGUROBI=true], [LIBGUROBI=false], [-lpthread -lm])

###############################################################################
# check if we can build executable that uses GLPK
###############################################################################
if test "$ac_cv_header_glpk_h" == no
then
        BUILD_GLPL_EXEC=no
	AC_MSG_WARN([GLPK header is missing!]);
fi

if test "$LIBGLPK" == false
then
        BUILD_GLPL_EXEC=no
	AC_MSG_WARN([CPLEX library 'libglpk' is missing!]);
fi
###############################################################################

###############################################################################
# check if we can build executable that uses CPLEX
###############################################################################
if test "$ac_cv_header_ilcplex_cplex_h" == no
then
        BUILD_CPLEX_EXEC=no
	AC_MSG_WARN([CPLEX header is missing!]);
fi

if test "$LIBPTHREAD" == false
then
        BUILD_CPLEX_EXEC=no
	AC_MSG_WARN([pthread library is NOT missing!]);
fi

if test "$LIBCPLEX" == false
then
        BUILD_CPLEX_EXEC=no
	AC_MSG_WARN([CPLEX library 'libcplex' is missing!]);
fi
###############################################################################

###############################################################################
# check if we can build executable that uses GUROBI
###############################################################################
if test "$ac_cv_header_gurobi_c_h" == no
then
        BUILD_GUROBI_EXEC=no
	AC_MSG_WARN([GUROBI header is missing!]);
fi

if test "$LIBPTHREAD" == false
then
        BUILD_GUROBI_EXEC=no
	AC_MSG_WARN([pthread library is NOT missing!]);
fi

if test "$LIBGUROBI" == false
then
        BUILD_GUROBI_EXEC=no
	AC_MSG_WARN([GUROBI library 'libgurobi75' is missing!]);
fi
###############################################################################

###############################################################################
# Print which of the executables we can build and which we cannot
###############################################################################
if test "$BUILD_CPLEX_EXEC" == no
then
	AC_MSG_WARN([Executable using CPLEX is NOT going to be built!]);
else
        AC_MSG_NOTICE([Executable using CPLEX is going to be built!]);
fi

if test "$BUILD_GUROBI_EXEC" == no
then
	AC_MSG_WARN([Executable using GUROBI is NOT going to be built!]);
else
        AC_MSG_NOTICE([Executable using GUROBI is going to be built!]);
fi

if test "$BUILD_GLPK_EXEC" == no
then
	AC_MSG_WARN([Executable using GLPK is NOT going to be built!]);
else
        AC_MSG_NOTICE([Executable using GLPK is going to be built!]);
fi
###############################################################################


###############################################################################
# create main target for Makefile
###############################################################################
if test "$BUILD_GLPK_EXEC" == yes
then
   if test "$BUILD_CPLEX_EXEC" == yes 
   then
       if test "$BUILD_GUROBI_EXEC" == yes 
       then
           # we build all three supported version: GPLK, CPLEX, GUROBI
           AC_SUBST(TARGET,["mhsCalculator_linprog_gurobi mhsCalculator_linprog_glpk mhsCalculator_linprog_cplex mhsCalculator_linprog_cplex_solpool"])
       else
           # we build: GPLK, CPLEX
           AC_SUBST(TARGET,["mhsCalculator_linprog_glpk mhsCalculator_linprog_cplex mhsCalculator_linprog_cplex_solpool"])
       fi
   else
       if test "$BUILD_GUROBI_EXEC" == yes 
       then
           # we build: GPLK, GUROBI
           AC_SUBST(TARGET,["mhsCalculator_linprog_gurobi mhsCalculator_linprog_glpk"])
       else
           # we build: GPLK
           AC_SUBST(TARGET,["mhsCalculator_linprog_glpk"])
       fi
   fi
else
   if test "$BUILD_CPLEX_EXEC" == yes 
   then
       if test "$BUILD_GUROBI_EXEC" == yes 
       then
           # we build: CPLEX, GUROBI
           AC_SUBST(TARGET,["mhsCalculator_linprog_gurobi mhsCalculator_linprog_cplex mhsCalculator_linprog_cplex_solpool"])
       else
           # we build: CPLEX
           AC_SUBST(TARGET,["mhsCalculator_linprog_cplex mhsCalculator_linprog_cplex_solpool"])
       fi
   else
       if test "$BUILD_GUROBI_EXEC" == yes 
       then
           # we build: GUROBI
           AC_SUBST(TARGET,["mhsCalculator_linprog_gurobi"])
       else
           # no supported linear programming solver found
           AC_MSG_WARN([No buildable target for linear programming version could be found -> just build Berge algorithm])
       fi
   fi
fi
###############################################################################

###############################################################################
# create Makefile
###############################################################################
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
###############################################################################
