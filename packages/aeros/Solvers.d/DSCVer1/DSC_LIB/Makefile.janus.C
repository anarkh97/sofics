#!
# Makefile for DSC-solver library
# written by Padma Raghavan 
#
#
# Use -DDBL_R_NUM to create lib with datatype double
# without this flag, lib created is datatype float
#
#
#
#
CSRC =	io_stats.c group_comm.c	 environ_dep_mpi.c	free_routines.c
COBJ =	io_stats.o group_comm.o	 environ_dep_mpi.o	free_routines.o
#
#
OSRC =	order_util_spd1.c	order_util_spd2.c  order_util_spd3.c order_algs.c
OOBJ =	order_util_spd1.o 	order_util_spd2.o  order_util_spd3.o order_algs.o
#
#
NSSRC =	s_factor.c		s_factor_util_spd.c\
	n_factor.c		n_factor_util_spd.c\
	d_n_factor_spd.c	n_work_blas3.c\
	n_factor_blas1.c	n_factor_blas2_3.c\
	n_solve.c		n_solve_util_spd.c\
	d_n_solve_spd.c		n_solve_algs.c
NSOBJ =	s_factor.o		s_factor_util_spd.o\
	n_factor.o		n_factor_util_spd.o\
	d_n_factor_spd.o	n_work_blas3.o\
	n_factor_blas1.o	n_factor_blas2_3.o\
	n_solve.o		n_solve_util_spd.o\
	d_n_solve_spd.o		n_solve_algs.o
#
#
DSCSRC	= ${NSSRC} ${CSRC} ${OSRC} interface.c
DSCOBJ	= ${NSOBJ} ${COBJ} ${OOBJ} interface.o
#
#rules
#
#
.c.o:
	$(CC) $(SWITCHES) $(CPPFLAGS) -c  $<
#
#
#
lib_dbl: 
	make -f Makefile.janus.C\
       "CC=cicc" "SWITCHES=-O3 -DDSC_LBLAS3 -DDSC_DBLAS2 -DF_NEEDS_UNDSC  -DDBL_R_NUM " dsc_lib_dbl
#
# -DDBL_R_NUM -- use doubles (if not specified, uses floats)
# -DF_NEEDS_UNDSC --  calls fortran functions with _
# -DDSC_EFF -- if not defined, does extra handshaking to fail better
#
#
lib_flt: 
	make -f Makefile.janus.C \
       "CC=cicc" "SWITCHES=-O3 -DDSC_LBLAS3 -DDSC_DBLAS2 -DF_NEEDS_UNDSC " dsc_lib_flt
#
#  if using C++ compiler add flag -DDSC_CPLUS
#
dsc_lib_dbl:${DSCSRC} ${DSCOBJ} 	
	xar  cr dsclibdbl.a ${DSCOBJ}
	make -f Makefile.janus.C clean
#
#
dsc_lib_flt:${DSCSRC} ${DSCOBJ} 	
	xar  cr dsclibflt.a ${DSCOBJ}
	make -f Makefile.janus.C clean
#
#
clean:	${DSCOBJ} 
	rm ${DSCOBJ} 
#
#
#

