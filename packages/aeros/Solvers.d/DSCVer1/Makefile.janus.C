#!
# Makefile for examples using solver library
# written by Padma Raghavan
#
#___________________________________________________
#
#
#
#
#
LDFLAGS 	= -lcsmath_cop -lf    -lmpi  -lm  
# this gives the best performance (at present)
#LDFLAGS 	= -lcsmath_r   -lf  -lmpi  -lm 
# janus specific -- non proc2 version of lib
#LDFLAGS 	= -lcsmath     -lf  -lmpi  -lm  
# janus specific -- non proc2 version of lib
#LDFLAGS 	= -lkmath        -lmpi  -lm 
# janus specific -- non proc2 version of lib
#
#
SRC = generate_test_mat.c
OBJ = generate_test_mat.o
#
#
#rules
#
#
.c.o:
	$(CC) $(SWITCHES) $(CPPFLAGS) -c  $<
#
#_________________________________________________________________________________________________
# Make  Example 1
#_________________________________________________________________________________________________
#
Example1_dbl:
	make -f Makefile.janus.C \
       "CC=cicc" "SWITCHES= -O3 -DDSC_LBLAS3 -DDSC_DBLAS2 -DF_NEEDS_UNDSC -DDBL_R_NUM" \
       "DSC_LIBRARY = ./DSC_LIB/dsclibdbl.a" example1
#
#
Example1_flt:
	make -f Makefile.janus.C \
       "CC=cicc" "SWITCHES= -O3 -DDSC_LBLAS3 -DDSC_DBLAS2  -DF_NEEDS_UNDSC " \
        "DSC_LIBRARY = ./DSC_LIB/dsclibflt.a" example1
#
example1: example1.c ${OBJ} ${DSC_LIBRARY}
	$(CC) $(SWITCHES)   -o Solve1 example1.c  \
	${OBJ} ${DSC_LIBRARY} $(LDFLAGS)

#  
#_________________________________________________________________________________________________
# Make  Example 2
#_________________________________________________________________________________________________
#
#
Example2_dbl:
	make -f Makefile.janus.C \
       "CC=cicc" "SWITCHES= -O3 -DDSC_LBLAS3 -DDSC_DBLAS2 -DF_NEEDS_UNDSC -DDBL_R_NUM" \
       "DSC_LIBRARY = ./DSC_LIB/dsclibdbl.a" example2
#
#
Example2_flt:
	make -f Makefile.janus.C \
       "CC=cicc" "SWITCHES= -O3 -DDSC_LBLAS3 -DDSC_DBLAS2 -DF_NEEDS_UNDSC " \
       "DSC_LIBRARY = ./DSC_LIB/dsclibflt.a" example2
#
#
example2: example2.c ${OBJ} ${DSC_LIBRARY}
	$(CC) $(SWITCHES) -o Solve2 example2.c  \
	${OBJ} ${DSC_LIBRARY} $(LDFLAGS)

#  
#
