#include	<mpi.h>
#include	<cstdio>
#include	<cstdlib>
#include	<cmath>
#ifdef  DBL_R_NUM
	#define		real_number_type  	double
#else
	#define		real_number_type	float
#endif

#define		DSC_MAX_NONZ_PER_ROW		500
/* change this if any column has more than 500 nonzeroes */
#define		DSC_BLAS_BLK_SIZE		32
/* change this if needed; a larger number may lead to better performance */
#define         DSC_MAX_MEM1                   100
#define         DSC_MAX_MEM2                   100 
#define		DSC_EMPTY			-1
#define		DSC_PARAM0			0
#define		DSC_PARAM1			1
#define		DSC_PARAM2			2
#define		DSC_MAX_REPS			12
#define		DSC_LINK_and_LINK		1
#define		DSC_LINK_and_NOT_LINK		2
#define		DSC_MAX_ARGS			10
#define         DSC_MAX_P_LEVEL     		6
#define		DSC_no_error			0
#define		DSC_argument_error		1
#define		DSC_derived_data_error		2
#define		DSC_fopen_error			4
#define		DSC_internal_error		10
#define		DSC_malloc_error	  	20	
#define		DSC_TWICE_MAX_PROCS	 	130	
#define		DSC_Solve_Gather_Type		80000
#define		DSC_Order_Gather_Type1		84000
#define		DSC_Order_Gather_Type2		85000
#define		DSC_Order_Gather_Type3		86000
#define		DSC_CONTINUE_TYPE		90000
#define		DSC_STAT_TYPE		        99000	
#define		DSC_STOP_TYPE			95000
#define		DSC_NO_ERROR			0
#define		DSC_stage_undefined		-1
#define		DSC_stage_clean			0
#define		DSC_stage_init			1
#define		DSC_stage_order			2
#define		DSC_stage_s_factor		3
#define		DSC_stage_nonz_input		4
#define		DSC_stage_n_factor		5	
#define		DSC_stage_n_factor_si		6	
#define		DSC_stage_f_solve		7	
#define		DSC_stage_b_solve		8	
#define		DSC_stage_f_solve_si		9	
#define		DSC_stage_b_solve_si		10	
#include	"dscstats.h"
#include	"dsctypes.h"
