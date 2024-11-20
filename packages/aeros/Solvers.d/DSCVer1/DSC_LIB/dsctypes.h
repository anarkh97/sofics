typedef struct sparseastruc {
		int	global_ns;
		int 	*a_index, *a_struc;
		int     *replication;
		int	*a_nonz_index;
		real_number_type	*a_nonz;
} SPARSE_A_STRUC;

typedef struct sparsea {
		int	*my_n_numbers, *new_n_numbers;
		int	*owned_by, *new_ns_numbers;
		int	*local_old_ns_numbers;
		int 	*a_index,	*a_size, *a_struc;
		real_number_type	*a_nonz,	*b;
} SPARSE_A;

typedef struct mfl {
	int 	**factor_struc,	*factor_struc_sizes, 
   	        *factor_index_list, *factor_nonz_sizes,	
	        **ptr_to_parent_factorstruc,
	        done_factor_columns,	
                done_solution_elements,
	        max_factor_columns;
       real_number_type	**factor_nonz,	*factor_b;   
} MF_L;  

typedef struct mfstack{
        int 	*index_list, *stack_nonz_sizes,
		*stack_map,
 		stack_ptr,	child_ptr,	now_ptr, 
		max_stack_columns,	
		max_mat_size_in_d_phase,
		max_mat_size;
	real_number_type	*LAST_FREE;
	real_number_type	**stack_nonz,	*stack_b;   
	int	LAST_MEM;
} MF_STACK;

typedef struct mftree{
	int 	*chain_index,	*tree_chains,
		*tree_child,	*tree_parent, *tree_sibling,
		*tree_local_column, *o_n_tree_chains,
		*processor_list, 
		*tree_start_procs, *tree_count_procs, *tree_end_procs,
		d_tree_size,	local_phase_root, tree_size;
	real_number_type	*b_tree_chains;
} MF_TREE;


typedef struct sparsetmp{
	int 	
		*tmp_local, *tmp_global,
		*tmp_list, tmp_list_size, max_tmp_list_size,
		*subs_list, 	subs_list_size, max_subs_list_size,
		*factor_struc_received, factor_struc_received_size, 
		max_factor_struc_received_size,
		*send_vector, send_vector_size, max_send_vector_size,
		max_b_size,
		*have_contrib,	*contrib_nonz_sizes,
		have_contrib_next;

	real_number_type	**contrib_nonz,		*contrib_b;
		
	real_number_type	 *col1,	*col2, *col3, *col_zero;
} SPARSE_TMP;
		

typedef struct sparsesizes{
	int	my_world_pid, my_world_P;

	int	*arg_int[DSC_MAX_ARGS],
		my_N,   N, g_N, P,  my_pid, my_NS,
		local_nonz;

	int	scheme_code, number_of_rhs;
        MPI_Group DSC_ALL_GROUP, DSC_SOLVER_GROUP;
	MPI_Comm  DSC_ALL_COMM,  DSC_SOLVER_COMM;
       /*---------------------------------------------------------
	   my_N      =  number of columns owned by this processor
                       \sum_{i=0} ^ {P-1 { my_N}  = g_N 
           my_NS     =  number of distintc sparsity structures 
		        giving rise to my_N columns
	   local_nonz:  nonzeroes over all my_N columns 
                       owned by this processor
           N         = number of columns in local domain + columns 
                      in all separators on path to root; N > my_N
	----------------------------------------------------------*/

		
} SPARSE_SIZES; 
typedef struct sparsestatus{
	int	stage_completed,
		stage_in_progress,
		processor_rank_all_group,
		processor_rank_solver_group,
		solve_pass,
		error_code,
		cont_or_stop;
	char	
		name_stage_completed[80],
		name_stage_in_progress[80],
		name_function_done_last[80],
		name_function_in_progress[80],
		error_message[80];
} SPARSE_STATUS;
		


