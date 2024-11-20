extern	int     DSC_Analyze (int , int *, int *, int *);

extern  void    DSC_Check_Global_Known_X(char *, int , int *, real_number_type *,
		real_number_type, int);

extern  void    DSC_Check_Local_Known_X(char *, real_number_type , int );

extern	void	DSC_Clean_Up(void);

extern 	double 	DSC_Clock0(void);

extern 	void   	DSC_Do_Stats(char *, int);

extern	void	DSC_Error_Display(void);

extern	void    DSC_Final_Free_All(void);

extern	int 	DSC_Input_Rhs_Global_Vec(real_number_type *, int);

extern  int     DSC_Input_Rhs_Local_Vec(real_number_type *, int );

extern  int     DSC_N_Fact (int , real_number_type *);

extern	int  	DSC_N_Solve(int);

extern	void    DSC_Open0(int, int *, MPI_Comm);

extern  int  	DSC_Order_A_Struct (int , int *, int *,  int *,
			int *, int *, int *, int *,  int **, int **, int **, int **);

extern  int     DSC_Order_S_N_Fact (int, int, int *, int *, int *, int *, real_number_type *);

extern	void 	DSC_Re_Init();

extern  int     DSC_Remove_Null_Structs(int *, int *, int *, int *);

extern  int     DSC_S_Fact (int *, int *, int );

extern  int     DSC_S_N_Fact (int, real_number_type *);

extern	void  	DSC_Set_Error(char *, int);

extern	int     DSC_Solve_Gather (int *, real_number_type *);

extern 	void 	DSC_Sync0(void);

