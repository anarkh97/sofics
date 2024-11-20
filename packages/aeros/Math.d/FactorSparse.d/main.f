************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*
*   main.f:
*
*     This is a driver that solves a sparse, symmetric, positive
*     semi-definite linear systems via sparse LDL' factorization.
*     The most efficient known algorithms and subroutines are
*     included: the multiple minimum degree subroutines come from
*     the most recent release of SPARSPAK; most of the symbolic
*     factorization subroutines and all of the numerical
*     factorization subroutines were developed at the Oak Ridge
*     National Laboratory.
*
*     ******************************************************************
*     ******************************************************************
*
*     This driver is not a practical instrument for solving a user's
*     sparse symmetric positive semi-definite linear system.  In
*     particular, the driver does not allocate and deallocate memory
*     in an efficient manner.  We emphasize, however, that the
*     subroutines themselves are designed to enable efficient memory
*     management, by either dynamic allocation or by "manual"
*     allocation (and de-allocation) from a single working array.
*
*     This driver and its comments
*
*     (1) indicate for each subroutine which variables and arrays
*         constitute its input and which constitute its output,
*
*     (2) point out when the contents of an array are no longer
*         needed,
*
*     (3) demonstrate appropriate error trapping, and
*
*     (4) give enough information to enable efficient memory
*         management:  prior to each subroutine call, the amount
*         of memory required for each array argument is known.
*
*     Note also that the subroutines used to create the graph of
*     the coeffcient matrix, to input the numerical values, and to
*     generate the right-hand side are not general-purpose tools for
*     these tasks.  They are extremely simple subroutines that the
*     user will need to replace with subroutines of his/her own.
*     The driver reads in the structure of the matrix from an input
*     file, which contains the lower adjacency structure stored in
*     a column-oriented format.
*
*     ******************************************************************
*     ******************************************************************
*
*     The amount of storage required to solve a linear system depends
*     on the following parameters, each of which is either read in from
*     a file or computed during the solution process.
*
*     (1) N      -- the number of rows (columns) in matrix A.
*
*     (2) NSUPER -- the number of supernodes (NSUPER <= N).
*
*     (3) NNZERO -- the number of nonzero entries in the lower
*                   triangular part of A, including entries on the
*                   main diagonal.
*
*     (4) NSUB   -- the number of row subscripts needed to
*                   represent the zero-nonzero structure of the
*                   lower triangular factor L.
*
*     (5) NNZL   -- the number of nonzero entries in L, including
*                   entries on the main diagonal.
*
*     (6) TMPSIZ -- the size of the floating point work array (TMPVEC)
*                   required by the factorization subroutine (BLKLDL).
*
*     With the exception of the IWORK(*) and RWORK(*) arrays, the
*     length of each array should be N+1, NSUPER+1, or one of the
*     preceding six parameters.  In this driver the length of each
*     arrays is declared at run time, which is, of course, before the
*     appropriate length is known.  The following parameters,
*     initialized in the code below, are upper bounds on the seven
*     variables listed above and are used to declare the length of
*     all arrays used by the program.
*
*     (1) NMAX   -- upper bound on N and NSUPER.
*     (2) NZMAX  -- upper bound on NNZERO.
*     (3) SUBMAX -- upper bound on NSUB.
*     (4) LNZMAX -- upper bound on NNZL.
*     (5) TMPMAX -- upper bound on TMPSIZ.
*     (6) IWMAX  -- upper bound on IWSIZ.
*     (7) RWMAX  -- upper bound on RWSIZ.
*
*     The amount of integer work space (IWORK) required by each
*     subroutine varies from one subroutine to another.  The maximum
*     required is 7*N+3.
*
************************************************************************
************************************************************************
*
*       --------------------------------------------------------
*       Parameters ...
*
*           NMAX   -- upper bound on N and NSUPER.
*           NZMAX  -- upper bound on NNZERO.
*           SUBMAX -- upper bound on NSUB.
*           LNZMAX -- upper bound on NNZL.
*           MAXDEF -- upper bound on the size of the null space.
*           TMPMAX -- upper bound on TMPSIZ.
*           IWMAX  -- upper bound on IWSIZ.
*           RWMAX  -- upper bound on RWSIZ.
*       --------------------------------------------------------
        INTEGER             NMAX  , NZMAX
        INTEGER             SUBMAX, LNZMAX
        INTEGER             MAXDEF
        INTEGER             TMPMAX, IWMAX , RWMAX
*
        PARAMETER       (
     &                      NMAX        =    15 000,
     &                      NZMAX       =   300 000,
     &                      LNZMAX      = 1 000 000,
     &                      SUBMAX      =   300 000,
     &                      MAXDEF      =        25,
     &                      TMPMAX      =    50 000,
     &                      IWMAX       =  7*NMAX+3,
     &                      RWMAX       =    TMPMAX
     &                  )
*
*
        DOUBLE PRECISION    ONE   , ZERO
*
        PARAMETER       (
     &                      ONE         = 1.0D0 ,
     &                      ZERO        = 0.0D0
     &                  )
*
************************************************************************
************************************************************************
*
*       ----------------
*       Declarations ...
*       ----------------
*
*       --------------------------------------------------
*       Variables related to the HARWELL-BOEING format ...
*       --------------------------------------------------
        CHARACTER           TITLE*72
        CHARACTER           KEY*8
        CHARACTER           TYPE*3
        CHARACTER           PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20
        INTEGER             N, NCOLS , NNZERO
        INTEGER             NRHS
        INTEGER             PTRCRD, INDCRD, VALCRD, RHSCRD, TOTCRD
*
        INTEGER             COLPTR(NMAX+1)
        INTEGER             ROWIDX(NZMAX)
        DOUBLE PRECISION    VALUES(NZMAX)
*
*       --------------------------------
*       The full adjacency structure ...
*       --------------------------------
        INTEGER             NADJ
        INTEGER             XADJ(NMAX+1), ADJ(2*NZMAX-NMAX)
*
*       -----------------------------------------------
*       The right-hand side vector and the solution ...
*       -----------------------------------------------
        DOUBLE PRECISION    RHS(NMAX), SOLN(NMAX)
        DOUBLE PRECISION    NEWRHS(NMAX)
*
*       ------------------------------------------
*       The permutation vector and its inverse ...
*       ------------------------------------------
        INTEGER             PERM(NMAX), INVP(NMAX)
*
*       -------------------------------------------------
*       Supernode partition and column nonzero counts for
*       the LDL' factorization ...
*       -------------------------------------------------
        INTEGER             NSUPER
        INTEGER             XSUPER(NMAX+1), SNODE(NMAX)
        INTEGER             COLCNT(NMAX)
*
*       ----------------------------------------------------------
*       Compressed (supernodal) representation of the zero-nonzero
*       structure of the LDL' factorization ...
*       ----------------------------------------------------------
        INTEGER             XLINDX(NMAX+1), LINDX(SUBMAX)
        INTEGER             NSUB
*
*       --------------------------------------------------
*       Data structure for the nonzero entries of the LDL'
*       factorization ...
*       --------------------------------------------------
        INTEGER             XLNZ(NMAX+1)
        DOUBLE PRECISION    LNZ(LNZMAX)
        INTEGER             NNZL
*
*       ----------------------------------------------
*       Rank deficiency information and null space ...
*       ----------------------------------------------
        INTEGER             DEFBLK, NDEF  , LBDEF
        INTEGER             DEF(NMAX)
        INTEGER             IPROW(NMAX), IPCOL(NMAX)
        DOUBLE PRECISION    NS(NMAX,MAXDEF)
*
*       ----------------------
*       Various work space ...
*       ----------------------
        INTEGER             IWSIZE
        INTEGER             IWORK(IWMAX)
        INTEGER             RWSIZE
        DOUBLE PRECISION    RWORK(RWMAX)
        INTEGER             TMPSIZ
        DOUBLE PRECISION    TMPVEC(TMPMAX)
*
*       -------------------
*       Other variables ...
*       -------------------
*
        INTEGER             DEBUG
*
        INTEGER             MATUNT, OUTUNT
        CHARACTER*60        MATFIL, OUTFIL
*
        INTEGER             MAXSUP, VALKEY
*
        DOUBLE PRECISION    ANORM , EPS   , TOL
        DOUBLE PRECISION    MXCOMP
        DOUBLE PRECISION    E     , ERROR , RESID
*
        INTEGER             I     , ICASE , IFLAG , LDNS
        DOUBLE PRECISION    T
*
        REAL                TIMBEG, TIMEND, TIME
*
        CHARACTER*40        ORDRNG(2)
        DATA  ORDRNG(1)     / 'NATURAL' /,
     &        ORDRNG(2)     / 'MULTIPLE MINIMUM DEGREE' /
*
*       ----------------------
*       External Functions ...
*       ----------------------
        REAL                GTIMER
        DOUBLE PRECISION    DLAMCH
        EXTERNAL            DLAMCH, GTIMER
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            BFINIT, BLKLDL, BLKNS , BLKSLV, CHKNS ,
     &                      GETADJ, GETNRM, GETRES, GETRHS, INPMAT,
     &                      INPNV , LSTATS, LTFREP, ORDMMD2, ORDNAT,
     &                      SFINIT, SYMFCT
*
*       -----------------------
*       Intrinsic Functions ...
*       -----------------------
        INTRINSIC           ABS
*
************************************************************************
************************************************************************
*
*       ----------------
*       Debugging level.
*       ----------------
        DEBUG = 0
*
*       ---------------------
*       Set I/O unit numbers.
*       ---------------------
        OUTUNT = 6
        MATUNT = 10
*
*       *****************************************
*       Read in problem data from standard input.
*       *****************************************
*
*       -------------------------------------------
*       Read name of output file and open the file.
*       -------------------------------------------
        READ *, OUTFIL
        OPEN  (UNIT=OUTUNT,FILE=OUTFIL)
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  '*** OUTPUT FILE: ', OUTFIL
*
*       -------------------------------------------
*       Read name of matrix file and open the file.
*       -------------------------------------------
        READ *, MATFIL
        OPEN (UNIT=MATUNT,FILE=MATFIL)
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  '*** MATRIX FILE: ', MATFIL
*
*       -------------------------------------------------------
*       Read flag for numerical values.
*       If VALKEY is 0, then read numerical values from MATFIL.
*       Otherwise, INPMAT will generate numerical values.
*       -------------------------------------------------------
        READ *, VALKEY
        WRITE (OUTUNT,*)  ' '
        IF  ( VALKEY .EQ. 0 )  THEN
            WRITE (OUTUNT,*)  '*** NUMERICAL VALUES IN MATRIX FILE'
        ELSE
            WRITE (OUTUNT,*)  '*** NUMERICAL VALUES TO BE GENERATED'
        END IF
*
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)
     &      '----------------------------------------------------------'
*
*       -----------------------------------------------
*       Read choice of ordering and check its validity.
*       -----------------------------------------------
        READ *, ICASE
        IF  ( ICASE .LE. 0  .OR.  ICASE .GE. 3 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,1)  ICASE
    1       FORMAT ( '***** ERROR *****' /
     &               '***** ORDERING OPTION ICASE = ', I5 /
     &               '***** SHOULD HAVE 1 <= ICASE <= 2' )
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,2)  ICASE, ORDRNG(ICASE)
    2   FORMAT ( '*** ORDERING OPTION: ', I2, ' - ', A40 )
*
*       ----------------------------
*       Read maximum supernode size.
*       ----------------------------
        READ *, MAXSUP
        IF  ( MAXSUP .LE. 0 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,3)  MAXSUP
    3       FORMAT ( '***** ERROR *****' /
     &               '***** NEGATIVE MAX SUP SIZE = ', I5, /
     &               '***** SHOULD HAVE 1 <= ICASE' )
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,4)  MAXSUP
    4   FORMAT ( '*** MAXIMUM SUPERNODE SIZE: ', I6 )
*
*       ---------------------------------------------------
*       Read size of last block (for deficiency detection).
*       ---------------------------------------------------
        READ *, DEFBLK
        IF  ( DEFBLK .LT. 0 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,5)  DEFBLK
    5       FORMAT ( '***** ERROR *****' /
     &               '***** NEGATIVE BLOCK SIZE= ', I5 )
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,6)  DEFBLK
    6   FORMAT ( '*** SIZE OF RANK-DEFICIENT BLOCK : ', I6 )
*
************************************************************************
************************************************************************
*
*       *************************************************
*       Read matrix file.
*       Matrix is assumed to be in HARWELL-BOEING format.
*       *************************************************
*
*       -----------------------
*       Get matrix description.
*       -----------------------
        TIMBEG = GTIMER()
        READ (MATUNT,7,END=100)
     &          TITLE , KEY   ,
     &          TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     &          TYPE  , N     , NCOLS , NNZERO, NRHS  ,
     &          PTRFMT, INDFMT, VALFMT, RHSFMT
    7   FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
        GO TO 200
  100   CONTINUE
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  '***** ERROR *****'
        WRITE (OUTUNT,*)  '***** END OF FILE READING MATRIX DESCRIPTION'
        GO TO 300
*
  200   CONTINUE
        WRITE (OUTUNT,8)  TITLE, KEY, N
    8   FORMAT (/ '*** ', A72 / '*** ' A8 / '*** N:',I10)
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)
     &      '----------------------------------------------------------'
*
*       ------------------------------
*       Is matrix of the correct type?
*       ------------------------------
*
        IF  ( TYPE .NE. 'RSA' .AND. TYPE .NE. 'PSA' )  THEN
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)
     &          '***** MATRIX TYPE ', TYPE, ' IS NOT RSA OR PSA'
            GO TO 300
        END IF
*
*       --------------------------------
*       Stop if the matrix is too large.
*       --------------------------------
        IF  ( N .GT. NMAX  .OR.  NNZERO .GT. NZMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,9)  TITLE
    9       FORMAT ( A72 )
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,10)
     &          'NUMBER OF EQUATIONS                   = ', N
   10       FORMAT ( A40, I10 )
            WRITE (OUTUNT,10)
     &          'NUMBER OF NONZEROS (INCLUDING DIAG.)  = ', NNZERO
            WRITE (OUTUNT,10)
     &          'NUMBER OF NONZEROS (EXCLUDING DIAG.)  = ', NNZERO-N
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** MATRIX IS TOO LARGE'
            IF  ( N .GT. NMAX )  THEN
                WRITE (OUTUNT,*)  '***** NUMBER OF EQUATIONS = ', N
                WRITE (OUTUNT,*)  '***** IS LARGER THAN NMAX = ', NMAX
            END IF
            IF  ( NNZERO .GT. NZMAX )  THEN
                WRITE (OUTUNT,*)
     &              '***** NUMBER OF OFF-DIAGONAL NONZEROS (IN A) = ',
     &              NNZERO
                WRITE (OUTUNT,*)
     &              '***** IS LARGER THAN NZMAX = ', NZMAX
            END IF
            GO TO 300
        END IF
*
*       -------------------------------------------------------
*       INPMAT ...  get input matrix.
*
*       Input:      MATUNT, VALKEY, PTRFMT, INDFMT, VALFMT,
*                   N, NNZERO
*       Output:     COLPTR, ROWIDX, VALUES, IFLAG
*       Work:       IWORK(N)
*
*       If VALKEY = 0, then values is read in from matrix file.
*       Otherwise, values is generated by INPMAT.
*       -------------------------------------------------------
        CALL  INPMAT ( MATUNT, VALKEY, PTRFMT, INDFMT, VALFMT,
     &                 N, NNZERO, COLPTR, ROWIDX, VALUES, IWORK,
     &                 IFLAG )
        IF  ( IFLAG .EQ. 1 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** ERROR READING POINTERS'
            GO TO 300
        ELSE IF  ( IFLAG .EQ. 2 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** ERROR READING ROW INDICES'
            GO TO 300
        ELSE IF  ( IFLAG .EQ. 3 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** ERROR READING NUMERICAL VALUES'
            GO TO 300
        END IF
*
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,9)  TITLE
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,10)
     &      'NUMBER OF EQUATIONS                   = ', N
        WRITE (OUTUNT,10)
     &      'NUMBER OF NONZEROS (INCLUDING DIAG.)  = ', NNZERO
        WRITE (OUTUNT,10)
     &      'NUMBER OF NONZEROS (EXCLUDING DIAG.)  = ', NNZERO - N
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR READING THE MATRIX           = ', TIME
   11   FORMAT ( A40, F10.3 )
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING INPMAT ***'
            PRINT *,  'N = ', N
            PRINT *,  'NNZERO = ', NNZERO
        END IF
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  'COLPTR ...'
            PRINT *,  (COLPTR(I),I=1,N+1)
            PRINT *,  'ROWIDX ...'
            PRINT *,  (ROWIDX(I),I=1,NNZERO)
        END IF
        IF  ( DEBUG .GT. 2 )  THEN
            PRINT *,  'VALUES ...'
            PRINT *,  (VALUES(I),I=1,NNZERO)
        END IF
*
*       ---------------------------------------------------------------
*       LTFREP ...  generate a full representation of the input matrix.
*
*       Input:      N, COLPTR, ROWIDX, VALUES
*       Output:     COLPTR, ROWIDX, VALUES
*       Work:       XLINDX, LINDX, LNZ, IWORK(N)
*       ---------------------------------------------------------------
        TIMBEG = GTIMER()
        CALL  LTFREP ( N, COLPTR, ROWIDX, VALUES, XLINDX, LINDX,
     &                 LNZ, IWORK )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR CONVERTING THE MATRIX        = ', TIME
        IF  ( DEBUG .GT. 2 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING LTFREP ***'
            PRINT *,  'COLPTR ...'
            PRINT *,  (COLPTR(I),I=1,N+1)
            PRINT *,  'ROWIDX ...'
            PRINT *,  (ROWIDX(I),I=1,NNZERO)
            PRINT *,  'VALUES ...'
            PRINT *,  (VALUES(I),I=1,NNZERO)
        END IF
*
        NADJ = COLPTR(N+1) - 1 - N
        IF  ( NADJ+N .GT. NZMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** NZMAX HAS TO BE AT LEAST ', NADJ+N
            GO TO 300
        END IF
*       ------------------------------------------------------
*       GETADJ ...  generate full adjacency from matrix input.
*
*       Input:      N, COLPTR, ROWIDX
*       Output:     XADJ, ADJ
*       ------------------------------------------------------
        TIMBEG = GTIMER()
        CALL  GETADJ ( N, COLPTR, ROWIDX, XADJ, ADJ )
        NADJ = XADJ(N+1)-1
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
*
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR GENERATING FULL ADJACENCY    = ', TIME
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING GETADJ ***'
            PRINT *,  'NADJ = ', NADJ
            PRINT *,  'XADJ ...'
            PRINT *,  (XADJ(I),I=1,N+1)
            PRINT *,  'ADJ ...'
            PRINT *,  (ADJ(I),I=1,XADJ(N+1)-1)
        END IF
*
************************************************************************
************************************************************************
*
*       *****************************************************
*       Compute true solution and associated right-hand side.
*       *****************************************************
*
*       -----------------------------
*       Construct the exact solution.
*       -----------------------------
        TIMBEG = GTIMER()
        T = ONE
        DO  I = 1, N
            SOLN(I) = T
            T       = - T
        END DO
*       -------------------------------------------------------
*       GETRHS ...  construct right hand side vector associated
*                   with solution in SOLN.
*
*       Input:      N, COLPTR, ROWIDX, VALUES, SOLN
*       Output:     RHS
*       -------------------------------------------------------
        CALL  GETRHS ( N, COLPTR, ROWIDX, VALUES, SOLN, RHS )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR CONSTRUCTING RHS             = ', TIME
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING GETRHS ***'
            PRINT *,  'SOLN ...'
            PRINT *,  (SOLN(I),I=1,N)
            PRINT *,  'RHS ...'
            PRINT *,  (RHS(I),I=1,N)
        END IF
*
*       ----------------------------------------------------
*       GETNRM ...  compute L-infinity norm of input matrix.
*
*       Input:      N, COLPTR, ROWIDX, VALUES
*       Output:     ANORM
*       Work:       RWORK(N)
*       ----------------------------------------------------
        TIMBEG = GTIMER()
        CALL  GETNRM ( N, COLPTR, ROWIDX, VALUES, RWORK, ANORM )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,12)
     &      'NORM OF INPUT MATRIX                  = ', ANORM
   12   FORMAT ( A40, 1PD20.10 )
        WRITE (OUTUNT,11)
     &      'TIME FOR COMPUTING NORM OF MATRIX     = ', TIME
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING GETNRM ***'
            PRINT *,  'NORM OF MATRIX = ', ANORM
        ENDIF
*
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)
     &      '----------------------------------------------------------'
*
************************************************************************
************************************************************************
*
*       *******************
*       Reorder the matrix.
*       *******************
*
        IF  ( ICASE .EQ. 1 )  THEN
*
*           ---------------------------------------------------
*           ORDNAT ...  natural ordering (i.e., compute and use
*                       the identity permutation.)
*
*           Input:      N
*           Output:     INVP, PERM
*           ---------------------------------------------------
            TIMBEG = GTIMER()
            CALL  ORDNAT ( N, PERM, INVP )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,11)
     &          'TIME FOR ORDERING                     = ', TIME
*
        ELSE IF  ( ICASE .EQ. 2 )  THEN
*
*           -------------------------------------------
*           Copy matrix structure from (XADJ,ADJ) to
*           (XLINDX,LINDX) (because matrix structure is
*           destroyed by the minimum degree ordering
*           subroutine).
*           -------------------------------------------
            TIMBEG = GTIMER()
            DO  I = 1, N+1
                XLINDX(I) = XADJ(I)
            END DO
            DO  I = 1, NADJ
                LINDX(I) = ADJ(I)
            END DO
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,11)
     &          'TIME FOR COPYING ADJACENCY STRUCT.    = ', TIME
*
*           ---------------------------------------------
*           ORDMMD2 ...  multiple minimum degree ordering.
*
*           Input:      N, XLINDX, LINDX, IWSIZE
*           Output:     PERM, INVP, NSUB, IFLAG
*           Work:       IWORK(4*N)
*           ---------------------------------------------
            IWSIZE = 4*N
            TIMBEG = GTIMER()
            CALL  ORDMMD2 ( N, XLINDX, LINDX, INVP, PERM,
     &                     IWSIZE, IWORK, NSUB, IFLAG )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            IF  ( IFLAG .EQ. 11 )  THEN
                WRITE (OUTUNT,*)  ' '
                WRITE (OUTUNT,*)  '***** ERROR *****'
                WRITE (OUTUNT,*)  '***** SIZE OF IWORK = ', IWSIZE
                WRITE (OUTUNT,*)  '***** IS LARGER THAN IWMAX = ',
     &                            IWMAX
                GO TO 300
            END IF
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,11)
     &          'TIME FOR ORDERING                     = ', TIME
            IF  ( DEBUG .GT. 1 )  THEN
                PRINT *,  ' '
                PRINT *,  'AFTER ORDERING ***'
                PRINT *,  'PERM ...'
                PRINT *,  (PERM(I),I=1,N)
                PRINT *,  'INVP ...'
                PRINT *,  (INVP(I),I=1,N)
            END IF
*
        END IF
*
************************************************************************
************************************************************************
*
*       ***********************
*       Symbolic factorization.
*       ***********************
*
*       -------------------------------------------------------------
*       SFINIT ...  symbolic factorization initialization, which
*                   computes supernode partition and storage
*                   requirements for symbolic factorization;
*                   new ordering is a postordering of the nodal
*                   elimination tree.
*
*       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, MAXSUP, DEFBLK,
*                   IWSIZE
*       Output:     PERM, INVP, COLCNT, NNZL, NSUB, NSUPER, XSUPER,
*                   SNODE, IFLAG
*       Work:       IWORK(7*N+3) ... the max required any subroutine.
*       -------------------------------------------------------------
        IWSIZE = 7*N + 3
        TIMBEG = GTIMER()
        CALL  SFINIT ( N, NADJ, XADJ, ADJ, PERM, INVP, MAXSUP, DEFBLK,
     &                 COLCNT, NNZL, NSUB, NSUPER, XSUPER, SNODE,
     &                 IWSIZE, IWORK, IFLAG )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. 21 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** SIZE OF IWORK = ', IWSIZE
            WRITE (OUTUNT,*)  '***** IS LARGER THAN IWMAX = ', IWMAX
            GO TO 300
        END IF
        IF  ( NNZL .GT. LNZMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** NUMBER OF NONZEROS IN L = ', NNZL
            WRITE (OUTUNT,*)  '***** IS LARGER THAN LNZMAX = ', LNZMAX
            GO TO 300
        END IF
        IF  ( NSUB .GT. SUBMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** NUMBER OF FACTOR SUBSCRIPTS = ',
     &                        NSUB
            WRITE (OUTUNT,*)  '***** IS LARGER THAN SUBMAX = ', SUBMAX
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,10)
     &      'NUMBER OF SUPERNODES                  = ', NSUPER
        WRITE (OUTUNT,10)
     &      'NUMBER OF COMPRESSED INDICES          = ', NSUB
        WRITE (OUTUNT,10)
     &      'NUMBER OF NONZEROS IN L               = ', NNZL
        WRITE (OUTUNT,11)
     &      'TIME FOR SYMBOLIC FACT. SETUP         = ', TIME
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING SFINIT ***'
            PRINT *,  'NSUPER = ', NSUPER
            PRINT *,  'NNZL = ', NNZL
            PRINT *,  'NSUB = ', NSUB
        END IF
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  'PERM ...'
            PRINT *,  (PERM(I),I=1,N)
            PRINT *,  'INVP ...'
            PRINT *,  (INVP(I),I=1,N)
            PRINT *,  'COLCNT ...'
            PRINT *,  (COLCNT(I),I=1,N)
            PRINT *,  'XSUPER ...'
            PRINT *,  (XSUPER(I),I=1,NSUPER+1)
            PRINT *,  'SNODE ...'
            PRINT *,  (SNODE(I),I=1,N)
        END IF
*
*       ------------------------------------------------------
*       SYMFCT ...  perform supernodal symbolic factorization.
*
*       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, COLCNT,
*                   NSUPER, XSUPER, SNODE , NSUB, IWSIZE
*       Output:     XLINDX, LINDX, XLNZ, IFLAG
*       WORK:       IWORK(NSUPER+2*N+1)
*
*       No longer needed: ADJ, XADJ, COLCNT
*       ------------------------------------------------------
        IWSIZE = NSUPER + 2*N + 1
        TIMBEG = GTIMER()
        CALL  SYMFCT ( N, NADJ, XADJ, ADJ, PERM, INVP, COLCNT,
     &                 NSUPER, XSUPER, SNODE, NSUB, XLINDX, LINDX,
     &                 XLNZ, IWSIZE, IWORK, IFLAG )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. 22 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** SIZE OF IWORK = ', IWSIZE
            WRITE (OUTUNT,*)  '***** IS LARGER THAN IWMAX = ', IWMAX
            GO TO 300
        END IF
        IF  ( IFLAG .EQ. 23 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)
     &          '***** INCONSISTENCY IN THE INPUT TO SYMFCT'
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR SYMBOLIC FACTORIZATION       = ', TIME
        IF  ( DEBUG .GT. 2 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING SYMFCT ***'
            PRINT *,  'XLINDX ...'
            PRINT *,  (XLINDX(I),I=1,NSUPER+1)
            PRINT *,  'LINDX ...'
            PRINT *,  (LINDX(I),I=1,XLINDX(NSUPER+1)-1)
            PRINT *,  'XLNZ ...'
            PRINT *,  (XLNZ(I),I=1,N+1)
        END IF
*
************************************************************************
************************************************************************
*
*       ***************************************************
*       Numerical input into data structure for sparse LDL'
*       factorization.
*       ***************************************************
*
*       --------------------------------------------------------
*       INPNV ...   input numerical values into data structures.
*
*       Input:      N, COLPTR, ROWIDX, VALUES, PERM, INVP,
*                   NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
*       Output:     LNZ
*       Work:       IWORK(N)
*       --------------------------------------------------------
        IWSIZE = N
        TIMBEG = GTIMER()
        CALL  INPNV ( N, COLPTR, ROWIDX, VALUES, PERM, INVP,
     &                NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
     &                LNZ, IWORK )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR NUMERICAL INPUT              = ', TIME
        IF  ( DEBUG .GT. 2 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING INPNV ***'
            PRINT *,  'LNZ ...'
            PRINT *,  (LNZ(I),I=1,XLNZ(N+1)-1)
        END IF
*
************************************************************************
************************************************************************
*
*       ************************
*       Numerical factorization.
*       ************************
*
*       ---------------------------------------------------
*       BFINIT ...  initialization for block factorization.
*
*       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX
*       Output:     TMPSIZ, RWSIZE
*       ---------------------------------------------------
        TIMBEG = GTIMER()
        CALL  BFINIT ( NSUPER, XSUPER, SNODE, XLINDX, LINDX,
     &                 TMPSIZ, RWSIZE )
        TMPSIZ = 2*TMPSIZ
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( TMPSIZ .GT. TMPMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** WORKING SPACE REQUIREMENT = ',
     &                        TMPSIZ
            WRITE (OUTUNT,*)  '***** IS LARGER THAN TMPMAX = ', TMPMAX
            GO TO 300
        END IF
        IF  ( RWSIZE .GT. RWMAX )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** WORKING SPACE REQUIREMENT = ',
     &                        RWSIZE
            WRITE (OUTUNT,*)  '***** IS LARGER THAN RWMAX = ', RWMAX
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR FACTORIZATION INIT.          = ', TIME
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING BFINIT ***'
            PRINT *,  'TMPSIZ = ', TMPSIZ
            PRINT *,  'RWSIZE = ', RWSIZE
        END IF
*
*       -------------------------------------------------------
*       BLKLDL ...  numerical factorization.
*
*       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX, XLNZ,
*                   LNZ, DEFBLK, TOL, TMPSIZ, IWSIZE, RWSIZE
*       Output:     LNZ, NDEF, LBDEF, DEF, IPROW, IPCOL, IFLAG
*       Work:       TMPVEC(TMPSIZ), IWORK(2*N+2*NSUPER),
*                   RWORK(RWSIZE)
*       -------------------------------------------------------
        EPS = DLAMCH ( 'EPS' )
        TOL = EPS*ANORM
        IWSIZE = 3*N + 2*NSUPER
        TIMBEG = GTIMER()
        CALL  BLKLDL ( NSUPER, XSUPER, SNODE, XLINDX, LINDX,
     &                 XLNZ, LNZ, DEFBLK, NDEF, LBDEF, DEF, TOL,
     &                 IPROW, IPCOL, TMPSIZ, TMPVEC, IWSIZE, IWORK,
     &                 RWSIZE, RWORK, IFLAG )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        IF  ( IFLAG .EQ. 31 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** INSUFFICIENT WORK SPACE IN TMPVEC'
            GO TO 300
        ELSEIF ( IFLAG .EQ. 32 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** INSUFFICIENT WORK SPACE IN IWORK'
            GO TO 300
        ELSEIF ( IFLAG .EQ. 33 )  THEN
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,*)  '***** ERROR *****'
            WRITE (OUTUNT,*)  '***** INSUFFICIENT WORK SPACE IN RWORK'
            GO TO 300
        END IF
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,12)
     &      'TOLERANCE FOR DETERMINING DEFICIENCY  = ', TOL
        WRITE (OUTUNT,10)
     &      'RANK DEFICIENCY                       = ', NDEF
        WRITE (OUTUNT,10)
     &      'RANK DEFICIENCY OF LAST PANEL         = ', LBDEF
        WRITE (OUTUNT,11)
     &      'TIME FOR NUMERICAL FACTORIZATION      = ', TIME
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING BLKLDL ***'
            PRINT *,  'TOL = ', TOL
            PRINT *,  'NDEF = ', NDEF
            PRINT *,  'LBDEF = ', LBDEF
            PRINT *,  'DEF ...'
            PRINT *,  (DEF(I),I=1,NDEF)
        END IF
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  'IPROW ...'
            PRINT *,  (IPROW(I),I=1,DEFBLK)
            PRINT *,  'IPCOL ...'
            PRINT *,  (IPCOL(I),I=1,DEFBLK)
        END IF
*
************************************************************************
************************************************************************
*
*       *******************
*       Compute null space.
*       *******************
*
        IF  ( NDEF .NE. 0  .AND.  NDEF .LE. MAXDEF )  THEN
*           -------------------------------------------------------
*           BLKNS ...   compute the null space.
*
*           Input:      NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ,
*                       DEFBLK, NDEF, LBDEF, DEF, IPCOL, INVP, LDNS
*           Output:     NS
*           Work:       RWORK(N)
*           -------------------------------------------------------
            LDNS = NMAX
            TIMBEG = GTIMER()
            CALL  BLKNS ( NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ,
     &                    DEFBLK, NDEF, LBDEF, DEF, IPCOL , INVP, NS,
     &                    LDNS, RWORK )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,11)
     &          'TIME FOR COMPUTING NULL SPACE         = ', TIME
*
*           -----------------------------------------------------
*           CHKNS ...   check validity of null space.
*
*           Input:      N, COLPTR, ROWIDX, VALUES, NDEF, NS, LDNS
*           Output:     MXCOMP
*           Work:       RWORK(N)
*           -----------------------------------------------------
            TIMBEG = GTIMER()
            CALL  CHKNS ( N, COLPTR, ROWIDX, VALUES, NDEF, NS, LDNS,
     &                    MXCOMP, RWORK )
            TIMEND = GTIMER()
            TIME = TIMEND - TIMBEG
            WRITE (OUTUNT,*)  ' '
            WRITE (OUTUNT,12)
     &          'MAX ENTRY IN MATRIX * NULL SPACE      = ', MXCOMP
            WRITE (OUTUNT,11)
     &          'TIME FOR CHECKING NULL SPACE          = ', TIME
            IF  ( DEBUG .GT. 0 )  THEN
                PRINT *,  ' '
                PRINT *,  'AFTER CALLING BLKNS AND CHKNS ***'
                PRINT *,  'MXCOMP = ', MXCOMP
            END IF
        END IF
*
************************************************************************
************************************************************************
*
*       ********************
*       Triangular solution.
*       ********************
*
*       --------------------------------------------------------------
*       BLKSLV ...  numerical solution.
*
*       Input:      NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ, DEFBLK,
*                   NDEF, LBDEF, DEF, IPROW, IPCOL, PERM, INVP, NEWRHS
*       Output:     SOLN
*       Work:       RWORK(N)
*       --------------------------------------------------------------
        DO  I = 1, N
            NEWRHS(I) = RHS(I)
            RWORK(I) = ZERO
        END DO
        TIMBEG = GTIMER()
        CALL  BLKSLV ( NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ,
     &                 DEFBLK, NDEF, LBDEF, DEF, IPROW, IPCOL, PERM,
     &                 INVP, NEWRHS, SOLN, RWORK )
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,11)
     &      'TIME FOR TRIANGULAR SOLUTIONS         = ', TIME
        IF  ( DEBUG .GT. 1 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING BLKSLV ***'
            PRINT *,  'SOLN ...'
            PRINT *,  (SOLN(I),I=1,N)
        END IF
*
*       ------------------------------------------------
*       GETRES ...  compute residual in solution.
*
*       Input:      N, COLPTR, ROWIDX, VALUES, RHS, SOLN
*       Output:     RESID
*       Work:       RWORK(N)
*       ------------------------------------------------
        TIMBEG = GTIMER()
        CALL  GETRES ( N, COLPTR, ROWIDX, VALUES, RHS, SOLN, RESID,
     &                 RWORK )
*       --------------------------------------------------
*       Compute error (meaningless if a is semi-definite).
*       --------------------------------------------------
        ERROR = ZERO
        T = ONE
        DO  I = 1, N
            E = ABS(SOLN(I) - T)
            T = - T
            IF  ( E .GT. ERROR )  ERROR = E
        END DO
        TIMEND = GTIMER()
        TIME = TIMEND - TIMBEG
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,12)
     &      'RESIDUAL                              = ', RESID
        WRITE (OUTUNT,12)
     &      'MAXIMUM RELATIVE ERROR                = ', ERROR
        WRITE (OUTUNT,11)
     &      'TIME FOR COMPUTING RESIDUAL AND ERROR = ', TIME
        IF  ( DEBUG .GT. 0 )  THEN
            PRINT *,  ' '
            PRINT *,  'AFTER CALLING GETRES ***'
            PRINT *,  'RESID = ', RESID
            PRINT *,  'ERROR = ', ERROR
        END IF
*
************************************************************************
************************************************************************
*
*       ------------------------------------------------------------
*       LSTATS ...  compute and print statistics.
*
*       Input:      NSUPER, XSUPER, XLINDX, XLNZ, TMPSIZ, RWSIZE,
*                   OUTUNT
*       ------------------------------------------------------------
        CALL  LSTATS ( NSUPER, XSUPER, XLINDX, XLNZ, TMPSIZ, RWSIZE,
     *                 OUTUNT )
*
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)
     &      '----------------------------------------------------------'
*
*       -------------------
*       NORMAL TERMINATION.
*       -------------------
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  '***** NORMAL TERMINATION *****'
        STOP
*
************************************************************************
************************************************************************
*
*       ---------------------
*       ABNORMAL TERMINATION.
*       ---------------------
  300   CONTINUE
        WRITE (OUTUNT,*)  ' '
        WRITE (OUTUNT,*)  '***** ABNORMAL TERMINATION *****'
        STOP
*
************************************************************************
************************************************************************
*
      END
